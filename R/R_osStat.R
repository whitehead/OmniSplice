library(aod)
library(data.table)
library(parallel)


pass_min_read <- function(d, min_cover, min_unspliced, do_ambi) {
  g1_succ <- sum(d$y[d$grp == "control"])
  g1_fail <- sum(d$fail[d$grp == "control"])
  g2_succ <- sum(d$y[d$grp == "treatment"])
  g2_fail <- sum(d$fail[d$grp == "treatment"])

  if ((d$ambigious[1] == "true") && (do_ambi == "true")) {
    return(FALSE)
  }
  if ((g1_fail + g1_succ) < min_cover || (g2_fail + g2_succ) < min_cover) {
    return(FALSE)
  }
  if (g1_fail < min_unspliced && g2_fail < min_unspliced) {
    return(FALSE)
  }
  TRUE
}

test_junction_betabin <- function(success, failure, group,
                                  group_levels = c("control", "treatment")) {
  res <- tryCatch({
    group <- factor(group, levels = group_levels)
    df <- data.frame(success = success, failure = failure, group = group)
    
    null_fit <- betabin(cbind(success, failure) ~ 1, ~ 1, data = df, warnings = FALSE)
    full_fit <- betabin(cbind(success, failure) ~ group, ~ 1, data = df, warnings = FALSE)
    
    lr_stat <- deviance(null_fit) - deviance(full_fit)
    p_value <- pchisq(lr_stat, df = 1, lower.tail = FALSE)
    
    coefs <- summary(full_fit)@Coef
    beta_group <- coefs[2, "Estimate"]
    se_group   <- coefs[2, "Std. Error"]
    
    odds_ratio <- exp(beta_group)
    ci_lower   <- exp(beta_group - 1.96 * se_group)
    ci_upper   <- exp(beta_group + 1.96 * se_group)
    
    phi_aod <- summary(full_fit)@Phi[1, "Estimate"]
    theta <- (1 - phi_aod) / phi_aod
    
    list(p_value = p_value, odds_ratio = odds_ratio,
         ci_lower = ci_lower, ci_upper = ci_upper,
         theta = theta, phi_aod = phi_aod)
  }, error = function(e) {
    list(p_value = NaN, odds_ratio = NaN,
         ci_lower = NaN, ci_upper = NaN,
         theta = NaN, phi_aod = NaN)
  })
  
  data.frame(p_value_betabin = res$p_value, odds_ratio_betabin = res$odds_ratio,
             ci_lower_betabin = res$ci_lower, ci_upper_betabin = res$ci_upper,
             theta_betabin = res$theta, phi_aod_betabin = res$phi_aod)
}


run_betabin_row <- function(d) {
  out <- test_junction_betabin(d$y, d$fail, d$grp,
                               group_levels = c("control", "treatment"))
  out$row_id <- d$row_id[1]
  out
}

run_ttest <- function(d) {
  res <- tryCatch({
    prop <- d$y / (d$y + d$fail)
    ctrl_vals <- prop[d$grp == "control"]
    trt_vals  <- prop[d$grp == "treatment"]
    tt <- t.test(ctrl_vals, trt_vals)   # Welch's by default (var.equal = FALSE)
    tt$p.value
  }, error = function(e) NaN)
  data.frame(row_id = d$row_id[1], p_value_ttest = res)
}

test_junction_fisher <- function(success, failure, group,
                                 group_levels = c("control", "treatment")) {
  res <- tryCatch({
    group <- factor(group, levels = group_levels)

    g1_succ <- sum(success[group == group_levels[1]])
    g1_fail <- sum(failure[group == group_levels[1]])
    g2_succ <- sum(success[group == group_levels[2]])
    g2_fail <- sum(failure[group == group_levels[2]])

    tab <- matrix(c(g1_succ, g1_fail, g2_succ, g2_fail), nrow = 2, byrow = TRUE)
    ft <- fisher.test(tab)

    list(p_value = ft$p.value, odds_ratio = unname(ft$estimate),
         ci_lower = ft$conf.int[1], ci_upper = ft$conf.int[2])
  }, error = function(e) {
    list(p_value = NaN, odds_ratio = NaN, ci_lower = NaN, ci_upper = NaN)
  })

  data.frame(p_value_fisher = res$p_value, odds_ratio_fisher = res$odds_ratio,
             ci_lower_fisher = res$ci_lower, ci_upper_fisher = res$ci_upper)
}

run_fisher_row <- function(d) {
  out <- test_junction_fisher(d$y, d$fail, d$grp,
                              group_levels = c("control", "treatment"))
  out$row_id <- d$row_id[1]
  out
}

test_junction_logit <- function(success, failure, group,
                                group_levels = c("control", "treatment")) {
  res <- tryCatch({
    group <- factor(group, levels = group_levels)
    df <- data.frame(success = success, failure = failure, group = group)

    null_fit <- glm(cbind(success, failure) ~ 1, data = df, family = binomial)
    full_fit <- glm(cbind(success, failure) ~ group, data = df, family = binomial)

    lr_stat <- deviance(null_fit) - deviance(full_fit)
    p_value <- pchisq(lr_stat, df = 1, lower.tail = FALSE)

    coefs <- summary(full_fit)$coefficients
    beta_group <- coefs[2, "Estimate"]
    se_group   <- coefs[2, "Std. Error"]

    odds_ratio <- exp(beta_group)
    ci_lower   <- exp(beta_group - 1.96 * se_group)
    ci_upper   <- exp(beta_group + 1.96 * se_group)

    list(p_value = p_value, odds_ratio = odds_ratio,
         ci_lower = ci_lower, ci_upper = ci_upper)
  }, error = function(e) {
    list(p_value = NaN, odds_ratio = NaN, ci_lower = NaN, ci_upper = NaN)
  })

  data.frame(p_value_logit = res$p_value, odds_ratio_logit = res$odds_ratio,
             ci_lower_logit = res$ci_lower, ci_upper_logit = res$ci_upper)
}

run_logit_row <- function(d) {
  out <- test_junction_logit(d$y, d$fail, d$grp,
                             group_levels = c("control", "treatment"))
  out$row_id <- d$row_id[1]
  out
}


reader <- function(dt) {
  
  # unique id per original row, so we can group back later
  
  # split the comma-lists into numeric vectors
  cs <- strsplit(dt$control_success, ",")
  cf <- strsplit(dt$control_failures, ",")
  ts <- strsplit(dt$treatment_success, ",")
  tf <- strsplit(dt$treatment_failures, ",")
  
  # build one long data.frame, row by row, but efficiently with rbindlist
  long_dt <- rbindlist(lapply(seq_len(nrow(dt)), function(i) {
    data.frame(
      row_id = dt$row_id[i],
      y      = as.numeric(c(cs[[i]], ts[[i]])),
      fail   = as.numeric(c(cf[[i]], tf[[i]])),
      grp    = c(rep("control", length(cs[[i]])), rep("treatment", length(ts[[i]]))),
      ambigious = dt$ambigious[i]
    )
  }))
}



args <- commandArgs(trailingOnly = TRUE)
# Check if arguments are provided
if (length(args) == 0) {
  stop("No arguments provided")
}

# Print the arguments
print(args)
# Example: Access individual arguments
myfile <- args[1]
outfile <- args[2]
thread <- args[3]
min_cover <- args[4]
min_unspliced <- args[5]
low_repl <- args[6]
do_ambi <- args[7]


all_lines <- readLines(myfile)
is_comment <- grepl("^#", all_lines)
comment_lines <- all_lines[is_comment]
n_comment <- length(comment_lines)

dt <- fread(myfile)
dt[, row_id := .I]  

long_dt = reader(dt)
split_list <- split(long_dt, long_dt$row_id)

filter_flags <- rbindlist(lapply(split_list, function(d) {
  data.frame(row_id = d$row_id[1],
             pass_filter = pass_min_read(d, min_cover, min_unspliced, do_ambi))
}))

pass_ids <- filter_flags[pass_filter == TRUE, row_id]
split_pass <- split_list[as.character(pass_ids)]


n_cores <- thread
all_ids <- data.table(row_id = dt$row_id)

if (low_repl == "true") {

  fisher_results <- rbindlist(mclapply(split_pass, run_fisher_row, mc.cores = n_cores))
  logit_results  <- rbindlist(mclapply(split_pass, run_logit_row,  mc.cores = n_cores))

  results <- merge(fisher_results, logit_results, by = "row_id")

  results_full <- merge(all_ids, results, by = "row_id", all.x = TRUE)
  results_full <- merge(results_full, filter_flags, by = "row_id")

  results_full[, padj_fisher := p.adjust(p_value_fisher, method = "BH")]  # NAs auto-excluded from adjustment
  results_full[, padj_logit  := p.adjust(p_value_logit,  method = "BH")]

  dt <- merge(dt, results_full, by = "row_id")
  dt <- dt[order(dt$padj_fisher), ]

} else {

  bb_results <- rbindlist(mclapply(split_pass, run_betabin_row, mc.cores = n_cores))
  t_results  <- rbindlist(mclapply(split_pass, run_ttest,       mc.cores = n_cores))

  results <- merge(bb_results, t_results, by = "row_id")

  results_full <- merge(all_ids, results, by = "row_id", all.x = TRUE)
  results_full <- merge(results_full, filter_flags, by = "row_id")

  results_full[, padj_betabin := p.adjust(p_value_betabin, method = "BH")]  # NAs auto-excluded from adjustment
  results_full[, padj_ttest   := p.adjust(p_value_ttest,   method = "BH")]

  dt <- merge(dt, results_full, by = "row_id")
  dt <- dt[order(dt$padj_betabin), ]

}

# keep all p-value-adjusted columns grouped at the end, gene_transcript_intron last of all
padj_cols  <- grep("^padj_", names(dt), value = TRUE)
other_cols <- setdiff(names(dt), c(padj_cols, "gene_transcript_intron"))
setcolorder(dt, c(other_cols, padj_cols, "gene_transcript_intron"))
dt[, row_id := NULL]
# write comments first (raw text, no quoting/escaping)
writeLines(comment_lines, outfile)

# then append the data table in tsv format
fwrite(dt, outfile, sep = "\t", append = TRUE, col.names = TRUE, na = "NA",quote = FALSE)
 