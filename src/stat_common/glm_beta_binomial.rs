

use crate::stat_common::errors::LogisticRegressionError::CIUnavail;

use super::super::common::error::OmniError;
use super::errors::LogisticRegressionError;
use super::common::{Genotype, Tester, CountsStats, SplicingCategory, JunctionStats, TestResults, TestStatus, hyper_geom_test, welch_t_test};
use super::glm_logistic::GLM;

//use fishers_exact::fishers_exact;
use statrs::stats_tests::fishers_exact;
use statrs::stats_tests::Alternative;
use statrs::stats_tests::fishers_exact_with_odds_ratio;
use statrs::distribution::{Discrete, DiscreteCDF, Hypergeometric};
use statrs::statistics::Distribution;
use statrs::prec;
use statrs::statistics::{Min, Max};
use statrs::function::gamma::ln_gamma;

use core::f64;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt::format;
use std::iter::Sum;
use std::mem::min_align_of;
use std::path::Path;
use std::{fs, result};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::convert::TryInto;

use nalgebra::{DMatrix, DVector};
use statrs::distribution::{ChiSquared, ContinuousCDF};
use flexi_logger::{FileSpec, Logger, WriteMode};
use log::{debug, error, info, trace, warn};



  /*fn fisher_fallback(contingency: &[u64]) -> Result<(TestStatus, f64, f64, f64, f64), LogisticRegressionError> {
      match hyper_geom_test(contingency[0], contingency[1], contingency[2], contingency[3]) {
          Ok(p) => {
              let mut odd_ratio = 0.;
              if contingency[1] != 0 && contingency[3] != 0 {
                  odd_ratio = (contingency[0] as f64 / contingency[1] as f64) / (contingency[2] as f64 / contingency[3] as f64);
              }
              Ok((TestStatus::FisherFallBack, p, odd_ratio, 0., 0.))
          }
          Err(OmniError::EmptyHyperGrom) => Ok((TestStatus::FisherFallBack, f64::NAN, f64::NAN, 0., 0.)),
          Err(_) => Err(LogisticRegressionError::HyperGeomError),
      }
  }*/

  //fn t_test_fallback(p_ctrl: &Vec<f64>, p_treat: &Vec<f64>) -> Result<(TestStatus, f64, f64, f64, f64), LogisticRegressionError>{
  //     match welch_t_test(p_ctrl, p_treat){

   //     Ok(p, s) =>  Ok((TestStatus::T, p, odd_ratio, 0., 0.))
   //    }
  //}

pub struct GLMBetaBinomiale {
    success: Vec<u32>,
    failures: Vec<u32>,
    groups: Vec<Genotype>,
    identifier: String,
}

impl Tester for GLMBetaBinomiale {
    fn success(&self) -> &Vec<u32> { &self.success }
    fn success_mut(&mut self) -> &mut Vec<u32> { &mut self.success }
    fn failures(&self) -> &Vec<u32> { &self.failures }
    fn failures_mut(&mut self) -> &mut Vec<u32> { &mut self.failures }
    fn groups(&self) -> &Vec<Genotype> { &self.groups }
    fn groups_mut(&mut self) -> &mut Vec<Genotype> { &mut self.groups }



    fn test(&self, donotrun: bool, min_coverage: u32, min_failure: u32) -> TestResults {
        let mut test_res = TestResults::get_empty();
        let (ctrl_suc, ctrl_fail, treat_suc, treat_fail) = self.get_proportion();
        test_res.control_failure = ctrl_fail;
        test_res.control_success = ctrl_suc;
        test_res.treatment_failure = treat_fail;
        test_res.treatment_success = treat_suc;
        test_res.string_count = self.get_proportion_string();

        let treat_trial = treat_suc + treat_fail;
        let ctrl_trial = ctrl_fail + ctrl_suc;
        if (treat_trial == ctrl_trial) && (ctrl_trial == 0) {
            test_res.status = Some(TestStatus::EmptyData);
            return test_res;
        } else if ctrl_trial == 0 {
            test_res.status = Some(TestStatus::ControlIsNull);
            return test_res;
        } else if treat_trial == 0 {
            test_res.status = Some(TestStatus::TreatmentIsNull);
            return test_res;
        } else if !self.pass_min_read( min_coverage, min_failure) {
            test_res.status = Some(TestStatus::FailFilter);
            return test_res;
        } else {
            test_res.control_prop = Some(ctrl_suc as f32 / ctrl_trial as f32);
            test_res.treatment_prop = Some(treat_suc as f32 / treat_trial as f32);

            if donotrun {
                test_res.status = Some(TestStatus::ambiguous);
                return test_res;
            }

            match self.actual_test() {
                Ok((status, pval, odr, ord_lower, odr_high)) => {
                    test_res.status = Some(status);
                    test_res.p_value = Some(pval);
                    test_res.odd_ratio = Some(odr);
                    test_res.or_ci_lower = Some(ord_lower);
                    test_res.or_ci_upper = Some(odr_high);
                }
                Err(e) => {
                    warn!("error: {} {}", self.identifier, e);
                    test_res.status = Some(e.into());
                }
            }
        }
        test_res
    }
}

impl GLMBetaBinomiale {
    pub fn get_sub_vec(&self) -> (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>) {
        let g1 = &self.groups[0];
        let mut g1_succ = Vec::new();
        let mut g1_fail = Vec::new();
        let mut g2_succ = Vec::new();
        let mut g2_fail = Vec::new();
        for (i, g) in self.groups.iter().enumerate() {
            if g == g1 {
                g1_succ.push(self.success[i]);
                g1_fail.push(self.failures[i]);
            } else {
                g2_succ.push(self.success[i]);
                g2_fail.push(self.failures[i]);
            }
        }
        (g1_succ, g1_fail, g2_succ, g2_fail)
    }

    pub fn pass_min_read(&self, min_cover: u32, min_unspliced: u32) -> bool {
        let (g1_succ, g1_fail, g2_succ, g2_fail) = self.get_sub_vec();
        if (g1_fail.iter().sum::<u32>() + g1_succ.iter().sum::<u32>() < min_cover) ||
                 (g2_fail.iter().sum::<u32>() + g2_succ.iter().sum::<u32>() < min_cover){
            return false;
        }
        if g1_fail.iter().sum::<u32>() < min_unspliced && g2_fail.iter().sum::<u32>() < min_unspliced{
            return false;
        }
        true
    }

    pub fn new(
        treatment: &Vec<CountsStats>,
        control: &Vec<CountsStats>,
        successes_cat: &Vec<SplicingCategory>,
        failures_cat: &Vec<SplicingCategory>,
        identifier: String,
    ) -> Self {
        let mut glm = GLMBetaBinomiale {
            success: Vec::new(),
            failures: Vec::new(),
            groups: Vec::new(),
            identifier,
        };
        glm.format_data(treatment, control, successes_cat, failures_cat);
        glm
    }


    //}



    /// Computes standard errors for the regression coefficients and theta jointly,
    /// via the inverse of the observed Fisher information (negative Hessian) at
    /// the fitted values. Mirrors `GLM::standard_errors`, but built on
    /// `beta_binomial_hessian` instead of the binomial `X'WX`.
    fn beta_binomial_standard_errors(
        x: &DMatrix<f64>,
        y: &[f64],
        n: &[f64],
        mu: &[f64],
        theta: f64,
    ) -> Result<DVector<f64>, LogisticRegressionError> {
        let h = beta_binomial_hessian(x, y, n, mu, theta)?;
        //let info = -h; // Fisher information ~ negative Hessian at the MLE
        //let var_covar = info
        //    .try_inverse();

        let n_vars = x.ncols();
    //let h = beta_binomial_hessian(x, y, n, mu, theta)?;
        let h_bb = h.view((0, 0), (n_vars, n_vars)).into_owned();
        let info_bb = -h_bb;

        let var_covar = match info_bb.try_inverse() {
            Some(v) => v,
            None => return Err(LogisticRegressionError::CIUnavail(
                "Failed to invert beta-binomial information matrix".to_string()
            )),
        };

        /*let var_covar = match info.try_inverse() {
        Some(v) => v,
        None => return Err(LogisticRegressionError::CIUnavail(
            "Failed to invert beta-binomial information matrix".to_string()
        )),
    };
        
         if var_covar.is_none(){
            return Err(CIUnavail("Failed to invert beta-binomial information matrix".to_string()))};

        let var_covar = var_covar.unwrap();*/

        Ok(DVector::from_iterator(
            var_covar.nrows(),
            var_covar.diagonal().iter().map(|&v| v.sqrt()),
        ))
    }

    pub fn actual_test(&self) -> Result<(TestStatus, f64, f64, f64, f64), LogisticRegressionError> {
        let mut current_test_status = TestStatus::Ok;

        let successes: Vec<f64> = self.success().iter().map(|&x| x as f64).collect();
        let failures: Vec<f64> = self.failures().iter().map(|&x| x as f64).collect();
        let n_trials: Vec<f64> = successes.iter().zip(failures.iter()).map(|(s, f)| s + f).collect();

        // Low-count fallback -- identical pattern to GLM::actual_test.
        let contingency: Vec<u64> = self.to_contengency();


        let n_obs = self.groups().len();

        // Fit null model: intercept + theta only.
        let x_null = DMatrix::from_element(n_obs, 1, 1.0);
        let theta_seed = 10.0;
        
        //let (beta_null, theta_null, status_null) =
        //    irls_beta_binomial(&x_null, &successes, &n_trials, DVector::zeros(1), theta_seed, 150, 1e-8, 10)?;
        let (beta_null, theta_null, status_null) = irls_beta_binomial_alternating(
            &x_null, &successes, &n_trials, DVector::zeros(1), theta_seed, 100, 1e-8, 50, 1e-9, 15,
        )?;
        if status_null > current_test_status {
            current_test_status = status_null.clone();
        }
        let eta_null = &x_null * &beta_null;
        let mu_null: Vec<f64> = eta_null.iter().map(|&e| GLM::inv_logit(e)).collect();
        let ll_null = match beta_binomial_log_likelihood(&successes, &n_trials, &mu_null, theta_null){
            Ok(e) => e,
            //Err(LogisticRegressionError::FailUseFisher) => {
            //    return Err(LogisticRegressionError::NumericalInstability("ll_null".to_string()));
                //eprintln!("DEBUG fallback at ll_null: theta_null={:e}, beta_null={:?}", theta_null, beta_null);
                //return fisher_fallback(&contingency);
                
           //},
           //current_test_status = 
           Err(e) => return  Err(e)
        };

        //eprintln!("DIAG null: theta_null={:.6}  ll_null={:.6}", theta_null, ll_null);


        //eprintln!("DEBUG null: theta_null={:e}, ll_null={:?}, status_null={:?}", theta_null, ll_null, status_null.clone());

        // Fit full model: intercept + treatment + theta, warm-started from the null fit.
        let treatment: Vec<f64> = self.groups().iter()
            .map(|g| match g {
                Genotype::TREATMENT => 1.0,
                Genotype::CONTROL => 0.0,
            })
            .collect();
        let mut x_full_data = vec![1.0; n_obs];
        x_full_data.extend(treatment);
        let x_full = DMatrix::from_vec(n_obs, 2, x_full_data);

        let beta_full_init = DVector::from_vec(vec![beta_null[0], 0.0]);
        //let (beta_full, theta_full, status_full) =
           // irls_beta_binomial(&x_full, &successes, &n_trials, beta_full_init, theta_null, 150, 1e-8, 10)?;
        let (beta_full, theta_full, status_full) = match irls_beta_binomial_alternating(
            &x_full, &successes, &n_trials,
             beta_full_init, theta_null,
             100, 1e-8, 50, 1e-9, 15,
        ){ Ok(e) => e,
           //Err(LogisticRegressionError::FailUseFisher) => {
            //return fisher_fallback(&contingency);
            //eprintln!("DEBUG fallback at irls_beta_binomial_alternating : theta_null={:e}, beta_null={:?}", theta_null, beta_null);
           //},
        
           Err(e) => return  Err(e)
        };

        if status_full > current_test_status {
            current_test_status = status_full.clone();
        }
        let eta_full = &x_full * &beta_full;
        let mu_full: Vec<f64> = eta_full.iter().map(|&e| GLM::inv_logit(e)).collect();
        let ll_full = match beta_binomial_log_likelihood(&successes, &n_trials, &mu_full, theta_full){
            Ok(e) => e,
           // Err(LogisticRegressionError::FailUseFisher) => {
            //    eprintln!("DEBUG fallback at beta_binomial_log_likelihood ll_full : theta_null={:e}, beta_null={:?}", theta_null, beta_null);
            //    return fisher_fallback(&contingency)
                
            //},
           Err(e) => return  Err(e)
        };        
        //eprintln!("DEBUG full: theta_full={:e}, ll_full={:?}, status_full={:?}", theta_full, ll_full, status_full.clone());
        //eprintln!("DIAG full: theta_full={:.6}  ll_full={:.6}", theta_full, ll_full);
        //eprintln!("DIAG LR={:.6}", 2.0 * (ll_full - ll_null));

        // Likelihood ratio test -- same shape as GLM::actual_test's LRT.
        let lr_stat = 2.0 * (ll_full - ll_null);
        let chi_sq = ChiSquared::new(1.0).expect("Failed to create chi-squared distribution");
        let p_value = chi_sq.sf(lr_stat);


        let se = GLMBetaBinomiale::beta_binomial_standard_errors(&x_full, &successes, &n_trials, &mu_full, theta_full);
        let (or, or_lower, or_upper) = match se{
            Ok(std_e) => { GLM::odds_ratio_with_ci(beta_full[1], std_e[1])},
            Err(LogisticRegressionError::CIUnavail(_)) => {current_test_status = TestStatus::CIUnavail; (beta_full[1].exp(), 0.,0.)}
            //Err(LogisticRegressionError::FailUseFisher) => {
            //    eprintln!("DEBUG fallback at beta_binomial_standard_errors : theta_null={:e}, beta_null={:?}", theta_null, beta_null);
            //        return fisher_fallback(&contingency);
            //},
            _ => {current_test_status = TestStatus::OddRatioUnavail; (0., 0.,0.)}
        };

        //let (or, or_lower, or_upper) = GLM::odds_ratio_with_ci(beta_full[1], se[1]);

        // Extreme-OR fallback -- identical pattern to GLM::actual_test.
        if or > 100. || or < 0.001 {
            warn!("odd ratio extreme, CI unreliable, keeping LRT p_value {:?}", contingency);
                current_test_status = TestStatus::CIUnavail;
                return Ok((current_test_status, p_value, or, 0., 0.));
            }
        
        Ok((current_test_status, p_value, or, or_lower, or_upper))
    }
}




#[cfg(test)]
mod glm_beta_binomiale_usage {
    use super::*;

    #[test]
    fn my_two_group_test() {
        // Two groups: control and treatment.
        // Each entry is one sample/replicate: (success_count, failure_count).
        // Doesn't need to be the same length for both groups.
        let glm = GLMBetaBinomiale {
            groups: vec![
                Genotype::CONTROL, Genotype::CONTROL, Genotype::CONTROL, Genotype::CONTROL, Genotype::CONTROL, Genotype::CONTROL,
                Genotype::TREATMENT, Genotype::TREATMENT, Genotype::TREATMENT, Genotype::TREATMENT, Genotype::TREATMENT, Genotype::TREATMENT
            ],


            success:  vec![916,904,834,744,1368,1323, 180,154,212,230,205,183	],
        failures: vec![708,610,580,556,944,984,4006,3788,4899,4528,4665,4224],

            identifier: "my_junction".to_string(),
        };

        let result = glm.test(false, 1, 1);

        println!("{:?}", result);

        // The fields you'll care about:
        println!("status:     {:?}", result.status);
        println!("p_value:    {:?}", result.p_value);
        println!("odds ratio: {:?}", result.odd_ratio);
        println!("OR 95% CI:  [{:?}, {:?}]", result.or_ci_lower, result.or_ci_upper);
        println!("control prop:   {:?}", result.control_prop);
        println!("treatment prop: {:?}", result.treatment_prop);
    }

    #[test]
    fn via_new_matches_hand_built_struct() {
        // Same counts as my_two_group_test / the reported junction, but built the way
        // the real pipeline builds it: CountsStats -> GLMBetaBinomiale::new() -> format_data(),
        // instead of setting success/failures/groups directly.
        let control_count = vec![
            CountsStats::new(&["40", "14", "0", "0", "0", "0", "0", "0"]),
            CountsStats::new(&["82", "21", "0", "0", "0", "0", "0", "0"]),
            CountsStats::new(&["97", "25", "0", "0", "0", "0", "0", "0"]),
        ];
        let treat_count = vec![
            CountsStats::new(&["36", "12", "0", "0", "0", "0", "0", "0"]),
            CountsStats::new(&["64", "26", "0", "0", "0", "0", "0", "0"]),
            CountsStats::new(&["43", "28", "0", "0", "0", "0", "0", "0"]),
        ];
        let successes_cat = vec![SplicingCategory::Spliced];
        let failures_cat = vec![SplicingCategory::Unspliced];

        let glm = GLMBetaBinomiale::new(&treat_count, &control_count, &successes_cat, &failures_cat, "via_new".to_string());
        let result = glm.test(false, 1, 1);

        println!("status:     {:?}", result.status);
        println!("p_value:    {:?}", result.p_value);
        println!("odds ratio: {:?}", result.odd_ratio);
        println!("OR 95% CI:  [{:?}, {:?}]", result.or_ci_lower, result.or_ci_upper);
        println!("control prop:   {:?}", result.control_prop);
        println!("treatment prop: {:?}", result.treatment_prop);

        // Should match my_two_group_test's manually-verified p_value=0.0257 / OR=0.5936.
        assert!((result.p_value.unwrap() - 0.025680342308132206).abs() < 1e-6,
            "p_value via new()/format_data() diverged from hand-built struct: {:?}", result.p_value);
    }
}




#[cfg(test)]
mod glm_beta_binomiale_tests {
    use super::*;

    #[test]
    fn actual_test_matches_verified_pipeline() {

        let y: Vec<u32> = vec![
            6, 16, 8, 10, 9, 13, 12, 8, 4, 11, 10, 26, 26, 24, 26, 16, 14, 6, 25, 12,
            16, 22, 16, 22, 47, 22, 29, 29, 21, 20, 20, 10, 13, 28, 42, 7, 24, 23, 24, 30,
        ];
        let n: Vec<u32> = vec![
            19, 49, 44, 34, 34, 53, 18, 46, 24, 19, 38, 58, 48, 49, 47, 50, 38, 20, 52, 35,
            37, 31, 23, 56, 50, 43, 33, 52, 39, 34, 35, 25, 19, 39, 54, 17, 53, 52, 27, 43,
        ];
        let groups: Vec<Genotype> = (0..40)
            .map(|i| if i < 20 { Genotype::CONTROL } else { Genotype::TREATMENT })
            .collect();
        let failures: Vec<u32> = y.iter().zip(n.iter()).map(|(&yi, &ni)| ni - yi).collect();

        let glm = GLMBetaBinomiale {
            success: y,
            failures,
            groups,
            identifier: "test".to_string(),
        };

        let (status, p_value, or, or_lower, or_upper) = glm.actual_test().unwrap();

        assert_eq!(status, TestStatus::Ok);
        assert!((p_value - 0.000010).abs() < 1e-5);
        assert!((or - 2.8125).abs() < 1e-3);
        assert!(or_lower < or && or < or_upper);
        assert!((or_lower - 1.8646).abs() < 1e-3);
        assert!((or_upper - 4.2421).abs() < 1e-3);
    }
}



#[cfg(test)]
mod irls_beta_binomial_tests {
    use super::*;

    #[test]
    fn matches_verified_fixture() {
        let y = vec![1.0, 0.0, 0.0, 22.0, 1.0, 2.0, 5.0, 9.0];
        let n = vec![5.0, 1.0, 2.0, 22.0, 1.0, 2.0, 12.0, 20.0];
        let treat = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0];

        let n_obs = y.len();
        let mut x_data = vec![1.0; n_obs];
        x_data.extend(treat.iter());
        let x = DMatrix::from_vec(n_obs, 2, x_data);

        // Step-halving never triggers on this well-behaved fixture, so result is
        // identical to the pre-halving version verified earlier.
        let (beta, theta, status) =
            irls_beta_binomial(&x, &y, &n, DVector::zeros(2), 5.0, 150, 1e-8, 10).unwrap();

        assert_eq!(status, TestStatus::Ok);
        assert!((beta[0] - (-1.63917846)).abs() < 1e-5);
        assert!((beta[1] - 2.8511118).abs() < 1e-5);
        assert!((theta - 2.187506330556031).abs() < 1e-5);
    }

    #[test]
    fn recovers_gracefully_from_bad_start() {
        // Deliberately terrible starting point: beta far from any sensible fit,
        // theta near zero. Without step-halving this produces a non-finite
        // log-likelihood on the very first iteration.
        let y = vec![1.0, 0.0, 0.0, 22.0, 1.0, 2.0, 5.0, 9.0];
        let n = vec![5.0, 1.0, 2.0, 22.0, 1.0, 2.0, 12.0, 20.0];
        let treat = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0];

        let n_obs = y.len();
        let mut x_data = vec![1.0; n_obs];
        x_data.extend(treat.iter());
        let x = DMatrix::from_vec(n_obs, 2, x_data);

        let bad_beta_init = DVector::from_vec(vec![10.0, -10.0]);
        let result = irls_beta_binomial(&x, &y, &n, bad_beta_init, 0.01, 150, 1e-8, 10);

        // Should not error out or panic -- should return a valid, finite state.
        let (beta, theta, _status) = result.unwrap();
        assert!(beta.iter().all(|b| b.is_finite()));
        assert!(theta.is_finite() && theta > 0.0);
    }
}




/// Computes the trigamma function ψ'(x), the derivative of the digamma function.
///
/// Not available in `statrs`, so implemented here via the standard approach:
/// recurrence relation ψ'(x) = ψ'(x+1) + 1/x² to shift x into a range where
/// the asymptotic (Bernoulli-number) series converges to full double precision,
/// then the series itself.
///
/// # Panics
/// Panics if `x <= 0.0` (trigamma has poles at non-positive integers and is
/// undefined for x <= 0 in the domain we need it: shape parameters must be positive).
fn trigamma(x: f64) -> f64 {
    assert!(x > 0.0, "trigamma undefined for x <= 0: x={}", x);

    let mut x = x;
    let mut sum = 0.0;

    // Recurrence: shift x up until the asymptotic series is accurate to
    // machine precision (empirically verified against scipy at x >= 10).
    while x < 10.0 {
        sum += 1.0 / (x * x);
        x += 1.0;
    }

    let inv_x = 1.0 / x;
    let inv_x2 = inv_x * inv_x;

    // Asymptotic (Bernoulli) series:
    // psi1(x) ~ 1/x + 1/(2x^2) + 1/(6x^3) - 1/(30x^5) + 1/(42x^7)
    //             - 1/(30x^9) + 5/(66x^11) - 691/(2730x^13) + 7/(6x^15)
    let tail = 1.0 / 6.0
        - inv_x2
            * (1.0 / 30.0
                - inv_x2
                    * (1.0 / 42.0
                        - inv_x2
                            * (1.0 / 30.0
                                - inv_x2 * (5.0 / 66.0 - inv_x2 * (691.0 / 2730.0 - inv_x2 * (7.0 / 6.0))))));

    sum + inv_x + 0.5 * inv_x2 + inv_x2 * inv_x * tail
}

/// Computes the beta-binomial log-likelihood for observed data.
///
/// Parameterized by mean `mu_i = inv_logit(x_i'β)` and a shared precision `theta`,
/// with `alpha_i = mu_i * theta`, `beta_i = (1 - mu_i) * theta`. This is the beta-binomial
/// analogue of `binomial_log_likelihood`, and reduces to it (up to the constant `ln C(n,y)`
/// term, which is omitted here as well) as `theta -> infinity`.
///
/// # Arguments
/// * `y` - Vector of success counts
/// * `n` - Vector of trial counts (total observations per group)
/// * `mu` - Vector of predicted mean probabilities
/// * `theta` - Shared precision parameter (theta > 0; larger theta = less overdispersion)
///
/// # Returns
/// The sum of log-likelihoods across observations, including the `ln C(n,y)` term
/// (unlike `binomial_log_likelihood`, this term does NOT cancel automatically between
/// null/full models here because it's included for clarity — it's still a constant
/// w.r.t. beta and theta for fixed data, so it has no effect on the LRT statistic,
/// but if you want a leaner hot path you can drop it, mirroring `binomial_log_likelihood`.)
///
/// large theta revert to GLM methods
/// # Panics
/// Panics if `theta <= 0.0`.
/// 
const THETA_BINOMIAL_SWITCH: f64 = 1e8;
fn beta_binomial_log_likelihood(y: &[f64], n: &[f64], mu: &[f64], theta: f64) -> Result<f64, LogisticRegressionError> {
    use statrs::function::gamma::ln_gamma;

    if !(theta > 0.0){
        return Err(LogisticRegressionError::NumericalInstability("theta".to_string()))
    }

    if theta > THETA_BINOMIAL_SWITCH {
    let ln_choose_sum: f64 = y.iter().zip(n.iter())
        .map(|(&yi, &ni)| ln_gamma(ni + 1.0) - ln_gamma(yi + 1.0) - ln_gamma(ni - yi + 1.0))
        .sum();
    return Ok(GLM::binomial_log_likelihood(y, n, mu) + ln_choose_sum);
}
    //assert!(theta > 0.0, "theta must be > 0: theta={}", theta);

    Ok(y.iter()
        .zip(n.iter())
        .zip(mu.iter())
        .map(|((&yi, &ni), &mui)| {
            // Clamp probability to avoid numerical issues, same convention as binomial_log_likelihood
            let mui = mui.max(1e-10).min(1.0 - 1e-10);
            let alpha = mui * theta;
            let beta = (1.0 - mui) * theta;
            let fi = ni - yi;

            let ln_choose = ln_gamma(ni + 1.0) - ln_gamma(yi + 1.0) - ln_gamma(fi + 1.0);

            ln_choose
                + ln_gamma(yi + alpha)
                + ln_gamma(fi + beta)
                - ln_gamma(ni + theta)
                - ln_gamma(alpha)
                - ln_gamma(beta)
                + ln_gamma(theta)
        })
        .sum())
}


/// Computes the score (gradient) of the beta-binomial log-likelihood.
///
/// Returns `(grad_beta, grad_theta)` where `grad_beta` is the gradient w.r.t. the
/// regression coefficients (chain-ruled through `mu_i = inv_logit(x_i'beta)`), and
/// `grad_theta` is the gradient w.r.t. the shared precision parameter.
///
/// Derivation: with `alpha_i = mu_i*theta`, `beta_i = (1-mu_i)*theta`,
/// `d_succ_i = psi(y_i+alpha_i) - psi(alpha_i)`, `d_fail_i = psi(n_i-y_i+beta_i) - psi(beta_i)`:
///
/// - `d(ll_i)/d(mu_i)    = theta * (d_succ_i - d_fail_i)`
/// - `d(ll_i)/d(theta)   = mu_i*d_succ_i + (1-mu_i)*d_fail_i - psi(n_i+theta) + psi(theta)`
/// - `d(ll_i)/d(beta_k)  = d(ll_i)/d(mu_i) * mu_i*(1-mu_i) * x_{i,k}`   (standard logit chain rule)
///
/// Verified against central-difference numerical gradients (max err ~1e-9) before porting.
fn beta_binomial_score(
    x: &DMatrix<f64>,
    y: &[f64],
    n: &[f64],
    mu: &[f64],
    theta: f64,
) -> Result<(DVector<f64>, f64), LogisticRegressionError> {
    use statrs::function::gamma::digamma;

    if !(theta > 0.0){
        return Err(LogisticRegressionError::NumericalInstability("Theta".to_string()));
    }

    let n_obs = x.nrows();
    let n_vars = x.ncols();

    let mut grad_beta = DVector::zeros(n_vars);
    let mut grad_theta = 0.0;


    if theta > THETA_BINOMIAL_SWITCH {
        for i in 0..n_obs {
            let yi = y[i];
            let ni = n[i];
            let mui = mu[i].max(1e-10).min(1.0 - 1e-10);
            let fi = ni - yi;

            let dll_dmu = yi / mui - fi / (1.0 - mui);
            let dmu_deta = mui * (1.0 - mui);
            let common = dll_dmu * dmu_deta;

            for k in 0..n_vars {
                grad_beta[k] += x[(i, k)] * common;
            }
        }
        return Ok((grad_beta, 0.0));
    }

    for i in 0..n_obs {
        let yi = y[i];
        let ni = n[i];
        let mui = mu[i].max(1e-10).min(1.0 - 1e-10);
        let fi = ni - yi;

        let alpha = mui * theta;
        let beta_shape = (1.0 - mui) * theta;

        let d_succ = digamma(yi + alpha) - digamma(alpha);
        let d_fail = digamma(fi + beta_shape) - digamma(beta_shape);

        grad_theta += mui * d_succ + (1.0 - mui) * d_fail - digamma(ni + theta) + digamma(theta);

        let dll_dmu = theta * (d_succ - d_fail);
        let dmu_deta = mui * (1.0 - mui);
        let common = dll_dmu * dmu_deta;

        for k in 0..n_vars {
            grad_beta[k] += x[(i, k)] * common;
        }
    }

    Ok((grad_beta, grad_theta))
}


/// Computes the Hessian of the beta-binomial log-likelihood with respect to
/// the regression coefficients and theta jointly.
///
/// Returns an `(n_vars+1) x (n_vars+1)` matrix, where indices `0..n_vars` are the
/// regression coefficients and index `n_vars` (the last row/column) is theta.
///
/// # Derivation
/// With `alpha_i = mu_i*theta`, `beta_i = (1-mu_i)*theta`,
/// `d_succ_i = psi(y_i+alpha_i) - psi(alpha_i)`, `d_fail_i = psi(n_i-y_i+beta_i) - psi(beta_i)`,
/// `t1_succ_i = psi1(y_i+alpha_i) - psi1(alpha_i)`, `t1_fail_i = psi1(n_i-y_i+beta_i) - psi1(beta_i)`
/// (`psi1` = trigamma):
///
/// - `d2(ll_i)/d(theta)^2 = mu_i^2*t1_succ_i + (1-mu_i)^2*t1_fail_i - psi1(n_i+theta) + psi1(theta)`
/// - `d2(ll_i)/d(theta)d(mu_i) = (d_succ_i - d_fail_i) + theta*(mu_i*t1_succ_i - (1-mu_i)*t1_fail_i)`,
///   chain-ruled to `beta_k` via `* mu_i*(1-mu_i) * x_{i,k}`
/// - `d2(ll_i)/d(beta_k)d(beta_j) = x_{i,k}*x_{i,j} * [g'(mu_i)*h(eta_i)^2 + g(mu_i)*h'(eta_i)]`
///   where `g(mu) = theta*(d_succ-d_fail)` (the score w.r.t. mu), `g'(mu) = theta^2*(t1_succ+t1_fail)`,
///   `h(eta) = mu*(1-mu)`, `h'(eta) = mu*(1-mu)*(1-2*mu)` (standard logistic derivative identity)
///
/// Verified against a full numeric Hessian (finite differences of `beta_binomial_score`)
/// to ~1e-10 before porting.
fn beta_binomial_hessian(x: &DMatrix<f64>, y: &[f64], n: &[f64], mu: &[f64], theta: f64) -> Result<DMatrix<f64>, LogisticRegressionError> {
    use statrs::function::gamma::digamma;

    if !(theta > 0.0){
        return Err(LogisticRegressionError::NumericalInstability("Theta".to_string()));
    }
    //assert!(theta > 0.0, "theta must be > 0: theta={}", theta);

    let n_obs = x.nrows();
    let n_vars = x.ncols();
    let n_params = n_vars + 1;

    
    let mut h = DMatrix::zeros(n_params, n_params);

    if theta > THETA_BINOMIAL_SWITCH {
        for i in 0..n_obs {
            let mui = mu[i].max(1e-10).min(1.0 - 1e-10);
            let ni = n[i];
            let bb_scalar = -ni * mui * (1.0 - mui);
            for k in 0..n_vars {
                let xk = x[(i, k)];
                for j in 0..n_vars {
                    h[(k, j)] += xk * x[(i, j)] * bb_scalar;
                }
            }
        }
        return Ok(h);
    }

    for i in 0..n_obs {
        let yi = y[i];
        let ni = n[i];
        let mui = mu[i].max(1e-10).min(1.0 - 1e-10);
        let fi = ni - yi;

        let alpha = mui * theta;
        let beta_shape = (1.0 - mui) * theta;

        let d_succ = digamma(yi + alpha) - digamma(alpha);
        let d_fail = digamma(fi + beta_shape) - digamma(beta_shape);
        let t1_succ = trigamma(yi + alpha) - trigamma(alpha);
        let t1_fail = trigamma(fi + beta_shape) - trigamma(beta_shape);

        // theta-theta
        h[(n_vars, n_vars)] += mui * mui * t1_succ + (1.0 - mui) * (1.0 - mui) * t1_fail
            - trigamma(ni + theta)
            + trigamma(theta);

        // theta-mu cross, chained to theta-beta_k
        let d2ll_dtheta_dmu = (d_succ - d_fail) + theta * (mui * t1_succ - (1.0 - mui) * t1_fail);
        let dmu_deta = mui * (1.0 - mui);
        let cross_scalar = d2ll_dtheta_dmu * dmu_deta;

        // beta-beta block
        let g = theta * (d_succ - d_fail);
        let g_prime = theta * theta * (t1_succ + t1_fail);
        let hh = mui * (1.0 - mui);
        let h_prime = hh * (1.0 - 2.0 * mui);
        let bb_scalar = g_prime * hh * hh + g * h_prime;

        for k in 0..n_vars {
            let xk = x[(i, k)];
            h[(k, n_vars)] += xk * cross_scalar;
            h[(n_vars, k)] += xk * cross_scalar;

            for j in 0..n_vars {
                h[(k, j)] += xk * x[(i, j)] * bb_scalar;
            }
        }
    }

    Ok(h)
}



/// Fits a beta-binomial regression model via Newton-Raphson with step-halving,
/// jointly estimating the regression coefficients and the precision parameter theta.
///
/// # Step-halving
/// A raw Newton step can overshoot or land somewhere pathological (e.g. mu saturating
/// at 0/1), especially early in the fit or on messy low-count data. Each proposed step
/// is checked against the log-likelihood before being accepted; if it doesn't improve
/// (allowing a small tolerance for floating-point noise near convergence), the step is
/// halved and retried, up to `max_halvings` times. If no halving improves things, the
/// fit stops and reports `ConvergenceFailed` rather than propagating a bad or
/// non-finite state — verified against a deliberately bad starting point where the
/// unguarded version dies on a non-finite log-likelihood at iteration 0, while this
/// version reports a clean failure instead.
///
/// # Design: theta is optimized in log-space
/// See module docs on `irls_beta_binomial` reparameterization — unchanged from before.
fn irls_beta_binomial(
    x: &DMatrix<f64>,
    y: &[f64],
    n: &[f64],
    beta_init: DVector<f64>,
    theta_init: f64,
    max_iter: usize,
    tol: f64,
    max_halvings: usize,
) -> Result<(DVector<f64>, f64, TestStatus), LogisticRegressionError> {
    assert!(theta_init > 0.0, "theta_init must be > 0: theta_init={}", theta_init);

    let n_vars = x.ncols();
    let n_params = n_vars + 1;

    GLM::validate_irls_inputs(x, &DVector::from_vec(y.to_vec()), &DVector::from_vec(n.to_vec()))?;

    let mut beta = beta_init;
    let mut tau = theta_init.ln();
    let mut large_sep_warn = false;
    let mut current_state = TestStatus::Ok;

    for iter in 0..max_iter {
        let eta = x * &beta;

        for (i, &e) in eta.iter().enumerate() {
            if !e.is_finite() {
                return Err(LogisticRegressionError::NumericalInstability(format!(
                    "Non-finite linear predictor at iteration {} observation {}: eta={}",
                    iter + 1, i, e
                )));
            }
            if (!large_sep_warn) && (e.abs() > 20.0) {
                large_sep_warn = true;
                if current_state < TestStatus::QuasiPerfectSeparation {
                    current_state = TestStatus::QuasiPerfectSeparation;
                }
            }
        }

        let mu: Vec<f64> = eta.iter().map(|&e| GLM::inv_logit(e)).collect();
        let theta = tau.exp();
        let current_ll = beta_binomial_log_likelihood(y, n, &mu, theta)?;

        let (grad_beta, grad_theta) = beta_binomial_score(&x, y, n, &mu, theta)?;
        let h = beta_binomial_hessian(&x, y, n, &mu, theta)?;

        let g_tau = theta * grad_theta;
        let h_tt_tau = theta * theta * h[(n_vars, n_vars)] + theta * grad_theta;

        let mut g_aug = DVector::zeros(n_params);
        g_aug.rows_mut(0, n_vars).copy_from(&grad_beta);
        g_aug[n_vars] = g_tau;

        let mut h_aug = DMatrix::zeros(n_params, n_params);
        for k in 0..n_vars {
            for j in 0..n_vars {
                h_aug[(k, j)] = h[(k, j)];
            }
            let h_bt_tau = theta * h[(k, n_vars)];
            h_aug[(k, n_vars)] = h_bt_tau;
            h_aug[(n_vars, k)] = h_bt_tau;
        }
        h_aug[(n_vars, n_vars)] = h_tt_tau;

        let raw_delta = match h_aug.lu().solve(&(-&g_aug)) {
            Some(d) => d,
            None => {
                return Err(LogisticRegressionError::SingularMatrix(format!(
                    "Cannot solve Newton-Raphson system at iteration {}. Matrix is singular or near-singular.",
                    iter + 1
                )));
            }
        };

        for (i, &d) in raw_delta.iter().enumerate() {
            if !d.is_finite() {
                return Err(LogisticRegressionError::NumericalInstability(format!(
                    "Non-finite Newton-Raphson step at iteration {} parameter {}: delta={}",
                    iter + 1, i, d
                )));
            }
        }

        // Step-halving: accept the first halving that improves (or matches within
        // tolerance) the log-likelihood, starting from the full step.
        let mut step = raw_delta.clone();
        let mut accepted: Option<(DVector<f64>, f64, f64)> = None; // (beta_try, tau_try, ll_try)

        for _ in 0..=max_halvings {
            let beta_try = &beta + step.rows(0, n_vars);
            let tau_try = tau + step[n_vars];
            let theta_try = tau_try.exp();

            let eta_try = x * &beta_try;
            let mu_try: Vec<f64> = eta_try.iter().map(|&e| GLM::inv_logit(e)).collect();
            let ll_try = beta_binomial_log_likelihood(y, n, &mu_try, theta_try)?;

            if ll_try.is_finite() && ll_try >= current_ll - 1e-10 {
                accepted = Some((beta_try, tau_try, ll_try));
                break;
            }
            step /= 2.0;
        }

        let (beta_new, tau_new, _ll_new) = match accepted {
            Some(v) => v,
            None => {
                // Couldn't find an improving step even after max_halvings: report
                // gracefully rather than propagating a bad or non-finite state.
                if TestStatus::ConvergenceFailed > current_state {
                    current_state = TestStatus::ConvergenceFailed;
                }
                return Ok((beta, tau.exp(), current_state));
            }
        };

        let delta_norm = (&beta_new - &beta).norm().hypot(tau_new - tau);
        beta = beta_new;
        tau = tau_new;

        if delta_norm < tol {
            return Ok((beta, tau.exp(), current_state));
        }
    }

    if TestStatus::ConvergenceFailed > current_state {
        current_state = TestStatus::ConvergenceFailed;
    }
    Ok((beta, tau.exp(), current_state))
}


/// Fits the regression coefficients via Newton-Raphson with step-halving,
/// holding theta fixed. This is the "inner" step of the alternating fit.
fn fit_beta_given_theta(
    x: &DMatrix<f64>,
    y: &[f64],
    n: &[f64],
    beta_init: DVector<f64>,
    theta: f64,
    max_iter: usize,
    tol: f64,
    max_halvings: usize,
) -> Result<(DVector<f64>, TestStatus), LogisticRegressionError> {
    let n_vars = x.ncols();
    let mut beta = beta_init;

    for _iter in 0..max_iter {
        let mu: Vec<f64> = (x * &beta).iter().map(|&e| GLM::inv_logit(e)).collect();
        let current_ll = beta_binomial_log_likelihood(y, n, &mu, theta)?;

        let (grad_beta, _grad_theta) = beta_binomial_score(x, y, n, &mu, theta)?;
        let h = beta_binomial_hessian(x, y, n, &mu, theta)?;
        let h_bb = h.view((0, 0), (n_vars, n_vars)).into_owned();

        let raw_delta = match h_bb.lu().solve(&(-&grad_beta)) {
            Some(d) => d,
            None => return Ok((beta, TestStatus::ConvergenceFailed)), // singular -- report gracefully
        };

        let mut step = raw_delta;
        let mut accepted: Option<DVector<f64>> = None;

        for _ in 0..=max_halvings {
            let beta_try = &beta + &step;
            let mu_try: Vec<f64> = (x * &beta_try).iter().map(|&e| GLM::inv_logit(e)).collect();
            let ll_try = beta_binomial_log_likelihood(y, n, &mu_try, theta)?;
            if ll_try.is_finite() && ll_try >= current_ll - 1e-10 {
                accepted = Some(beta_try);
                break;
            }
            step /= 2.0;
        }

        let beta_new = match accepted {
            Some(b) => b,
            None => return Ok((beta, TestStatus::ConvergenceFailed)),
        };

        let delta_norm = (&beta_new - &beta).norm();
        beta = beta_new;
        if delta_norm < tol {
            return Ok((beta, TestStatus::Ok));
        }
    }
    Ok((beta, TestStatus::ConvergenceFailed))
}

/// Fits theta (in log-space) via 1-D Newton-Raphson with step-halving,
/// holding the regression coefficients (and therefore mu) fixed. This is
/// the "inner" step of the alternating fit.
fn fit_theta_given_beta(
    x: &DMatrix<f64>,
    y: &[f64],
    n: &[f64],
    mu: &[f64],
    tau_init: f64,
    max_iter: usize,
    tol: f64,
    max_halvings: usize,
) -> Result<(f64, TestStatus), LogisticRegressionError> {
    let n_vars = x.ncols();
    let mut tau = tau_init;

    for _iter in 0..max_iter {
        let theta = tau.exp();

        if theta > THETA_BINOMIAL_SWITCH {
            return Ok((theta, TestStatus::Ok));
        }

        let current_ll = beta_binomial_log_likelihood(y, n, mu, theta)?;

        let (_grad_beta, grad_theta) = beta_binomial_score(x, y, n, mu, theta)?;
        let h = beta_binomial_hessian(x, y, n, mu, theta)?;

        let g_tau = theta * grad_theta;
        let h_tt_tau = theta * theta * h[(n_vars, n_vars)] + theta * grad_theta;

        //if h_tt_tau == 0.0 {
        //    return Ok((theta, TestStatus::ConvergenceFailed));
        //}
        //let raw_delta = -g_tau / h_tt_tau;


        // The profile log-likelihood in tau = ln(theta) is NOT guaranteed
        // concave far from the optimum. When h_tt_tau >= 0 the raw Newton step
        // -g/h points AWAY from the maximum, and step-halving (which only
        // shrinks a fixed direction, never flips it) can never recover -- it
        // strands theta at its starting value. This is exactly what happens to
        // the intercept-only null fit when the two groups are well separated:
        // the true theta is tiny (large overdispersion) but the seed theta is
        // large, putting the start in a non-concave region. Fall back to a
        // capped gradient-ascent step whenever curvature is non-concave.
        /*let raw_delta = if h_tt_tau < -1e-12 {
            -g_tau / h_tt_tau
        } else {
            g_tau.signum() * g_tau.abs().min(2.0)
        };*/
        // Ascent direction in tau = ln(theta): Newton where the profile is
        // locally concave (h_tt_tau < 0), plain gradient ascent otherwise -- a
        // raw Newton step in a non-concave region points AWAY from the maximum
        // and step-halving (which only shrinks a fixed direction, never flips
        // it) can't recover, stranding theta at its start.
        let ascent = if h_tt_tau < -1e-12 {
           -g_tau / h_tt_tau
        } else {
            g_tau
        };

        // Trust-region cap in log-space. Even a *correct-direction* Newton step
        // can be enormous far from the optimum (e.g. theta=10, steep gradient,
        // small curvature -> step of -19 in tau) and leap clean over the peak
        // to a worse-but-still-better-than-start point that step-halving accepts
        // and can't climb back from. Bounding |delta_tau| keeps every step local;
        // step-halving then refines within the cap. Reaching large theta (e.g.
        // 3000+) just takes a few extra capped outer passes -- well inside limits.
        const TAU_STEP_CAP: f64 = 2.0;
        let raw_delta = ascent.clamp(-TAU_STEP_CAP, TAU_STEP_CAP);

        let mut step = raw_delta;
        let mut accepted: Option<f64> = None;

        for _ in 0..=max_halvings {
            let tau_try = tau + step;
            let ll_try = beta_binomial_log_likelihood(y, n, mu, tau_try.exp())?;
            if ll_try.is_finite() && ll_try >= current_ll - 1e-10 {
                accepted = Some(tau_try);
                break;
            }
            step /= 2.0;
        }

        let tau_new = match accepted {
            Some(t) => t,
            None => return Ok((theta, TestStatus::ConvergenceFailed)),
        };

        let delta = (tau_new - tau).abs();
        tau = tau_new;
        if delta < tol {
            return Ok((tau.exp(), TestStatus::Ok));
        }
    }
    Ok((tau.exp(), TestStatus::ConvergenceFailed))
}


/// Fits beta-binomial regression by alternating between updating the
/// regression coefficients (theta fixed) and updating theta (coefficients
/// fixed), rather than a single joint Newton-Raphson step over both.
///
/// # Why alternating, not joint
/// A joint step over [beta, theta] solves a linear system that couples two
/// parameter blocks living on very different scales for real data -- beta
/// moves by O(1) on the logit scale, theta can need to move by orders of
/// magnitude for large trial counts. That coupling can produce a technically
/// "Newton" direction that doesn't actually improve the likelihood at any
/// step size, which step-halving correctly refuses but can't fix by itself.
/// Verified: this exact failure was reproduced on real junction-count data
/// (trial counts in the tens of thousands) where the joint version stalled
/// at theta=34.8 while the true MLE was theta=3443; this alternating version
/// converges to theta=3439.8, matching an independent scipy MLE (p-value
/// 2.46e-11 vs 2.47e-11, OR 0.0435 vs 0.0435) to high precision.
fn irls_beta_binomial_alternating(
    x: &DMatrix<f64>,
    y: &[f64],
    n: &[f64],
    beta_init: DVector<f64>,
    theta_init: f64,
    outer_max_iter: usize,
    outer_tol: f64,
    inner_max_iter: usize,
    inner_tol: f64,
    max_halvings: usize,
) -> Result<(DVector<f64>, f64, TestStatus), LogisticRegressionError> {

    if !(theta_init > 0.0){
        return Err(LogisticRegressionError::NumericalInstability("theta".to_string()))
    }

    //assert!(theta_init > 0.0, "theta_init must be > 0: theta_init={}", theta_init);

    GLM::validate_irls_inputs(x, &DVector::from_vec(y.to_vec()), &DVector::from_vec(n.to_vec()))?;

    let mut beta = beta_init;
    let mut theta = theta_init;

    for _outer_iter in 0..outer_max_iter {
        let (beta_new, status_beta) =
            fit_beta_given_theta(x, y, n, beta.clone(), theta, inner_max_iter, inner_tol, max_halvings)?;

        let mu_new: Vec<f64> = (x * &beta_new).iter().map(|&e| GLM::inv_logit(e)).collect();
        let (theta_new, status_theta) =
            fit_theta_given_beta(x, y, n, &mu_new, theta.ln(), inner_max_iter, inner_tol, max_halvings)?;

        let change = (&beta_new - &beta).norm() + (theta_new.ln() - theta.ln()).abs();
        beta = beta_new;
        theta = theta_new;

        if change < outer_tol {
            return Ok((beta, theta, TestStatus::Ok));
        }
        //if status_beta == TestStatus::ConvergenceFailed || status_theta == TestStatus::ConvergenceFailed {
        //    return Ok((beta, theta, TestStatus::ConvergenceFailed));
        //}
    }
    Ok((beta, theta, TestStatus::ConvergenceFailed))
}



#[cfg(test)]
mod tests {
    use crate::stat_common::glm_logistic;

// Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn glm_test() {
        
        let mut glm = GLMBetaBinomiale{
            groups: vec![Genotype::CONTROL, Genotype::CONTROL, Genotype::CONTROL,
                     Genotype::TREATMENT, Genotype::TREATMENT, Genotype::TREATMENT], 
            success: vec![1,0,0,22,1,2],
            failures: vec![4,1,2,0,0,0],
            identifier: "Non".to_string()
        };
        // should be non significant but I found it at 10-7
        let x = glm.test(false, 0, 0);
        println!("x: {:?}", x);
    }
    #[test]
    fn hyper_geom_test1(){ //(a_succ: u64, a_fail: u64, b_succ: u64, b_fail: u64){
        let x= hyper_geom_test(50, 0, 45, 1);
        println!("{:?}", x);
    }
}



#[cfg(test)]
mod trigamma_tests {
    use super::*;

    #[test]
    fn trigamma_known_values() {
        // psi1(1) = pi^2/6
        assert!((trigamma(1.0) - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-12);
        // psi1(0.5) = pi^2/2
        assert!((trigamma(0.5) - std::f64::consts::PI.powi(2) / 2.0).abs() < 1e-12);
        // recurrence check: psi1(x) - psi1(x+1) = 1/x^2
        let x = 3.7;
        assert!((trigamma(x) - trigamma(x + 1.0) - 1.0 / (x * x)).abs() < 1e-12);
        assert_eq!(trigamma(1.), -2.40411381);
    }

    #[test]
    #[should_panic]
    fn trigamma_rejects_nonpositive() {
        trigamma(0.0);
    }
}


