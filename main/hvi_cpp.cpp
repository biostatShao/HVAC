// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List vi_bayes_paper(
    const arma::vec& mu_afr, const arma::vec& sigma_afr,
    const arma::vec& mu_eur, const arma::vec& sigma_eur,
    int max_iter = 100, double tol = 1e-6,
    double sigma0_sq = 1e6, double a = 0.001, double b = 0.001
) {
  
  const int N = mu_afr.n_elem;
  const vec sigma_afr_sq = square(sigma_afr);
  const vec sigma_eur_sq = square(sigma_eur);
  
  vec m = mu_afr;
  vec s2 = sigma_afr_sq;
  vec n = mu_eur;
  vec t2 = sigma_eur_sq;
  
  double mu_q = 0.5 * (mean(m) + mean(n));
  double sigma_q_sq = 1e6;
  double tau_sq = 1.0;
  
  double elbo_prev = -datum::inf;
  bool converged = false;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    for (int i = 0; i < N; ++i) {
      double sigma_inv = 1.0 / sigma_afr_sq[i];
      double tau_inv = 1.0 / tau_sq;
      m[i] = (sigma_inv * mu_afr[i] + tau_inv * mu_q) / (sigma_inv + tau_inv);
      s2[i] = 1.0 / (sigma_inv + tau_inv);
    }
    
    for (int i = 0; i < N; ++i) {
      double sigma_inv = 1.0 / sigma_eur_sq[i];
      double tau_inv = 1.0 / tau_sq;
      n[i] = (sigma_inv * mu_eur[i] + tau_inv * mu_q) / (sigma_inv + tau_inv);
      t2[i] = 1.0 / (sigma_inv + tau_inv);
    }
    
    double sum_m = sum(m);
    double sum_n = sum(n);
    mu_q = (sum_m + sum_n) / (2 * N + 1e-10);
    sigma_q_sq = 1.0 / (1.0 / sigma0_sq + 2 * N / tau_sq);
    
    double numerator = sum(s2 + square(m - mu_q)) + sum(t2 + square(n - mu_q)) + 2 * b;
    tau_sq = numerator / (2 * N + 2 * a + 2);
    
    double elbo = 0.0;
    
    for (int i = 0; i < N; ++i) {
      elbo += -0.5 * log(2 * M_PI * sigma_afr_sq[i]) 
      - 0.5 * (pow(mu_afr[i] - m[i], 2) + s2[i]) / sigma_afr_sq[i];
      elbo += -0.5 * log(2 * M_PI * sigma_eur_sq[i]) 
        - 0.5 * (pow(mu_eur[i] - n[i], 2) + t2[i]) / sigma_eur_sq[i];
    }
    
    double tau_inv = 1.0 / tau_sq;
    for (int i = 0; i < N; ++i) {
      elbo += -0.5 * log(2 * M_PI * tau_sq) 
      - 0.5 * tau_inv * (pow(m[i] - mu_q, 2) + s2[i]);
      elbo += -0.5 * log(2 * M_PI * tau_sq) 
        - 0.5 * tau_inv * (pow(n[i] - mu_q, 2) + t2[i]);
    }
    
    elbo += -0.5 * log(2 * M_PI * sigma0_sq) 
      - 0.5 * (pow(mu_q, 2) + sigma_q_sq) / sigma0_sq;
    elbo += a * log(b) - lgamma(a) 
      - (a + 1) * log(tau_sq) - b * tau_inv;
    
    for (int i = 0; i < N; ++i) {
      elbo += 0.5 * log(2 * M_PI * s2[i]) + 0.5;
      elbo += 0.5 * log(2 * M_PI * t2[i]) + 0.5;
    }

    if (iter > 0 && std::abs(elbo - elbo_prev) < tol) {
      converged = true;
      break;
    }
    elbo_prev = elbo;
  }
  
  return List::create(
    Named("mu_global") = mu_q,
    Named("tau_sq") = tau_sq,
    Named("beta_afr") = m,
    Named("beta_afr_var") = s2,
    Named("elbo") = elbo_prev
  );
}
