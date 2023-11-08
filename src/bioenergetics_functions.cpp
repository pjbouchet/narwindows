// R and 'Eigen' integration using 'Rcpp'. 
// 'Eigen' is a C++ template library 
// for linear algebra: matrices, vectors, 
// numerical solvers and related algorithms.

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>
#include <RcppEigen.h>
#include <Rcpp.h> // For lists
#include <iostream>
#include <cstdlib>
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// // [[Rcpp::export]]
// double binom(double x, double mu, double sigma){
//  return R::pnorm5(x,mu,sigma,1,0);
// }


// 
// float startage(int cohort){
//   if(cohort == 1){
//     return R::runif(1/365, 1);
//   } else if(cohort >= 2 && cohort <= 3){
//     return R::runif(1, 9);
//   } else {
//     return R::runif(9, 69);
//   }
// }

// Vector form
// Rcpp::NumericVector startage(int n, int cohort){
//   if(cohort == 1){
//     return Rcpp::runif(n, 1/365, 1);
//   } else if(cohort >= 2 && cohort <= 3){
//     return Rcpp::runif(n, 1, 9);
//   } else {
//     return Rcpp::runif(n, 9, 69);
//   }
// }
// 
// double age2length_cpp(float age){  // Age in years
//   double a, b, c; // Parameters of the Gompertz growth curves (see Fortune et al. 2020)
//     if(age <= 0.79){
//       a = R::rnorm(1067.353, 20.479);
//       b = R::rnorm(-0.923, 0.058);
//       c = R::rnorm(-3.075, 0.315);
//     } else {
//       a = R::rnorm(1360.675, 19.501); 
//       b = R::rnorm(-0.361, 0.023); 
//       c = R::rnorm(-0.166, 0.026); 
//     }
//     return a * exp(b*exp(c*age));
// }
// 

// // [[Rcpp::export]]
// Eigen::MatrixXd agL_test(double age, int n = 1){
//   
//   double a, b, c; // Parameters of the Gompertz growth curves        
//   Eigen::MatrixXd out(n,3);
//   
//   for(int i = 0; i < n; i++){
//     
//     if(age <= 0.79){
//       a = R::rnorm(1067.19, 19.67);
//       b = R::rnorm(-0.93, 0.08);
//       c = R::rnorm(-3.11, 0.28);
//     } else {
//       a = R::rnorm(1362.755, 22.88); 
//       b = R::rnorm(-0.37, 0.03); 
//       c = R::rnorm(-0.18, 0.03); 
//     }
//     
//     out(i,0) = a;
//     out(i,1) = b;
//     out(i,2) = c;
//     
//   }
//   return(out);
// }
// 
// // [[Rcpp::export]]
// Rcpp::NumericVector random_multivariate_normal_test(const Eigen::MatrixXd mu, const Eigen::MatrixXd Sigma){
//   
//   Rcpp::NumericVector out(2);
//   int P = mu.rows(), i = 0;
//   Eigen::MatrixXd y(Eigen::MatrixXd(P, 1).setZero()); // Create a column matrix of zeroes
//   Eigen::MatrixXd z(Eigen::MatrixXd(P, 1).setZero());
//   
//   for(i = 0 ; i < P ; i++) z(i, 0) = Rf_rnorm(0, 1); // To get original R api, use Rf_*
//   
//   y = mu + Sigma.llt().matrixL() * z;
//   
//   out[0] = y(0,0);
//   out[1] = y(1,0);
//   return out;
// }
// 
// // [[Rcpp::export]]
// Eigen::MatrixXd mL_test(int n = 1){
//   
//   Eigen::MatrixXd out = Eigen::MatrixXd(n,2);
//   
//   // Multivariate parameters of the mass-at-length function
//   // See growth_curves.R script
//   Eigen::MatrixXd mu = Eigen::MatrixXd(2,1);
//   Eigen::MatrixXd sigma = Eigen::MatrixXd(2,2);
//   
//   // Means
//   mu(0,0) = -4.834189;
//   mu(1,0) = 2.984353;
//   
//   // Variance-covariance matrix
//   sigma(0,0) = 0.21128304;
//   sigma(0,1) = -0.07515154;
//   sigma(1,0) = -0.07515154;
//   sigma(1,1) = 0.02686353;
//   
//   for (int i = 0; i < n; i++) {
//     Rcpp::NumericVector deviates = random_multivariate_normal_test(mu, sigma);
//     out(i,0) = deviates[0];
//     out(i,1) = deviates[1];
//   }
//   return(out);
// }

// // [[Rcpp::export]]
// Rcpp::NumericVector age2length_v(Rcpp::NumericVector age){  // Age in years
//   int n = age.size();
//   Rcpp::NumericVector L (n);
//   double a, b, c; // Parameters of the Gompertz growth curves (see Fortune et al. 2020)
//   for(int i = 0; i < n; ++i){
//     if(age[i] > 0.79){
//       a = 1360.68;
//       b = -0.36;
//       c = -0.16;
//     } else {
//       a = 1067.35;
//       b = -0.923;
//       c = -3.08;
//     }
//     L[i] = (a * std::exp(b*std::exp(c*age[i])));
//   }
//   return L;
// }
// 
// double sumC(Rcpp::NumericVector x) {
//   int n = x.size();
//   double total = 0;
//   for(int i = 0; i < n; ++i) {
//     total += x[i];
//   }
//   return total;
// }
// 

// 
// float length2mass_cpp(float L, float age){  // Body length in cm
//   float a, b; // Intercept and slope of the logarithmic mass-at-length relationship (see Fortune et al. 2020)
//   if(age <= 0.79){
//     a = R::rnorm(-5.091821, 0.2578327);
//     b = R::rnorm(3.077823, 0.08325852);
//   } else {
//     a = R::rnorm(-5.096379, 0.2592405);
//     b = R::rnorm(3.079408, 0.08360103);
//   }
//   return pow(10, a + b*log10(L));
// }

// // [[Rcpp::export]]
// Rcpp::NumericVector start_percfat(int n){ // Body condition (% fat mass)
//  return Rcpp::rbeta(n, 6, 20);
// }
// 
// 
// // [[Rcpp::export]]
// Rcpp::NumericVector age2mass(Rcpp::NumericVector age){  // Body length in cm
//   double a, b; // Intercept and slope of the logarithmic mass-at-length relationship (see Fortune et al. 2020)
//   int n = age.size();
// 
//   Rcpp::NumericVector L = age2length_v(age);
//   Rcpp::NumericVector mass (n);
// 
//   for(int i = 0; i < n; ++i){
// 
//     if(age[i] > 0.79){
//       a = -5.0964;
//       b = 3.0794;
//     } else {
//       a = -5.0918;
//       b = 3.0778;
//     }
//     mass[i] = pow(10, a + b*log10(L[i]));
//   }
//   return mass;
// }
// 
// 
// // [[Rcpp::export]]
// Rcpp::NumericVector construct_array(Rcpp::NumericVector a, Rcpp::IntegerVector dim) {
//   a.attr("dim") = dim;
//   return a;
// }


// 
// // [[Rcpp::export]]
// Rcpp::NumericVector first_five(Rcpp::NumericVector x) {
//   Rcpp::IntegerVector idx = Rcpp::IntegerVector::create(0, 1, 2, 3, 4);
//   return x[idx];
// }

// Rcpp::NumericVector construct_array2(Rcpp::NumericVector v, int nparam, int nsim) {
//   Rcpp::NumericVector zeroes (nparam * nsim * 364);
//   Rcpp::IntegerVector dim {nparam, nsim, 365};
//   int n = zeroes.size();
//   for (int i=0; i < n; ++i) {
//     v.push_back(zeroes[i] );
//   }
//   v.attr("dim") = dim;
//   return v;
// }
  