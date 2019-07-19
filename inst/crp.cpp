//#########################################################
//#                                                       #
//#            crp method : rcpp code                     #
//#                                                       #
//#########################################################
#define ARMA_DONT_PRINT_ERRORS
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp/Benchmark/Timer.h>
#include <string>
#include <iostream>
#include <vector>
#include <cstring>
// #include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


//#############################################################################
//############################ basic functions  ###############################
//#############################################################################

// Function seq
arma::vec sequence_le(int a,int b,int le){
 arma::vec res = arma::zeros<arma::vec>(le);
 if(le == 1){
  res(0) = a;
  return(res);
 }else{
  double a_tmp = a;
  double b_tmp = b;
  double by = (b_tmp-a_tmp)/(le-1) ;

  for(int i=0 ; i<le ; i++){
   res(i) = round(a_tmp + i*by);
  }
  res(le-1) = b;

  return res;
 }
}

// Sample in 0:n-1.
int sample_cpp (int n){
 double u = R::runif(0,1) * n;
 int res = trunc(u) ;
 if(res == n) res = n-1;
 return res;
}

// Vector sample in 0:n-1.
arma::vec sample_cpp (int nbre, int n){
 arma::vec res = arma::zeros<arma::vec>(nbre);
 for (int i = 0 ; i < nbre ; ++i) {
  res(i) = sample_cpp(n);
 }
 return res;
}

// Weighted sample in 0:n-1.
int sample_weight (arma::vec proba){
 if(sum(proba) == 0)   proba = ones<arma::vec>(proba.n_rows) ;
 arma::vec proba_cum = cumsum(proba)/sum(proba) ;

 unsigned ret = 0;
 double u = R::runif(0,1);
 while (ret <= proba_cum.n_rows and u > proba_cum(ret)) {
  ret++;
 }
 return ret;
}

// Vector sample in 0:n-1.
arma::vec sample_no_replace_cpp (int nbre, int n){
 arma::vec res(nbre);
 arma::vec probs = arma::ones<arma::vec>(n) ;
 for (int i = 0 ; i < nbre ; ++i) {
  res(i) = sample_weight(probs);
  probs(res(i)) = 0 ;
 }
 return res;
}

// Vector weighted sample in 0:n-1.
arma::vec sample_weight (int n, arma::vec & proba){
 arma::vec res = arma::zeros<arma::vec>(n);
 for (unsigned i = 0 ; i < n ; ++i) {
  res(i) = sample_weight(proba);
 }
 return res;
}

//
int is_nonnegative_num(int val) {
 if (val > 0) return(1) ; else return(0) ;
}

//
arma::vec is_nonnegative_vec(arma::vec & v){
 int N = v.size();
 arma::vec res(N) ;

 for(int i=0 ; i<N ; ++i){
  res(i) = is_nonnegative_num(v(i)) ;
 }

 return(res) ;
}

//
arma::vec which_nonnegative(arma::vec & v){
 arma::vec v_nonneg = is_nonnegative_vec(v);
 double size = sum(v_nonneg) ;
 arma::vec res(size);

 double n = v.size();
 int j = 0;
 for(int i=0 ; i<n ; ++i){
  if(v_nonneg(i) == 1){
   res(j) = i ;
   j += 1 ;
  }
  // Rcpp::Rcout << res <<  std::endl;
 }
 return(res) ;
}

//
arma::vec which_empty_table(arma::vec & v){
 arma::vec v_zero = is_nonnegative_vec(v) * (-1) + 1;
 int size = sum(v_zero) ;
 arma::vec res(size);

 // if(size != 0){
 // arma::vec res_tmp(size);
 // res = res_tmp;

 int n = v.size();
 int j = 0;
 for(int i=0 ; i<n ; ++i){
  if(v_zero(i) == 1){
   res(j) = i ;
   j += 1 ;
  }
  // Rcpp::Rcout << res <<  std::endl;
 }
 // }

 return(res) ;
}


//
arma::mat expand_mat (arma::mat & m, int n){
 arma::mat res =  arma::zeros<arma::mat>(m.n_rows,n) ;

 res.submat(0,0,
            m.n_rows-1,m.n_cols-1) = m ;
 return(res) ;
}

//
arma::cube expand_cube (arma::cube & cub, int n){
 arma::cube res =  arma::zeros<arma::cube>(cub.n_rows,cub.n_cols,n) ;

 res.subcube(0,0,0,
             cub.n_rows-1,cub.n_cols-1,cub.n_slices-1) = cub;

 return(res) ;
}

//
arma::mat reduce (arma::mat & m){
 int old_size = m.n_cols ;
 arma::mat res;

 int index = old_size;
 for(int k=old_size-1 ; k>=0 ; --k){
  if(sum(m.col(k))!=0){
   index = k;
   break ;
  }
 }
 if(index < old_size){
  res = m.cols(0,index);
 }else{
  res = m ;
 }

 return(res) ;
}

//
arma::cube reduce_cube (arma::cube & cub){
 int old_size = cub.n_slices ;
 arma::cube res;
 arma::mat m_tmp ;
 arma::rowvec v_tmp ;

 int index = old_size;
 for(int k=old_size-1 ; k>=0 ; --k){
  m_tmp = cub.slice(k);
  v_tmp = sum(m_tmp) ;
  if(sum(v_tmp)!=0){
   index = k;
   break ;
  }
 }
 if(index < old_size){
  res = cub.slices(0,index);
 }else{
  res = cub ;
 }

 return(res) ;
}

//#############################################################################
//############################ crp functions  #################################
//#############################################################################

// simulate a crp : assign
// [[Rcpp::export]]
arma::vec rcrp_assign_cpp(double alpha, int N){
 arma::vec crp = arma::zeros<arma::vec>(N);
 arma::vec probs = arma::zeros<arma::vec>(N);

 int new_customer=1;
 int K = 1;

 probs(0) = alpha ;
 crp(0)   = new_customer ;
 probs(1) = 1.0 ;


 for(int i=1 ; i<N ; ++i){
  new_customer = sample_weight(probs) ;

  if(new_customer == 0){
   K += 1 ;
   new_customer = K ;
  }
  crp(i) = new_customer ;
  probs(new_customer) += 1 ;
 }

 return(crp) ;
}

// simulate a crp : size
// [[Rcpp::export]]
arma::vec rcrp_size_cpp(double alpha, int N){
 arma::vec crp = arma::zeros<arma::vec>(N);
 arma::vec probs = arma::zeros<arma::vec>(N);

 int new_customer=1;
 int K = 1;

 probs(0) = alpha ;
 crp(0)   = new_customer ;
 probs(1) = 1.0 ;


 for(int i=1 ; i<N ; ++i){
  new_customer = sample_weight(probs) ;

  if(new_customer == 0){
   K += 1 ;
   new_customer = K ;
  }
  crp(i) = new_customer ;
  probs(new_customer) += 1 ;
 }

 arma::vec crp_size = arma::zeros<arma::vec>(K);
 for(int i=0 ; i<N ; ++i){
  crp_size(crp(i)-1) += 1 ;
 }

 return(crp_size) ;
}

// convert an assigned crp to a crp size
// [[Rcpp::export]]
arma::vec size_to_assign(arma::vec & crp_size){
 int K = crp_size.size() ;
 int N = sum(crp_size);
 arma::vec crp_assign(N);
 arma::vec crp_size_tmp = crp_size;

 arma::vec which_nonneg = which_nonnegative(crp_size_tmp) ;
 arma::vec probs = arma::zeros<arma::vec>(K) ;

 crp_assign(0) = which_nonneg(0);
 crp_size_tmp(which_nonneg(0)) -= 1;

 if(crp_size_tmp(which_nonneg(0)) != 0)
  probs(which_nonneg(0)) = 1 ;
 probs(which_nonneg(1)) = 1 ;

 int N_tmp = which_nonneg.size() ;
 int tmp = 1;

 for(int i=1 ; i<N ; ++i){
  crp_assign(i) = sample_weight(probs) ;
  crp_size_tmp(crp_assign(i)) -= 1 ;

  if(crp_size_tmp(crp_assign(i)) == 0){
   probs(crp_assign(i)) = 0;
  }
  if(crp_size_tmp(crp_assign(i)) == crp_size(crp_assign(i))-1){
   tmp += 1 ;
   if(tmp < N_tmp)
    probs(which_nonneg(tmp)) = 1 ;
  }
 }

 return(crp_assign+1);
}

// convert an assigned crp to a crp size
// [[Rcpp::export]]
arma::vec assign_to_size(arma::vec & crp_assign){
 int K = crp_assign.max();
 int N = crp_assign.size();
 arma::vec crp_size = arma::zeros<arma::vec>(K);

 for(int i=0 ; i<N ; ++i){
  crp_size(crp_assign(i)-1) += 1 ;
 }

 return(crp_size) ;
}

// [[Rcpp::export]]
arma::vec subcrp_size_cpp(arma::vec & crp_size,  double prob) {
 arma::vec v(crp_size.size());

 for(int k=0 ; k<v.size() ; ++k){
  v(k) = R::rbinom(crp_size(k),prob) ;
 }

 // std::transform( crp_size.begin(), crp_size.end(), v.begin(),
 //                 [=](double size){ return R::rbinom(size, prob); });
 return(v);
}

// [[Rcpp::export]]
arma::vec augm_crp_assign_cpp(arma::vec & crp_assign, double alpha, int N_tilde) {
 int N = crp_assign.size() ;
 if(N_tilde > N){
  arma::vec crp(N_tilde);
  arma::vec crp_size = assign_to_size(crp_assign);
  int K = crp_size.size();

  arma::vec probs = arma::zeros<arma::vec>(N_tilde+1);

  crp.subvec(0,N-1) = crp_assign;
  probs(0) = alpha ;
  probs.subvec(1,K) = crp_size;

  arma::vec empty_tables = which_empty_table(crp_size);
  // Rcpp::Rcout << empty_tables <<  std::endl;

  int new_customer;
  int index_new_add;

  for(int i=N ; i<N_tilde ; ++i){
   new_customer = sample_weight(probs) ;

   if(new_customer == 0){
    if(empty_tables.size() != 0){
     index_new_add = sample_cpp(empty_tables.size()+1) ;
     // Rcpp::Rcout << "index_new_add = " << index_new_add <<  std::endl;
     if(index_new_add == 0){
      K += 1 ;
      new_customer = K;
     }else{
      new_customer = empty_tables(index_new_add-1)+1;
      empty_tables.shed_row(index_new_add-1) ;
      // Rcpp::Rcout << "new_customer = " << new_customer <<  std::endl;
     }
    }else{
     // Rcpp::Rcout << "new table "<<  std::endl;
     K += 1 ;
     new_customer = K ;
    }
   }

   // Rcpp::Rcout << "new_customer = " << new_customer <<  std::endl;

   crp(i) = new_customer ;
   probs(new_customer) += 1 ;
   // Rcpp::Rcout << "probs = " << probs <<  std::endl;

  }

  return(crp) ;
 }else{
  return(crp_assign);
 }
}

// [[Rcpp::export]]
arma::vec augm_crp_size_cpp(arma::vec & crp_size, double alpha, int N_tilde) {
 arma::vec crp_assign = size_to_assign(crp_size) ;
 int N = crp_assign.size() ;

 if(N_tilde > N){
  arma::vec crp(N_tilde);
  int K = crp_size.size();

  arma::vec probs = arma::zeros<arma::vec>(N_tilde+1);

  crp.subvec(0,N-1) = crp_assign;
  probs(0) = alpha ;
  probs.subvec(1,K) = crp_size;

  arma::vec empty_tables = which_empty_table(crp_size);
  // Rcpp::Rcout << empty_tables <<  std::endl;

  int new_customer;
  int index_new_add;

  for(int i=N ; i<N_tilde ; ++i){
   new_customer = sample_weight(probs) ;

   if(new_customer == 0){
    if(empty_tables.size() != 0){
     index_new_add = sample_cpp(empty_tables.size()+1) ;
     // Rcpp::Rcout << "index_new_add = " << index_new_add <<  std::endl;
     if(index_new_add == 0){
      K += 1 ;
      new_customer = K;
     }else{
      new_customer = empty_tables(index_new_add-1)+1;
      empty_tables.shed_row(index_new_add-1) ;
      // Rcpp::Rcout << "new_customer = " << new_customer <<  std::endl;
     }
    }else{
     // Rcpp::Rcout << "new table "<<  std::endl;
     K += 1 ;
     new_customer = K ;
    }
   }

   // Rcpp::Rcout << "new_customer = " << new_customer <<  std::endl;

   crp(i) = new_customer ;
   probs(new_customer) += 1 ;
   // Rcpp::Rcout << "probs = " << probs <<  std::endl;

  }

  return(crp) ;
 }else{
  return(crp_assign);
 }
}

//#############################################################################
//####################### Functions for Sampler  ##############################
//#############################################################################

// [[Rcpp::export]]
List rprior_cpp (List & data, List & hyperparam){

 double lambda1 = hyperparam["lambda1"];
 double lambda2 = hyperparam["lambda2"];
 int kappa = hyperparam["kappa"];
 double pi = hyperparam["pi"];
 double eta1 = hyperparam["eta1"];
 double eta2 = hyperparam["eta2"];

 arma::vec c = data["c"];
 int N = data["N"];

 double alpha = R::rgamma(lambda1,1./lambda2) ;
 int r = R::rnbinom( kappa, 1-pi ) ;
 double p = R::rbeta( eta1, eta2) ;
 int N_tilde = N + r ;

 arma::vec c_tilde = augm_crp_size_cpp(c,alpha,N_tilde) ;

 arma::vec c_tilde_size = assign_to_size(c_tilde) ;
 c_tilde_size = reverse(sort(c_tilde_size));
 int K_tilde = c_tilde_size.size() ;

 return  List::create(_["alpha"]=alpha,
                      _["c_tilde"]=c_tilde_size,
                      _["c_tilde_assign"]=c_tilde,
                      _["K_tilde"]=K_tilde,
                      _["p"]=p,
                      _["r"]=r
 );
}

// [[Rcpp::export]]
List compute_starting_point_cpp (List & data, List & hyperparam){
 List res;

 res = rprior_cpp(data,hyperparam);

 return(res);
}

// [[Rcpp::export]]
double log_dlkh_cpp (List & theta, List & data, double temperature) {
 double res;

 int N = data["N"];
 int K = data["K"];
 arma::vec c = data["c"];

 // double alpha = theta["alpha"];
 double p = theta["p"];
 double r = theta["r"];
 arma::vec c_tilde = theta["c_tilde"] ;

 double choose_tmp = 0;
 for( int k=0 ; k<K ; ++k){
  choose_tmp += R::lchoose(c_tilde(k),c(k)) ;
  // Rcpp::Rcout << "choose_tmp = "<< choose_tmp << std::endl;
 }

 res = temperature * (
  N * log(p) + r * log(1-p) + choose_tmp
 ) ;

 return(res);
}

// [[Rcpp::export]]
arma::mat log_dtarget_cpp (List & theta, List & data, List & hyperparam){
 arma::mat res(1,6);

 int N = data["N"];
 int K = data["K"];
 arma::vec c = data["c"];

 double p = theta["p"];
 double r = theta["r"];
 double alpha = theta["alpha"];
 arma::vec c_tilde = theta["c_tilde"] ;
 int K_tilde = theta["K_tilde"] ;

 int kappa = hyperparam["kappa"] ;
 double pi = hyperparam["pi"] ;
 double lambda1 = hyperparam["lambda1"] ;
 double lambda2 = hyperparam["lambda2"] ;
 double eta1 = hyperparam["eta1"] ;
 double eta2 = hyperparam["eta2"] ;

 double choose_tmp = 0;
 for( int k=0 ; k<K ; ++k){
  choose_tmp += R::lchoose(c_tilde(k),c(k)) ;
  // Rcpp::Rcout << "choose_tmp = "<< choose_tmp << std::endl;
 }

 double lgamma_tmp = 0 ;
 for( int k=0 ; k<K_tilde ; ++k){
  lgamma_tmp += lgamma(c_tilde(k)) ;
 }

 res(0,2) = N * log(p) + r * log(1-p) + choose_tmp;
 res(0,4) = K_tilde * log(alpha) + lgamma_tmp + lgamma(alpha) - // c_tilde prior
  lgamma(alpha + r + N) +
  R::lchoose(kappa + r - 1 , r) + kappa * log(1-pi) + r * log(pi) + // r prior
  lambda1 * log(lambda2) - lgamma(lambda1) + (lambda1 -1) * log(alpha) - // alpha prior
  lambda2 * alpha -
  R::lbeta(eta1,eta2) + (eta1-1) * log(p) + (eta2-1) *log(1-p) ; // p prior

 res(0,0) = res(0,2) + res(0,4);
 res(0,1) = std::exp(res(0,0));
 res(0,3) = std::exp(res(0,2));
 res(0,5) = std::exp(res(0,4));

 return(res);
}

// [[Rcpp::export]]
arma::mat log_dtarget_sample_cpp (List & trace, List & data, List & hyperparam){
 arma::vec alpha = trace["alpha"] ;
 int n = alpha.size() ;

 arma::mat res(n,6) ;

 int N = data["N"];
 int K = data["K"];
 arma::vec c = data["c"];

 arma::vec p = trace["p"];
 arma::vec r = trace["r"];
 arma::mat c_tilde = trace["c_tilde"] ;
 arma::vec K_tilde = trace["K_tilde"] ;

 int kappa = hyperparam["kappa"] ;
 double pi = hyperparam["pi"] ;
 double lambda1 = hyperparam["lambda1"] ;
 double lambda2 = hyperparam["lambda2"] ;
 double eta1 = hyperparam["eta1"] ;
 double eta2 = hyperparam["eta2"] ;

 for( int i=0 ; i<n ; ++i ){
  double choose_tmp = 0;
  for( int k=0 ; k<K ; ++k){
   choose_tmp += R::lchoose(c_tilde(i,k),c(k)) ;
   // Rcpp::Rcout << "choose_tmp = "<< choose_tmp << std::endl;
  }

  double lgamma_tmp = 0 ;
  for( int k=0 ; k<K_tilde(i) ; ++k){
   lgamma_tmp += lgamma(c_tilde(i,k)) ;
  }

  res(i,2) = N * log(p(i)) + r(i) * log(1-p(i)) + choose_tmp;
  res(i,4) = K_tilde(i) * log(alpha(i)) + lgamma_tmp + lgamma(alpha(i)) - // c_tilde prior
   lgamma(alpha(i) + r(i) + N) +
   R::lchoose(kappa + r(i) - 1 , r(i)) + kappa * log(1-pi) + r(i) * log(pi) + // r prior
   lambda1 * log(lambda2) - lgamma(lambda1) + (lambda1 -1) * log(alpha(i)) - // alpha prior
   lambda2 * alpha(i) -
   R::lbeta(eta1,eta2) + (eta1-1) * log(p(i)) + (eta2-1) *log(1-p(i)) ; // p prior

  res(i,0) = res(i,2) + res(i,4);
  res(i,1) = std::exp(res(i,0));
  res(i,3) = std::exp(res(i,2));
  res(i,5) = std::exp(res(i,4));
 }

 return(res) ;
}

// [[Rcpp::export]]
arma::mat log_dtarget_sample_cpp_PT (List & trace, List & data, List & hyperparam,
                                     int which_chain ){
 arma::mat alpha_mat = trace["alpha"] ;
 arma::vec alpha = alpha_mat.col(which_chain-1) ;
 int n = alpha.size() ;

 arma::mat res(n,6) ;

 int N = data["N"];
 int K = data["K"];
 arma::vec c = data["c"];

 arma::mat p_mat = trace["p"];
 arma::mat r_mat = trace["r"];
 arma::mat K_tilde_mat = trace["K_tilde"] ;
 arma::cube c_tilde_cube = trace["c_tilde"] ;

 arma::vec p = p_mat.col(which_chain-1) ;
 arma::vec r = r_mat.col(which_chain-1) ;
 arma::vec K_tilde = K_tilde_mat.col(which_chain-1) ;
 arma::mat c_tilde = c_tilde_cube.col(which_chain-1) ;

 int kappa = hyperparam["kappa"] ;
 double pi = hyperparam["pi"] ;
 double lambda1 = hyperparam["lambda1"] ;
 double lambda2 = hyperparam["lambda2"] ;
 double eta1 = hyperparam["eta1"] ;
 double eta2 = hyperparam["eta2"] ;

 for( int i=0 ; i<n ; ++i ){
  double choose_tmp = 0;
  for( int k=0 ; k<K ; ++k){
   choose_tmp += R::lchoose(c_tilde(i,k),c(k)) ;
   // Rcpp::Rcout << "choose_tmp = "<< choose_tmp << std::endl;
  }

  double lgamma_tmp = 0 ;
  for( int k=0 ; k<K_tilde(i) ; ++k){
   lgamma_tmp += lgamma(c_tilde(i,k)) ;
  }

  res(i,2) = N * log(p(i)) + r(i) * log(1-p(i)) + choose_tmp;
  res(i,4) = K_tilde(i) * log(alpha(i)) + lgamma_tmp + lgamma(alpha(i)) - // c_tilde prior
   lgamma(alpha(i) + r(i) + N) +
   R::lchoose(kappa + r(i) - 1 , r(i)) + kappa * log(1-pi) + r(i) * log(pi) + // r prior
   lambda1 * log(lambda2) - lgamma(lambda1) + (lambda1 -1) * log(alpha(i)) - // alpha prior
   lambda2 * alpha(i) -
   R::lbeta(eta1,eta2) + (eta1-1) * log(p(i)) + (eta2-1) *log(1-p(i)) ; // p prior

  res(i,0) = res(i,2) + res(i,4);
  res(i,1) = std::exp(res(i,0));
  res(i,3) = std::exp(res(i,2));
  res(i,5) = std::exp(res(i,4));
 }

 return(res) ;
}

// [[Rcpp::export]]
List MH_proposal_cpp (List & theta, List & data, List & hyperparam, double temperature){

 List theta_star = rprior_cpp(data,hyperparam);
 double ratio_log_target = log_dlkh_cpp(theta_star, data, temperature) -
  log_dlkh_cpp(theta, data, temperature) ;

 double acc_ratio=1;
 if(ratio_log_target < 0){
  acc_ratio= exp(ratio_log_target);
 }

 return  List::create(_["theta_star"]=theta_star,
                      _["acc_ratio"]=acc_ratio
 );
}

//#############################################################################
//############################## Sampler  #####################################
//#############################################################################

// [[Rcpp::export]]
List MH_cpp (List & param, List & data, List & hyperparam, bool verbose) {
 // Initialization
 // int burnin;
 // if( param.containsElementNamed("burnin") ){
 //  burnin = param["burnin"];
 // }else{
 //  burnin = 0 ;
 // }

 std::string target ;
 if( param.containsElementNamed("target") ){
  target = as<std::string>(param["target"]) ;
 }else{
  target = "posterior" ;
 }

 double temperature ;
 if( param.containsElementNamed("temperature") ){
  temperature = param["temperature"] ;
 }else{
  temperature = 1 ;
 }

 if ( target == "prior") {
  temperature = 0 ;
 }

 int n_iter ;
 if( param.containsElementNamed("n_iter") ){
  n_iter = param["n_iter"];
 }else{
  stop("The option 'n_iter' is missing in the 'param' input.\n");
 }

 int nbre_print;
 if( param.containsElementNamed("nbre_print") ){
  nbre_print = param["nbre_print"];
 }else{
  nbre_print = 100 ;
 }
 if(nbre_print>=n_iter) nbre_print = n_iter-1;


 int N = data["N"] ;
 double pi = hyperparam["pi"] ;
 List res_proposal ;
 List theta_star ;
 double u;


 // Compute a starting point
 List theta ;
 if( param.containsElementNamed("starting_point") ){
  theta = param["starting_point"] ;
 }else{
  theta = rprior_cpp(data,hyperparam) ;
 }

 // Initialize the trace
 arma::vec acc_ratio( n_iter ) ;
 arma::vec accepted = arma::zeros<arma::vec>( n_iter ) ;

 arma::vec r( n_iter ) ;
 arma::vec p( n_iter) ;
 arma::vec K_tilde( n_iter ) ;
 arma::vec alpha( n_iter ) ;

 arma::mat c_tilde_assign = arma::zeros<arma::mat>( n_iter , floor(N/(1-pi)) ) ;
 arma::mat c_tilde_size = arma::zeros<arma::mat>( n_iter , floor(N/(1-pi)) ) ;
 arma::vec c_tilde_assign_tmp;
 arma::vec c_tilde_size_tmp;

 // Compute what to verbose
 arma::mat n_iters;
 Timer timer;

 if(verbose){
  int tmp = ceil(n_iter / (nbre_print+1)) ;
  int tmp2 = n_iter % tmp ;

  if(tmp2 == 0){
   n_iters =  arma::mat(nbre_print+1,2);
  }else{
   n_iters =  arma::mat(nbre_print+2,2);
  }

  for(int k=0 ; k<n_iters.n_rows  ; ++k){
   n_iters(k,0) = k*tmp ;
   n_iters(k,1) = (k+1)*tmp-1 ;
  }
  if(tmp2 != 0){
   n_iters(nbre_print+1,0) = (nbre_print+1)*tmp ;
   n_iters(nbre_print+1,1) = n_iter-1 ;
  }
 }else{
  n_iters = arma::mat(1,2);
  n_iters(0,0) = 0;
  n_iters(0,1) = n_iter -1;
 }

 arma::vec percents = sequence_le(0,100,nbre_print+1) ;
 std::string to_print = "" ;
 std::string to_add = "=" ;
 CharacterVector string_tmp ;
 arma::vec when_to_print = sequence_le(0,n_iters.n_rows-1,20) ;
 int count = 0;

 NumericVector timer_vec ;
 double start_time ;
 double end_time ;
 double remain_time ;
 double time_mean = 0;

 // Start the loop
 for(int j=0 ; j<n_iters.n_rows ; ++j){
  // Stuff to print
  if(verbose){
   if(j == when_to_print(count)){
    string_tmp  = wrap(to_print + to_add) ;
    to_print = string_tmp(0);
    count += 1;
   }
   if(j!=0){
    timer_vec = NumericVector(timer);
    start_time = timer_vec(2*j-2) ;
    end_time = timer_vec(2*j-1) ;

    time_mean = (time_mean * (j-1) + (end_time - start_time)) / j ;
    remain_time = time_mean * (n_iters.n_rows-j) / 1000000000 ;

    if(j<nbre_print+1){
     std::stringstream strs;
     if(remain_time  > 60){
      if(remain_time > 3600){
       strs << "[" << to_print << " "<<  percents(j) << "%]  " <<
        "(remaining: " << remain_time  / 3600<< " hours)" ;
      }
      strs << "[" << to_print << " "<<  percents(j) << "%]  " <<
       "(remaining: " << remain_time  / 60 << " mins)" ;
     }else{
      strs << "[" << to_print << " "<<  percents(j) << "%]  " <<
       "(remaining: " << remain_time  << " seconds)" ;
     }
     std::string temp_str = strs.str();
     char const* char_type = temp_str.c_str();

     REprintf("\r");
     REprintf("%s", char_type);
     REprintf("\r");

     R_FlushConsole();
     R_CheckUserInterrupt();
    }
   }
  }

  // The MH iteration
  timer.step("step");
  for(int i=n_iters(j,0) ; i<=n_iters(j,1) ; ++i){
   res_proposal = MH_proposal_cpp(theta, data, hyperparam, temperature) ;
   theta_star = res_proposal["theta_star"] ;
   acc_ratio(i) = res_proposal["acc_ratio"] ;

   // Accept/reject
   u = R::runif(0,1) ;
   if(u < acc_ratio(i)){
    theta = theta_star ;
    accepted(i)  = 1 ;
   }

   // Update traces ;
   r(i) = theta["r"];
   p(i) = theta["p"];
   K_tilde(i) = theta["K_tilde"];
   alpha(i) = theta["alpha"];


   c_tilde_assign_tmp = as<arma::vec>(theta["c_tilde_assign"]) ;
   if(c_tilde_assign_tmp.size() > c_tilde_assign.n_cols){
    c_tilde_assign = expand_mat(c_tilde_assign,c_tilde_assign_tmp.size()) ;
   }
   c_tilde_assign.row(i).subvec(0,c_tilde_assign_tmp.size()-1) =
    trans(c_tilde_assign_tmp) ;

   c_tilde_size_tmp = as<arma::vec>(theta["c_tilde"]) ;
   if(c_tilde_size_tmp.size() > c_tilde_size.n_cols){
    c_tilde_size = expand_mat(c_tilde_size,c_tilde_size_tmp.size()) ;
   }
   c_tilde_size.row(i).subvec(0,c_tilde_size_tmp.size()-1) =
    trans(c_tilde_size_tmp) ;

  }
  timer.step("step");
 }

 // Rearrange c_tilde traces
 c_tilde_assign = reduce(c_tilde_assign);
 c_tilde_size = reduce(c_tilde_size);

 // Return the result
 return  List::create(_["r"]=r,
                      _["p"]=p,
                      _["K_tilde"]=K_tilde,
                      _["alpha"]=alpha,
                      _["c_tilde_assign"]=c_tilde_assign,
                      _["c_tilde"]=c_tilde_size,
                      _["accepted"]=accepted,
                      _["acc_ratio"]=acc_ratio
 );
}



// [[Rcpp::export]]
List PT_cpp (List & param, List & data, List & hyperparam, bool verbose) {
 // Initialization param
 int n_iter ;
 if( param.containsElementNamed("n_iter") ){
  n_iter = param["n_iter"];
 }else{
  stop("The option 'n_iter' is missing in the 'param' input.\n");
 }

 int burnin;
 if( param.containsElementNamed("burnin") ){
  burnin = param["burnin"];
 }else{
  burnin = 1 ;
 }
 if(burnin < 1) burnin = 1 ;

 int thin;
 if( param.containsElementNamed("thin") ){
  thin = param["thin"];
 }else{
  thin = 1 ;
 }
 if(thin < 1) thin = 1 ;

 arma::vec temperatures;
 if( param.containsElementNamed("temperatures") ){
  temperatures = as<arma::vec>(param["temperatures"]) ;
 }else{
  stop("The temperatures should be specified.") ;
 }
 int n_MH = temperatures.size() ;

 double lambda;
 if( param.containsElementNamed("lambda") ){
  lambda = param["lambda"];
 }else{
  lambda = n_MH ;
 }

 int nbre_print;
 if( param.containsElementNamed("nbre_print") ){
  nbre_print = param["nbre_print"];
 }else{
  nbre_print = 100 ;
 }
 if(nbre_print>=n_iter) nbre_print = n_iter-1;

 // Initialization data and param objects
 int N = data["N"] ;
 double pi = hyperparam["pi"] ;
 List res_proposal ;
 List theta_star ;
 List theta_tmp ;
 double u;

 // Initialization of objects related to the swap moves
 int n_swap;
 arma::vec js(2);
 double acc_ratio_tmp;

 // Initialization computational timer
 Timer timer_tot;
 Timer timer;
 timer_tot.step("all");

 // Compute a starting point
 List MH_res(n_MH) ;
 List theta(n_MH) ;
 if( param.containsElementNamed("starting_point") ){
  for(int j=0 ; j<n_MH ; ++j){
   theta[j] = param["starting_point"] ;
  }
 }else{
  for(int j=0 ; j<n_MH ; ++j){
   theta[j] = rprior_cpp(data,hyperparam) ;
  }
 }


 // Initialize the trace
 arma::mat acc_ratio(n_iter,n_MH) ;
 arma::mat accepted = arma::zeros<arma::mat>(n_iter,n_MH) ;

 arma::mat r(n_iter,n_MH) ;
 arma::mat p(n_iter,n_MH) ;
 arma::mat K_tilde(n_iter,n_MH) ;
 arma::mat alpha(n_iter,n_MH) ;

 arma::cube c_tilde_assign = arma::zeros<arma::cube>( n_iter , n_MH, floor(N/(1-pi)) ) ;
 arma::cube c_tilde_size = arma::zeros<arma::cube>( n_iter  , n_MH, floor(N/(1-pi)) ) ;
 arma::vec c_tilde_assign_tmp;
 arma::vec c_tilde_size_tmp;


 // Compute what to verbose
 arma::mat n_iters;

 if(verbose){
  int tmp = ceil(n_iter / (nbre_print+1)) ;
  int tmp2 = n_iter % tmp ;

  if(tmp2 == 0){
   n_iters =  arma::mat(nbre_print+1,2);
  }else{
   n_iters =  arma::mat(nbre_print+2,2);
  }


  for(int k=0 ; k<n_iters.n_rows  ; ++k){
   n_iters(k,0) = k*tmp ;
   n_iters(k,1) = (k+1)*tmp-1 ;
  }
  if(tmp2 != 0){
   n_iters(nbre_print+1,0) = (nbre_print+1)*tmp ;
   n_iters(nbre_print+1,1) = n_iter-1 ;
  }
 }else{
  n_iters = arma::mat(1,2);
  n_iters(0,0) = 0;
  n_iters(0,1) = n_iter -1;
 }

 arma::vec percents = sequence_le(0,100,nbre_print+1) ;
 std::string to_print = "" ;
 std::string to_add = "=" ;
 CharacterVector string_tmp ;
 arma::vec when_to_print = sequence_le(0,n_iters.n_rows-1,20) ;
 int count = 0;

 NumericVector timer_vec ;
 double start_time ;
 double end_time ;
 double remain_time ;
 double time_mean = 0;

 // Start the burnin loop
 if(verbose) Rcpp::Rcout << "Burnin " ;
 for(int i=0 ; i<burnin ; ++i){
  // The MH iteration
  // For each parallel chains
  for(int l=0 ; l<n_MH ; ++l){
   theta_tmp = theta[l] ;

   res_proposal = MH_proposal_cpp(theta_tmp, data, hyperparam, temperatures(l)) ;
   theta_star = res_proposal["theta_star"] ;
   acc_ratio_tmp = res_proposal["acc_ratio"] ;

   // Accept/reject
   u = R::runif(0,1) ;
   if(u < acc_ratio_tmp){
    theta_tmp = theta_star ;
    theta[l] = theta_tmp ;
   }
  }
 }
 if(verbose) Rcpp::Rcout << "done." << std::endl ;

 // Start the loop
 for(int j=0 ; j<n_iters.n_rows ; ++j){
  // Stuff to print
  if(verbose){
   if(j == when_to_print(count)){
    string_tmp  = wrap(to_print + to_add) ;
    to_print = string_tmp(0);
    count += 1;
   }
   if(j!=0){
    timer_vec = NumericVector(timer);
    start_time = timer_vec(2*j-2) ;
    end_time = timer_vec(2*j-1) ;

    time_mean = (time_mean * (j-1) + (end_time - start_time)) / j ;
    remain_time = time_mean * (n_iters.n_rows-j) / 1000000000 ;

    if(j<nbre_print+1){
     std::stringstream strs;
     if(remain_time  > 60){
      if(remain_time > 3600){
       strs << "[" << to_print << " "<<  percents(j) << "%]  " <<
        "(remaining: " << remain_time  / 3600<< " hours)" ;
      }
      strs << "[" << to_print << " "<<  percents(j) << "%]  " <<
       "(remaining: " << remain_time  / 60 << " mins)" ;
     }else{
      strs << "[" << to_print << " "<<  percents(j) << "%]  " <<
       "(remaining: " << remain_time  << " seconds)" ;
     }
     std::string temp_str = strs.str();
     char const* char_type = temp_str.c_str();

     REprintf("\r");
     REprintf("%s", char_type);
     REprintf("\r");

     R_FlushConsole();
     R_CheckUserInterrupt();
    }
   }
  }

  // The MH iteration
  timer.step("step");
  for(int i=n_iters(j,0) ; i<=n_iters(j,1) ; ++i){
   // Thinned iterations
   for(int t=0 ; t<thin-1 ; ++t){
    // For each parallel chains
    for(int l=0 ; l<n_MH ; ++l){
     theta_tmp = theta[l] ;

     res_proposal = MH_proposal_cpp(theta_tmp, data, hyperparam, temperatures(l)) ;
     theta_star = res_proposal["theta_star"] ;
     acc_ratio_tmp = res_proposal["acc_ratio"] ;

     // Accept/reject
     u = R::runif(0,1) ;
     if(u < acc_ratio_tmp){
      theta_tmp = theta_star ;
      theta[l] = theta_tmp ;
     }
    }

    // Swap moves
    n_swap = R::rpois(lambda) ;
    if(n_swap == 0) n_swap = 1 ;
    for(int i_swap=0 ; i_swap<n_swap ; ++i_swap){
     js = sample_no_replace_cpp(2,n_MH);

     theta_tmp = theta[js(1)] ;
     acc_ratio_tmp = log_dlkh_cpp(theta_tmp,data,
                                  temperatures(js(0)) - temperatures(js(1)) ) ;
     theta_tmp = theta[js(0)] ;
     acc_ratio_tmp -= log_dlkh_cpp(theta_tmp,data,
                                   temperatures(js(0)) - temperatures(js(1)) ) ;
     acc_ratio_tmp = exp(acc_ratio_tmp) ;

     // Accept/reject
     u = R::runif(0,1) ;
     if(u < acc_ratio_tmp){
      theta_tmp = theta[js(0)]  ;
      theta[js(0)]  = theta[js(1)]  ;
      theta[js(1)]  = theta_tmp  ;
     }
    }
   }

   // For each parallel chains
   for(int l=0 ; l<n_MH ; ++l){
    theta_tmp = theta[l] ;

    res_proposal = MH_proposal_cpp(theta_tmp, data, hyperparam, temperatures(l)) ;
    theta_star = res_proposal["theta_star"] ;
    acc_ratio(i,l) = res_proposal["acc_ratio"] ;

    // Accept/reject
    u = R::runif(0,1) ;
    if(u < acc_ratio(i,l)){
     theta_tmp = theta_star ;
     theta[l] = theta_tmp ;
     accepted(i,l) = 1 ;
    }
   }

   // Swap moves
   n_swap = R::rpois(lambda) ;
   if(n_swap == 0) n_swap = 1 ;
   for(int i_swap=0 ; i_swap<n_swap ; ++i_swap){
    js = sample_no_replace_cpp(2,n_MH);

    theta_tmp = theta[js(1)] ;
    acc_ratio_tmp = log_dlkh_cpp(theta_tmp,data,
                                 temperatures(js(0)) - temperatures(js(1)) ) ;
    theta_tmp = theta[js(0)] ;
    acc_ratio_tmp -= log_dlkh_cpp(theta_tmp,data,
                                  temperatures(js(0)) - temperatures(js(1)) ) ;
    acc_ratio_tmp = exp(acc_ratio_tmp) ;

    // Accept/reject
    u = R::runif(0,1) ;
    if(u < acc_ratio_tmp){
     theta_tmp = theta[js(0)]  ;
     theta[js(0)]  = theta[js(1)]  ;
     theta[js(1)]  = theta_tmp  ;
    }
   }

   // Update traces ;
   for(int l=0 ; l<n_MH ; ++l){
    theta_tmp = theta[l] ;

    r(i,l) = theta_tmp["r"];
    p(i,l) = theta_tmp["p"];
    K_tilde(i,l) = theta_tmp["K_tilde"];
    alpha(i,l) = theta_tmp["alpha"];;

    c_tilde_assign_tmp = as<arma::vec>(theta_tmp["c_tilde_assign"]) ;
    if(c_tilde_assign_tmp.size() > c_tilde_assign.n_slices){
     c_tilde_assign = expand_cube(c_tilde_assign,c_tilde_assign_tmp.size()) ;
    }
    c_tilde_assign.subcube(i,l,0,
                           i,l,c_tilde_assign_tmp.size()-1) =
                            trans(c_tilde_assign_tmp) ;


    c_tilde_size_tmp = as<arma::vec>(theta_tmp["c_tilde"]) ;
    if(c_tilde_size_tmp.size() > c_tilde_size.n_slices){
     c_tilde_size = expand_cube(c_tilde_size,c_tilde_size_tmp.size()) ;
    }
    c_tilde_size.subcube(i,l,0,
                         i,l,c_tilde_size_tmp.size()-1) =
                          trans(c_tilde_size_tmp) ;
   }
  }
  timer.step("step");
 }

 // Rearrange c_tilde traces
 c_tilde_assign = reduce_cube(c_tilde_assign);
 c_tilde_size = reduce_cube(c_tilde_size);

 // Computational time
 timer_tot.step("all");
 NumericVector timer_vec_tot ;
 timer_vec_tot = NumericVector(timer_tot);
 double comput_time = (timer_vec_tot(1) - timer_vec_tot(0))/ 1000000000;

 // Return the result
 return  List::create(_["trace"]=List::create(
  _["r"]=r,
  _["p"]=p,
  _["K_tilde"]=K_tilde,
  _["alpha"]=alpha,
  _["c_tilde_assign"]=c_tilde_assign,
  _["c_tilde"]=c_tilde_size,
  _["accepted"]=accepted,
  _["acc_ratio"]=acc_ratio),
  _["comput_time"]=comput_time
 );
}



// [[Rcpp::export]]
List PT_bs_cpp (List & param, List & datas, List & hyperparam, bool verbose) {
 // Initialization param
 int n_iter ;
 if( param.containsElementNamed("n_iter") ){
  n_iter = param["n_iter"];
 }else{
  stop("The option 'n_iter' is missing in the 'param' input.\n");
 }

 int burnin;
 if( param.containsElementNamed("burnin") ){
  burnin = param["burnin"];
 }else{
  burnin = 1 ;
 }
 if(burnin < 1) burnin = 1 ;

 int thin;
 if( param.containsElementNamed("thin") ){
  thin = param["thin"];
 }else{
  thin = 1 ;
 }
 if(thin < 1) thin = 1 ;

 arma::vec temperatures;
 if( param.containsElementNamed("temperatures") ){
  temperatures = as<arma::vec>(param["temperatures"]) ;
 }else{
  stop("The temperatures should be specified.") ;
 }
 int n_MH = temperatures.size() ;

 double lambda;
 if( param.containsElementNamed("lambda") ){
  lambda = param["lambda"];
 }else{
  lambda = n_MH ;
 }

 int nbre_print;
 if( param.containsElementNamed("nbre_print") ){
  nbre_print = param["nbre_print"];
 }else{
  nbre_print = 100 ;
 }
 if(nbre_print>=n_iter) nbre_print = n_iter-1;

 // Initialization data and param objects
 int bs = datas.size() ;
 List data = datas[0] ;
 int N = data["N"] ;
 double pi = hyperparam["pi"] ;
 List res_proposal ;
 List theta_star ;
 List theta_tmp ;
 double u;

 // Initialization of objects related to the swap moves
 int n_swap;
 arma::vec js(2);
 double acc_ratio_tmp;

 // Initialization computational timer
 Timer timer_tot;
 Timer timer;
 timer_tot.step("all");

 // Compute a starting point
 List MH_res(n_MH) ;
 List theta(n_MH) ;
 for(int j=0 ; j<n_MH ; ++j){
  if( param.containsElementNamed("starting_point") ){
   theta[j] = param["starting_point"] ;
  }else{
   theta[j] = rprior_cpp(data,hyperparam) ;
  }
 }


 // Initialize the trace
 arma::mat acc_ratio(bs*n_iter,n_MH) ;
 arma::mat accepted = arma::zeros<arma::mat>(bs*n_iter,n_MH) ;

 arma::mat r(bs*n_iter,n_MH) ;
 arma::mat p(bs*n_iter,n_MH) ;
 arma::mat K_tilde(bs*n_iter,n_MH) ;
 arma::mat alpha(bs*n_iter,n_MH) ;

 arma::cube c_tilde_assign = arma::zeros<arma::cube>( bs*n_iter , n_MH, floor(N/(1-pi)) ) ;
 arma::cube c_tilde_size = arma::zeros<arma::cube>( bs*n_iter  , n_MH, floor(N/(1-pi)) ) ;
 arma::vec c_tilde_assign_tmp;
 arma::vec c_tilde_size_tmp;


 // Compute what to verbose
 arma::mat n_iters;

 if(verbose){
  int tmp = ceil(n_iter / (nbre_print+1)) ;
  int tmp2 = n_iter % tmp ;

  if(tmp2 == 0){
   n_iters =  arma::mat(nbre_print+1,2);
  }else{
   n_iters =  arma::mat(nbre_print+2,2);
  }


  for(int k=0 ; k<n_iters.n_rows  ; ++k){
   n_iters(k,0) = k*tmp ;
   n_iters(k,1) = (k+1)*tmp-1 ;
  }
  if(tmp2 != 0){
   n_iters(nbre_print+1,0) = (nbre_print+1)*tmp ;
   n_iters(nbre_print+1,1) = n_iter-1 ;
  }
 }else{
  n_iters = arma::mat(1,2);
  n_iters(0,0) = 0;
  n_iters(0,1) = n_iter -1;
 }

 arma::vec percents = sequence_le(0,100,nbre_print+1) ;
 std::string to_print = "" ;
 std::string to_add = "=" ;
 CharacterVector string_tmp ;
 arma::vec when_to_print = sequence_le(0,n_iters.n_rows-1,20) ;
 int count = 0;

 NumericVector timer_vec ;
 double start_time ;
 double end_time ;
 double remain_time ;
 double time_mean = 0;

 for(int b=0 ; b<bs ; ++b){
  data = datas[b] ;

  to_print = "" ;
  to_add = "=" ;
  count = 0;
  // Start the burnin loop
  // if(verbose) Rcpp::Rcout << "Burnin." << std::endl ;
  for(int i=0 ; i<burnin ; ++i){
   // The MH iteration
   // For each parallel chains
   for(int l=0 ; l<n_MH ; ++l){
    theta_tmp = theta[l] ;

    res_proposal = MH_proposal_cpp(theta_tmp, data, hyperparam, temperatures(l)) ;
    theta_star = res_proposal["theta_star"] ;
    acc_ratio_tmp = res_proposal["acc_ratio"] ;

    // Accept/reject
    u = R::runif(0,1) ;
    if(u < acc_ratio_tmp){
     theta_tmp = theta_star ;
     theta[l] = theta_tmp ;
    }
   }
  }

  // Start the loop
  for(int j=0 ; j<n_iters.n_rows ; ++j){
   // Stuff to print
   if(verbose){
    if(j == when_to_print(count)){
     string_tmp  = wrap(to_print + to_add) ;
     to_print = string_tmp(0);
     count += 1;
    }
    if(j!=0){
     timer_vec = NumericVector(timer);
     start_time = timer_vec(2*j-2) ;
     end_time = timer_vec(2*j-1) ;

     time_mean = (time_mean * (j-1) + (end_time - start_time)) / j ;
     remain_time = (time_mean * (n_iters.n_rows-j) +
      (bs-b-1) * time_mean * n_iters.n_rows)  / 1000000000 ;

     if(j<nbre_print+1){
      std::stringstream strs;
      if(remain_time  > 60){
       if(remain_time > 3600){
        strs << "Dataset "<< b+1 << ": [" << to_print << " "<<  percents(j) << "%]  " <<
         "(remaining: " << remain_time  / 3600<< " hours)                            " ;
       }
       strs << "Dataset "<< b+1 << ": [" << to_print << " "<<  percents(j) << "%]  " <<
        "(remaining: " << remain_time  / 60 << " mins)                               " ;
      }else{
       strs << "Dataset "<< b+1 << ": [" << to_print << " "<<  percents(j) << "%]  " <<
        "(remaining: " << remain_time  << " seconds)                                 " ;
      }
      std::string temp_str = strs.str();
      char const* char_type = temp_str.c_str();

      REprintf("\r");
      REprintf("%s", char_type);
      REprintf("\r");

      R_FlushConsole();
      R_CheckUserInterrupt();
     }
    }
   }

   // The MH iteration
   timer.step("step");
   for(int i=n_iters(j,0) ; i<=n_iters(j,1) ; ++i){
    // Thinned iterations
    for(int t=0 ; t<thin-1 ; ++t){ // voir si avec un -1 au thin ca marche quand thin = 1
     // For each parallel chains
     for(int l=0 ; l<n_MH ; ++l){
      theta_tmp = theta[l] ;

      res_proposal = MH_proposal_cpp(theta_tmp, data, hyperparam, temperatures(l)) ;
      theta_star = res_proposal["theta_star"] ;
      acc_ratio_tmp = res_proposal["acc_ratio"] ;

      // Accept/reject
      u = R::runif(0,1) ;
      if(u < acc_ratio_tmp){
       theta_tmp = theta_star ;
       theta[l] = theta_tmp ;
      }
     }

     // Swap moves
     n_swap = R::rpois(lambda) ;
     if(n_swap == 0) n_swap = 1 ;
     for(int i_swap=0 ; i_swap<n_swap ; ++i_swap){
      js = sample_no_replace_cpp(2,n_MH);

      theta_tmp = theta[js(1)] ;
      acc_ratio_tmp = log_dlkh_cpp(theta_tmp,data,
                                   temperatures(js(0)) - temperatures(js(1)) ) ;
      theta_tmp = theta[js(0)] ;
      acc_ratio_tmp -= log_dlkh_cpp(theta_tmp,data,
                                    temperatures(js(0)) - temperatures(js(1)) ) ;
      acc_ratio_tmp = exp(acc_ratio_tmp) ;

      // Accept/reject
      u = R::runif(0,1) ;
      if(u < acc_ratio_tmp){
       theta_tmp = theta[js(0)]  ;
       theta[js(0)]  = theta[js(1)]  ;
       theta[js(1)]  = theta_tmp  ;
      }
     }
    }

    // For each parallel chains
    for(int l=0 ; l<n_MH ; ++l){
     theta_tmp = theta[l] ;

     res_proposal = MH_proposal_cpp(theta_tmp, data, hyperparam, temperatures(l)) ;
     theta_star = res_proposal["theta_star"] ;
     acc_ratio(i+b*n_iter,l) = res_proposal["acc_ratio"] ;

     // Accept/reject
     u = R::runif(0,1) ;
     if(u < acc_ratio(i+b*n_iter,l)){
      theta_tmp = theta_star ;
      theta[l] = theta_tmp ;
      accepted(i+b*n_iter,l) = 1 ;
     }
    }

    // Swap moves
    n_swap = R::rpois(lambda) ;
    if(n_swap == 0) n_swap = 1 ;
    for(int i_swap=0 ; i_swap<n_swap ; ++i_swap){
     js = sample_no_replace_cpp(2,n_MH);

     theta_tmp = theta[js(1)] ;
     acc_ratio_tmp = log_dlkh_cpp(theta_tmp,data,
                                  temperatures(js(0)) - temperatures(js(1)) ) ;
     theta_tmp = theta[js(0)] ;
     acc_ratio_tmp -= log_dlkh_cpp(theta_tmp,data,
                                   temperatures(js(0)) - temperatures(js(1)) ) ;
     acc_ratio_tmp = exp(acc_ratio_tmp) ;

     // Accept/reject
     u = R::runif(0,1) ;
     if(u < acc_ratio_tmp){
      theta_tmp = theta[js(0)]  ;
      theta[js(0)]  = theta[js(1)]  ;
      theta[js(1)]  = theta_tmp  ;
     }
    }

    // Update traces ;
    for(int l=0 ; l<n_MH ; ++l){
     theta_tmp = theta[l] ;

     r(i+b*n_iter,l) = theta_tmp["r"];
     p(i+b*n_iter,l) = theta_tmp["p"];
     K_tilde(i+b*n_iter,l) = theta_tmp["K_tilde"];
     alpha(i+b*n_iter,l) = theta_tmp["alpha"];;

     c_tilde_assign_tmp = as<arma::vec>(theta_tmp["c_tilde_assign"]) ;
     if(c_tilde_assign_tmp.size() > c_tilde_assign.n_slices){
      c_tilde_assign = expand_cube(c_tilde_assign,c_tilde_assign_tmp.size()) ;
     }
     c_tilde_assign.subcube(i+b*n_iter,l,0,
                            i+b*n_iter,l,c_tilde_assign_tmp.size()-1) =
                             trans(c_tilde_assign_tmp) ;


     c_tilde_size_tmp = as<arma::vec>(theta_tmp["c_tilde"]) ;
     if(c_tilde_size_tmp.size() > c_tilde_size.n_slices){
      c_tilde_size = expand_cube(c_tilde_size,c_tilde_size_tmp.size()) ;
     }
     c_tilde_size.subcube(i+b*n_iter,l,0,
                          i+b*n_iter,l,c_tilde_size_tmp.size()-1) =
                           trans(c_tilde_size_tmp) ;
    }
   }
   timer.step("step");
  }
 }

 // Rearrange c_tilde traces
 c_tilde_assign = reduce_cube(c_tilde_assign);
 c_tilde_size = reduce_cube(c_tilde_size);

 // Computational time
 timer_tot.step("all");
 NumericVector timer_vec_tot ;
 timer_vec_tot = NumericVector(timer_tot);
 double comput_time = (timer_vec_tot(1) - timer_vec_tot(0))/ 1000000000;

 // Return the result
 return  List::create(_["trace"]=List::create(
  _["r"]=r,
  _["p"]=p,
  _["K_tilde"]=K_tilde,
  _["alpha"]=alpha,
  _["c_tilde_assign"]=c_tilde_assign,
  _["c_tilde"]=c_tilde_size,
  _["accepted"]=accepted,
  _["acc_ratio"]=acc_ratio),
  _["comput_time"]=comput_time
 );
}
