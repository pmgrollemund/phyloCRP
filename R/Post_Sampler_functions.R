################################# ----
#' compute_summaries
################################# ----
#' @description Compute post sampling summaries to analyze the resulting Markov chain.
#' @details The summaries are (for each model parameters): the posterior mean,
#' the posterior variance, the difference between posterior mean and true value
#' @return a numerical vector.
#' @param posterior_sample a list, the trace of the sampler algorithm.
#' @param truth a list, the true values of the model parameters (see \code{\link{truth1}}).
#' @export
#' @examples
#' res_PT <- PT(param1,data1,hyperparam1)
#' compute_summaries(res_PT$trace,truth1)
compute_summaries <- function(posterior_sample,truth,chain=1){
 alpha_mean <- mean(posterior_sample$alpha[,chain])
 alpha_var  <- var(posterior_sample$alpha[,chain])
 p_mean <- mean(posterior_sample$p[,chain])
 p_var  <- var(posterior_sample$p[,chain])
 r_mean <- mean(posterior_sample$r[,chain])
 r_var  <- var(posterior_sample$r[,chain])
 K_mean <- mean(posterior_sample$K_tilde[,chain])
 K_var  <- var(posterior_sample$K_tilde[,chain])

 if(ncol(posterior_sample$c_tilde[,chain,]) < 10){
  add <- 10 - ncol(posterior_sample$c_tilde[,chain,])
  c_tilde_mean <- c(apply(posterior_sample$c_tilde[,chain,],2,mean),
                    rep(0,add))
  c_tilde_var <- c(apply(posterior_sample$c_tilde[,chain,],2,var),
                   rep(0,add))
 }else{
  c_tilde_mean <- apply(posterior_sample$c_tilde[,chain,1:10],2,mean)
  c_tilde_var <- apply(posterior_sample$c_tilde[,chain,1:10],2,var)
 }

 alpha_truth   <- truth$alpha
 p_truth       <- truth$p
 r_truth       <- truth$r
 K_truth       <- truth$K_tilde
 c_tilde_truth <- truth$c_tilde

 alpha_error <- alpha_truth - alpha_mean
 p_error     <- p_truth     - p_mean
 r_error     <- r_truth     - r_mean
 K_error     <- K_truth     - K_mean
 if(length(c_tilde_truth) == length(c_tilde_mean)){
  c_tilde_error <- sqrt( sum(c_tilde_truth - c_tilde_mean)^2 )
 }
 if(length(c_tilde_truth) < length(c_tilde_mean)){
  c_tilde_truth <- c(c_tilde_truth,rep(0,length(c_tilde_mean)-length(c_tilde_truth)))
  c_tilde_error <- sqrt( sum(c_tilde_truth - c_tilde_mean)^2 )
 }
 if(length(c_tilde_truth) > length(c_tilde_mean)){
  c_tilde_mean_tmp <- c(c_tilde_mean,rep(0,-length(c_tilde_mean)+length(c_tilde_truth)))
  c_tilde_error <- sqrt( sum(c_tilde_truth - c_tilde_mean_tmp)^2 )
 }

 acc_ratio_mean <- mean(posterior_sample$acc_ratio[,chain],na.rm=T)

 res <- c(alpha_mean,alpha_var,alpha_error,
          p_mean, p_var, p_error,
          r_mean, r_var, r_error,
          K_mean, K_var, K_error,
          c_tilde_mean, c_tilde_var, c_tilde_error,
          acc_ratio_mean)
 names(res) <- c("alpha_mean","alpha_var","alpha_error",
                 "p_mean", "p_var", "p_error",
                 "r_mean", "r_var", "r_error",
                 "K_mean", "K_var", "K_error",
                 rep("c_tilde_mean",length(c_tilde_mean)),
                 rep("c_tilde_var",length(c_tilde_var)),
                 rep("c_tilde_error",length(c_tilde_error)),
                 "acc_ratio_mean")
 return(res)
}

################################# ----
#' compute_CI
################################# ----
#' @description Compute the credible interval for each model parameters.
#' @return A list containing the credible intervals.
#' @param posterior_sample a list, the trace of the sampler algorithm.
#' @param prob a nonnegative value, the coverage probability of the credible interval.
#' @export
#' @examples
#' res_PT <- PT(param1,data1,hyperparam1)
#' compute_CI(res_PT$trace)
compute_CI <- function(posterior_sample,prob=0.95,chain=1){
 N <- length(posterior_sample$alpha[,chain])
 index_quantiles <- floor(c((1-prob)/2,1-((1-prob)/2))*N)

 ### alpha
 alpha_tmp <- sort(posterior_sample$alpha[,chain])
 CI_alpha <- alpha_tmp[index_quantiles]

 ### p
 p_tmp <- sort(posterior_sample$p[,chain])
 CI_p <- p_tmp[index_quantiles]

 ### r
 r_tmp <- sort(posterior_sample$r[,chain])
 CI_r <- r_tmp[index_quantiles]

 ### r
 K_tmp <- sort(posterior_sample$K_tilde[,chain])
 CI_K <- K_tmp[index_quantiles]

 ### c_tilde
 CI_c_tilde <- matrix(0,ncol=ncol(posterior_sample$c_tilde[,chain,]),nrow=2)
 for(i in 1:ncol(posterior_sample$c_tilde[,chain,])){
  c_tilde_tmp <- sort(posterior_sample$c_tilde[,chain,i])
  CI_c_tilde[,i] <- c_tilde_tmp[index_quantiles]
 }

 ### return
 return(list(
  alpha = CI_alpha,
  p     = CI_p,
  r     = CI_r,
  K     = CI_K,
  c_tilde = CI_c_tilde ))
}

################################# ----
#' CI_covering
################################# ----
#' @description Check if the true parameter values belong to the credible intervals.
#' @return A list of boolean indicating if the credible intervals covers the
#' true values.
#' @param posterior_sample a list, the trace of the sampler algorithm.
#' @param truth a list, the true values of the model parameters (see \code{\link{truth1}}).
#' @export
#' @examples
#' res_PT <- PT(param1,data1,hyperparam1)
#' CI_covering(res_PT$trace,truth1)
CI_covering <- function(posterior_sample,truth){
 CI <- compute_CI(posterior_sample)

 res_alpha <- truth$alpha %between% CI$alpha
 res_p <- truth$p %between% CI$p
 res_r <- truth$r %between% CI$r
 res_K <- truth$K %between% CI$K
 res_c_tilde <- NULL
 for(i in 1:min(length(truth$c_tilde),ncol(CI$c_tilde))){
  res_c_tilde <- c(res_c_tilde,
                   truth$c_tilde[i] %between% CI$c_tilde[,i])
 }

 return(list(
  alpha = res_alpha,
  p = res_p,
  r = res_r,
  K = res_K,
  c_tilde = res_c_tilde
 ))
}
