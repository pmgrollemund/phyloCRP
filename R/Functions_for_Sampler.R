################################# ----
#' rprior
################################# ----
#' @description Sample from the hierarchical prior distribution#'
#' @return a list containing the sampled values: \code{alpha}, \code{c_tilde},
#' \code{K_tilde}, \code{p} and \code{r}.
#' @param data see the example object \code{data1}
#' @param hyperparam see the example object \code{hyperparam1}
#' @details Prefer to use \code{rprior_cpp}.
#' @export
#' @examples
#' rprior(data1,hyperparam1)
rprior <- function(data,hyperparam){
 alpha <- stats::rgamma(1,shape=hyperparam$lambda1,rate=hyperparam$lambda2)
 r     <- stats::rnbinom(1,size=hyperparam$kappa,prob=1-hyperparam$pi)
 p     <- stats::rbeta(1,hyperparam$eta1,hyperparam$eta2)
 N_tilde <- data$N + r

 c_tilde <- rev(sort(augm_rcrp(data$c,alpha,N_tilde)))
 K_tilde <- length(c_tilde)

 theta <- list(
  alpha = alpha,
  c_tilde = c_tilde,
  K_tilde = K_tilde,
  p = p,
  r=r
 )

 return(theta)
}

################################# ----
#' compute_starting_point
################################# ----
#' @description Compute a starting point for the posterior sampler.
#' @return a list containing the sampled values: \code{alpha}, \code{c_tilde},
#' \code{K_tilde}, \code{p} and \code{r}.
#' @param data see the example object \code{data1}
#' @param hyperparam see the example object \code{hyperparam1}
#' @details Prefer to use \code{compute_starting_point_cpp}.
#' @export
#' @examples
#' compute_starting_point(data1,hyperparam1)
compute_starting_point <- function(data,hyperparam){
 alpha <- stats::rgamma(1,shape=hyperparam$lambda1,rate=hyperparam$lambda2)
 p     <- stats::rbeta(1,hyperparam$eta1,hyperparam$eta2)
 r     <- floor(data$N * (1-p) / p)
 d     <- hyperparam$d

 N_tilde <- data$N + r
 c_tilde <- data$c

 K_add   <- floor(data$K * (1-p) / p)
 c_tilde <- c(c_tilde,rep(1,K_add))

 if(sum(c_tilde) < N_tilde){
  while( sum(c_tilde) != N_tilde ){
   probs <- c(c_tilde)
   choice <- sample(1:length(c_tilde),1,prob=probs)
   c_tilde[choice] <- c_tilde[choice] + 1
  }
 }
 if(sum(c_tilde) > N_tilde){
  while( sum(c_tilde) != N_tilde ){
   K_tilde <- length(c_tilde)
   c_tmp <- c(data$c,rep(0,K_tilde - data$K))

   probs <- (c_tilde - 1) * (c_tilde > c_tmp)
   choice <- sample(1:length(c_tilde),1,prob=probs)
   c_tilde[choice] <- c_tilde[choice] - 1
   c_tilde <- c_tilde[c_tilde > 0]
  }
 }

 c_tilde <- rev(sort(c_tilde))
 K <- length(c_tilde)
 t <- rep(1:K,times=c_tilde)

 theta <- list(
  alpha = alpha,
  K_tilde =K,
  c_tilde = c_tilde,
  p = p,
  r=r
 )

 return(theta)
}

################################# ----
#' MH_proposal
################################# ----
#' @description Compute a proposal state for the Metropolis-Hastings algorithm.
#' @return a list containing:
#' \describe{
#'   \item{theta_star}{a list containing the proposed values for each parameter}
#'   \item{acc_ratio}{the acceptance probability.}
#' }
#' @param theta a list, the current state of the Markov Chain.
#' @param data see the example object \code{data1}
#' @param hyperparam see the example object \code{hyperparam1}
#' @param temperature a nonnegative value, the temperature parameter required to
#' compute the targer distribution.
#' @details Prefer to use \code{MH_proposal_cpp}.
#' @export
#' @examples
#' theta <- compute_starting_point(data1,hyperparam1)
#' MH_proposal(theta,data1,hyperparam1,1)
#' MH_proposal(theta,data1,hyperparam1,0.2)
MH_proposal  <- function(theta,data,hyperparam,temperature){
 # Initialization
 theta_star <- rprior(data,hyperparam)

 log_target <- log_dlkh(theta_star,data,hyperparam, temperature=temperature) -
  log_dlkh(theta,data,hyperparam,temperature=temperature)

 ### Compute acceptance ratio and return it
 acc_ratio <- min(1,exp(log_target))

 ################ Return the result ----
 return(list(theta_star = theta_star,
             acc_ratio  = acc_ratio)
 )
}

################################# ----
#' log_dlkh
################################# ----
#' @description Compute the tempered log likelihood.
#' @return a numerical value.
#' @param theta a list of parameter values.
#' @param data see the example object \code{data1}
#' @param hyperparam see the example object \code{hyperparam1}
#' @param temperature a nonnegative value, the temperature parameter required to
#' compute the targer distribution.
#' @details Prefer to use \code{log_dlkh_cpp}.
#' @export
#' @examples
#' theta <- compute_starting_point(data1,hyperparam1)
#' log_dlkh(theta,data1,hyperparam1,temperature=1)
#' log_dlkh(theta,data1,hyperparam1,temperature=0.5)
log_dlkh <- function(theta,data,hyperparam,temperature=1){
 c_tmp <- c( data$c , rep(0,theta$K_tilde - data$K) )

 res <- temperature *(
  data$N * log(theta$p) + theta$r * log(1-theta$p) + sum(lchoose(theta$c_tilde,c_tmp))
 )

 return(res)
}

################################# ----
#' log_dtarget
################################# ----
#' @description For one \code{theta} value, compute the (log) likelihood
#' posterior density, prior density and the likelihood, according to the target
#'  distribution (with \code{temperature = 1}).
#' @return a data.frame containing the log posterior density, the posterior density,
#'  the log likelihood, the likelihood, the log prior density and the prior density.
#' @param theta a list containing the \code{theta} values.
#' @param data see the example object \code{data1}
#' @param hyperparam see the example object \code{hyperparam1}
#' @details Prefer to use \code{log_dtarget_cpp}.
#' @export
#' @examples
#' theta <- compute_starting_point(data1,hyperparam1)
#' log_dtarget(theta,data1,hyperparam1)
log_dtarget <- function(theta,data,hyperparam){
 c_tmp <- c( data$c , rep(0,theta$K_tilde - data$K) )

 res <- matrix(NA,nrow=1,ncol=6)
 colnames(res) <- c("log_posterior","posterior",
                    "log_lkh","lkh",
                    "log_prior","prior")

 res[1,3] <- data$N * log(theta$p) + theta$r * log(1-theta$p) + sum(lchoose(theta$c_tilde,c_tmp))
 res[1,5] <-
  theta$K_tilde * log(theta$alpha) + sum(lgamma(theta$c_tilde)) + # c_tilde prior
  lgamma(theta$alpha) - lgamma(theta$alpha + theta$r + data$N) +

  lchoose(hyperparam$kappa + theta$r - 1, theta$r) + # r prior
  hyperparam$kappa * log(1-hyperparam$pi) + theta$r * log(hyperparam$pi) +

  hyperparam$lambda1 * log(hyperparam$lambda2) - lgamma(hyperparam$lambda1) +  # alpha prior
  (hyperparam$lambda1-1) * log(theta$alpha) -
  hyperparam$lambda2 * theta$alpha -

  lbeta(hyperparam$eta1,hyperparam$eta2) + (hyperparam$eta1-1) * log(theta$p) + # p prior
  (hyperparam$eta2-1) * log(1-theta$p)

 res[1,1] <- res[1,3] + res[1,5]
 res[1,c(2,4,6)] <- exp(res[1,c(1,3,5)])

 return(as.data.frame(res))
}


################################# ----
#' log_dtarget_sample
################################# ----
#' @description For a sample of \code{theta}, compute the (log) likelihood
#' posterior density, prior density and the likelihood, according to the target
#'  distribution (with \code{temperature = 1}).
#' @return a data.frame containing the log posterior density, the posterior density,
#'  the log likelihood, the likelihood, the log prior density and the prior density.
#' @param trace a list of parameter values.
#' @param data see the example object \code{data1}
#' @param hyperparam see the example object \code{hyperparam1}
#' @details Prefer to use \code{log_dtarget_sample_cpp}.
#' @export
#' @examples
#' res_MH <- MH(param1,data1,hyperparam1,verbose=FALSE)
#' log_dtarget_sample(res_MH,data1,hyperparam1)
log_dtarget_sample <- function(trace,data,hyperparam){
 N <- length(trace$alpha)
 res <- matrix(NA,nrow=N,ncol=6)
 colnames(res) <- c("log_posterior","posterior",
                    "log_lkh","lkh",
                    "log_prior","prior")

 for(i in 1:N){
  c_tmp <- c( data$c , rep(0,trace$K_tilde[i] - data$K) )

  c_tilde <- trace$c_tilde[i,]
  c_tilde <- c_tilde[c_tilde!=0]

  res[i,3] <- data$N * log(trace$p[i]) + trace$r[i] * log(1-trace$p[i]) + sum(lchoose(c_tilde,c_tmp))
  res[i,5] <-
   trace$K_tilde[i] * log(trace$alpha[i]) + sum(lgamma(c_tilde)) + lgamma(trace$alpha[i]) - # c_tilde prior
   lgamma(trace$alpha[i] + trace$r[i] + data$N)  +

   lchoose(hyperparam$kappa + trace$r[i] - 1, trace$r[i]) + # r prior
   hyperparam$kappa * log(1-hyperparam$pi) + trace$r[i] * log(hyperparam$pi) +

   hyperparam$lambda1 * log(hyperparam$lambda2) - lgamma(hyperparam$lambda1) +  # alpha prior
   (hyperparam$lambda1-1) * log(trace$alpha[i]) -
   hyperparam$lambda2 * trace$alpha[i] -

   lbeta(hyperparam$eta1,hyperparam$eta2) + (hyperparam$eta1-1) * log(trace$p[i]) + # p prior
   (hyperparam$eta2-1) * log(1-trace$p[i])
 }

 res[,1] <- res[,3] + res[,5]
 res[,c(2,4,6)] <- exp(res[,c(1,3,5)])

 return(as.data.frame(res))
}
