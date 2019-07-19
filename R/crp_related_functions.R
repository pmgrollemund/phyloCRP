################################# ----
#' rcrp
################################# ----
#' @description Sample from the Chinese Restaurant Process.
#' @return A vector of table sizes.
#' @param alpha a nonnegative value, the CRP parameter.
#' @param N a nonnegative, the number of customers.
#' @details Prefer to use \code{rcrp_size_cpp} (or \code{rcrp_assign_cpp}).
#' @export
#' @examples
#' rcrp(5,10)
rcrp <- function(alpha,N){
 crp <- as.vector(1)
 K <- 1

 for(k in 2:N){
  probs <- c( crp , alpha )/( alpha + k - 1)
  choice <- sample(1:(K+1),1,prob=probs)

  if(choice == K+1){
   crp <- c(crp, 1)
   K <- K+1
  }else{
   crp[choice] <- crp[choice] + 1
  }
 }
 return(crp)
}

################################# ----
#' sub_rcrp
################################# ----
#' @description Subsample from a CRP with a Binomial model.
#' @return A vector of table sizes.
#' @param crp a numerical vector, the output of \code{\link{rcrp}}.
#' @param p a numerical value belonging to (0,1), the sampling proportion.
#' @details Prefer to use \code{subcrp_size_cpp}.
#' @export
#' @examples
#' crp <- rcrp(5,15)
#' print(crp)
#' print(sub_rcrp(crp,0.6))
#' print(sub_rcrp(crp,0.6))
sub_rcrp  <- function(crp,p){
 K <- length(crp)
 # N_tilde <- sum(crp)

 for(k in 1:K){
  crp[k] <- rbinom(1,crp[k],p)
 }

 K   <- sum(crp != 0)
 crp <- rev(sort(crp))
 crp <- crp[1:K]

 return(crp)
}

################################# ----
#' augm_rcrp
################################# ----
#' @description Augment a CRP.
#' @return A vector of table sizes.
#' @param crp a numerical vector, the output of \code{\link{rcrp}}.
#' @param alpha a nonnegative value, the CRP parameter (could be different
#' from the value used to generate \code{crp}).
#' @param N_tilde a nonnegative integer, the number of additional customers.
#' @details Prefer to use \code{augm_crp_size_cpp} (or \code{augm_crp_assign_cpp}).
#' @export
#' @examples
#' crp <- rcrp(5,15)
#' print(crp)
#' print(augm_rcrp(crp,5,0.6))
#' print(augm_rcrp(crp,5,0.6))
augm_rcrp <- function(crp,alpha,N_tilde){
 K <- length(crp)
 N <- sum(crp)

 if(N_tilde > N){
  for(k in (N+1):(N_tilde)){
   probs <- c( crp , alpha )/( alpha + k - 1)
   # probs <- c( crp , alpha )/( alpha + sum(crp))
   choice <- sample(1:(K+1),1,prob=probs)

   if(choice == K+1){
    crp <- c(crp, 1)
    K <- K+1
   }else{
    crp[choice] <- crp[choice] + 1
   }
  }
 }

 return(crp)
}

################################# ----
#' alpha_estimator
################################# ----
#' @description Compute the Maximum Likelihood estimate of \code{alpha} from a
#' \code{crp}.
#' @return A list containing
#' \describe{
#' \item{alpha_hat}{the ML estimate,}
#' \item{elln_alpha}{a useful function to derive the estimate.}
#' }
#' @param crp a numerical vector, the output of \code{\link{rcrp}}.
#' @param alpha_lim a two-vector, the range of the potential values for the estimate.
#' @param nbre_alpha an interger, the number of the potential values for the estimate.
#' @export
#' @examples
#' crp <- rcrp(5,15)
#' print(crp)
#' alpha_ML <- alpha_estimator(crp,c(0.1,10))
#' alpha_ML$alpha_hat
#'
#' plot(seq(0.1,10,le=1e4),alpha_ML$elln_alpha,type="l",xlab="",ylab="")
#' abline(h=length(crp),col="gray")
#' abline(v=alpha_ML$alpha_hat,col="gray")
alpha_estimator <- function(crp,alpha_lim,nbre_alpha=1e4){
 N <- sum(crp)
 K <- length(crp)

 alphas      <- seq(alpha_lim[1],
                    alpha_lim[2],
                    le=nbre_alpha)
 elln_alphas <- alphas*(digamma(alphas+N)-digamma(alphas))

 dist_from_K <- abs(elln_alphas - K)
 alpha_hat   <- alphas[which.min(dist_from_K)]

 return(list(alpha_hat  = alpha_hat,
             elln_alpha = elln_alphas))
}
