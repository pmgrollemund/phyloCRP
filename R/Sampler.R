################################# ----
#' Metropolis-Hastings
################################# ----
#' @description A Metropolis-Hastings algorithm to sample from the prior or the
#' posterior.
#' @details Should only be used for obtaining a prior sample: \code{param$target = "prior"}
#' @return a list containing the \code{theta} values and other useful quantities.
#' @param param a list containing options for the sampler algorithm:
#' \describe{
#'   \item{n_iter}{the size of the ouput posterior sample}
#'   \item{temperature}{(optional) the temperatures used to define the chain invariant distributions}
#'   \item{lambda}{(optional) related to the random number of swap moves per iteration}
#'   \item{thin}{(optional) keep only one iteration every \code{thin}-iterations.}
#'   \item{burnin}{(optional) an integer, drop the first iterations.}
#'   \item{target}{(optional) a string, "posterior" (default) or "prior".}
#' }
#' @param data see the example object \code{data1}
#' @param hyperparam see the example object \code{hyperparam1}
#' @param verbose Write stuff.
#' @details Prefer to use \code{MH_cpp}.
#' @export
#' @examples
#' param1$target <- "prior"
#' res_MH <- MH(param1,data1,hyperparam1)
MH <- function(param,data,hyperparam,verbose=T){
 ## Initialization
 n_iter <- param$n_iter
 burnin <- param$burnin
 target <- param$target
 if(is.null(burnin)) burnin <- floor(n_iter/10)
 if(is.null(target)) target <- "posterior"

 ## Options
 temperature <- param$temperature
 if(is.null(temperature)) temperature <- 1
 if(target=="prior")      temperature <- 0

 ## Compute a starting point
 theta <- param$starting_point
 if(is.null(theta)) theta <- compute_starting_point(data,hyperparam)

 ## Initialize 'trace'
 trace <- list()
 trace$acc_ratio <- rep(NA,n_iter)
 trace$accepted  <- rep(NA,n_iter)

 trace$r <- rep(NA,n_iter)
 trace$p <- rep(NA,n_iter)
 trace$K_tilde <- rep(NA,n_iter)
 trace$alpha <- rep(NA,n_iter)
 trace$c_tilde <- matrix(0,nrow=n_iter,ncol= floor(data$N/(1-hyperparam$pi)) )

 ## Start the loop
 i_verbose <- floor(seq(1,n_iter,le=11))
 if(verbose) cat("MH: \n")
 for(i in 1:n_iter){
  ## Display the progress
  if(verbose & i %in% i_verbose)
   cat("\t",(which(i_verbose == i)-1) * 10,"%\n",sep="")

  ## Compute the proposal and the acceptance ratio
  res_proposal <- MH_proposal(theta,data,hyperparam,temperature=temperature)
  theta_star   <- res_proposal$theta_star
  acc_ratio    <- res_proposal$acc_ratio

  ## Accept/Reject the proposal
  u <- runif(1,0,1)
  accepted <- u <= acc_ratio
  if(accepted){
   theta <- theta_star
  }

  ## Update the trace
  trace$acc_ratio[i] <- acc_ratio
  trace$accepted[i] <- accepted

  trace$r[i] <- theta$r
  trace$p[i] <- theta$p
  trace$K_tilde[i] <- theta$K_tilde
  trace$alpha[i] <- theta$alpha
  if( length(theta$c_tilde) > ncol(trace$c_tilde) ){
   trace$c_tilde <- expand_mat(trace$c_tilde,length(theta$c_tilde))
  }
  trace$c_tilde[i,1:length(theta$c_tilde)] <- theta$c_tilde
 }

 trace$c_tilde <- reduce_mat(trace$c_tilde)

 ## Return the trace
 return(trace)
}


################################# ----
#' Parallele Tempering
################################# ----
#' @description A Parallel Tempering algorithm to sample from the target distribution.
#' @return a list containing the computational time and the list of the \code{MH}
#' output for each parallel chains (see \code{\link{MH}}).
#' @param param a list containing options for the sampler algorithm:
#' \describe{
#'   \item{n_iter}{the size of the ouput posterior sample}
#'   \item{temperatures}{(optional) the temperatures used to define the chain invariant distributions}
#'   \item{lambda}{(optional) related to the random number of swap moves per iteration}
#'   \item{thin}{(optional) keep only one iteration every \code{thin}-iterations.}
#'   \item{burnin}{(optional) an integer, drop the first iterations.}
#'   \item{target}{(optional) a string, "posterior" (default) or "prior".}
#' }
#' @param data see the example object \code{data1}
#' @param hyperparam see the example object \code{hyperparam1}
#' @param verbose Write stuff.
#' @details Prefer to use \code{PT_cpp}.
#' @export
#' @examples
#' res_PT <- PT(param1,data1,hyperparam1)
PT <- function(param,data,hyperparam,verbose=T){
 # Initialization
 n_iter <- param$n_iter
 temperature <- param$temperatures
 N <- length(temperature)
 MH_res <- list() ; length(MH_res) <- N
 lambda <- param$lambda
 if(is.null(lambda)) lambda <- N
 thin <- param$thin
 if(is.null(thin)) thin <- 1
 burnin <- param$burnin
 if(is.null(burnin)) burnin <- 0

 t_start <- Sys.time()
 ## Initialize theta and trace
 theta <- list() ; length(theta) <- N ;
 for(j in 1:N){
  if(is.null(param$starting_point)){
   theta[[j]] <- compute_starting_point(data,hyperparam)
  }else{
   theta[[j]] <- param$starting_point
  }
 }

 acc_ratio <- rep(NA,N)
 accepted <- rep(NA,N)

 trace <- list()
 trace$acc_ratio <- matrix(NA,n_iter,N)
 trace$accepted  <- matrix(NA,n_iter,N)

 trace$r <- matrix(NA,n_iter,N)
 trace$p <- matrix(NA,n_iter,N)
 trace$K_tilde <- matrix(NA,n_iter,N)
 trace$alpha <- matrix(NA,n_iter,N)
 trace$c_tilde <- array(0, dim=c(n_iter,N,floor(data$N/(1-hyperparam$pi))))

 ## Start the loop
 count_thin <- thin
 i_verbose <- floor(seq(burnin+1,burnin+n_iter*thin,le=11))
 if(verbose) cat("PT: \n")
 if(burnin > 0) cat("\tBurnin.  \n")
 for(i in 1:(n_iter*thin+burnin)){
  ## Display the progress
  if(verbose & i %in% i_verbose)
   cat("\t",(which(i_verbose == i)-1) * 10,"%\n",sep="")

  for(j in 1:N){
   ## Compute the proposal and the acceptance ratio
   res_proposal <- MH_proposal(theta[[j]],data,hyperparam,temperature=temperature[j])
   theta_star   <- res_proposal$theta_star
   acc_ratio[j]  <- res_proposal$acc_ratio

   ## Accept/Reject the proposal
   u <- runif(1,0,1)
   accepted[j] <- u <= acc_ratio[j]
   if(accepted[j]){
    theta[[j]] <- theta_star
   }
  }

  K <- rpois(1,lambda)
  if(K == 0) K <- 1
  for(k in 1:K){
   js <- sample(N,2)
   acc_ratio_tmp <- log_dlkh(theta[[js[2]]],data,hyperparam,temperature[js[1]] - temperature[js[2]]) -
    log_dlkh(theta[[js[1]]],data,hyperparam,temperature[js[1]] - temperature[js[2]])

   acc_ratio_tmp <- exp(acc_ratio_tmp)

   u <- runif(1,0,1)
   accepted_tmp <- u <= acc_ratio_tmp
   if(accepted_tmp){
    theta[[js[1]]] -> theta_tmp
    theta[[js[1]]] <- theta[[js[2]]]
    theta[[js[2]]] <- theta_tmp
   }
  }

  if(i > burnin){
   count_thin <- count_thin - 1
   if(count_thin == 0){
    for(j in 1:N){
     ## Update the trace
     trace$acc_ratio[i/thin-burnin,j] <- acc_ratio[j]
     trace$accepted[i/thin-burnin,j] <- accepted[j]

     trace$r[i/thin-burnin,j] <- theta[[j]]$r
     trace$p[i/thin-burnin,j] <- theta[[j]]$p
     trace$K_tilde[i/thin-burnin,j] <- theta[[j]]$K_tilde
     trace$alpha[i/thin-burnin,j] <- theta[[j]]$alpha
     if( length(theta$c_tilde) > dim(trace$c_tilde)[3] ){
      trace$c_tilde <- expand_cube(trace$c_tilde,length(theta$c_tilde))
     }
     trace$c_tilde[i/thin-burnin ,j,1:length(theta[[j]]$c_tilde)] <- theta[[j]]$c_tilde
    }
    count_thin <- thin
   }

  }
 }
 trace$c_tilde <- reduce_cube(trace$c_tilde)
 t_end <- Sys.time()

 comput_time <- t_end-t_start
 ## Return the trace
 return(list(trace=trace,
             comput_time=comput_time))
}


################################# ----
#' Parallele Tempering (with bootstrap)
################################# ----
#' @description A Parallel Tempering algorithm to sample from the target distribution,
#' with bootstrap input data.
#' @return a list  containing the result (for each
#'  bootstrap chain) of the \code{MH} output for each parallel chains
#'  (see \code{\link{MH}}).
#' @param param a list containing options for the sampler algorithm:
#' \describe{
#'   \item{n_iter}{the size of the ouput posterior sample}
#'   \item{temperatures}{(optional) the temperatures used to define the chain invariant distributions}
#'   \item{lambda}{(optional) related to the random number of swap moves per iteration}
#'   \item{thin}{(optional) keep only one iteration every \code{thin}-iterations.}
#'   \item{burnin}{(optional) an integer, drop the first iterations.}
#'   \item{target}{(optional) a string, "posterior" (default) or "prior".}
#' }
#' @param data a list of data, each item should be as the example object \code{data1}
#' @param hyperparam see the example object \code{hyperparam1}
#' @param verbose Write stuff.
#' @details Prefer to use \code{PT_bs_cpp}.
#' @export
#' @examples
#' datas <- list(data1,data1,data1)
#' res_PT <- PT_bs(param1,datas,hyperparam1)
PT_bs <- function(param,datas,hyperparam,verbose=T){
 # Initialization
 n_iter <- param$n_iter
 temperature <- param$temperatures
 N <- length(temperature)
 MH_res <- list() ; length(MH_res) <- N
 lambda <- param$lambda
 if(is.null(lambda)) lambda <- N
 burnin <- param$burnin
 thin <- param$thin
 if(is.null(thin)) thin <- 1
 if(is.null(burnin)) burnin <- max(floor(n_iter/10),1)

 t_start <- Sys.time()
 ## Initialize theta and trace
 theta <- list() ; length(theta) <- N ;
 for(j in 1:N){
  if(is.null(param$starting_point)){
   theta[[j]] <- compute_starting_point(datas[[1]],hyperparam)
  }else{
   theta[[j]] <- param$starting_point
  }
 }

 acc_ratio <- rep(NA,N)
 accepted <- rep(NA,N)

 bs <- length(datas)
 trace <- list()
 trace$acc_ratio <- matrix(NA,bs*n_iter,N)
 trace$accepted  <- matrix(NA,bs*n_iter,N)

 trace$r <- matrix(NA,bs*n_iter,N)
 trace$p <- matrix(NA,bs*n_iter,N)
 trace$K_tilde <- matrix(NA,bs*n_iter,N)
 trace$alpha <- matrix(NA,bs*n_iter,N)
 trace$c_tilde <- array(0, dim=c(bs*n_iter,N,floor(datas[[1]]$N/(1-hyperparam$pi))))

 for(b in 1:bs){
  data <- datas[[b]]

  ## Start the loop
  count_thin <- thin
  i_verbose <- floor(seq(1,n_iter,le=11))
  if(verbose) cat("PT: \n")
  for(i in 1:(burnin+n_iter)){
   ## Display the progress
   if(verbose & i %in% i_verbose)
    cat("\t",(which(i_verbose == i)-1) * 10,"%\n",sep="")

   for(j in 1:N){
    ## Compute the proposal and the acceptance ratio
    res_proposal <- MH_proposal(theta[[j]],data,hyperparam,temperature=temperature[j])
    theta_star   <- res_proposal$theta_star
    acc_ratio[j] <- res_proposal$acc_ratio

    ## Accept/Reject the proposal
    u <- runif(1,0,1)
    accepted <- u <= acc_ratio[j]
    if(accepted){
     theta[[j]] <- theta_star
    }
   }

   K <- rpois(1,lambda)
   if(K == 0) K <- 1
   for(k in 1:K){
    js <- sample(N,2)
    acc_ratio_tmp <- log_dlkh(theta[[js[2]]],data,hyperparam,temperature[js[1]] - temperature[js[2]]) -
     log_dlkh(theta[[js[1]]],data,hyperparam,temperature[js[1]] - temperature[js[2]])

    acc_ratio_tmp <- exp(acc_ratio_tmp)

    u <- runif(1,0,1)
    accepted_tmp <- u <= acc_ratio_tmp
    if(accepted_tmp){
     theta[[js[1]]] -> theta_tmp
     theta[[js[1]]] <- theta[[js[2]]]
     theta[[js[2]]] <- theta_tmp
    }
   }

   if(i > burnin){
    count_thin <- count_thin - 1
    if(count_thin == 0){
     for(j in 1:N){
      ## Update the trace
      trace$acc_ratio[i/thin-burnin + (b-1)*n_iter,j] <- acc_ratio[j]
      trace$accepted[i/thin-burnin + (b-1)*n_iter,j] <- accepted[j]

      trace$r[i/thin-burnin + (b-1)*n_iter,j] <- theta[[j]]$r
      trace$p[i/thin-burnin + (b-1)*n_iter,j] <- theta[[j]]$p
      trace$K_tilde[i/thin-burnin + (b-1)*n_iter,j] <- theta[[j]]$K_tilde
      trace$alpha[i/thin-burnin + (b-1)*n_iter,j] <- theta[[j]]$alpha
      if( length(theta$c_tilde) > dim(trace$c_tilde)[3] ){
       trace$c_tilde <- expand_cube(trace$c_tilde,length(theta$c_tilde))
      }
      trace$c_tilde[i/thin-burnin + (b-1)*n_iter,j,
                    1:length(theta[[j]]$c_tilde)] <- theta[[j]]$c_tilde
     }
     count_thin <- thin
    }
   }
  }
 }
 trace$c_tilde <- reduce_cube(trace$c_tilde)
 t_end <- Sys.time()

 comput_time <- t_end-t_start
 ## Return the trace
 return(list(trace=trace,
             comput_time=comput_time))
}
