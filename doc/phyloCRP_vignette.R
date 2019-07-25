## ----setup, include=FALSE-------------------------------------------------------------------------
options(width = 100)
knitr::opts_chunk$set(echo = TRUE, dev = c('png','pdf'), fig.width= 6, 
                      fig.height=4,
                      collapse = TRUE,
                      comment = "#>"
)

## ---- echo=FALSE----------------------------------------------------------------------------------
# Remove all the object from memory.
rm(list=ls())

# Define paths
path_Data <- "../inst/extdata/"

## ----libs ,warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------
# Packages
library(coda)
library(ggplot2)
library(phyloCRP)
library(Rcpp)
library(RcppArmadillo)

if(!exists("MH_cpp")){
 # Source functions
 sourceCpp(file.path("..","inst","crp.cpp"))
}

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
N_tilde <- 50
p_truth <- 0.7
alpha_truth <- alpha <- 10

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
# Simulate the true subtree sizes 
c_tilde <- rev(sort(rcrp_size_cpp(alpha,N_tilde)))
cat("c_tilde =",c_tilde," \n")

# Simulate the observed subtree sizes
c <- subcrp_size_cpp(c_tilde,p_truth)
c <- rev(sort(c[c!=0]))
cat("c =",c )

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
data <- list(
 c = c, 
 N = sum(c), 
 K = length(c)
)

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
truth <- list(
 alpha   = alpha_truth, 
 p       = p_truth, 
 K_tilde = length(c_tilde),
 c_tilde = c_tilde,
 # t = rep(1:length(c_tilde),times=c_tilde),
 r       = N_tilde - data$N 
)

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
### Empirical Bayes ----
alpha_hat <- alpha_estimator(data$c,alpha_lim=c(1,50))$alpha_hat
alpha_hat

### Prior information ----
p0 <- p_truth
r0 <- data$N*(1-p0)/p0

### Hyperparam ----
hyperparam <- list(
 lambda1  = alpha_hat^2,
 lambda2  = alpha_hat,
 eta1 = data$N+1,
 eta2 = data$N*(1-p0)/p0+1,
 pi    = 1-p0,         
 kappa = data$N
)

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
### Options ----
param <- list(n_iter = 1e4)

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
param_prior <- param
param_prior$target <- "prior"

prior_sample <- MH_cpp(param_prior,data,hyperparam,TRUE)

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
### Posterior sampler ----
param$temperatures <- c(1,0.8,0.6,0.4)
posterior_sample <- PT_cpp(param,data,hyperparam,TRUE)

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
### The object structure ----
str(prior_sample)

### Traces ----
par(mfrow=c(2,2))
plot(prior_sample$alpha,type="l",main="alpha",xlab="",ylab="")
plot(prior_sample$p,type="l",main="p",xlab="",ylab="")
plot(prior_sample$K_tilde,type="l",main="K_tilde",xlab="",ylab="")
matplot(prior_sample$c_tilde,type="l",main="c_tilde",xlab="",ylab="")

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
### The object structure ----
str(posterior_sample)

### Traces ----
par(mfrow=c(2,2))
plot(posterior_sample$trace$alpha[,1],type="l",main="alpha",xlab="",ylab="")
plot(posterior_sample$trace$p[,1],type="l",main="p",xlab="",ylab="")
plot(posterior_sample$trace$K_tilde[,1],type="l",main="K_tilde",xlab="",ylab="")
matplot(posterior_sample$trace$c_tilde[,1,],type="l",main="c_tilde",xlab="",ylab="")

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
threshold <- 0.01

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
mcmc_alpha <- as.mcmc(posterior_sample$trace$alpha[,1])

geweke_z_alpha    <- geweke.diag(mcmc_alpha)
geweke_z_alpha    <- 1-pnorm(abs(geweke_z_alpha$z),0,1)

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
if(geweke_z_alpha > threshold) 
 cat("The p-value of the Geweke test is :", geweke_z_alpha, ".",
     "Then the null hypothesis is not rejected: the means of the first part and the last part",
     "of the Markov chain are not significatively different, which could be an",
     "evidence in support of the convergence of the chain. ") else
 cat("The p-value of the Geweke test is :", geweke_z_alpha, ".",
     "Then the null hypothesis is rejected: the means of the first part and the last part",
     "of the Markov chain are significatively different, which could be an",
     "evidence in support of the divergence of the chain. ")

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
heildel_sta_alpha <- heidel.diag(mcmc_alpha)[3]

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
if(heildel_sta_alpha > threshold) 
 cat("The p-value of the Heidelberger and Welch's test is :", heildel_sta_alpha, ".",
     "Then the null hypothesis is not rejected: the sample values come from a", 
     "stationary distribution. ") else
 cat("The p-value of the Heidelberger and Welch's test is :", heildel_sta_alpha, ".",
     "Then the null hypothesis is rejected: the sample values do not come from a", 
     "stationary distribution. ")

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
summaries <- compute_summaries(posterior_sample$trace,truth)

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
cat("The difference between the true alpha and the estimate is: ",summaries["alpha_error"])

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
cat("The difference between the true p and the estimate is: ",summaries["p_error"])

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
cat("The difference between the true r and the estimate is: ",summaries["r_error"])

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
cat("The difference between the true K_tilde and the estimate is: ",summaries["K_error"])

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
CIs <- compute_CI(posterior_sample$trace)

summaries2 <- CI_covering(posterior_sample$trace,truth)
summaries2$c_tilde <- as.logical(prod(summaries2$c_tilde))

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
if(summaries2$alpha)
 cat("The true value of alpha (",truth$alpha,") is in the credible interval: ",
     "[",CIs$alpha[1],",",CIs$alpha[2],"]") else 
  cat("The true value of alpha (",truth$alpha,") is not in the credible interval: ",
      "[",CIs$alpha[1],",",CIs$alpha[2],"]")

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
if(summaries2$p)
 cat("The true value of p (",truth$p,") is in the credible interval: ",
     "[",CIs$p[1],",",CIs$p[2],"]") else 
  cat("The true value of p (",truth$p,") is not in the credible interval: ",
      "[",CIs$p[1],",",CIs$p[2],"]")

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
if(summaries2$r)
 cat("The true value of r (",truth$r,") is in the credible interval: ",
     "[",CIs$r[1],",",CIs$r[2],"]") else 
  cat("The true value of r (",truth$r,") is not in the credible interval: ",
      "[",CIs$r[1],",",CIs$r[2],"]")

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
if(summaries2$K)
 cat("The true value of K (",truth$K,") is in the credible interval: ",
     "[",CIs$K[1],",",CIs$K[2],"]") else 
  cat("The true value of K (",truth$K,") is not in the credible interval: ",
      "[",CIs$K[1],",",CIs$K[2],"]")

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
xlim <- c(1, 
          max(ncol(CIs$c_tilde),length(truth$c_tilde)))
ylim <- c(0,
          max(CIs$c_tilde[,1],truth$c_tilde[1]))

matplot(t(CIs$c_tilde),type="o",lty=2,pch=16,col="gray",xlim=xlim,ylim=ylim)
lines(truth$c_tilde,type="o")

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
if(summaries2$c_tilde)
 cat("In this case, we see that the true subtree sizes are in the credible interval.") else 
  cat("In this case, we see that the true subtree sizes are not in the credible interval,",
      " since at least one true subtree size is not in the marginal credible interval.")

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
### CI posterior predictive ----
burnin <- 10
q1 <- 0.025 ; q2 <- 1-q1
in_CI <- NULL
for(k in (burnin+2):length(posterior_sample$trace$alpha[,1])){
 c_tmp <- posterior_sample$trace$c_tilde[k,1,]
 c_tmp <- c_tmp[c_tmp != 0]
 
 p_tmp <- posterior_sample$trace$p[k,1]
 
 c_data_tmp <- data$c
 c_data_tmp <- c(c_data_tmp,rep(0,length(c_tmp)-length(c_data_tmp)))
 
 in_CI_tmp <- TRUE
 for(j in 1:length(c_tmp)){
  in_CI_tmp <- in_CI_tmp & (phyloCRP:::"%between%"(c_data_tmp[j], 
                                                   qbinom(p_tmp,c_tmp[j],prob=c(q1,q2))) )
 }
 in_CI <- c(in_CI,in_CI_tmp)
}

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
cat("For all the iterations, the observed subtree sizes are: \n",
    sum(in_CI), "times in the credible interval of the posterior predictive distribution, and \n ",
    sum(1-in_CI), "times outside the credible interval of the posterior predictive distribution.")

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
### CI posterior predictive ----
q1 <- 0.025 ; q2 <- 1-q1

indexes <- which(names(summaries) == "c_tilde_mean")
c_tmp <- round(summaries[indexes],0)

p_tmp <- summaries["p_mean"]

c_data_tmp <- data$c
if(length(c_tmp)-length(c_data_tmp) > 0)
 c_data_tmp <- c(c_data_tmp,rep(0,length(c_tmp)-length(c_data_tmp)))

in_CI_tmp <- TRUE
CI_values <- NULL
for(j in 1:length(c_tmp)){
 CI_values <- cbind(CI_values,
                    qbinom(p_tmp,c_tmp[j],prob=c(q1,q2)))
 in_CI_tmp <- in_CI_tmp & (phyloCRP:::"%between%"(c_data_tmp[j], 
                                                   qbinom(p_tmp,c_tmp[j],prob=c(q1,q2))))
}


## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
if(in_CI_tmp)
 cat("With the estimate of p and c_tilde, the observed subtree sizes are ",
    "in the credible interval of the posterior predictive distribution.") else
 cat("With the estimate of p and c_tilde, the observed subtree sizes are not ",
    "in the credible interval of the posterior predictive distribution.") 

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
ESS <- c(
 effectiveSize(posterior_sample$trace$alpha[,1]),
 effectiveSize(posterior_sample$trace$K_tilde[,1]),
 effectiveSize(posterior_sample$trace$p[,1]),
 effectiveSize(posterior_sample$trace$r[,1])
)
names(ESS) <- c("alpha","K_tilde","p","r")
ESS

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
densities <- log_dtarget_sample_cpp_PT(posterior_sample$trace,data,hyperparam,1)
truth_density <- log_dtarget(truth,data,hyperparam)

head(densities)
densities <- as.data.frame(densities)
colnames(densities) <- c("log_posterior","posterior",
                    "log_lkh","lkh",
                    "log_prior","prior")

plot(densities$log_posterior,type="l",main="log posterior density") ; abline(h=truth_density$log_posterior,col=2)

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
alpha_est <- mean(posterior_sample$trace$alpha[,1])
p_est <- mean(posterior_sample$trace$p[,1])
r_est <- median(posterior_sample$trace$r[,1])
K_est <- median(posterior_sample$trace$K_tilde[,1])

c_est <- NULL
for(j in 1:ncol(posterior_sample$trace$c_tilde[,1,])){
 c_est <- c(c_est,
            median(posterior_sample$trace$c_tilde[,1,j]))
}

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
cat("The estimate of alpha is:",alpha_est,"\n")
cat("The estimate of p is:",p_est,"\n")
cat("The estimate of r is:",r_est,"\n")
cat("The estimate of K_tilde is:",K_est,"\n")
cat("The estimate of c_tilde is:",c_est[c_est!=0],"\n")

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
CIs <- compute_CI(posterior_sample$trace)

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
posterior_sample$trace$N_tilde <- apply(posterior_sample$trace$c_tilde[,1,],1,sum)
PIs <- posterior_sample$trace$K[,1] / posterior_sample$trace$N_tilde

PI_estimate <- mean(PIs)
PI_CI <- quantile(sort(PIs),c(0.025,0.975))

plot(PIs,type="l",ylab="Importance rate",xlab="Iteration")

## ----warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------------
cat("The estimate of the proportion of introductions is ",PI_estimate,
    " (",PI_CI[1],"-",PI_CI[2],").")

