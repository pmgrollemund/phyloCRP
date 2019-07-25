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
path_RDS   <- paste(path_Data,"FakeData/",sep="")
path_tree  <- paste(path_Data,"Tree/",sep="")
path_plot  <- paste(path_Data,"Colored_tree/",sep="")
path_clade <- paste(path_Data,"Clades/",sep="")

## ----libs ,warning=FALSE,message=FALSE,echo=FALSE-------------------------------------------------
# Packages
library(big.phylo)
library(coda)
library(data.table)
library(ggplot2)
library(memisc)
library(stringr)

library(phyloCRP)
library(Rcpp)
library(RcppArmadillo)

if(!exists("MH_cpp")){
 # Source functions
 sourceCpp(file.path("..","inst","crp.cpp"))
}

# Options 
name_HXB2_file      <- "./LANL/HXB2.fasta"
name_HXB2           <- c("B.FR.1983.HXB2-LAI-IIIB-BRU.K03455.M",
                         "B.FR.1983.HXB2_copy_LAI.AF033819.M")
bs.n <- 2
subtypes <- c("01_AE","A1")
min_occur <- 10 # the threshold to determine the most represented subtypes 

biannual <-
 c("1980-1982","1983-1984","1985-1986","1987-1988","1989-1990","1991-1992",
   "1993-1994","1995-1996","1997-1998","1999-2000","2001-2002","2003-2004",
   "2005-2006","2007-2008","2009-2010","2011-2012","2013-2014","2015-2016",
   "2017-2018")
quinquennium <- 
 c("1980-1984","1984-1988","1989-1993","1994-1998","1999-2003","2004-2008",
   "2009-2013","2014-2018")

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
person <- readRDS(paste(path_RDS,"person_modified.rds",sep=""))
sequences_meta <- readRDS(paste(path_RDS,"sequences_meta_modified.rds",sep=""))

diagnosed_n_sequenced <- person$newnum %in% sequences_meta$newnum
p0 <- as.numeric( prop.table(table(diagnosed_n_sequenced))["TRUE"] )

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
trait <- person$main_transm2 == "IDU"
 
diagnosed_n_sequenced <- person$newnum[trait] %in% sequences_meta$newnum
p0_IDU <- as.numeric( prop.table(table(diagnosed_n_sequenced))["TRUE"] )

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
load(paste(path_clade,"Clades_save.RData",sep=""))
c <- rev(sort( unlist( lapply(clades,function(clade) length(clade$tip.label)) ) ))
data <- list(
 c = as.numeric(c), 
 N = sum(c), 
 K = length(c)
)

plot(data$c,type="o",pch=15,cex=0.7,xlab="k",ylab="c_k")

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
datas <- list() ; length(datas) <- bs.n+1

datas[[1]] <- data
for(b in 1:bs.n){
load(paste(path_clade,"Clades_save_",sprintf("%03d",b-1),".RData",sep=""))
c <- rev(sort( unlist( lapply(clades,function(clade) length(clade$tip.label)) ) ))
data_tmp <- list(
 c = as.numeric(c), 
 N = sum(c), 
 K = length(c)
)
data_tmp$c <- data_tmp$c[data_tmp$c != 0]
datas[[b+1]] <- data_tmp
}

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
r0 <- data$N*(1-p0)/p0
K0 <- data$K / p0
alpha_estimate <- alpha_estimator(data$c,c(0.1,data$N),nbre_alpha=1e4)
alpha0 <- alpha_estimate$alpha_hat

hyperparam <- list(
 lambda1  = alpha0^2,
 lambda2  = alpha0,
 eta1 = data$N+1,
 eta2 = data$N*(1-p0)/p0+1,
 pi    = 1-p0,         
 kappa = data$N
)

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
param <- list(n_iter = 1e3, 
              burnin=1e2,
              lambda=5,
              thin=5,
              temperatures = rev(seq(0,1,le=5+1)[-1]))
param_prior <- param
param_prior$target <- "prior"

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
prior_sample <- MH(param_prior,data,hyperparam)
prior_sample_cpp <- MH_cpp(param_prior,data,hyperparam,TRUE)

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
posterior_sample <- PT(param,data,hyperparam)
posterior_sample_cpp <- PT_cpp(param,data,hyperparam,TRUE)

## ---- echo=FALSE, eval=TRUE-----------------------------------------------------------------------
time_R <- as.numeric(posterior_sample$comput_time)
if(attr(posterior_sample$comput_time,which="units") == "mins"){
 time_R <- as.numeric(posterior_sample$comput_time) * 60
}
if(attr(posterior_sample$comput_time,which="units") == "hours"){
 time_R <- as.numeric(posterior_sample$comput_time) * 3600
}
cat("The Rcpp code is ",as.numeric(time_R  / posterior_sample_cpp$comput_time), 
    "times faster. \n")

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
posterior_sample_bs <- PT_bs(param,datas,hyperparam)
posterior_sample_bs_cpp <- PT_bs_cpp(param,datas,hyperparam,TRUE)

## ---- echo=FALSE, eval=TRUE-----------------------------------------------------------------------
time_R <- as.numeric(posterior_sample_bs$comput_time)
if(attr(posterior_sample_bs$comput_time,which="units") == "mins"){
 time_R <- as.numeric(posterior_sample_bs$comput_time) * 60
}
if(attr(posterior_sample_bs$comput_time,which="units") == "hours"){
 time_R <- as.numeric(posterior_sample_bs$comput_time) * 3600
}
cat("The Rcpp code is ",as.numeric(time_R  / posterior_sample_bs_cpp$comput_time), 
    "times faster. \n")

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
plot(posterior_sample_cpp$trace$alpha[,1],type="l")

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
matplot(posterior_sample_cpp$trace$alpha,type="l",lty=1)
legend(0, max(posterior_sample_cpp$trace$alpha), legend=paste("Temp:", param$temperatures),
       col=1:length(param$temperatures), lty=1, cex=0.8)


plot(density(posterior_sample_cpp$trace$alpha[,1]),main="",xlab="")
for(i in 2:length(param$temperatures))
 lines(density(posterior_sample_cpp$trace$alpha[,i]),col=i)
legend(min(posterior_sample_cpp$trace$alpha[,1]), 
       max(density(posterior_sample_cpp$trace$alpha[,1])$y), legend=paste("Temp:", param$temperatures),
       col=1:length(param$temperatures), lty=1, cex=0.8)

## ----warning=FALSE,message=FALSE,echo=TRUE, eval=TRUE---------------------------------------------
threshold <- 0.01

## ----warning=FALSE,message=FALSE,echo=TRUE, eval=TRUE---------------------------------------------
mcmc_alpha <- as.mcmc(posterior_sample_cpp$trace$alpha[,1])

geweke_z_alpha    <- geweke.diag(mcmc_alpha)
geweke_z_alpha    <- 1-pnorm(abs(geweke_z_alpha$z),0,1)

## ----warning=FALSE,message=FALSE,echo=FALSE, eval=TRUE--------------------------------------------
if(geweke_z_alpha > threshold) 
 cat("The p-value of the Geweke test is :", geweke_z_alpha, ".",
     "Then the null hypothesis is not rejected: the means of the first part and the last part",
     "of the Markov chain are not significatively different, which could be an",
     "evidence in support of the convergence of the chain. ") else
 cat("The p-value of the Geweke test is :", geweke_z_alpha, ".",
     "Then the null hypothesis is rejected: the means of the first part and the last part",
     "of the Markov chain are significatively different, which could be an",
     "evidence in support of the divergence of the chain. ")

## ----warning=FALSE,message=FALSE,echo=TRUE, eval=TRUE---------------------------------------------
heildel_sta_alpha <- heidel.diag(mcmc_alpha)[3]

## ----warning=FALSE,message=FALSE,echo=FALSE, eval=TRUE--------------------------------------------
if(heildel_sta_alpha > threshold) 
 cat("The p-value of the Heidelberger and Welch's test is :", heildel_sta_alpha, ".",
     "Then the null hypothesis is not rejected: the sample values come from a", 
     "stationary distribution. ") else
 cat("The p-value of the Heidelberger and Welch's test is :", heildel_sta_alpha, ".",
     "Then the null hypothesis is rejected: the sample values do not come from a", 
     "stationary distribution. ")

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
plot(posterior_sample_bs_cpp$trace$alpha[,1],type="l")

## ---- echo=TRUE, eval=TRUE------------------------------------------------------------------------
matplot(posterior_sample_bs_cpp$trace$alpha,type="l",lty=1)
legend(0, max(posterior_sample_bs_cpp$trace$alpha), legend=paste("Temp:", param$temperatures),
       col=1:length(param$temperatures), lty=1, cex=0.8)


plot(density(posterior_sample_bs_cpp$trace$alpha[,1]),main="",xlab="")
for(i in 2:length(param$temperatures))
 lines(density(posterior_sample_bs_cpp$trace$alpha[,i]),col=i)
legend(mean(posterior_sample_bs_cpp$trace$alpha[,1]), 
       max(density(posterior_sample_bs_cpp$trace$alpha[,1])$y), legend=paste("Temp:", param$temperatures),
       col=1:length(param$temperatures), lty=1, cex=0.8)

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
threshold <- 0.01

## ----warning=FALSE,message=FALSE,echo=TRUE--------------------------------------------------------
mcmc_alpha <- as.mcmc(posterior_sample_bs_cpp$trace$alpha[,1])

geweke_z_alpha    <- geweke.diag(mcmc_alpha)
geweke_z_alpha    <- 1-pnorm(abs(geweke_z_alpha$z),0,1)

## ----warning=FALSE,message=FALSE,echo=FALSE, eval=TRUE--------------------------------------------
if(geweke_z_alpha > threshold) 
 cat("The p-value of the Geweke test is :", geweke_z_alpha, ".",
     "Then the null hypothesis is not rejected: the means of the first part and the last part",
     "of the Markov chain are not significatively different, which could be an",
     "evidence in support of the convergence of the chain. ") else
 cat("The p-value of the Geweke test is :", geweke_z_alpha, ".",
     "Then the null hypothesis is rejected: the means of the first part and the last part",
     "of the Markov chain are significatively different, which could be an",
     "evidence in support of the divergence of the chain. ")

## ----warning=FALSE,message=FALSE,echo=TRUE, eval=TRUE---------------------------------------------
heildel_sta_alpha <- heidel.diag(mcmc_alpha)[3]

## ----warning=FALSE,message=FALSE,echo=FALSE, eval=TRUE--------------------------------------------
if(heildel_sta_alpha > threshold) 
 cat("The p-value of the Heidelberger and Welch's test is :", heildel_sta_alpha, ".",
     "Then the null hypothesis is not rejected: the sample values come from a", 
     "stationary distribution. ") else
 cat("The p-value of the Heidelberger and Welch's test is :", heildel_sta_alpha, ".",
     "Then the null hypothesis is rejected: the sample values do not come from a", 
     "stationary distribution. ")

## ----warning=FALSE,message=FALSE,echo=FALSE, eval=TRUE--------------------------------------------
posterior_sample_cpp$trace$N_tilde <- apply(posterior_sample_cpp$trace$c_tilde[,1,],1,sum)

PIs <- posterior_sample_cpp$trace$K_tilde[,1] / posterior_sample_cpp$trace$N_tilde
plot(PIs,type="l",ylab="Proportion of introductions",xlab="Iteration")

## ----warning=FALSE,message=FALSE,echo=TRUE, eval=TRUE---------------------------------------------
PI_estimate <- mean(PIs[-(1:param$burnin)])
CI_PI <- quantile(PIs,c(0.025,0.975))

cat("Estimate: ", round(PI_estimate,3), " (",round(CI_PI[1],3),"-",
    round(CI_PI[2],3),")",sep="")

## ----warning=FALSE,message=FALSE,echo=FALSE, eval=TRUE--------------------------------------------
posterior_sample_bs_cpp$trace$N_tilde <- apply(posterior_sample_bs_cpp$trace$c_tilde[,1,],1,sum)

PIs <- posterior_sample_bs_cpp$trace$K_tilde[,1] / posterior_sample_bs_cpp$trace$N_tilde
plot(PIs,type="l",ylab="Proportion of introductions",xlab="Iteration")

## ----warning=FALSE,message=FALSE,echo=TRUE, eval=TRUE---------------------------------------------
PI_estimate <- mean(PIs[-(1:param$burnin)])
CI_PI <- quantile(PIs,c(0.025,0.975))

cat("Estimate: ", round(PI_estimate,3), " (",round(CI_PI[1],3),"-",
    round(CI_PI[2],3),")",sep="")

