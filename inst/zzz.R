.onLoad <- function(libname, pkgname){
 library(Rcpp)
 library(RcppArmadillo)

 path_source <- paste(.libPaths()[1],"/phyloCRP/crp.cpp",sep="")

  sourceCpp(path_source)
}
