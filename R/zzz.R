.onLoad <- function(libname, pkgname){
 paths_source <- c(
  paste(.libPaths()[1],"/phyloCRP/crp.cpp",sep=""),
  system.file(package="phyloCRP","crp.cpp"),
  file.path("inst","crp.cpp")
 )

 index <- which(file.exists(paths_source))
 if(length(index)>0){
  Rcpp::sourceCpp(paths_source[index[1]])
 }
}
