################################# ----
#' between
################################# ----
#' @description Check if a number belong to a given interval.
#' @return a logical value.
#' @param value a numerical value.
#' @param interval a numerical vector: (lower,upper).
#' @examples
#' 1 %between% c(0,2)
#' 2 %between% c(0,2)
#' 3 %between% c(0,2)
"%between%" <- function(x,I){
 !(x < I[1] || x > I[2])
}

################################# ----
#' reduce_mat
################################# ----
#' @description Internal function to reduce matrix dimensions when it is necessary.
#' @return a matrix.
#' @param mat a matrix.
reduce_mat <- function(mat){
 col_sum <- apply(mat,2,sum )
 index <- which(col_sum ==0)[1]
 res <- mat[,1:(index-1)]

 return(res)
}

################################# ----
#' reduce_cube
################################# ----
#' @description Internal function to reduce cube dimensions when it is necessary.
#' @return a cube (3D array).
#' @param cube a cube (3D array).
reduce_cube <- function(cube){
 slice_sum <- apply(cube,3,sum)
 index <- which(slice_sum ==0)[1]
 res <- cube[,,1:(index-1)]

 return(res)
}

################################# ----
#' expand_mat
################################# ----
#' @description Internal function to expand matrix dimensions when it is necessary.
#' @return a matrix.
#' @param mat a matrix.
expand_mat <- function(mat,n){
 n_cols <- ncol(mat)
 if(n_cols>n){
  return(mat)
 }else{
  res <- cbind(mat,matrix(0,nrow(mat),n-n_cols))
  return(res)
 }
}

################################# ----
#' expand_cube
################################# ----
#' @description Internal function to expand cube dimensions when it is necessary.
#' @return a cube (3D array).
#' @param cube a cube (3D array).
expand_cube <- function(cube,n){
 n_slices <- dim(cube)[3]
 if(n_slices>n){
  return(cube)
 }else{
  res <- array(0,c(dim(cube)[1:2],n))
  res[,,1:n_slices] <- cube

  return(res)
 }
}






