#' a list of objects
#'
#' A data object for CRP model
#' @format a list of data
#' \describe{
#'   \item{c}{the table sizes}
#'   \item{N}{the number of customers}
#'   \item{K}{the number of tables}
#' }
"data1"

#'
#'
#' A list of the true parameter values
#' @format a list of the true parameter values
#' \describe{
#'   \item{alpha}{the value used to generate the true complete table sizes}
#'   \item{p}{the true sampling proportion}
#'   \item{K_tilde}{the true number of tables}
#'   \item{c_tilde}{the true complete table sizes}
#'   \item{r}{the true number of missing customers}
#' }
"truth1"

#'
#'
#' A hyperparam object for the prior specification
#' @format a hyperparam object (list)
#' \describe{
#'   \item{lambda1}{shape in the Gamma prior of \code{alpha}}
#'   \item{lambda2}{rate in the Gamma prior of \code{alpha}}
#'   \item{eta1}{shape in the Beta prior of \code{p}}
#'   \item{eta2}{shape in the Beta prior of \code{p}}
#'   \item{pi}{probability of success in the Negative Binomial distribution of \code{r}}
#'   \item{kappa}{number of failures in the Negative Binomial distribution of \code{r}}
#' }
"hyperparam1"

#'
#'
#' A param object for the posterior sampler
#' @format a parma object (list)
#' \describe{
#'   \item{n_iter}{the size of the ouput posterior sample}
#'   \item{temperature}{the temperatures used to define the chain invariant distributions}
#'   \item{lambda}{related to the random number of swap moves per iteration}
#'   \item{thin}{keep only one iteration every \code{thin}-iterations.}
#' }
"param1"

