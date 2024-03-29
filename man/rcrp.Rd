% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crp_related_functions.R
\name{rcrp}
\alias{rcrp}
\title{rcrp}
\usage{
rcrp(alpha, N)
}
\arguments{
\item{alpha}{a nonnegative value, the CRP parameter.}

\item{N}{a nonnegative, the number of customers.}
}
\value{
A vector of table sizes.
}
\description{
Sample from the Chinese Restaurant Process.
}
\details{
Prefer to use \code{rcrp_size_cpp} (or \code{rcrp_assign_cpp}).
}
\examples{
rcrp(5,10)
}
