% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_for_Sampler.R
\name{MH_proposal}
\alias{MH_proposal}
\title{MH_proposal}
\usage{
MH_proposal(theta, data, hyperparam, temperature)
}
\arguments{
\item{theta}{a list, the current state of the Markov Chain.}

\item{data}{see the example object \code{data1}}

\item{hyperparam}{see the example object \code{hyperparam1}}

\item{temperature}{a nonnegative value, the temperature parameter required to
compute the targer distribution.}
}
\value{
a list containing:
\describe{
  \item{theta_star}{a list containing the proposed values for each parameter}
  \item{acc_ratio}{the acceptance probability.}
}
}
\description{
Compute a proposal state for the Metropolis-Hastings algorithm.
}
\details{
Prefer to use \code{MH_proposal_cpp}.
}
\examples{
theta <- compute_starting_point(data1,hyperparam1)
MH_proposal(theta,data1,hyperparam1,1)
MH_proposal(theta,data1,hyperparam1,0.2)
}
