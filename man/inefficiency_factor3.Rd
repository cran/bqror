% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODRII.R
\name{inefficiency_factor3}
\alias{inefficiency_factor3}
\title{Inefficiency Factor for Ordinal Models
with 3 outcomes}
\usage{
inefficiency_factor3(beta_draws, nlags = 2, sigma_draws)
}
\arguments{
\item{beta_draws}{Gibbs draw of coeffcients of dimension \eqn{(k x iter)}.}

\item{nlags}{scalar variable with default = 2.}

\item{sigma_draws}{Gibbs draw of scale factor.}
}
\value{
Returns a list with components
\itemize{
\item{\code{inefficiency_beta}: }{a vector with inefficiency facor for each \eqn{\beta}.}
\item{\code{inefficiency_sigma}: }{a vector with inefficiency factor for each \eqn{\sigma}.}
}
}
\description{
This function calculates the inefficiency factor from the MCMC draws
of \eqn{(\beta, \sigma)} for an ordinal model with 3 outcomes. The
inefficiency factor is calculated using the batch-means method.
}
\details{
Calculates the inefficiency factor of \eqn{(\beta, \sigma)} using the batch-means
method.
}
\examples{
set.seed(101)
data("data25j3")
x <- data25j3$x
y <- data25j3$y
p <- 0.25
ans <- quan_reg3(y, x, mc = 50, p)
beta_draws <- ans$beta_draws
sigma_draws <- ans$sigma_draws

inefficiency <- inefficiency_factor3(beta_draws, 2, sigma_draws)

# inefficiency_beta
#   1.322590
#   1.287309
#   1.139322
# inefficiency_sigma
#   1.392045

}
\references{
Greenberg, E. (2012). “Introduction to Bayesian Econometrics.”
 Cambridge University Press, Cambridge.
}
\seealso{
pracma
}
