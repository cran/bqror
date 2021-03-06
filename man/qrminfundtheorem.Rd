% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODRI.R
\name{qrminfundtheorem}
\alias{qrminfundtheorem}
\title{Minimize the negative of log-likelihood}
\usage{
qrminfundtheorem(
  deltain,
  y,
  x,
  beta,
  cri0,
  cri1,
  stepsize,
  maxiter,
  h,
  dh,
  sw,
  p
)
}
\arguments{
\item{deltain}{initialization of cut-points.}

\item{y}{dependent variable i.e. ordinal outcome values.}

\item{x}{covariate matrix of dimension \eqn{(n x k)} including a column of ones.}

\item{beta}{column vector of coeffcients of dimension \eqn{(k x 1)}.}

\item{cri0}{initial criterion, \eqn{cri0 = 1}.}

\item{cri1}{criterion lies between (0.001 to 0.0001).}

\item{stepsize}{learning rate lies between (0.1, 1).}

\item{maxiter}{maximum number of iteration.}

\item{h}{change in value of each \eqn{\delta}, holding other \eqn{\delta}
constant for first derivatives.}

\item{dh}{change in each value of \eqn{\delta}, holding other \eqn{\delta} constant
for second derivaties.}

\item{sw}{iteration to switch from BHHH to inv(-H) algorithm.}

\item{p}{quantile level or skewness parameter, p in (0,1).}
}
\value{
Returns a list with components
\itemize{
\item{\code{dmin}: }{a vector with cutpoints that minimize the log-likelihood function.}
\item{\code{sumlogl}: }{a scalar with sum of log-likelihood values.}
\item{\code{logl}: }{a vector with log-likelihood values.}
\item{\code{G}: }{a gradient vector, \eqn{(n x k)} matrix with i-th row as the score
for the i-th unit.}
\item{\code{H}: }{represents Hessian matrix.}
}
}
\description{
This function minimizes the negative of the log-likelihood for an
ordinal quantile model with respect to the cut-points \eqn{\delta} using the
Fundamental Theorem of Calculus.
}
\details{
First derivative from first principle
\deqn{dy/dx=[f(x+h)-f(x-h)]/2h}

Second derivative from First principle

\deqn{f'(x-h)=(f(x)-f(x-h))/h}

\deqn{f''(x)= [{(f(x+h)-f(x))/h} - (f(x)-f(x-h))/h]/h}

      \deqn{= [(f(x+h)+f(x-h)-2 f(x))]/h^2}

cross partial derivatives

\deqn{f(x) = [f(x+dh,y)-f(x-dh,y)]/2dh}

\deqn{f(x,y)=[{(f(x+dh,y+dh) - f(x+dh,y-dh))/2dh} - {(f(x-dh,y+dh) -
f(x-dh,y-dh))/2dh}]/2dh}

\deqn{= 0.25* [{(f(x+dh,y+dh)-f(x+dh,y-dh))} -{(f(x-dh,y+dh)-f(x-dh,y-dh))}]/dh2}
}
\examples{
set.seed(101)
deltain <- c(-0.9026915, -2.2488833)
data("data25j4")
x <- data25j4$x
y <- data25j4$y
p <- 0.25
beta <- c(-1.429465, 1.135585, 2.107666)
cri0     <- 1
cri1     <- 0.001
stepsize <- 1
maxiter  <- 10
h        <- 0.002
dh       <- 0.0002
sw       <- 20
ans <- qrminfundtheorem(deltain, y, x, beta, cri0, cri1, stepsize, maxiter, h, dh, sw, p)

# deltamin
#   0.2674061 -0.6412074
# negsum
#   247.9525
# logl
#   -2.30530839
#   -1.60437267
#   -0.52085599
#   -0.93506872
#   -0.91064423
#   -0.49535299
#   -1.53635828
#   -1.36311002
#   -0.35753865
#   -0.55554991.. soon
# G
#   0.84555485  0.00000000
#   0.84555485  0.00000000
#   0.00000000  0.00000000
#   -0.32664119 -0.13166332
#   -0.32664119 -0.13166332
#   -0.32664119 -0.13166332
#   0.93042126  0.00000000
#   -0.32664119 -0.13166332
#   -0.32664119 -0.13166332
#   0.00000000  0.00000000
#   -0.32664119 -0.13166332.. soon
# H
#   -47.266464  -2.379509
#   -2.379509 -13.830474
# checkoutput
#   0    0    0    0    0    0    0 ... soon

}
\references{
Rahman, M. A. (2016). “Bayesian
Quantile Regression for Ordinal Models.”
Bayesian Analysis, 11(1): 1-24.
}
\seealso{
differential calculus, functional maximization,
\link[MASS]{ginv}, \link[pracma]{mldivide}
}
