#' Bayesian Quantile Regression for Ordinal Models
#'
#' @description
#' This package serves the following 3 purposes for Ordinal
#' Models under bayesian analysis:
#' \itemize{
#' \item{Package provides an estimation technique for
#' Bayesian quantile regression in ordinal models. Two algorithms are considered}
#' \itemize{
#' \item{one for an ordinal
#' model with three outcomes.}
#' \item{second for an ordinal
#' model with more than three outcomes.}
#' }
#' \item{Package provides model
#' performance criteria's.}
#' \item{It also provides trace plots for Markov chain Monte Carlo (MCMC) draws.}
#' }
#'
#' @details
#' \deqn{Package: bqror}
#' \deqn{Type: Package}
#' \deqn{Version: 1.1.0}
#' \deqn{License: GPL (>=2)}
#'
#' Package \strong{bqror} provides the following functions:
#'
#' \itemize{
#' \item{For an Ordinal Model with three outcomes:}
#' }
#' \code{\link[bqror]{quantreg_or2}}, \code{\link[bqror]{drawlatent_or2}},
#' \code{\link[bqror]{drawbeta_or2}}, \code{\link[bqror]{drawsigma_or2}},
#' \code{\link[bqror]{drawnu_or2}}, \code{\link[bqror]{deviance_or2}},
#' \code{\link[bqror]{negLoglikelihood}}, \code{\link[bqror]{rndald}},
#' \code{\link[bqror]{traceplot_or2}}, \code{\link[bqror]{infactor_or2}},
#' \code{\link[bqror]{logmargLikelihood_or2}}
#'
#' \itemize{
#' \item{For an Ordinal Model with more than three outcomes:}
#' }
#' \code{\link[bqror]{quantreg_or1}}, \code{\link[bqror]{qrminfundtheorem}},
#' \code{\link[bqror]{qrnegloglikensum}}, \code{\link[bqror]{drawbeta_or1}},
#' \code{\link[bqror]{draww_or1}}, \code{\link[bqror]{drawlatent_or1}},
#' \code{\link[bqror]{drawdelta_or1}}, \code{\link[bqror]{deviance_or1}},
#' \code{\link[bqror]{alcdfstd}}, \code{\link[bqror]{alcdf}},
#' \code{\link[bqror]{traceplot_or1}}, \code{\link[bqror]{infactor_or1}},
#' \code{\link[bqror]{covariateEffect_or1}}, \code{\link[bqror]{logmargLikelihood_or1}}
#'
#'
#' @author
#' Prof. Mohammad Arshad Rahman
#'
#' Prajual Maheshwari <prajual1391@gmail.com>
#'
#' @references
#' Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24, <doi:10.1214/15-BA939>.
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin B. P. and Linde A. (2002).
#' “Bayesian Measures of Model Complexity and Fit.” Journal of the
#' Royal Statistical Society B, Part 4: 583-639,  <doi:10.1111/1467-9868.00353>.
#'
#' Greenberg, E. (2012). “Introduction to Bayesian Econometrics.”
#' Cambridge University Press, Cambridge, <doi:10.1017/CBO9781139058414>.
#'
#' @seealso \link[GIGrvg]{rgig}, \link[MASS]{mvrnorm}, \link[MASS]{ginv},
#' \link[truncnorm]{rtruncnorm}, \link[NPflow]{mvnpdf},
#' \link[invgamma]{rinvgamma}, \link[pracma]{mldivide},
#' \link[pracma]{rand}, \link[stats]{qnorm},
#' \link[stats]{rexp}, \link[stats]{rnorm},
#' \link[pracma]{std}, \link[stats]{sd}, \link[stats]{acf},
#' \link[pracma]{Reshape}, \link[tcltk]{setTkProgressBar},
#' \link[tcltk]{tkProgressBar}, \link[invgamma]{dinvgamma}
#'
#' @docType package
#' @name bqror
NULL


