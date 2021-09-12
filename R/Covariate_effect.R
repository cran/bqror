#' Covariate Effect for Bayesian Quantile Regression for Ordinal Model
#' with more than 3 outcomes
#'
#' This function estimates the change in probability of different ordinal
#' outcomes due to change in an independent variable, marginalized over
#' the parameters and values of other covariates
#'
#' @param model     outcome of the ODR I (quantreg_or1) model.
#' @param y         observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param modX      matrix x with suitable modification to an independent variable including a column of ones with or without column names.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Function estimates the covariate effect of the change in an independent
#' variable. The effect is marginalized over other parameters and values of other covariates.
#' This function can only be used for ordinal models with more than 3 outcomes. Function uses
#' the output of ODR I model to compute the covariate effect.
#'
#' @return Returns a list with components:
#' \itemize{
#' \item{\code{avgDiffProb}: }{a vector with change in predicted
#' probabilities for each outcome category.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Jeliazkov, I., Graves, J., and Kutzbach, M. (2008). “Fitting and Comparison of Models
#' for Multivariate Ordinal Outcomes.” Advances in Econometrics: Bayesian Econo-
#' metrics, 23: 115–156.
#'
#' @importFrom "stats" "sd"
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' k <- dim(x)[2]
#' J <- dim(as.array(unique(y)))[1]
#' D0 <- 0.25*diag(J - 2)
#' output <- quantreg_or1(y = y,x = x, B0 = 10*diag(k), D0 = D0,
#' mcmc = 50, p = 0.25, tune = 1)
#' modX <- x
#' modX[,3] <- modX[,3] + 0.02
#' res <- covariateEffect_or1(output, y, x, modX, p = 0.25)
#'
#' # avgDiffProb
#' #   -0.007821158 -0.001108564 -0.008400506  0.009509070
#'
#' @export
covariateEffect_or1 <- function(model, y, x, modX, p) {
    cols <- colnames(x)
    cols1 <- colnames(modX)
    names(modX) <- NULL
    names(y) <- NULL
    names(x) <- NULL
    x <- as.matrix(x)
    modX <- as.matrix(modX)
    y <- as.matrix(y)
    J <- dim(as.array(unique(y)))[1]
    if ( J <= 3 ){
        stop("This function is only available for models with more than 3 outcome
                variables.")
    }
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(modX))){
        stop("each entry in modX must be numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    N <- dim(model$beta)[2]
    m <- (N)/(1.25)
    burn <- 0.25 * m
    n <- dim(x)[1]
    k <- dim(x)[2]
    betaBurnt <- model$beta[, (burn + 1):N]
    deltaBurnt <- model$delta[, (burn + 1):N]
    expdeltaBurnt <- exp(deltaBurnt)
    gammacpCov <- array(0, dim = c(J-1, m))
    for (j in 2:(J-1)) {
            gammacpCov[j, ] <- sum(expdeltaBurnt[1:(j - 1), 1])
    }
    mu <- 0
    sigma <- 1
    oldProb <- array(0, dim = c(n, m, J))
    newProb <- array(0, dim = c(n, m, J))
    oldComp <- array(0, dim = c(n, m, (J-1)))
    newComp <- array(0, dim = c(n, m, (J-1)))
    for (j in 1:(J-1)) {
        for (i in 1:m) {
            for (k in 1:n) {
                oldComp[k, i, j] <- alcdf((gammacpCov[j, i] - (x[k, ] %*% betaBurnt[, i])), mu, sigma, p)
                newComp[k, i, j] <- alcdf((gammacpCov[j, i] - (modX[k, ] %*% betaBurnt[, i])), mu, sigma, p)
            }
            if (j == 1) {
                oldProb[, i, j] <- oldComp[, i, j]
                newProb[, i, j] <- newComp[, i, j]
            }
            else {
                oldProb[, i, j] <- oldComp[, i, j] - oldProb[, i, (j-1)]
                newProb[, i, j] <- newComp[, i, j] - newProb[, i, (j-1)]
            }
        }
    }
    oldProb[, , J] = 1 - oldComp[, , (J-1)]
    newProb[, , J] = 1 - newComp[, , (J-1)]
    diffProb <- newProb - oldProb
    avgDiffProb = (colMeans(diffProb, dims = 2))
    result <- list("avgDiffProb" = avgDiffProb)

    return(result)
}
