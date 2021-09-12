#' Bayesian Quantile Regression for Ordinal Model
#' with 3 outcomes
#'
#' This function estimates Bayesian quantile regression for ordinal model with
#' 3 outcomes and reports the posterior mean, posterior standard deviation, and 95
#' percent credible intervals of \eqn{(\beta, \sigma)}.
#'
#' @param y                 observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x                 covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param b0                prior mean for normal distribution to sample \eqn{\beta}. Default is 0.
#' @param B0                prior variance for normal distribution to sample \eqn{\beta}
#' @param n0                prior for shape parameter to sample \eqn{\sigma} from inverse gamma distribution, default is 5.
#' @param d0                prior for scale parameter to sample \eqn{\sigma} from inverse gamma distribution, default is 8.
#' @param gamma             one and only cut-point other than 0.
#' @param mcmc              number of MCMC iterations, post burn-in.
#' @param p                 quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Function implements the Bayesian quantile regression for
#' ordinal model with 3 outcomes using a Gibbs sampling
#' procedure.
#'
#' Function initializes prior and then iteratively
#' samples \eqn{\beta}, \eqn{\sigma} and latent variable z.
#' Burn-in is taken as \eqn{0.25*mcmc} and \eqn{nsim = burn}-\eqn{in + mcmc}.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{postMeanbeta}: }{a vector with mean of sampled
#'  \eqn{\beta} for each covariate.}
#' \item{\code{postMeansigma}: }{a vector with mean of sampled
#'  \eqn{\sigma}.}
#' \item{\code{postStdbeta}: }{a vector with standard deviation
#'  of sampled \eqn{\beta} for each covariate.}
#'  \item{\code{postStdsigma}: }{a vector with standard deviation
#'  of sampled \eqn{\sigma}.}
#'  \item{\code{allQuantDIC}: }{results of the DIC criteria.}
#'  \item{\code{logMargLikelihood}: }{a scalar value for log marginal likelihood.}
#'  \item{\code{beta}: }{a matrix with all sampled values for \eqn{\beta}.}
#'  \item{\code{sigma}: }{a matrix with all sampled values for \eqn{\sigma}.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Yu, K. and Moyeed, R. A. (2001). “Bayesian Quantile Regression.” Statistics and
#' Probability Letters, 54(4): 437–447.
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @importFrom "tcltk" "tkProgressBar" "setTkProgressBar"
#' @importFrom "stats" "sd"
#' @importFrom "stats" "quantile"
#' @importFrom "pracma" "inv"
#' @seealso tcltk, \link[stats]{rnorm}, \link[stats]{qnorm},
#' Gibbs sampling
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' k <- dim(x)[2]
#' output <- quantreg_or2(y = y, x = x, B0 = 10*diag(k),
#' mcmc = 50, p = 0.25)
#'
#' # Number of burn-in draws : 12.5
#' # Number of retained draws : 50
#' # Summary of estimated beta :
#'
#' #            Post Mean Post Std Upper Credible Lower Credible
#' #    beta_0   -4.0802   0.6525        -3.0715        -5.4562
#' #    beta_1    6.0618   0.7255         7.3667         5.0018
#' #    beta_2    5.6469   0.8592         7.4881         4.3384
#' #    sigma     1.2251   0.1796         1.7111         1.0553
#'
#' # Model Marginal Likelihood: -508.04
#' # Model DIC: 779.46
#'
#' @export
quantreg_or2 <- function(y, x, b0 = 0, B0 , n0 = 5, d0 = 8, gamma = 3, mcmc = 15000, p) {
    cols <- colnames(x)
    names(x) <- NULL
    names(y) <- NULL
    x <- as.matrix(x)
    y <- as.matrix(y)
    if ( dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( length(mcmc) != 1){
        stop("parameter mcmc must be scalar")
    }
    if ( !is.numeric(mcmc)){
        stop("parameter mcmc must be a numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if ( any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if ( length(b0) != 1){
        stop("parameter b0 must be scalar")
    }
    if ( !all(is.numeric(b0))){
        stop("parameter b0 must be numeric")
    }
    if ( length(n0) != 1){
        stop("parameter n0 must be scalar")
    }
    if ( !all(is.numeric(n0))){
        stop("parameter n0 must be numeric")
    }
    if ( length(d0) != 1){
        stop("parameter d0 must be scalar")
    }
    if ( !all(is.numeric(d0))){
        stop("parameter d0 must be numeric")
    }
    J <- dim(as.array(unique(y)))[1]
    if ( J > 3 ){
        stop("This function is for 3 outcome
                variables. Please correctly specify the inputs
             to use quantreg_or2")
    }
    n <- dim(x)[1]
    k <- dim(x)[2]
    if ((dim(B0)[1] != (k)) | (dim(B0)[2] != (k))){
        stop("B0 is the prior variance to sample beta
             must have dimension kxk")
    }
    burn <- 0.25 * mcmc
    nsim <- burn + mcmc

    b0 <- array(rep(b0, k), dim = c(k, 1))
    invB0 <- inv(B0)
    invB0b0 <- invB0 %*% b0

    beta <- array(0, dim = c(k, nsim))
    sigma <- array(0, dim = c(1, nsim))
    btildeStore <- array(0, dim = c(k, nsim))
    BtildeStore <- array(0, dim = c(k, k, nsim))

    beta[, 1] <- array(rep(0,k), dim = c(k, 1))
    sigma[1] <- 2
    nu <- 5 * rep(1, n)
    gammacp <- array(c(-Inf, 0, gamma, Inf), dim = c(1, J+1))
    lambda <- 0.5
    theta <- (1 - 2 * p) / (p * (1 - p))
    tau <- sqrt(2 / (p * (1 - p)))
    tau2 <- tau^2

    z <- array( (rnorm(n, mean = 0, sd = 1)), dim = c(n, 1))

    pb <- tkProgressBar(title = "Simulation in Progress",
                        min = 0, max = nsim, width = 300)

    for (i in 2:nsim) {
        betadraw <- drawbeta_or2(z, x, sigma[(i - 1)], nu, tau2, theta, invB0, invB0b0)
        beta[, i] <- betadraw$beta
        btildeStore[, i] <- betadraw$btilde
        BtildeStore[, , i] <- betadraw$Btilde

        sigmadraw <- drawsigma_or2(z, x, beta[, i], nu, tau2, theta, n0, d0)
        sigma[i] <- sigmadraw$sigma

        nu <- drawnu_or2(z, x, beta[, i], sigma[i], tau2, theta, lambda)

        z <- drawlatent_or2(y, x, beta[, i], sigma[i], nu, theta, tau2, gammacp)

        setTkProgressBar(pb, i,
                         label = paste( round( (i / nsim) * 100, 0), "% done"))
    }
    close(pb)

    postMeanbeta <- rowMeans(beta[, (burn + 1):nsim])
    postStdbeta <- apply(beta[, (burn + 1):nsim], 1, sd)
    postMeansigma <- mean(sigma[(burn + 1):nsim])
    postStdsigma <- std(sigma[(burn + 1):nsim])

    allQuantDIC <- deviance_or2(y, x, gammacp, p,
                            postMeanbeta, postStdbeta,
                            postMeansigma, postStdsigma,
                            beta, sigma, burn, nsim)

    logMargLikelihood <- logmargLikelihood_or2(y, x, b0, B0,
                                                n0, d0, postMeanbeta,
                                                postMeansigma, btildeStore,
                                                BtildeStore, gamma, p)

    postMeanbeta <- array(postMeanbeta, dim = c(k, 1))
    postStdbeta <- array(postStdbeta, dim = c(k, 1))
    postMeansigma <- array(postMeansigma)
    postStdsigma <- array(postStdsigma)

    upperCrediblebeta <- array(apply(beta[ ,(burn + 1):nsim], 1, quantile, c(0.975)), dim = c(k, 1))
    lowerCrediblebeta <- array(apply(beta[ ,(burn + 1):nsim], 1, quantile, c(0.025)), dim = c(k, 1))
    upperCrediblesigma <- quantile(sigma[(burn + 1):nsim], c(0.975))
    lowerCrediblesigma <- quantile(sigma[(burn + 1):nsim], c(0.025))

    allQuantbeta <- cbind(postMeanbeta, postStdbeta, upperCrediblebeta, lowerCrediblebeta)
    allQuantsigma <- cbind(postMeansigma, postStdsigma, upperCrediblesigma, lowerCrediblesigma)
    summary <- rbind(allQuantbeta, allQuantsigma)
    name <- list('Post Mean', 'Post Std', 'Upper Credible', 'Lower Credible')
    dimnames(summary)[[2]] <- name
    dimnames(summary)[[1]] <- letters[1:(k+1)]
    j <- 1
    if (is.null(cols)) {
        rownames(summary)[j] <- c('Intercept')
        for (i in paste0("beta_",1:k-1)) {
            rownames(summary)[j] = i
            j = j + 1
        }
    }
    else {
        for (i in cols) {
            rownames(summary)[j] = i
            j = j + 1
        }
    }
    rownames(summary)[j] <- 'sigma'

    print(noquote(paste0('Number of burn-in draws : ', burn)))
    print(noquote(paste0('Number of retained draws : ', mcmc)))
    print(noquote('Summary of estimated beta : '))
    cat("\n")
    print(round(summary, 4))
    cat("\n")
    print(noquote(paste0('Model Marginal Likelihood: ', round(logMargLikelihood, 2))))
    print(noquote(paste0('Model DIC: ', round(allQuantDIC$DIC, 2))))

    result <- list("postMeanbeta" = postMeanbeta,
                   "postStdbeta" = postStdbeta,
                   "postMeansigma" = postMeansigma,
                   "postStdsigma" = postStdsigma,
                   "allQuantDIC" = allQuantDIC,
                   "logMargLikelihood" = logMargLikelihood,
                   "beta" = beta,
                   "sigma" = sigma)
    return(result)
}
#' Samples the latent variable z for an Ordinal Model
#' with 3 outcomes
#'
#' This function samples the latent variable z from a truncated
#' normal distribution for an ordinal model with 3 outcomes.
#'
#' @param y         observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      column vector of coefficients of dimension \eqn{(k x 1)}.
#' @param sigma     scale factor, a scalar value.
#' @param nu        modified scale factor, row vector.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param gammacp   row vector of cut-points including -Inf and Inf.
#'
#' @details
#' Function samples the latent variable z from a truncated normal
#' distribution.
#'
#' @return Returns a column vector of values for latent variable z.
#'
#' @references Albert, J. and Chib, S. (1993). “Bayesian Analysis of Binary and Polychotomous
#' Response Data.” Journal of the American Statistical Association, 88(422): 669–679.
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @seealso Gibbs sampling, truncated normal distribution,
#' \link[truncnorm]{rtruncnorm}
#' @importFrom "truncnorm" "rtruncnorm"
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' beta <- c(1.810504, 1.850332, 6.181163)
#' sigma <- 0.9684741
#' nu <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
#' theta <- 2.6667
#' tau2 <- 10.6667
#' gammacp <- c(-Inf, 0, 3, Inf)
#' output <- drawlatent_or2(y, x, beta, sigma, nu,
#' theta, tau2, gammacp)
#'
#' # output
#' #   1.257096 10.46297 4.138694
#' #   28.06432 4.179275 19.21582
#' #   11.17549 13.79059 28.3650 .. soon
#'
#' @export
drawlatent_or2 <- function(y, x, beta, sigma, nu, theta, tau2, gammacp) {
    if ( dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    if ( !all(is.numeric(nu))){
        stop("each entry in nu must be numeric")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("parameter tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("parameter theta must be numeric")
    }
    n <- dim(y)[1]
    z <- array(0, dim = c(n, 1))
    for (i in 1:n) {
        meancomp <- (x[i, ] %*% beta) + (theta * nu[i])
        std <- sqrt(tau2 * sigma * nu[i])
        temp <- y[i]
        a1 <- gammacp[temp]
        b1 <- gammacp[temp + 1]
        z[i, 1] <- rtruncnorm(n = 1, a = a1, b = b1, mean = meancomp, sd = std)
    }
    return(z)
}
#' Samples \eqn{\beta} for an Ordinal Model
#' with 3 outcomes
#'
#' This function samples \eqn{\beta} from its conditional
#' posterior distribution for an ordinal model with 3
#' outcomes.
#'
#' @param z         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param sigma     scale factor, a scalar value.
#' @param nu        modified scale factor, row vector.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param invB0     inverse of prior covariance matrix of normal distribution.
#' @param invB0b0   prior mean pre-multiplied by invB0.
#'
#' @details
#' Function samples a vector of \eqn{\beta} from a multivariate normal distribution.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{beta}: }{Returns a column vector of \eqn{\beta}
#' from a multivariate normal distribution.}
#' \item{\code{Btilde}: }{the variance parameter for the normal
#'  distribution.}
#' \item{\code{btilde}: }{the mean parameter for the
#' normal distribution.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988).
#' “The New S Language. Wadsworth & Brooks/Cole.”
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @importFrom "MASS" "mvrnorm"
#' @importFrom "pracma" "inv"
#' @seealso Gibbs sampling, normal distribution
#' , \link[GIGrvg]{rgig}, \link[pracma]{inv}
#' @examples
#' set.seed(101)
#' z <- c(21.01744, 33.54702, 33.09195, -3.677646,
#'  21.06553, 1.490476, 0.9618205, -6.743081, 21.02186, 0.6950479)
#' x <- matrix(c(
#'      1, -0.3010490, 0.8012506,
#'      1,  1.2764036, 0.4658184,
#'      1,  0.6595495, 1.7563655,
#'      1, -1.5024607, -0.8251381,
#'      1, -0.9733585, 0.2980610,
#'      1, -0.2869895, -1.0130274,
#'      1,  0.3101613, -1.6260663,
#'      1, -0.7736152, -1.4987616,
#'      1,  0.9961420, 1.2965952,
#'      1, -1.1372480, 1.7537353),
#'      nrow = 10, ncol = 3, byrow = TRUE)
#' sigma <- 1.809417
#' nu <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
#' tau2 <- 10.6667
#' theta <- 2.6667
#' invB0 <- matrix(c(
#'      1, 0, 0,
#'      0, 1, 0,
#'      0, 0, 1),
#'      nrow = 3, ncol = 3, byrow = TRUE)
#' invB0b0 <- c(0, 0, 0)
#'
#' output <- drawbeta_or2(z, x, sigma, nu, tau2, theta, invB0, invB0b0)
#'
#' # output$beta
#' #   -0.74441 1.364846 0.7159231
#'
#' @export
drawbeta_or2 <- function(z, x, sigma, nu, tau2, theta, invB0, invB0b0) {
    if ( !all(is.numeric(z))){
        stop("each entry in z must be numeric")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    if ( !all(is.numeric(nu))){
        stop("each entry in nu must be numeric")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("parameter tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("parameter theta must be numeric")
    }
    if ( !all(is.numeric(invB0))){
        stop("each entry in invB0 must be numeric")
    }
    if ( !all(is.numeric(invB0b0))){
        stop("each entry in invB0b0 must be numeric")
    }
    n <- dim(x)[1]
    k <- dim(x)[2]
    meancomp <- array(0, dim = c(n, k))
    varcomp <- array(0, dim = c(k, k, n))
    q <- array(0, dim = c(1, k))
    eye <- diag(k)
    for (i in 1:n) {
        meancomp[i, ] <- (x[i, ] * (z[i] - (theta * nu[i])) ) / (tau2 * sigma * nu[i])
        varcomp[, , i] <- ( x[i, ] %*% t(x[i, ])) / (tau2 * sigma * nu[i])
    }
    Btilde <- inv(invB0 + rowSums(varcomp, dims = 2))
    btilde <- Btilde %*% (invB0b0 + colSums(meancomp))
    L <- t(chol(Btilde))
    beta <- btilde + L %*%  (mvrnorm(n = 1, mu = q, Sigma = eye))

    betaReturns <- list("beta" = beta,
                   "Btilde" = Btilde,
                   "btilde" = btilde)

    return(betaReturns)
}
#' Samples the \eqn{\sigma} for an Ordinal Model
#' with 3 outcomes
#'
#' This function samples the \eqn{\sigma} from an inverse-gamma distribution
#' for an ordinal model with 3 outcomes.
#'
#' @param z         Gibbs draw of latent response variable, a column vector.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      Gibbs draw of coefficients of dimension \eqn{(k x 1)}.
#' @param nu        modified scale factor, row vector.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param n0        prior hyper-parameter for \eqn{\sigma}.
#' @param d0        prior hyper-parameter for \eqn{\sigma}.
#'
#' @details
#' Function samples the \eqn{\sigma} from an inverse
#' gamma distribution.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{sigma}: }{column vector of the \eqn{\sigma}
#' from an inverse gamma distribution.}
#' \item{\code{dtilde}: }{the scale parameter for the inverse
#' gamma distribution.}
#' }
#'
#' @importFrom "stats" "rgamma"
#'
#' @references Albert, J. and Chib, S. (1993). “Bayesian Analysis of Binary and Polychotomous
#' Response Data.” Journal of the American Statistical Association, 88(422): 669–679.
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @seealso \link[stats]{rgamma}, Gibbs sampling
#' @examples
#' set.seed(101)
#' z <- c(21.01744, 33.54702, 33.09195, -3.677646,
#'  21.06553, 1.490476, 0.9618205, -6.743081, 21.02186, 0.6950479)
#' x <- matrix(c(
#'      1, -0.3010490, 0.8012506,
#'      1,  1.2764036, 0.4658184,
#'      1,  0.6595495, 1.7563655,
#'      1, -1.5024607, -0.8251381,
#'      1, -0.9733585, 0.2980610,
#'      1, -0.2869895, -1.0130274,
#'      1,  0.3101613, -1.6260663,
#'      1, -0.7736152, -1.4987616,
#'      1,  0.9961420, 1.2965952,
#'      1, -1.1372480, 1.7537353),
#'      nrow = 10, ncol = 3, byrow = TRUE)
#' beta <- c(-0.74441, 1.364846, 0.7159231)
#' nu <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
#' tau2 <- 10.6667
#' theta <- 2.6667
#' n0 <- 5
#' d0 <- 8
#' output <- drawsigma_or2(z, x, beta, nu, tau2, theta, n0, d0)
#'
#' # output$sigma
#' #   3.749524
#'
#' @export
drawsigma_or2 <- function(z, x, beta, nu, tau2, theta, n0, d0) {
    if ( !all(is.numeric(z))){
        stop("each entry in z must be numeric")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( !all(is.numeric(nu))){
        stop("each entry in nu must be numeric")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("parameter tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("parameter theta must be numeric")
    }
    if ( length(n0) != 1){
        stop("parameter n0 must be scalar")
    }
    if ( !all(is.numeric(n0))){
        stop("parameter n0 must be numeric")
    }
    if ( length(d0) != 1){
        stop("parameter d0 must be scalar")
    }
    if ( !all(is.numeric(d0))){
        stop("parameter d0 must be numeric")
    }
    n <- dim(x)[1]
    ntilde <- n0 + (3 * n)
    temp <- array(0, dim = c(n, 1))
    for (i in 1:n) {
        temp[i, 1] <- (( z[i] - x[i, ] %*% beta - theta * nu[i] )^2) / (tau2 * nu[i])
    }
    dtilde <- sum(temp) + d0 + (2 * sum(nu))
    sigma <- 1/rgamma(n = 1, shape = (ntilde / 2), scale = (2 / dtilde))

    sigmaReturns <- list("sigma" = sigma,
                   "dtilde" = dtilde)
    return(sigmaReturns)
}

#' Samples the scale factor \eqn{\nu} for an Ordinal Model
#' with 3 outcomes
#'
#' This function samples the \eqn{\nu} from a generalized inverse Gaussian (GIG)
#' distribution for an ordinal model with 3 outcomes.
#'
#' @param z         Gibbs draw of latent response variable, a column vector.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      Gibbs draw of coefficients of dimension \eqn{(k x 1)}.
#' @param sigma     scale factor, a scalar.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param lambda    index parameter of GIG distribution which is equal to 0.5
#'
#' @details
#' Function samples the \eqn{\nu} from a GIG
#' distribution.
#'
#' @return Returns a row vector of the \eqn{\nu}
#' from GIG distribution.
#'
#' @references  Rahman, M. A. (2016), “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1), 1-24.
#'
#' Dagpunar, J. S. (1989). "An Easily Implemented Generalised Inverse Gaussian Generator."
#' Communication Statistics Simulation, 18: 703-710.
#'
#' @importFrom "GIGrvg" "rgig"
#' @seealso GIGrvg, Gibbs sampling, \link[GIGrvg]{rgig}
#' @examples
#' set.seed(101)
#' z <- c(21.01744, 33.54702, 33.09195, -3.677646,
#'  21.06553, 1.490476, 0.9618205, -6.743081, 21.02186, 0.6950479)
#' x <- matrix(c(
#'      1, -0.3010490, 0.8012506,
#'      1,  1.2764036, 0.4658184,
#'      1,  0.6595495, 1.7563655,
#'      1, -1.5024607, -0.8251381,
#'      1, -0.9733585, 0.2980610,
#'      1, -0.2869895, -1.0130274,
#'      1,  0.3101613, -1.6260663,
#'      1, -0.7736152, -1.4987616,
#'      1, 0.9961420, 1.2965952,
#'      1, -1.1372480, 1.7537353),
#'      nrow = 10, ncol = 3, byrow = TRUE)
#' beta <- c(-0.74441, 1.364846, 0.7159231)
#' sigma <- 3.749524
#' tau2 <- 10.6667
#' theta <- 2.6667
#' lambda <- 0.5
#' output <- drawnu_or2(z, x, beta, sigma, tau2, theta, lambda)
#'
#' # output
#' #   5.177456 4.042261 8.950365
#' #   1.578122 6.968687 1.031987
#' #   4.13306 0.4681557 5.109653
#' #   0.1725333
#'
#' @export
drawnu_or2 <- function(z, x, beta, sigma, tau2, theta, lambda) {
    if ( !all(is.numeric(z))){
        stop("each entry in z must be numeric")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("parameter tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("parameter theta must be numeric")
    }
    if ( length(lambda) != 1){
        stop("parameter lambda must be scalar")
    }
    if ( !all(is.numeric(lambda))){
        stop("parameter lambda must be numeric")
    }
    n <- dim(x)[1]
    tildegamma2 <- ( (theta ^ 2) / (tau2 * sigma)) + (2 / sigma)
    tildedelta2 <- array(0, dim = c(n, 1))
    nu <- array(0, dim = c(n, 1))
    for (i in 1:n) {
        tildedelta2[i, 1] <- ( (z[i] - x[i, ] %*%  beta)^2) / (tau2 * sigma)
        nu[i, 1] <- rgig(n = 1, lambda = lambda,
                         chi = tildedelta2[i, 1],
                         psi = tildegamma2)
    }
    return(nu)
}

#' Deviance Information Criteria for Ordinal Model
#' with 3 outcomes
#'
#' Function for computing the Deviance information criteria for ordinal
#' model with 3 outcomes.
#'
#' @param y              observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x              covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param gammacp        row vector of cut-points including -Inf and Inf.
#' @param p              quantile level or skewness parameter, p in (0,1).
#' @param postMeanbeta   mean value of \eqn{\beta} obtained from MCMC draws.
#' @param postStdbeta    standard deviation of \eqn{\beta} obtained from MCMC draws.
#' @param postMeansigma  mean value of \eqn{\sigma} obtained from MCMC draws.
#' @param postStdsigma   standard deviation of \eqn{\sigma} obtained from MCMC draws.
#' @param beta           MCMC draw of coefficients, dimension is \eqn{(k x nsim)}.
#' @param sigma          MCMC draw of scale factor, dimension is \eqn{(nsim x 1)}.
#' @param burn           number of discarded MCMC iterations.
#' @param nsim           total number of MCMC iterations including the burn-in.
#'
#' @details
#' The deviance is -2*(log likelihood) and has an important role in
#' statistical model comparison because of its relation with Kullback-Leibler
#' information criteria.
#'
#' @return Returns a list with components
#' \deqn{DIC = 2*avgdeviance - devpostmean}
#' \deqn{pd = avgdeviance - devpostmean}
#' \deqn{devpostmean = -2*(logLikelihood)}.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin B. P. and Linde A. (2002).
#' “Bayesian Measures of Model Complexity and Fit.” Journal of the
#' Royal Statistical Society B, Part 4: 583-639.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B.
#' “Bayesian Data Analysis.” 2nd Edition, Chapman and Hall.
#'
#' @seealso  decision criteria
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' k <- dim(x)[2]
#' output <- quantreg_or2(y = y, x = x, B0 = 10*diag(k),
#' mcmc = 50, p = 0.25)
#' gammacp <- c(-Inf, 0, 3, Inf)
#' postMeanbeta <- output$postMeanbeta
#' postStdbeta <- output$postStdbeta
#' postMeansigma <- output$postMeansigma
#' postStdsigma <- output$postStdsigma
#' beta <- output$beta
#' sigma <- output$sigma
#' mcmc = 50
#' burn <- 10
#' nsim <- burn + mcmc
#' deviance <- deviance_or2(y, x, gammacp, p = 0.25, postMeanbeta, postStdbeta,
#' postMeansigma, postStdsigma, beta, sigma, burn, nsim)
#'
#' # DIC
#' #   779.4566
#' # pd
#' #   5.233374
#' # devpostmean
#' #   768.9898
#'
#' @export
deviance_or2 <- function(y, x, gammacp, p, postMeanbeta, postStdbeta,
                      postMeansigma, postStdsigma,
                      beta, sigma, burn, nsim) {
    cols <- colnames(x)
    names(x) <- NULL
    names(y) <- NULL
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if ( !all(is.numeric(postMeanbeta))){
        stop("each entry in postMeanbeta must be numeric")
    }
    if ( !all(is.numeric(postStdbeta))){
        stop("each entry in postStdbeta must be numeric")
    }
    if ( !all(is.numeric(postMeansigma))){
        stop("each entry in postMeansigma must be numeric")
    }
    if ( !all(is.numeric(postStdsigma))){
        stop("each entry in postStdsigma must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( !all(is.numeric(sigma))){
        stop("each entry in sigma must be numeric")
    }
    if ( length(burn) != 1){
        stop("parameter burn must be scalar")
    }
    if ( length(nsim) != 1){
        stop("parameter nsim must be scalar")
    }
    k <- dim(x)[2]
    devpostmean <- array(0, dim = c(1))
    DIC <- array(0, dim = c(1))
    pd <- array(0, dim = c(1))
    devpostmean <- 2 * negLoglikelihood(y, x, gammacp, postMeanbeta, postMeansigma, p)

    postBurnin <- dim(beta[, (burn + 1):nsim])[2]
    Deviance <- array(0, dim = c(1, postBurnin))
    for (i in 1:postBurnin) {
        Deviance[1, i] <- 2 * negLoglikelihood(y, x, gammacp,
                                                   beta[ ,(burn + i)],
                                                   sigma[ ,(burn + i)],
                                                   p)
    }
    avgDeviance <- mean(Deviance)
    DIC <- (2 * avgDeviance) - devpostmean
    pd <- avgDeviance - devpostmean
    result <- list("DIC" = DIC,
                   "pd" = pd,
                   "devpostmean" = devpostmean)
    return(result)
}
#' NegLoglikelihood function for Ordinal Model with 3 outcomes
#'
#' This function computes the negative of the log-likelihood for quantile
#' ordinal model with 3 outcomes where the error is assumed to follow
#' an asymmetric Laplace distribution.
#'
#' @param y         observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param gammacp   row vector of cutpoints including -Inf and Inf.
#' @param beta      column vector of coefficients of dimension \eqn{(k x 1)}.
#' @param sigma     scale factor, a scalar.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Computes the negative of the log-likelihood for quantile
#' ordinal model with 3 outcomes where the error is assumed to follow
#' an asymmetric Laplace distribution.
#'
#' @return Returns the negative log-likelihood value.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' @seealso likelihood maximization
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' p <- 0.25
#' gammacp <- c(-Inf, 0, 3, Inf)
#' beta <- c(1.810504, 1.850332, 6.18116)
#' sigma <- 0.9684741
#' output <- negLoglikelihood(y, x, gammacp, beta, sigma, p)
#'
#' # output
#' #   791.5215
#'
#' @export
negLoglikelihood <- function(y, x, gammacp, beta, sigma, p) {
    cols <- colnames(x)
    names(x) <- NULL
    names(y) <- NULL
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    J <- dim(unique(y))[1]
    n <- dim(y)[1]
    lnpdf <- array(0, dim = c(n, 1))
    mu <- x %*% beta
    for (i in 1:n) {
        meanf <- mu[i]
        if (y[i] == 1) {
            lnpdf[i] <- log(alcdf(0, meanf, sigma, p))
        }
        else if (y[i] == J) {
            lnpdf[i] <- log(1 - alcdf(gammacp[J], meanf, sigma, p))
        }
        else {
            w <- (alcdf(gammacp[J], meanf, sigma, p) -
                      alcdf(gammacp[(J - 1)], meanf, sigma, p))
            lnpdf[i] <- log(w)
        }
    }
    negsuminpdf <- -sum(lnpdf)
    return(negsuminpdf)
}
#' Generates random numbers from an Asymmetric Laplace Distribution
#'
#' This function generates a vector of random numbers from an asymmetric
#' Laplace distribution with quantile p.
#'
#' @param sigma  scale factor, a scalar.
#' @param p      quantile or skewness parameter, p in (0,1).
#' @param n      number of observations
#'
#' @details
#' Generates a vector of random numbers from an asymmetric Laplace distribution,
#' as a mixture of normal–exponential distributions.
#'
#' @return Returns a vector \eqn{(n x 1)} of random numbers using an AL(0, \eqn{\sigma}, p)
#'
#' @references
#'  Kozumi, H. and Kobayashi, G. (2011). “Gibbs Sampling Methods for Bayesian Quantile Regression.”
#'  Journal of Statistical Computation and Simulation, 81(11): 1565–1578.
#'
#'  Koenker, R. and Machado, J. (1999). “Goodness of Fit and Related
#'  Inference Processes for Quantile Regression.”, Journal of
#'   American Statistics Association, 94(3): 1296-1309.
#'
#'  Keming Yu and Jin Zhang (2005). “A Three-Parameter Asymmetric
#'  Laplace Distribution.” Communications in Statistics - Theory and Methods: 1867-1879.
#'
#' @importFrom "stats" "rnorm" "rexp"
#' @seealso asymmetric Laplace distribution
#' @examples
#' set.seed(101)
#' sigma <- 2.503306
#' p <- 0.25
#' n <- 1
#' output <- rndald(sigma, p, n)
#'
#' # output
#' #   1.07328
#'
#' @export
rndald <- function(sigma, p, n){
    if ( any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if ( n != floor(n)){
        stop("parameter n must be an integer")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    u <- rnorm(n = n, mean = 0, sd = 1)
    w <- rexp(n = n, rate = 1)
    theta <- (1 - 2 * p) / (p * (1 - p))
    tau <- sqrt(2 / (p * (1 - p)))
    eps <- sigma * (theta * w + tau * sqrt(w) * u)
    return(eps)
}

#' Trace Plots for Ordinal Model with 3 outcomes
#'
#' This function generates trace plots of
#' MCMC samples for \eqn{(\beta ,\sigma)} in the quantile
#' regression model with 3 outcomes.
#'
#' @param beta      Gibbs draw of \eqn{\beta} vector of dimension \eqn{(k x nsim)}.
#' @param sigma     Gibbs draw of scale parameter, \eqn{\sigma}.
#' @param burn      number of discarded MCMC iterations.
#'
#' @details
#' Trace plot is a visual depiction of the values generated from the Markov chain
#' versus the iteration number.
#'
#' @return Returns trace plots for each element of \eqn{\beta} and \eqn{\sigma}.
#'
#' @importFrom "graphics" "plot"
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' @seealso traces in MCMC simulations
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' k <- dim(x)[2]
#' output <- quantreg_or2(y = y,x = x, B0 = 10*diag(k),
#' mcmc = 50, p = 0.25)
#' mcmc <- 50
#' beta <- output$beta
#' sigma <- output$sigma
#' traceplot_or2(beta, sigma, round(0.25*mcmc))
#'
#' @export
traceplot_or2 <- function(beta, sigma, burn) {
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( !all(is.numeric(sigma))){
        stop("each entry in sigma must be numeric")
    }

    nsim <- dim(beta)[2]
    k <- dim(beta)[1]
    for (t in 1:k) {
        plot(beta[t ,(burn + 1):nsim], xlab = '', ylab = '',type = "l", cex.main = 1.5, main = expression(paste(Trace~plot~beta)))
    }
    plot(sigma[, (burn + 1):nsim], xlab = '', ylab = '', type = "l", cex.main = 1.5, main = expression(paste(Trace~plot~sigma)))
}
#' Inefficiency Factor for Ordinal Model
#' with 3 outcomes
#'
#' This function calculates the inefficiency factor from the MCMC draws
#' of \eqn{(\beta, \sigma)} for an ordinal model with 3 outcomes. The
#' inefficiency factor is calculated using the batch-means method.
#'
#' @param x                         covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param beta                      Gibbs draw of coefficients of dimension \eqn{(k x nsim)}.
#' @param sigma                     Gibbs draw of scale factor.
#' @param autocorrelationCutoff     Cut-off to identify the number of lags.
#'
#' @details
#' Calculates the inefficiency factor of \eqn{(\beta, \sigma)} using the batch-means
#' method.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{inefficiencyBeta}: }{a vector with inefficiency factor for each \eqn{\beta}.}
#' \item{\code{inefficiencySigma}: }{a vector with inefficiency factor for each \eqn{\sigma}.}
#' }
#'
#' @importFrom "pracma" "Reshape" "std"
#' @importFrom "stats" "acf"
#'
#' @references Greenberg, E. (2012). “Introduction to Bayesian Econometrics.”
#'  Cambridge University Press, Cambridge.
#'
#' @seealso pracma, \link[stats]{acf}
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' k <- dim(x)[2]
#' output <- quantreg_or2(y = y, x = x, B0 = 10*diag(k),
#' mcmc = 50, p = 0.25)
#' beta <- output$beta
#' sigma <- output$sigma
#'
#' inefficiency <- infactor_or2(x, beta, sigma, 0.5)
#'
#' # Summary of Inefficiency Factor:
#' #            Inefficiency
#' # beta_0       1.0977
#' # beta_1       1.5840
#' # beta_2       1.5383
#' # sigma        2.8946
#'
#' @export
infactor_or2 <- function(x, beta, sigma, autocorrelationCutoff = 0.05) {
    cols <- colnames(x)
    names(x) <- NULL
    x <- as.matrix(x)
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    n <- dim(beta)[2]
    k <- dim(beta)[1]
    inefficiencyBeta <- array(0, dim = c(k, 1))
    for (i in 1:k) {
        autocorrelation <- acf(beta[i,], plot = FALSE)
        nlags <- min(which(autocorrelation$acf <= autocorrelationCutoff))
        nbatch <- floor(n / nlags)
        nuse <- nbatch * nlags
        b <- beta[i, 1:nuse]
        xbatch <- Reshape(b, nlags, nbatch)
        mxbatch <- colMeans(xbatch)
        varxbatch <- sum( (t(mxbatch) - mean(b)) *
                              (t(mxbatch) - mean(b))) / (nbatch - 1)
        nse <- sqrt(varxbatch / (nbatch))
        rne <- (std(b, 1) / sqrt( nuse )) / nse
        inefficiencyBeta[i, 1] <- 1 / rne
    }
    if ( !all(is.numeric(sigma))){
        stop("each entry in sigma must be numeric")
    }
    inefficiencySigma <- array(0, dim = c(1))
    autocorrelation <- acf(c(sigma), plot = FALSE)
    nlags <- min(which(autocorrelation$acf <= autocorrelationCutoff))
    nbatch2 <- floor(n / nlags)
    nuse2 <- nbatch2 * nlags
    b2 <- sigma[1:nuse2]
    xbatch2 <- Reshape(b2, nlags, nbatch2)
    mxbatch2 <- colMeans(xbatch2)
    varxbatch2 <- sum( (t(mxbatch2) - mean(b2)) *
                               (t(mxbatch2) - mean(b2))) / (nbatch2 - 1)
    nse2 <- sqrt(varxbatch2 / (nbatch2))
    rne2 <- (std(b2, 1) / sqrt( nuse2 )) / nse2
    inefficiencySigma <- 1 / rne2

    inefficiencyRes <- rbind(inefficiencyBeta, inefficiencySigma)
    name <- list('Inefficiency')
    dimnames(inefficiencyRes)[[2]] <- name
    dimnames(inefficiencyRes)[[1]] <- letters[1:(k+1)]
    j <- 1
    if (is.null(cols)) {
        rownames(inefficiencyRes)[j] <- c('Intercept')
        for (i in paste0("beta_",1:k-1)) {
            rownames(inefficiencyRes)[j] = i
            j = j + 1
        }
    }
    else {
        for (i in cols) {
            rownames(inefficiencyRes)[j] = i
            j = j + 1
        }
    }
    rownames(inefficiencyRes)[j] <- 'sigma'

    print(noquote('Summary of Inefficiency Factor: '))
    cat("\n")
    print(round(inefficiencyRes, 4))

    result <- list("inefficiencyBeta" = inefficiencyBeta,
                   "inefficiencySigma" = inefficiencySigma)

    return(result)
}
#' Marginal Likelihood for Ordinal Model
#' with 3 outcomes
#'
#' This function estimates the marginal likelihood for ordinal model with
#' 3 outcomes
#'
#' @param y                 observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x                 covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param b0                prior mean for normal distribution to sample \eqn{\beta}.
#' @param B0                prior variance for normal distribution to sample \eqn{\beta}
#' @param n0                prior for shape parameter to sample \eqn{\sigma} from inverse gamma distribution.
#' @param d0                prior for scale parameter to sample \eqn{\sigma} from inverse gamma distribution.
#' @param postMeanbeta      a vector with mean of sampled \eqn{\beta} for each covariate..
#' @param postMeansigma     a vector with mean of sampled \eqn{\sigma}.
#' @param btildeStore       a storage matrix for posterior mean of \eqn{\beta}.
#' @param BtildeStore       a storage matrix for posterior variance of \eqn{\beta}.
#' @param gamma             one and only cut-point other than 0.
#' @param p                 quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Function implements the marginal likelihood for
#' ordinal model with 3 outcomes using a Gibbs sampling
#' procedure.
#'
#' @return Returns a scalar for marginal likelihood
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' @importFrom "stats" "sd" "dnorm"
#' @importFrom "invgamma" "dinvgamma"
#' @importFrom "pracma" "inv"
#' @importFrom "NPflow" "mvnpdf"
#' @seealso \link[invgamma]{dinvgamma}, \link[NPflow]{mvnpdf}, \link[stats]{dnorm},
#' Gibbs sampling
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' k <- dim(x)[2]
#' output <- quantreg_or2(y = y, x = x, B0 = 10*diag(k),
#' mcmc = 50, p = 0.25)
#' # output$logMargLikelihood
#' #   -508.0386
#'
#' @export
logmargLikelihood_or2 <- function(y, x, b0, B0, n0, d0, postMeanbeta, postMeansigma, btildeStore, BtildeStore, gamma, p) {
    cols <- colnames(x)
    names(x) <- NULL
    names(y) <- NULL
    x <- as.matrix(x)
    y <- as.matrix(y)
    if ( dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if ( any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if ( !all(is.numeric(b0))){
        stop("each entry in b0 must be numeric")
    }
    if ( length(n0) != 1){
        stop("parameter n0 must be scalar")
    }
    if ( !all(is.numeric(n0))){
        stop("parameter n0 must be numeric")
    }
    if ( length(d0) != 1){
        stop("parameter d0 must be scalar")
    }
    if ( !all(is.numeric(d0))){
        stop("parameter d0 must be numeric")
    }
    J <- dim(as.array(unique(y)))[1]
    if ( J > 3 ){
        stop("This function is for 3 outcome
                variables. Please correctly specify the inputs
             to use quantreg_or2")
    }
    n <- dim(x)[1]
    k <- dim(x)[2]
    nsim <- dim(btildeStore)[2]
    burn <- (0.25 * nsim) / (1.25)
    nu <- 5 * rep(1, n)
    ntilde <- n0 + (3 * n)
    gammacp <- array(c(-Inf, 0, gamma, Inf), dim = c(1, J+1))
    lambda <- 0.5
    theta <- (1 - 2 * p) / (p * (1 - p))
    tau <- sqrt(2 / (p * (1 - p)))
    tau2 <- tau^2
    sigmaRedrun <- array(0, dim = c(1, nsim))
    dtildeStoreRedrun <- array(0, dim = c(1, nsim))
    z <- array( (rnorm(n, mean = 0, sd = 1)), dim = c(n, 1))
    b0 <- array(rep(b0, k), dim = c(k, 1))
    j <- 1
    postOrdbetaStore <- array(0, dim=c((nsim-burn),1))
    postOrdsigmaStore <- array(0, dim=c((nsim-burn),1))

    pb <- tkProgressBar(title = "Reduced Run in Progress",
                        min = 0, max = nsim, width = 300)
    for (i in 1:nsim) {
        sigmaStoreRedrun <- drawsigma_or2(z, x, postMeanbeta, nu, tau2, theta, n0, d0)
        sigmaRedrun[i] <- sigmaStoreRedrun$sigma
        dtildeStoreRedrun[i] <- sigmaStoreRedrun$dtilde

        nu <- drawnu_or2(z, x, postMeanbeta, sigmaRedrun[i], tau2, theta, lambda)

        z <- drawlatent_or2(y, x, postMeanbeta, sigmaRedrun[i], nu, theta, tau2, gammacp)

        setTkProgressBar(pb, i,
                         label = paste(round( (i / nsim) * 100, 0), "% done"))
    }
    close(pb)

    sigmaStar <- mean(sigmaRedrun[(burn + 1):nsim])

    for (i in (burn+1):(nsim)) {
        postOrdbetaStore[j] <- mvnpdf(x = matrix(postMeanbeta), mean = btildeStore[, i], varcovM = BtildeStore[, , i], Log = FALSE)
        postOrdsigmaStore[j] <- (dinvgamma(sigmaStar, shape = (ntilde / 2), scale = (2 / dtildeStoreRedrun[i])))
        j <- j  + 1
    }
    postOrdbeta <- mean(postOrdbetaStore)
    postOrdsigma <- mean(postOrdsigmaStore)

    priorContbeta <- mvnpdf(matrix(postMeanbeta), mean = b0, varcovM = B0, Log = FALSE)
    priorContsigma <- dinvgamma(postMeansigma, shape = (n0 / 2), scale = (2 / d0))

    logLikeCont <- -1 * negLoglikelihood(y, x, gammacp, postMeanbeta, postMeansigma, p)
    logPriorCont <- log(priorContbeta*priorContsigma)
    logPosteriorCont <- log(postOrdbeta*postOrdsigma)

    logMargLikelihood <- logLikeCont + logPriorCont - logPosteriorCont
    return(logMargLikelihood)
}
