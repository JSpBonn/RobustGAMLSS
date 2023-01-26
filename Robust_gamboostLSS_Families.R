library(mboost)
library(gamboostLSS)

library(gamlss.dist) # needed for qGamma (quantile of gamma distribution)
library(EnvStats)    # needed for our offset calculation for gamma distribution

############################################################################################################
# functions, which are hidden in the Boosting package, but publicly available:
############################################################################################################

# function for weighted sd
weighted.sd <- function(x, w, ...) {
  if (missing(w))
    w <- rep(1, length(x))
  m <- weighted.mean(x, w, ...)
  var <- weighted.mean((x - m)^2, w, ...) * sum(w) / (sum(w) - 1)
  return(sqrt(var))
}

############################################################################################################

## helper function that stabilizes the negative gradient if requested by the user
stabilize_ngradient <- function(ngr, w = 1, stabilization) {
  ## set which to MAD if gamboostLSS_stab_ngrad = TRUE and which == "none"
  if (stabilization == "none" && getOption("gamboostLSS_stab_ngrad"))
    stabilization <- "MAD"
  ## stabilization using the mean absolute deviation (MAD)
  if (stabilization == "MAD") {
    div <- weighted.median(abs(ngr - weighted.median(ngr, w = w, na.rm = TRUE)),
                           w = w, na.rm = TRUE)
    div <- ifelse(div < 0.0001, 0.0001, div)
    ngr <- ngr / div
  }
  if (stabilization == "L2") {
    div <- sqrt(weighted.mean(ngr^2, w =w,  na.rm = TRUE))
    div <- ifelse(div < 1e-04, 1e-04, div)
    div <- ifelse(div > 1e+04, 1e+04, div)
    ngr <- ngr / div
  }
  ngr
}

############################################################################################################

check_stabilization <- function(stabilization = c("none", "MAD", "L2")) {
  stabilization <- match.arg(stabilization)
  ## check if old stabilization interface is used and issue a warning
  if (getOption("gamboostLSS_stab_ngrad")) {
    warning("Usage of ", sQuote("options(gamboostLSS_stab_ngrad = TRUE)")," is deprecated.\n", "Use argument ", sQuote("stabilization")," in the fitting family. See ?Families for details.")
    if (stabilization == "none")
      warning(sQuote("stabilization"), " is set to ", dQuote("MAD"))
  }
  stabilization
}

qNormal <- function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  qnorm(p = p, mean = mu, sd = sigma, lower.tail = lower.tail, log.p = log.p)
}


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

# Robust gamboostLSS:

############################################################################################################
############################################################################################################


# Gaussian
# Robust GaussianLSS


# determination of data driven robustness constant "c" for quantiles "tau":
# use only outcome values with weight = 1
c_generate_Gaussian <- function(outcome , tau = 0.05 ){
  log_like_c <- dnorm(x = outcome, mean = mean(outcome), sd =sd(outcome), log = TRUE)
  quant_c <- quantile(log_like_c,probs = tau)
  bl <- max(exp(-quant_c)-1,0.000001)       # lower boundary to secure an adequate value of C
  c_tau <- max(log(bl),0.25,na.rm = F) # lower boundary to secure an adequate value of C
  return(c_tau)
}

############################################################################################################

robust_GaussianMu <- function (mu = NULL, sigma = NULL, stabilization,rob=rob) 
{
  loss <- function(sigma, y, f) -rho(z=dnorm(x = y, mean = f, sd = sigma, log = TRUE),rob=rob)  # changes within the loss function
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, sigma = sigma))
  }
  ngradient <- function(y, f, w = 1) {
    ngr <-     rho_ab(z=-0.5*log(2*pi)-log(sigma)-0.5*(y-f)^2/sigma^2 ,rob=rob)*(1/sigma^2) * (y - f) 
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    return(ngr)
  }
  offset <- function(y, w) {
    if (!is.null(mu)) {
      RET <- mu
    }
    else {
      RET <- weighted.mean(y, w = w, na.rm = TRUE)
    }
    return(RET)
  }
  rho_ab <- function(z, rob=rob){          # log-logistic function embedded 
    exp(rob+z)/(exp(rob+z)+1)}             # log-logistic function embedded 
  rho <- function(z, rob=rob){             # log-logistic function embedded 
    log((1+exp(z+rob))/(1+exp(rob)))}      # log-logistic function embedded 
  
  mboost::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) f, 
                 offset = offset, name = "robust Normal distribution: mu(id link)") # "robust"
}

############################################################################################################

robust_GaussianSigma <- function (mu = NULL, sigma = NULL, stabilization,rob=rob) 
{
  loss <- function(y, f, mu) -rho(z=dnorm(x = y, mean = mu, sd = exp(f), log = TRUE),rob=rob)   # changes within the loss function
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu))
  }
  ngradient <- function(y, f, w = 1) {
    ngr <-  rho_ab(z=-0.5*log(2*pi)-f-0.5*(y-mu)^2/exp(2*f) ,rob=rob)*(-1 + exp(-2 * f) * ((y - mu)^2) )
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    return(ngr)
  }
  offset <- function(y, w) {
    if (!is.null(sigma)) {
      RET <- log(sigma)
    }
    else {
      RET <- log(weighted.sd(y,w=w, na.rm = TRUE)) #log(sd(y, na.rm = TRUE)) # eigentlich: 
    }
    return(RET)
  }
  rho_ab <- function(z, rob=rob){          # log-logistic function embedded 
    exp(rob+z)/(exp(rob+z)+1)}             # log-logistic function embedded 
  rho <- function(z, rob=rob){             # log-logistic function embedded 
    log((1+exp(z+rob))/(1+exp(rob)))}      # log-logistic function embedded 
  
  mboost::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f), 
                 offset = offset, name = "robust Normal distribution: sigma (log link)") ## "robust"
}

############################################################################################################

robust_GaussianLSS <- function (mu = NULL, sigma = NULL, stabilization = c("none", "MAD", "L2"),rob=2){
  if ((!is.null(sigma) && sigma <= 0)) 
    stop(sQuote("sigma"), " must be greater than zero")
  stabilization <- check_stabilization(stabilization)
  Families(mu = robust_GaussianMu(mu = mu, sigma = sigma, stabilization = stabilization,rob=rob),
           sigma = robust_GaussianSigma(mu = mu, sigma = sigma, stabilization = stabilization,rob=rob),
           qfun = qNormal, name = "robust Gaussian") # 
}

############################################################################################################
############################################################################################################
############################################################################################################

# GAMMA 
# robust GammaLSS:

# determination of data driven robustness constant "c" for quantiles "tau":
# use only outcome values with weight = 1
c_generate_Gamma <- function(outcome , tau = 0.05 ){
  log_like_c <- dgamma(outcome,scale = egamma(outcome)$parameters[2],shape=egamma(outcome)$parameters[1],log = TRUE)
  quant_c <- quantile(log_like_c,probs = tau)
  bl <- max(exp(-quant_c)-1,0.000001)       # lower boundary to secure an adequate value of C
  c_tau <- max(log(bl),0.25,na.rm = F) # lower boundary to secure an adequate value of C
  return(c_tau)
}

############################################################################################################

## we use almost the same parameterization as GA but with sigma = sqrt(1/sigma)
qGamma <- function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  ## require gamlss.dist
  if (!requireNamespace("gamlss.dist", quietly = TRUE))
    stop("Please install package 'gamlss.dist' for using qGamma.")
  gamlss.dist::qGA(p = p, mu = mu, sigma = sqrt(1/sigma), lower.tail = lower.tail, log.p = log.p)
}

############################################################################################################

# robust GammaMu
robust_GammaMu <- function (mu = NULL, sigma = NULL, stabilization,rob=rob) 
{
  loss <- function(sigma, y, f) {
    cal1 <- lgamma(sigma) + sigma * y * exp(-f) - sigma * log(y) -sigma * log(sigma) + sigma * f + log(y) 
    return(-rho(z=-cal1,rob=rob))  # changes within the loss function
  }
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, sigma = sigma))
  }
  ngradient <- function(y, f, w = 1) {
    n_cal1 <- lgamma(sigma) + sigma * y * exp(-f) - sigma * log(y) -sigma * log(sigma) + sigma * f + log(y)
    ngr <- rho_ab(z=-n_cal1,rob=rob) *  (sigma * y * exp(-f) - sigma) # changes within the gradient
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    return(ngr)
  }
  offset <- function(y, w) {
    if (!is.null(mu)) {
      RET <- log(mu)
    }
    else {
      RET <- log(egammaAlt(y)$parameters[1])
    }
    return(RET)
  }
  check_y <- function(y) {
    if (!is.numeric(y) || !is.null(dim(y))) 
      stop("response is not a numeric vector but ", 
           sQuote("GammaLSS()"))
    if (any(y < 0)) 
      stop("response is not positive but ", sQuote("GammaLSS()"))
    y
  }
  rho_ab <- function(z, rob=rob){          # log-logistic function embedded 
    exp(rob+z)/(exp(rob+z)+1)}             # log-logistic function embedded 
  rho <- function(z, rob=rob){             # log-logistic function embedded 
    log((1+exp(z+rob))/(1+exp(rob)))}      # log-logistic function embedded 
  
  mboost::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f), 
                 offset = offset, check_y = check_y, name = "robust Gamma distribution: mu(log link)")
}

############################################################################################################

# robust GammaSigma
robust_GammaSigma <- function (mu = NULL, sigma = NULL, stabilization,rob=rob) 
{
  loss <- function(mu, y, f) {
    cal2 <-   lgamma(exp(f)) + (exp(f) * y )/mu - exp(f) * log(y) - f * exp(f) + exp(f) * log(mu) + log(y)
    return(-rho(z=-cal2,rob=rob)) # changes within the loss function
  }
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu))
  }
  ngradient <- function(y, f, w = 1) {
    n_cal2 <-( lgamma(exp(f)) + (exp(f) * y )/mu - exp(f) * log(y) - f * exp(f) + exp(f) * log(mu) + log(y)) 
    ngr <- rho_ab(-n_cal2,rob=rob)*(-digamma(exp(f))*exp(f) + (f+1)*exp(f) - log(mu)*exp(f) + exp(f)*log(y) - (y*exp(f))/mu) # changes within the gradient
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    return(ngr)
  }
  offset <- function(y, w) {
    if (!is.null(sigma)) {
      RET <- log(sigma)
    }
    else {
      RET <-  log(1/egammaAlt(y)$parameters[2]^2)
    }
    return(RET)
  }
  check_y <- function(y) {
    if (!is.numeric(y) || !is.null(dim(y))) 
      stop("response is not a numeric vector but ", 
           sQuote("GammaLSS()"))
    if (any(y < 0)) 
      stop("response is not positive but ", sQuote("GammaLSS()"))
    y
  }
  rho_ab <- function(z, rob=rob){          # log-logistic function embedded 
    exp(rob+z)/(exp(rob+z)+1)}             # log-logistic function embedded 
  rho <- function(z, rob=rob){             # log-logistic function embedded 
    log((1+exp(z+rob))/(1+exp(rob)))}      # log-logistic function embedded 
  
  mboost::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f), 
                 offset = offset, check_y = check_y, name = "robust Gamma distribution: sigma(log link)")
}

############################################################################################################

robust_GammaLSS <- function (mu = NULL, sigma = NULL, stabilization = c("none","MAD", "L2"),rob=2) 
{
  if ((!is.null(sigma) && sigma <= 0)) 
    stop(sQuote("sigma"), " must be greater than zero")
  if ((!is.null(mu) && mu <= 0)) 
    stop(sQuote("mu"), " must be greater than zero")
  stabilization <- check_stabilization(stabilization)
  Families(mu = robust_GammaMu(mu = mu, sigma = sigma, stabilization = stabilization,rob=rob), 
           sigma = robust_GammaSigma(mu = mu, sigma = sigma, stabilization = stabilization,rob=rob), 
           qfun = qGamma, name = "robust Gamma")
}

############################################################################################################
############################################################################################################
