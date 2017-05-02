
#' @title Integration via the tanh-sinh procedure
#' @description Evaluates an integral with the tanh-sinh procedure.
#'
#' @param f
#' @param lower,upper the bounds of integration; infinite allowed
#' @param ... additional arguments passed to \code{f}
#' @param error error estimate requested
#'
#' @return
#' @export
#'
#' @examples
#' # approximate pi/4
#' f <- function(x) sqrt(1-x^2)
#' result <- tanhsinh(f, 0, 1, error=1e-8)
#' result$value - pi/4
#' # Gauss integral
#' f <- function(x) exp(-x*x)
#' result <- tanhsinh(f, 0, Inf)
#' result$value - sqrt(pi)/2
#' result <- tanhsinh(f, -Inf, Inf)
#' result$value - sqrt(pi)
#' # singularity at 1
#' f <- function(x) exp(x)/sqrt(1-x)
#' result <- tanhsinh(f, 0, 1)
#' result$value - exp(1)*sqrt(pi)*(2*pnorm(sqrt(2))-1)
#' # Owen T-function at h=a=0.5
#' f <- function(x) exp (-0.25*(1+x*x)/2) /(1+x*x)
#' result <- tanhsinh(f, 0, 0.5)
#' result$value - 0.40519384189197524877
#' # Owen Q-function
#' f <- function(x, nu, delta, t){
#'   pnorm(t*x/sqrt(nu) - delta) *
#'     exp((nu-1)*log(x) - x*x/2 - ((nu/2)-1) * log(2) - lgamma(nu/2))
#' }
#' result <- tanhsinh(f, 0, 5, nu=10, delta=2, t=3)
#' result$value - 0.77383547873740988537
#' # Bessel K function
#' f <- function(x) x*besselK(x,0)^4
#' result <- tanhsinh(f, 0, Inf)
#' result$value - 7*1.2020569031595942854/8
tanhsinh <- function(f, lower, upper, ..., error=1e-16){
  if(lower > upper){
    stop("`lower`>`upper` is not allowed")
  }
  if(lower==-Inf && upper<Inf){
    return(
      tanhsinh(function(x) f(-x), lower=-upper, upper=Inf, ..., error=error)
    )
  }
  if(upper==Inf){
    if(lower>-1){
      g <- function(y) f(y, ...)
      ff <- function(x){
        g(x/(1-x))/(1-x)/(1-x)
      }
      upper <- 1
      lower <- lower/(1+lower)
    }else{
      g <- function(y) f(y, ...)
      ff <- function(x){
        (1+tan(x)*tan(x))*g(tan(x))
      }
      upper <- pi/2
      lower <- atan(lower)
    }
  }else{
    ff <- function(x) f(x, ...)
  }
  c <- (upper-lower)/2; d <- (lower+upper)/2
  error <- error/c
  integral <- ff(d)*halfpi
  abcissas <- Abcissas[[1L]]
  weights <- Weights[[1L]]
  for(i in 1L:3L){
    integral <- integral +
      weights[i]*(ff(c*abcissas[i] + d) + ff(-c*abcissas[i] + d))
  }
  currentDelta <- Inf
  h <- 1
  out <- NULL
  for(l in 2L:9L){
    abcissas <- Abcissas[[l]]
    weights <- Weights[[l]]
    h <- h/2
    newContribution <- 0
    for(i in 1L:layerSizes[l]){
      newContribution <- newContribution +
        weights[i]*(ff(c*abcissas[i] + d) + ff(-c*abcissas[i] + d));
    }
    newContribution <- newContribution*h
    previousDelta <- currentDelta
    currentDelta <- abs(integral/2 - newContribution)
    integral <- integral/2 + newContribution
    if(currentDelta == 0){
      out <- rbind(out,
                   c(c*integral, 0, nEvaluations[l]))
      break
    }
    r <- log(currentDelta)/log(previousDelta)
    if(r > 1.9 && r < 2.1){
      errorEstimate <- currentDelta*currentDelta
      cat("BIDULE")
    }else{
      errorEstimate <- currentDelta
    }
    out <- rbind(out,
                 c(c*integral, c*errorEstimate, nEvaluations[l]))
    if(errorEstimate < 0.1*error) break
  }
  # errorEstimate <- errorEstimate*c
  # evaluations <- nEvaluations[l]
  # out <- c*integral
  # attr(out, "error") <- errorEstimate
  # attr(out, "evaluations") <- evaluations
  out <- setNames(as.data.frame(out),
                  c("value", "error", "evaluations"))
  return(out)
}
