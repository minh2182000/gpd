#' Hypothesis test whether a censored data was generated from GPD distribution
#' 
#' \code{Htest.gpd} tests the Null Hypothesis that a right-censored data set follows GPD with the given parameters. Paremeters are estimated with mle.gpd if not given
#' 
#' @param time Time to failure or termination
#' @param censor Observation/censor indicator (1=observed, 0=right-censored), default to be all observations
#' @param k.hat optional, estimate of k, estimated by mle.gpd if either k.hat or a.hat is not given
#' @param a.hat optional, estimate of alpha, estimated by mle.gpd if either k.hat or a.hat is not given
#'
#' @return returns a list of the following elements:
#' @return \code{Null.Hypo}   States the Null Hypothesis
#' @return \code{test.stat}   Test statistic (Y-squared)
#' @return \code{df}          Degrees of freedom
#' @return \code{p.value}     p-value of the test
#' 
#' @details Methodology was described in:
#' Pham, M. H., Tsokos, C. P., Choi, B. (2018). Maximum Likelihood Estimation for the Generalized Pareto Distribution and Goodness-Of-Fit Test with Censored Data. Journal of Modern Applied Statistical Methods, 17(2).
#' 
#' @examples data(melanoma, package = "boot")
#' @examples status = ifelse(melanoma$status == 1, 1, 0)
#' @examples Htest.gpd(melanoma$time, status)
#' 
#' 
#' @import MASS
#' @import Matrix
#' @export

Htest.gpd = function(time, censor, k.hat = NULL, a.hat = NULL){
  x = time; delta = censor
  if (is.null(k.hat) | is.null(a.hat)){
    warning("parameters not given, estimated with mle.gpd")
    pars = mle.gpd(time, censor)
    k.hat = pars[1] ; a.hat = pars[2]
  }
  delta = delta[order(x)]; x = x[order(x)]
  if(length(x) != length(delta)){
      stop("vectors don't have equal lengths")
  }else{n = length(x)}
  
  # subinterval selection
  K = min(10, ceiling(n/10))
  E = vector("numeric", K)
  E[K] = sum(Clamda(x, k.hat,a.hat)); for(j in 1:(K-1)){E[j] = j/K*E[K]}
  
  # subinterval boundaries:
  a = vector("numeric", length = K + 1)
  a[1] = 0; a[K+1] = max(x)
  b = rep(NA, n)
  for (i in 1:n){
      b[i] = (n-i)*Clamda(x[i], k.hat, a.hat) + sum(Clamda(x[1:i], k.hat, a.hat))
  }
  for (j in 1:(max(1,K-1))){
      ii = match(TRUE, E[j]<=b[-1] & E[j]>=b[-n])
      a[j+1] = Clamda.inv((E[j] - sum(Clamda(x[1:(ii-1)],k.hat,a.hat)))/(n-ii+1), k.hat, a.hat)
  }
  
  # U vector-------------------------
  U = vector("numeric", length = K)
  for (j in 1:K){
      U[j] = sum(delta * (x > a[j]) * (x <= a[j+1]) )
  }
  # e vector --------------------------
  e = rep(E[K]/K, length = K)
  
  # Fisher information matrix ------------------------------
  I.hat = matrix(nrow = 2, ncol = 2)
  I.hat[1,1] = 1/n * sum( delta * (x/(a.hat - k.hat*x))^2 )
  I.hat[1,2] = I.hat[2,1] = 1/n * sum( delta * (-x/(a.hat - k.hat*x)^2) )
  I.hat[2,2] = 1/n * sum( delta / (a.hat - k.hat*x)^2 )
  
  # C-hat matrix ----------------------------------
  C.hat = matrix(nrow = 2, ncol = K)
  for (j in 1:K){
      C.hat[1,j] = 1/n * sum( (delta * x/(a.hat - k.hat*x)) * (x<=a[j+1] & x>a[j]) )
      C.hat[2,j] = 1/n * sum( (-delta / a.hat - k.hat*x) * (x<=a[j+1] & x>a[j]) )
  }
  
  # A-hat matrix ----------------------------------
  A.hat = diag(U, K, K)
  
  # Z-statistics vector
  Z = (U-e)/sqrt(n)
  
  # test statistics Y2:
  W = C.hat%*%solve(A.hat)%*%Z
  G = I.hat - C.hat%*%solve(A.hat)%*%t(C.hat)
  Y2 = sum((U - e)^2/U) + t(W)%*%ginv(G)%*%W
  
  # degree of freedom
  SIGMA = A.hat - t(C.hat)%*%solve(I.hat)%*%C.hat
  df = as.numeric(rankMatrix(ginv(SIGMA)))
  p.value = 1 - pchisq(Y2, df)
  
  # report result ---------------------------------
  NullH = paste("\n Null Hypothesis: the data follows GPD with k = ",k.hat, ", alpha = ",a.hat)
  return(list(Null.Hypo = NullH, test.stat = Y2, df = df, p.value = p.value))
}
