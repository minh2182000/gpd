#' Histogram with fit of GPD density function
#'
#' creates a plot of the data set's histogram and the fitted GPD density function
#' @import  graphics
#' @param time Time to failure or termination
#' @param censor Observation/censor indicator (1=observed, 0=right-censored), default to be all observations
#' @param k.hat optional, estimate of k, estimated by mle.gpd if either k.hat or a.hat is not given
#' @param a.hat optional, estimate of alpha, estimated by mle.gpd if either k.hat or a.hat is not given
#' @param show.all (default to False) whether to draw a historam with all data points or only observations. Default to show only observed points.
#' @param ... further arguments passed to hist()
#'
#' @examples data(melanoma, package = "boot")
#' @examples status = ifelse(melanoma$status == 1, 1, 0)
#' @examples hist.gpd(melanoma$time, status)
#'
#' @export


hist.gpd=function(time, censor = rep(1, times = length(time)), k.hat = NULL, a.hat = NULL, show.all = FALSE, ...){
  if (is.null(k.hat) | is.null(a.hat)){
    warning("parameters not given, estimated with mle.gpd")
    pars = mle.gpd(time, censor)
    k.hat = pars[1] ; a.hat = pars[2]
  }

  if (!show.all){
    obs.time = time[censor == 1]
    hist(obs.time, prob=T, xlim=range(time), ...)
  } else{
    hist(time, prob=T, xlim=range(time), ...)
  }

  curve(pdf.gpd(x,k=k.hat,alpha=a.hat),col=2,add=T)
}
