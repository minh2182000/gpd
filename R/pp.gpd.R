

#------------------- P-P plot ----------------------------------------------------- -
pp.gpd=function(time){
  path = getwd()
  source(paste(path, '/supporting scripts/CDF and PDF.R', sep = ""),
         local = TRUE)
  
  if (is.null(k.hat) | is.null(a.hat)){
    warning("parameters not given, estimated with mle.gpd")
    pars = mle.gpd(time, censor)
    k.hat = pars[1] ; a.hat = pars[2]
  }
  
  
  censor = rep(1, times = length(time))
  p=mle.gpd(time,censor)
  if (p[1]=="no MLE estimation"){cat("no MLE estimation")} else{
    time1=sort(time,decreasing=F)
    n=length(time1)
    
    theoretical=cdf.gpd(time1,k=p[1],alpha=p[2])
    empirical=1:n
    for (i in 1:n){
      empirical[i] = (i-0.5)/n
    }
    plot(x=theoretical,y=empirical,col="blue")
    curve((x),add=T,col=2,lwd=2)
  }
}
