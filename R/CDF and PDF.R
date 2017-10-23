#PDF
pdf.gpd=function(x,k,alpha) {
  if (!k==0) {
    1/alpha*(1-k*x/alpha)^(1/k-1)
  } else {
    1/alpha*e^(-x/alpha)
  }
}
#CDF
cdf.gpd=function(x,k,alpha) {
  if (!k==0) {
    1 - (1 - k*x/alpha)^(1/k)
  } else {
    1 - exp(-x/alpha)
  }
}

# hazard
lamda = function(x, k, alpha){
    if (k == 0) {
        1/alpha
    } else{
        1/(alpha - k*x)
    }
}

#cumulative hazard
Clamda = function(x, k, alpha){
    if (k==0){
        x/alpha
    } else {
        1/k*log(alpha/(alpha - k*x))
    }
}

# inverse of cumulative hazard
Clamda.inv = function(h, k, alpha){
  if (k==0){
     alpha*h
  } else {
    alpha/k*(1 - exp(-k*h))
  }
}