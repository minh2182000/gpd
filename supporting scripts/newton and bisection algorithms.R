# Newton-Raphson Procedure for h(theta): n_r(initial) -------------------- ----
n_r=function(initial){
  s0=initial
  repeat{
    s1=s0 - h(s0)/h_prime(s0)
    if (s1<thetaL|s1>thetaU){result="diverge"; break} else {
      if (abs(h(s1))>=abs(h(s0))) {result="diverge"; break} else {
        if (abs(h(s1))<epsilon|s1<thetaL|s1>thetaU){result=s1; break} else {s0=s1}
      }
    }
  }
  return(result)
}						 
# bisection method for h(theta)----------------------------------------------- ----
bisection=function(b1,b2){
  ## condition: b2>b1 ##
  repeat {
    b0=(b1+b2)/2
    if (h(b1)*h(b0)<0) {b1=b1; b2=b0} else {b1=b0; b2=b2}
    if (b2-b1<epsilon) {break}
    #final result
  }
  return((b1+b2)/2)
}
