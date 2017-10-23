#--------- Log likelihood function ------------------------------------------------- ----
    L = function(k,alpha){
      with(parent.frame(),
      sum(-delta*log(alpha)+(1/k-delta)*log(1-k/alpha*x))
      )
    }
#-------- k estimation from theta0 ------------------------------------------------- ----
    # note: (theta=k/alpha)
    # r = censoring weight, defined in mle function's environment
    k = function(theta0){
      if (!exists('r', envir = parent.frame()))){stop("this function is being run outside of the mother's function 'mle'")}
      r = get('r', envir = parent.frame()))
      with(parent.frame(),
        return(-1/r*sum(log(1-theta0*x)))
      )
    }
#------- alpha estimation from theta0 ---------------------------------------------- ----
    alpha = function(theta0){
      k(theta0)/theta0
    }
#------- function h(theta)--------------------------------------------------------- ----
    h = function (theta){
      with(parent.frame(),
        1/n*sum(log(1-theta*x))*1/r*sum((1-theta*x[1:r])^-1) + 1/n*sum((1-theta*x)^-1)-1
      )
    }
#------- function h'(theta) ------------------------------------------------------- ----
    h_prime = function(theta){
      with(parent.frame(),
        1/n*sum(-x/(1-theta*x))*1/r*sum((1-theta*x[1:r])^-1) +
        1/n*sum(log(1-theta*x))*1/r*sum(x[1:r]*(1-theta*x[1:r])^-2) +
        1/n*sum(x*(1-theta*x)^-2)
      )
    }
