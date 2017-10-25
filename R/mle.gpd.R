#' Parameter estimation for GPD
#' 
#' \code{mle.gpd} estimates parameters for Generalized Pareto	Distribution with right-censored data using Maximum Likelihood
#' 
#' @param Time Time to failure or termination
#' @param censor Observation/censor indicator (1=observed, 0=right-censored), default to be all observations
#' 
#' @return returns a vector of size 2. The first element is the estimate of k, the second is the estimate of alpha
#' @details 
#' Generalized Pareto distribution:			
#' CDF F(x,k,a) = 1 - (1-kx/a)^(1/k) if k!=0	or 1 - e^(-x/a) if k=0						
#' @details Algorithm was described in:
#' Pham, M. H., Tsokos, C. P., Choi, B. (2018). Maximum Likelihood Estimation for the Generalized Pareto Distribution and Goodness-Of-Fit Test with Censored Data. Journal of Modern Applied Statistical Methods, 17(2).
#'
#' @examples data(melanoma, package = "boot")
#' @examples status = ifelse(melanoma$status == 1, 1, 0)
#' @examples mle.gpd(melanoma$time, status)
#' @export

mle.gpd=function(Time, censor = rep(1, Times = length(Time))){
  pdf.gpd = supp$pdf.gpd; n_r = supp$n_r; lamda = supp$lamda; L = supp$L; k = supp$k; h_prime = supp$h_prime; h = supp$h; Clamda.inv = supp$Clamda.inv; Clamda = supp$Clamda; cdf.gpd = supp$cdf.gpd; bisection = supp$bisection; alpha = supp$alpha
  environment(alpha) = environment(bisection) = environment(cdf.gpd) = environment(Clamda) = environment(Clamda.inv) = environment(h) = environment(h_prime) = environment(k) = environment(L) = environment(lamda) = environment(n_r) = environment(pdf.gpd) = environment()
  x=Time
  delta=censor
	n=length(Time)
	xbar = mean(Time)
	options(digits=22)
		# re-order x and delta for better expression of h(theta) (see the math)---------- ----
				order = order(delta,decreasing=T)
				delta = delta[order]
				x = x[order]
				r = max(which(delta==1))

		# Theta estimation ---------------------------------------------------- ----
  			 #1
				epsilon = 10^-6/xbar
		       #2
				thetaL = 2*(min(x)-xbar)/(min(x))^2
				thetaU = 1/max(x)-epsilon
			 #3
				xbar_r = sum(x[1:r])/r
				lim_ddh = sum(1/n*x[1:n]^2) - 2*xbar*xbar_r
		       #4
				if (lim_ddh>0){
					epsilon.save=epsilon	   #save epsilon whne it is changed
				#A (thetaL,-epsilon)
					theta0_1 = n_r(thetaL)
					if (theta0_1=="diverge"){
						if (epsilon<10^-12){epsilon=10^-12}
									   #change epsilon to avoid bug
						theta0_1=bisection(thetaL,-epsilon)	#run bisection when NR fails
									}
					
				#B (epsilon,thetaU)
				 	theta0_2 = n_r(thetaU)
					if (theta0_2=="diverge"){
						if (epsilon<10^-12){epsilon=10^-12}
									   #change epsilon to avoid bug
						theta0_2=bisection(epsilon,thetaU)
									}
					epsilon = epsilon.save     #change epsilon again for further calculation
				#solution list:
					solutions=c(theta0_1,theta0_2)
						  }
		      #5
				if (lim_ddh<0){
				#A (thetaL,-epsilon)
					theta0_1=n_r(thetaL)
			
					if (!theta0_1=="diverge"){
				#B#C		#search for second 0
						if (h_prime(theta0_1)>0) {
							theta0_2=bisection(theta0_1,-epsilon)
				#D
										 } else {
							theta0_2=bisection(thetaL,theta0_1)
											  }
									         solution_L=c(theta0_1,theta0_2)
									 } else {solution_L=NULL}
				#E (epsilon,thetaU)
					theta0_3=n_r(thetaU)
							
					if (!theta0_3=="diverge"){
				#F#G		#search for second 0

						if (h_prime(theta0_3)>0) {
							theta0_4=bisection(epsilon,theta0_3)

										 } else {
				#H
							theta0_4=bisection(theta0_3,thetaU)
			
											  }
										   solution_U=c(theta0_3,theta0_4)
									 } else {solution_U=NULL}
					#solution list:
					solutions=c(solution_L,solution_U)
						  }
					#delete "diverge" in solutions list
					solutions=solutions[which(!solutions=="diverge")]
					solutions=as.numeric(solutions)
		if (length(solutions)==0){stop("no MLE estimation")} else {	
  			#6
				k_list=1 ; alpha_list=1
				for (i in 1:length(solutions)){k_list[i]=k(solutions[i])}
				for (i in 1:length(solutions)){alpha_list[i]=alpha(solutions[i])}
				L_list=1:length(solutions)
				for (i in 1:length(solutions)){
					L_list[i]=L(k_list[i],alpha_list[i])
						}
		      #7
				max_position=which(L_list==max(L_list))[1]
				local_max=c(k_list[max_position],alpha_list[max_position])
			    #8
				if(L(k=local_max[1],alpha=local_max[2])>-n*log(max(x))) {mle.parameters=local_max} else {mle.parameters=c(1,max(x))}
			names(mle.parameters)=c("k","alpha")
			options(digits=7)
			result = local_max
			names(result) = c("k", "alpha")
												  }
		return(result)
										}						
