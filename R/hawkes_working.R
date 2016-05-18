sayHello <- function(){
   print('hello')
}

loglik  <- function(params, arrivals) {
	n = length(arrivals)
	alpha_i <- params[1]
	beta_i	<- params[2]
	mu_i	<- params[3]
	#n <- 12
	term_1 <- -mu_i * arrivals[n]
	term_2 <- sum(alpha_i/beta_i * (exp( -beta_i * (arrivals[n] - arrivals)) - 1))
	Ai <- c(0,sapply(2:n, function(z) {
			sum(exp(-beta_i * (arrivals[z] - arrivals[1:(z-1)])))
		}))	
	term_3 <- sum(log(mu_i + alpha_i * Ai))
	return(-term_1- term_2 - term_3)	
}

loglikelihood <- function(params, arrivals, n){
	alpha_i <- params[1]
	beta_i <- params[2]
	mu_i <- params[3]

	term_1 <-(mu_i / beta_i) * (exp(-beta_i * arrivals[n]) - 1)
	term_2 <- sum(alpha_i/beta_i * (exp( -beta_i * (arrivals[n] - arrivals)) - 1))
	Ai <- c(0,sapply(2:n, function(z) {
			sum(exp(-beta_i * (arrivals[z] - arrivals[1:(z-1)])))
		}))	
	term_3 <- sum(log(mu_i * exp(-beta_i * arrivals[n]) + alpha_i * Ai))
	return(-term_1- term_2 - term_3)
}

neg.loglik <- function(params, data, opt=TRUE) {
  mu <- params[1]
  alpha <- params[2]
  beta <- params[3]
  t <- sort(data)
  r <- rep(0,length(t))
  for(i in 2:length(t)) {
    r[i] <- exp(-beta*(t[i]-t[i-1]))*(1+r[i-1])
  }
  loglik <- -tail(t,1)*mu
  loglik <- loglik+alpha/beta*sum(exp(-beta*(tail(t,1)-t))-1)
  loglik <- loglik+sum(log(mu+alpha*r))
  if(!opt) {
    return(list(negloglik=-loglik, mu=mu, alpha=alpha, beta=beta, t=t,
                r=r))
  }
  else {
    return(-loglik)
  }
}

# insert your values for (mu, alpha, beta) in par
# insert your times for data
#opt <- optim(par=c(1,2,0.1), fn=neg.loglik, data=c(0, 15, 36, 52, 104, 125, 154, 180, 232, 295, 717, 999))
#paste( c("alpha", "beta", "mu"), round(opt$par,2), sep="=")

#opt <- nlm(neg.loglik, c(1,2,0.1), hessian=TRUE, data=c(0, 15, 36, 52, 104, 125, 154, 180, 232, 295, 717, 999))
#print(opt)
#paste( c("alpha", "beta", "mu"), round(opt$par,2), sep="=")

#case1_solution1 <- optim(c(1,2,0.1), loglik, arrivals = c(0, 15, 36, 52, 104, 125, 154, 180, 232, 295, 717, 999), n = 12)
##paste( c("alpha", "beta", "mu"), round(case1_solution1$par,2), sep="=")

#case2_solution1 <- optim(c(0.001,0.0001,0.001), loglik, arrivals = c(0, 15, 36, 52, 104, 125, 154, 180, 232, 295, 717, 999), n = 12)
#paste( c("alpha", "beta", "mu"), round(case2_solution1$par,2), sep="=")

#case1_solution2 <- nlm(loglik, c(1, 2, 0.1), hessian = TRUE, arrivals = c(0, 19, 47, 218, 268, 346), n = 6)
#case1_solution2 <- nlm(loglik, c(1, 2, 0.1), hessian = TRUE, arrivals = c(0, 15, 36, 52, 104, 125, 154, 180, 232, 295, 717, 999), n = 12)
#case1_solution2 <- nlm(loglikelihood, c(1, 6, 0.1), hessian = TRUE, arrivals = c(0, 15, 36, 52, 104, 125, 154, 180, 232, 295, 717, 999, 1160, 1659, 1934, 3721, 4135, 4176, 9533, 9932), n=20)
#paste( c("alpha", "beta", "mu"), round(case1_solution2$par,2), sep="=")
#print(case1_solution2)
options(warn=-1)
myArgs <- commandArgs(trailingOnly = TRUE)
#myArgs[1]
myArrivals <- as.integer(unlist(strsplit(myArgs[1],",")))
myInitail <- as.double(unlist(strsplit(myArgs[2],",")))
#myVector
#print(length(myVector))
case1_solution3 <- nlm(loglik, myInitail, hessian = TRUE, arrivals = myArrivals)
#case1_solution3
#print(case1_solution3$estimate)
cat(case1_solution3$estimate)
#loglik(c())
#print(myArgs[1])
#
#args <- commandArgs(TRUE)
#eval(parse(text=args))
#a
#b