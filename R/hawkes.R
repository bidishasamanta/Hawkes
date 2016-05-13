sayHello <- function(){
   print('hello')
}

loglik  <- function(params, arrivals, n) {
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
case1_solution1 <- optim(c(1,2,0.1), loglik, arrivals = c(0, 15, 36, 52, 104, 125, 154, 180, 232, 295, 717, 999), n = 12)
paste( c("alpha", "beta", "mu"), round(case1_solution1$par,2), sep="=")

case2_solution1 <- optim(c(0.001,0.0001,0.001), loglik, arrivals = c(0, 15, 36, 52, 104, 125, 154, 180, 232, 295, 717, 999), n = 12)
paste( c("alpha", "beta", "mu"), round(case2_solution1$par,2), sep="=")

case1_solution2 <- nlm(loglik, c(1, 2, 0.1), hessian = TRUE, arrivals = c(0, 19, 47, 218, 268, 346), n = 6)
paste( c("alpha", "beta", "mu"), round(case1_solution2$par,2), sep="=")
