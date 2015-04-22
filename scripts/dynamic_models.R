
# Examples from dynamic crop models book
# Chapter 4

exp_model <- function(a, Y0, duration = 40, dt = 1){
	K = duration/dt + 1
	Y = rep(NA, K)
	time = 0
	k = 1

	while (time < duration){
		dY = a*Y[K]*dt
		Y[k + 1] + dY
		time = time + dt
		k = k + 1
	} #end while loop

	return(round(as.data.frame(cbind(time = seq(0, duration, by = dt), Y)), 10))

} #end exp_model

#main
a <- 0.1
Y0 <- 1.0 
duration <- 40
dt <- 2

x <- seq(0, duration, by = 0.1)
y <- Y0*exp(a*x)

#take a look at a plot
plot(x, y, type = "l", lwd = 3)
lines(exp_model(a, Y0, duration, dt), lty = 2, lwd = 3, col = "black")
legend("topleft",legend=c("Analytical Solution","Numerical
Solution"),lwd=c(2,2),lty=c(1,2),col=c("black","black"),
cex=0.75)

# improved euler method
exp_model2 <- function(a, Y0, duration = 40, dt = 1){
	K = duration/dt + 1  # number of time steps
	Y = rep(NA, K)      
	Y[1] = Y0
	time = 0
	k = 1

	while(time < duration){
		dY = a*Y[k]*dt # change in state variable time t
		YP = Y[k] + dY # predict value of state variable at time t+dt
		dYP = a*YP*dt  # estimate change of state variable at t+dt
		Y[k+1] = Y[k] + (1/2)*(dY+dYP) #value of state variable
		time = time + dt # update time
		k = k + 1
	} # end while loop

	return(round(as.data.frame(cbind(time = seq(0, duration, by = dt), Y)), 10))

} #end exp_model2


#main
a <- 0.1
Y0 <- 1.0 
duration <- 40
dt <- 2

x <- seq(0, duration, by = 0.1)
y <- Y0*exp(a*x)

# take a look at a plot
# approximates much closer to the real values

plot(x, y, type = "l", lwd = 3)
lines(exp_model2(a, Y0, duration, dt), lty = 2, lwd = 3, col = "black")
legend("topleft",legend=c("Analytical Solution","Numerical
Solution"),lwd=c(2,2),lty=c(1,2),col=c("black","black"),
cex=0.75)






