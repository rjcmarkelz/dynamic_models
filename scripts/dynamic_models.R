
# Examples from dynamic crop models book, some code directly from book
# some code modified by Cody Markelz
# Chapter 4

################
# euler integration method
################
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

################
# improved euler method
################

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


################
# multiple state variables Euler
################

# this has a lot of state variables to keep track of
pop_age_model <- function(rb = 3.5, mE = 0.017, rE = 0.172, m1 = 0.060,
                          r12 = 0.217, m2 = 0.032, r23 = 0.313, m3 = 0.022, 
                          r34 = 0.222, m4 = 0.020, r4P = 0.135, mP = 0.020, 
                          rPA = 0.099, mA = 0.027, iA = 0, duration = 100,
                          dt = 1){
	# create matrix of state variables, one per column, nrows for length of sim
	V = matrix(NA, ncol = 7, nrow = duration/dt+1, 
		dimnames = list(NULL, c("E", "L1", "L2", "L3", "L4", "P", "A")))

	# initialize state variables
	V[1, ] <- c(5, 0, 0, 0, 0, 0, 0)

	# simulation loop
	for(k in 1:(duration/dt)){
		# rates of change of state variables
		dE = (rb*V[k, "A"] - rE*V[k, "E"] - mE*V[k, "E"])*dt       # egg stage
		dL1 = (rE*V[k, "E"] - r12*V[k, "L1"] - m1*V[k, "L1"])*dt   # larvae1
		dL2 = (r12*V[k, "L1"] - r23*V[k, "L2"] - m2*V[k, "L2"])*dt # larvae2
		dL3 = (r23*V[k, "L2"] - r34*V[k, "L3"] - m3*V[k, "L3"])*dt # larvae3
		dL4 = (r34*V[k, "L3"] - r4P*V[k, "L4"] - m4*V[k, "L4"])*dt # larvae4
		dP = (r4P*V[k, "L4"] - rPA*V[k, "P"] - mP*V[k, "P"])*dt    # pupae
		dA = (rPA*V[k, "P"] - mA*V[k, "A"])*dt                     # adult

		dV = c(dE, dL1, dL2, dL3, dL4, dP, dA)
		# update variables
		V[k + 1, ] <- V[k, ] + dV
	} 
	# end simulation loop

	return(round(as.data.frame(cbind(time = (1: (duration/dt + 1))*dt - dt, V)),
		10))
}
# end pop_age_model


# non-neat book code for plotting, do not feel like re-writting
sim <- pop_age_model(rb=3.5,mE=0.017,rE=0.172,m1=0.060,
r12=0.217,m2=0.032,r23=0.313,m3=0.022,r34=0.222,m4=0.020,
r4P=0.135,mP=0.020,rPA=0.099,mA=0.027,iA=0,duration=40,dt=0.01)

sim[sim$time==0,]
sim[sim$time==10,]
sim[sim$time==40,]

graph.param=data.frame("V"=c("E","L1","L2","L3","L4","P","A"),
"lty"=c(2,2,2,2,3,4,1), "lwd"=c(2,1,1,1,1,2,3))
plot(c(0,max(sim$time)), c(0,max(sim[,-1])),type="n",xlab="time
(day)",ylab="population density")
null=sapply(c("E","L1","L2","L3","L4","P","A"),function(v)
lines(sim$time,sim[,v],
lty=graph.param[graph.param$V==v,"lty"],lwd=graph.param[graph.param$V==v,"lwd"]))
legend("topright", legend=graph.param$V, lty=graph.param$lty,
lwd = graph.param$lwd, cex=0.75)



################
# multiple state variables Runge-Kutta
################
pop_age_model_rk <- function(rb = 3.5, mE = 0.017, rE = 0.172, m1 = 0.060,
                          r12 = 0.217, m2 = 0.032, r23 = 0.313, m3 = 0.022, 
                          r34 = 0.222, m4 = 0.020, r4P = 0.135, mP = 0.020, 
                          rPA = 0.099, mA = 0.027, iA = 0, duration,
                          dt, method){

	# state variables starting values
	E0 = 5
	L10 = 0 
	L20 = 0
	L30 = 0
	L40 = 0
	P0  = 0
	A0  = 0

	# define ODEs
	pred_prey_ode <- function(time, state, pars){
		with(as.list(c(state, pars)), {
			dE  = (rb*A - rE*E - mE*E)*dt
			dL1 = (rE*E - r12*L1 - m1*L1)*dt
			dL2 = (r12*L1 - r23*L2- m2*L2)*dt
			dL3 = (r23*L1 - r34*L3- m3*L3)*dt
			dL4 = (r34*L3 - r4P*L4 - m4*L4)*dt
			dP  = (r4P*L4 - rPA*P - mP*P)*dt
			dA  = (rPA*P - mA*A+iA)*dt
            
            # return list
            return(list(c(dE, dL1, dL2, dL3, dL4, dP, dA)))
		})
	}

	# set the initial states in the ode and a call to ode
	sim = ode(y = c(E = E0, L1 = L10, L2 = L20, L3 = L30, L4 = L40, P = P0, A = A0),
	        times = seq(0, duration, by = dt), func = pred_prey_ode, 
	        parms = c(rb, mE, rE, m1, r12, m2, r23, m3, r34, m4, r4P, mP, rPA, mA, iA),
	        method = rkMethod(method))

	return(as.data.frame(sim))

}

#load deSolve
library(deSolve)

sim2 <- pop_age_model_rk(rb = 3.5, mE = 0.017, rE = 0.172, m1 = 0.060,
                         r12 = 0.217, m2 = 0.032, r23 = 0.313, m3 = 0.022, 
                          r34 = 0.222, m4 = 0.020, r4P = 0.135, mP = 0.020, 
                          rPA = 0.099, mA = 0.027, iA = 0, duration = 40,
                          dt = 1, method = "rk4")
sim2


sim2[sim2$time==0,]
sim2[sim2$time==10,]
sim2[sim2$time==40,]








#compressed plotting code from book
graph.param = data.frame("V"=c("E","L1","L2","L3","L4","P","A"),
"lty"=c(2,2,2,2,3,4,1), "lwd"=c(2,1,1,1,1,2,3))
plot(c(0,max(sim2$time)), c(0,max(sim2[,-1])),type="n",xlab="time
(day)",ylab="population density")
null=sapply(c("E","L1","L2","L3","L4","P","A"),function(v)
lines(sim2$time,sim2[,v],
lty=graph.param[graph.param$V==v,"lty"],lwd=graph.param[graph.param$V==v,"lwd"]))
legend("topright", legend=graph.param$V, lty=graph.param$lty,
lwd=graph.param$lwd, cex=0.75)



# simple maize model from chapter 1 using difference equations

maize <- function (Tbase, RUE, K, alpha, LAImax, TTM, TTL, weather, sdate, 
	ldate) {

    TT <- rep(NA, ldate)
    B <- rep(NA, ldate)
    LAI <- rep(NA, ldate)
    CumInt <- rep(NA, ldate)
    
    TT[sdate] <- 0
    B[sdate] <- 1
    LAI[sdate] <- 0.01
    CumInt[sdate] = 0

    for (day in sdate:(ldate - 1)) {

        dTT <- max((weather$Tmin[day] + weather$Tmax[day])/2 - Tbase, 0)

        if (TT[day] <= TTM) {
            dB <- RUE * (1 - exp(-K * LAI[day])) * weather$I[day]
        } else {
            dB <- 0
        }

        if (TT[day] <= TTL) {
            dLAI <- alpha * dTT * LAI[day] * max(LAImax - LAI[day], 0)
        } else {
            dLAI <- 0
        }

        TT[day + 1] = TT[day] + dTT
        B[day + 1]  B[day] + dB
        LAI[day + 1] <- LAI[day] + dLAI
        CumInt[day + 1] = CumInt[day] + weather$I[day]* (1 - exp(-K * LAI[day]))
    }

    return(data.frame(day = sdate:ldate, TT = TT[sdate:ldate], 
        LAI = LAI[sdate:ldate], B = B[sdate:ldate],
        CumInt = CumInt[sdate:ldate]))
}




	V = matrix(NA, ncol = 7, nrow = 50/2+1, 
		dimnames = list(NULL, c("E", "L1", "L2", "L3", "L4", "P", "A")))






