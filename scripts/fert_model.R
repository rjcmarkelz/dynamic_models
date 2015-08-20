
library(deSolve)

# fertilzer model
# Thornley and Johnson, 1990 page 467
# something funky going on with LAI need to debug

parameters <- c(Mc = 30, Mn = 62, J = 5*10^6,
    h = 43,000, Pc = 1, sigN = 3000,
    Kc = 0.05, Kn = 0.005, dr = 0.2, rhos = 1500, mu = 150, fc = 0.45, fn = 0.03, 
	f1 = 0.7, Num = 25, Ep = 2.5, Y = 0.75, alphan = 0.5, gammash = 0.1, 
	gammar = 0.01, Bn = 3*10^6)
state <- c(Wc = 0.015, Wn = 0.004, Wsh = 0.2, Wr = 0.2, L = 0.8, Ns = 33.33*10^-6)

times <- seq(0, 10, by = .1)
growth3 <- function(t, state, parameters){
	with(as.list(c(state, parameters)),{
		Wg <- Wsh + Wr
		C  <- Wc/Wg
		N  <- Wn/Wg
		Ws <- (Mc/12)*Wc + (Mn/14)*Wn
		fsh <- Wsh/Wg
		fr  <- Wr/Wg
		Wsht <- Wsh + fsh*Ws
		Wrt  <- Wr + fr*Ws
		Io <- J/h
		P  <- 43200*(12/44)*Pc
		Un   <- (sigN*Wr*Ns)/(1 + (Kc/C)*(1 + N/Kn))
		lamsh <- (fr*N/(N + fn)) /( (fsh*C)/(C + fc) + (fr*N)/(N + fn) )
		lamr <- (fsh*C/(C + fc)) /( (fsh*C)/(C + fc) + (fr*N)/(N + fn) )
		Gsh <- mu*C*N*lamsh*Wsh
		Gr <- mu*C*N*lamr*Wr
		Nu <- Num*(1 - Ep*C)

		# differential equations
		dWc <- P - fc*(Gsh + Gr)/Y - alphan*Un
		dWn <- Un - fn*(Gsh + Gr)
		dWsh <- Gsh - gammash*Wsh
		dWr  <- Gr - gammar*Wr
		dL   <- Nu*f1*Gsh - gammash*L # something wrong with the LAI 
		dNs  <- Bn - (Un/dr*rhos)

		#return rate of change
		list(c(dWc, dWn, dWsh, dWr, dL, dNs))
	})
}
out <- ode(y = state, times = times, func = growth3, parms = parameters)


head(out)
plot(out)

growth3 <- function(t, state, parameters){
	with(as.list(c(state, parameters)),{
		Pc <- (Pm/k)*log((alpha*k*Io + Pm(1-m))/(alpha*k*Io*exp(-k*L)))
		Wg <- Wsh + Wr
		C  <- Wc/Wg
		N  <- Wn/Wg
		Ws <- (Mc/12)*Wc + (Mn/14)*Wn
		fsh <- Wsh/Wg
		fr  <- Wr/Wg
		Wsht <- Wsh + fsh*Ws
		Wrt  <- Wr + fr*Ws
		Io <- J/h
		P  <- 43200*(12/44)*Pc
		Un   <- (sigN*Wr*Ns)/(1 + (Kc/C)*(1 + N/Kn))
		lamsh <- (fr*N/(N + fn)) /( (fsh*C)/(C + fc) + (fr*N)/(N + fn) )
		lamr <- (fsh*C/(C + fc)) /( (fsh*C)/(C + fc) + (fr*N)/(N + fn) )
		Gsh <- mu*C*N*lamsh*Wsh
		Gr <- mu*C*N*lamr*Wr
		Nu <- Num*(1 - Ep*C)

		# differential equations
		dWc <- P - fc*(Gsh + Gr)/Y - alphan*Un
		dWn <- Un - fn*(Gsh + Gr)
		dWsh <- Gsh - gammash*Wsh
		dWr  <- Gr - gammar*Wr
		dL   <- Nu*f1*Gsh - gammash*L # something wrong with the LAI 
		dNs  <- Bn - (Un/dr*rhos)

		#return rate of change
		list(c(dWc, dWn, dWsh, dWr, dL, dNs))
	})
}
out <- ode(y = state, times = times, func = growth3, parms = parameters)






?ode