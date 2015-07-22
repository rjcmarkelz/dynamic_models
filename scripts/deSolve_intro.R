library(deSolve)

#test model
parameters <- c(a = -8/3, b = -10, c = 28)
state <- c(X = 1, Y = 1, Z = 1)
times <- seq(0, 100, by = 0.01)

Lorenz <- function(t, state, parameters){
	with(as.list(c(state, parameters)),{
		# rate of change
		dX <- a*X + Y*Z
		dY <- b* (Y - Z)
		dZ <- -X*Y + c*Y - Z

		#return rate of change
		list(c(dX, dY, dZ))
	})
}

out <- ode(y = state, times = times, func = Lorenz, parms = parameters)
head(out)




#test model
parameters <- c(k = 0.11, m = 0.01, g = 0.1, n = 1)
state <- c(w = 0.0001)
times <- seq(0, 1000, by = 0.01)
growth <- function(t, state, parameters){
	with(as.list(c(state, parameters)),{
		# rate of change
		dw <- ((k*w/(1 + g*w^n))) - m*w

		#return rate of change
		list(c(dw))
	})
}
out <- ode(y = state, times = times, func = growth, parms = parameters)

head(out)
plot(out)


# two plant growth model

parameters <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.7, q = 2)
state <- c(w1 = 0.20, w2 = 0.2)
times <- seq(0, 1000, by = 0.01)
growth2 <- function(t, state, parameters){
	with(as.list(c(state, parameters)),{
		# rate of change
		dw1 <- (w1^q/(w1^q + C*w2^q))*(k*(w1+w2)/(1 + g*(w1+w2)^n)) - m*w1
		dw2 <- (C*(w2^q/(w1^q + C*w2^q)))*(k*(w1+w2)/(1 + g*(w1+w2)^n)) - m*w2

		#return rate of change
		list(c(dw1, dw2))
	})
}


out <- ode(y = state, times = times, func = growth2, parms = parameters)

head(out)
plot(out)

# two plant growth model nested
w <- w1 + w2
parameters <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.2, b = 25)
state <- c(w1 = 1, w2 = 1)
times <- seq(0, 1000, by = 0.01)

growth3 <- function(t, state, parameters){
	with(as.list(c(state, parameters)),{
		if ((1 - (b / (w1^(2/3) + w2^(2/3)))) >= 0)
			theta <- (1 - (b / (w1^(2/3) + w2^(2/3))))
		else 
			theta <- 0
		dw1 <- (w1/(w1 + C*theta*w2))*(k*w1/(1 + g*(w1^n))) - m*w1
		dw2 <- (w2/(w2 + (theta*w1)/C))*(k*w2/(1 + g*(w2^n))) - m*w2

		#return rate of change
		list(c(dw1, dw2))
	})
}
out <- ode(y = state, times = times, func = growth3, parms = parameters)

head(out)
plot(out)