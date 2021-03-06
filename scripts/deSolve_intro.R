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

parameters <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.99, q = 2)
state <- c(w1 = 0.20, w2 = 0.2)
times <- seq(0, 200, by = 0.1)
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
parameters1 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 1, b = 25)
parameters2 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 1, b = 10)
parameters3 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.5, b = 25)
parameters4 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.2, b = 10)
parameters5 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.2, b = 25)
parameters6 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.9, b = 25)
state <- c(w1 = 1, w2 = 1)
times <- seq(0, 1000, by = 0.1)

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
equal <- ode(y = state, times = times, func = growth3, parms = parameters1)
plot(equal)

equalclose <- ode(y = state, times = times, func = growth3, parms = parameters2)
plot(equalclose)

stable <- ode(y = state, times = times, func = growth3, parms = parameters3)
plot(stable)

unstable <- ode(y = state, times = times, func = growth3, parms = parameters4)
plot(unstable)

stable2 <- ode(y = state, times = times, func = growth3, parms = parameters5)
plot(stable2)

stable3 <- ode(y = state, times = times, func = growth3, parms = parameters6)
plot(stable3)



head(out)
plot(out)
head(equal)
library(ggplot2)
library(reshape2)

head(equal)
equal <- as.data.frame(equal)
head(equal)
equalmelt <- melt(equal, id.vars = "time")
head(equalmelt)

equalplot <- ggplot(equalmelt) + 
  geom_line(aes(x = time, y = value, color = variable), size = 3) +
  scale_colour_manual(values=c("black", "red"),
  	name ="Plant", labels=c("Plant 1", "Plant 2")) +
  xlab("Time") + ylab("Biomass") +
  theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16),
           axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16))
equalplot

head(equalclose)
equalclose <- as.data.frame(equalclose)
head(equalclose)
equalclosemelt <- melt(equalclose, id.vars = "time")
head(equalclosemelt)

equalcloseplot <- ggplot(equalclosemelt) + 
  geom_line(aes(x = time, y = value, color = variable), size = 3) +
  scale_colour_manual(values=c("black", "red"),
  	name ="Plant", labels=c("Plant 1", "Plant 2")) +
  xlab("Time") + ylab("Biomass") +
  theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16),
           axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16))
equalcloseplot


head(stable)
stable <- as.data.frame(stable)
head(stable)
stablemelt <- melt(stable, id.vars = "time")
head(stablemelt)

stableplot <- ggplot(stablemelt) + 
  geom_line(aes(x = time, y = value, color = variable), size = 3) +
  scale_colour_manual(values=c("black", "red"),
  	name ="Plant", labels=c("Plant 1", "Plant 2")) +
  xlab("Time") + ylab("Biomass") +
  theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16),
           axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16))
stableplot

head(unstable)
unstable <- as.data.frame(unstable)
head(unstable)
unstablemelt <- melt(unstable, id.vars = "time")
head(unstablemelt)

unstableplot <- ggplot(unstablemelt) + 
  geom_line(aes(x = time, y = value, color = variable), size = 3) +
  scale_colour_manual(values=c("black", "red"),
  	name ="Plant", labels=c("Plant 1", "Plant 2")) +
  xlab("Time") + ylab("Biomass") +
  theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16),
           axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16))
unstableplot


head(stable2)
stable2 <- as.data.frame(stable2)
head(stable2)
stable2melt <- melt(stable2, id.vars = "time")
head(stable2melt)

stable2plot <- ggplot(stable2melt) + 
  geom_line(aes(x = time, y = value, color = variable), size = 3) +
  scale_colour_manual(values=c("black", "red"),
  	name ="Plant", labels=c("Plant 1", "Plant 2")) +
  xlab("Time") + ylab("Biomass") +
  theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16),
           axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16))
stable2plot

head(stable3)
stable3 <- as.data.frame(stable3)
head(stable3)
stable3melt <- melt(stable3, id.vars = "time")
head(stable3melt)

stable3plot <- ggplot(stable3melt) + 
  geom_line(aes(x = time, y = value, color = variable), size = 3) +
  scale_colour_manual(values=c("black", "red"),
  	name ="Plant", labels=c("Plant 1", "Plant 2")) +
  xlab("Time") + ylab("Biomass") +
  theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16),
           axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16))
stable3plot

#######
#September
#######

parameters <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 1, q = 2)
parameters_n <- c(k = 0.11, m = 0.01, g = 0.1, n = 0.9, C = 1, q = 2)


parameters1 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 1, b = 25)
parameters2 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 1, b = 10)
parameters3 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.5, b = 25)
parameters4 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.2, b = 10)
parameters5 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.2, b = 25)
parameters6 <- c(k = 0.11, m = 0.01, g = 0.1, n = 1, C = 0.9, b = 25)



state <- c(w1 = 0.20, w2 = 0.20)
times <- seq(0, 400, by = 1)
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


out <- ode(y = state, times = times, func = growth2, parms = parameters_n)
head(out)
plot(out)

growth5 <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    if ((1 - (b / (w1^(2/3) + w2^(2/3)))) >= 0)
      theta <- (1 - (b / (w1^(2/3) + w2^(2/3))))
    else 
      theta <- 0
    dw1 <- (w1/(w1 + C1*theta*w2))*(k*w1/(1 + g1*(w1^n1))) - m1*w1
    dw2 <- (w2/(w2 + (theta*w1)/C1))*(k*w2/(1 + g2*(w2^n2))) - m2*w2

    #return rate of change
    list(c(dw1, dw2))
  })
}

# compare this to HN non-crowded total biomass for these RILs
parameters_1 <- c(k = 0.11, m1 = 0.01, m2 = 0.01, g1 = 0.1, g2 = 0.1, n1 = 1, n2 = 1, C1 = 1, b = 50)
state <- c(w1 = 0.1, w2 = 0.1)
times <- seq(0, 1000, by = 1)
out_1 <- ode(y = state, times = times, func = growth5, parms = parameters_1)
head(out_1)
plot(out_1)

parameters_2 <- c(k = 0.11, m1 = 0.01, m2 = 0.01, g1 = 0.1, g2 = 0.1, n1 = 1, n2 = 1, C1 = 1, b = 20)
state <- c(w1 = 0.1, w2 = 0.1)
times <- seq(0, 1000, by = 1)
out_2 <- ode(y = state, times = times, func = growth5, parms = parameters_2)
head(out_2)
plot(out_2)

parameters_3 <- c(k = 0.11, m1 = 0.01, m2 = 0.01, g1 = 0.1, g2 = 0.1, n1 = 1, n2 = 1, C1 = 0.94, b = 20)
state <- c(w1 = 0.1, w2 = 0.1)
times <- seq(0, 1000, by = 1)
out_3 <- ode(y = state, times = times, func = growth5, parms = parameters_3)
head(out_3)
plot(out_3)
out_3[3]
parameters_4 <- c(k = 0.11, m1 = 0.01, m2 = 0.01, g1 = 0.1, g2 = 0.1, n1 = 1, n2 = 1, C1 = 0.94, b = 12)
state <- c(w1 = 0.1, w2 = 0.1)
times <- seq(0, 1000, by = 1)
out_4 <- ode(y = state, times = times, func = growth5, parms = parameters_4)
head(out_4)
plot(out_4)

library(reshape2)
library(ggplot2)
out_1[,c(1,3)]
out_all <- as.data.frame(cbind(out_1[,c(1,3)], out_2[,3], out_3[,3], out_4[,3]))
head(out_all)
colnames(out_all) <- paste(c("time", "UN_HN", "UN_LN", "CR_HN", "CR_LN"))
head(out_all)

out_allmelt <- melt(out_all, id.vars = "time")
head(out_allmelt)

out_allplot <- ggplot(out_allmelt) + 
  geom_line(aes(x = time, y = value, color = variable), size = 3) +
  scale_colour_manual(values=c("black", "red", "blue", "green"),
    name ="Treatment", labels=c("UN_HN", "UN_LN", "CR_HN", "CR_LN")) +
  xlab("Time") + ylab("Biomass Potential") +
  theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16),
           axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16))
out_allplot




