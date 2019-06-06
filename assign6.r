
#---------------bite.mod()------------------
#Returns the predicted intake rates for the bite mass model

bite.mod <- function (Rmax, hbite, Sobs) {
    Ipred <- (Rmax * Sobs)/((Rmax * hbite) + Sobs)
    return(Ipred)
}

#---------------dens.mod()--------------------
#Returns the predicted intake rates for the density model

dens.mod <- function (Vmax, hdens, a, Dobs, Sobs) {
    Ipred <- (Vmax*a*sqrt(Dobs)*Sobs)/(1 + (a*hdens*Vmax*sqrt(Dobs)))
    return(Ipred)
}

#---------------bite.nll()--------------------
#Returns the negative log likelihood for the bite model with error from a normal distribution.

bite.nll <- function (p, Sobs, Iobs) {
    Ipred <- bite.mod(Rmax=p[1], hbite=p[2], Sobs=Sobs)
    nll <- -sum(dnorm(Iobs, mean=Ipred, sd=p[3], log=TRUE))
    return(nll)
}

#---------------dens.nll()--------------------
#Returns the negative log likelihood for the density model with error from a normal distribution.
dens.nll <- function (p, a, Sobs, Dobs, Iobs) {
    Ipred <-dens.mod(Vmax=p[1], hdens=p[2], a=a, Dobs=Dobs, Sobs=Sobs)
    nll <- -sum(dnorm(Iobs, mean=Ipred, sd=p[3], log=TRUE))
    return(nll)
}

#-------MAIN PROGRAM----------------------------------------------------------------------------

#-------Visualize the data by plotting

setwd("C:/Users/Cheng/Desktop/QuantBIO")
plant.data <- read.csv("student02_ass6.csv",header=TRUE)   

Sobs <- plant.data$S
Dobs <- plant.data$D
Iobs <- plant.data$I

plot(Sobs, Dobs)  #shows combinations of plant densities and sizes used in response surface design 
plot(Sobs, Iobs)  #As plant size increases, intake rate increases; looks nonlinear
plot(Dobs, Iobs)  #Plant density seems to have little effect on intake rate; many data points around 0 density

#-------Use gridsearch()

source('C:/Users/Cheng/Downloads/gridsearch.r')

#Narrow down parameter ranges for Bite model:

Rmax_range <- seq(0, 80, length.out=100)  
hbite_range <- exp(seq(log(1e-4), log(0.5), length.out=100))
sdbite_range <- seq(0, 60, length.out=100)

bitepvecs <- list(Rmax_range,hbite_range,sdbite_range)

bitegrid <- gridsearch(bitepvecs,bite.nll,mon=1,Sobs,Iobs)
bitegrid$par       #the best set of parameters found (corresponds to minimum nll)
bitegrid$value     #the nll corresponding to best parameters
head(bitegrid$profile)   #A matrix of the parameter combinations tried, plus corresponding nll values

par(mfrow=c(1,3))
plot(bitegrid$profile[,1], bitegrid$profile[,4], xlim=c(30,60), ylim=c(400,440), xlab="Rmax (rate of processing)", ylab="nll",main="NLL profile for Rmax") # 35<Rmax<40
plot(bitegrid$profile[,2], bitegrid$profile[,4], xlim=c(0,0.04), ylim=c(400,440), xlab="h (biting time)", ylab="nll", main="NLL profile for h") # 0<h<0.02
plot(bitegrid$profile[,3], bitegrid$profile[,4], xlim=c(4,15), ylim=c(410,450), xlab="sd (std dev of I)", ylab="nll", main="NLL profile for sd") # 8<sdbite<9


#Narrow down parameter ranges for Density model:

Vmax_range <- seq(0, 1000, length.out=100)
hdens_range <- exp(seq(log(1e-4), 0.5, length.out=100))
sddens_range <- seq(0, 50, length.out=100)

denspvecs <- list(Vmax_range,hdens_range,sddens_range)

densgrid <- gridsearch(denspvecs,dens.nll,mon=1,a=a,Sobs=Sobs,Dobs=Dobs,Iobs=Iobs)
densgrid$par     
densgrid$value     
head(densgrid$profile)   

#Negative log likelihood profile for Vmax values:
par(mfrow=c(1,3))
plot(densgrid$profile[,1],densgrid$profile[,4],xlim=c(0,1000),ylim=c(450,550),xlab="Vmax",ylab="nll",main="NLL profile for Vmax") # 0<Vmax<infinity
plot(densgrid$profile[,2],densgrid$profile[,4],xlim=c(0.05,0.2),ylim=c(450,600),xlab="hdens",ylab="nll",main="NLL profile for hdens") # 0.05<hdens<0.1
plot(densgrid$profile[,3],densgrid$profile[,4],xlim=c(10,25),ylim=c(450,600),xlab="sddens",ylab="nll",main="NLL profile for sddens") # 14<sddens<18


#-------Fit the bite model to the data

#Initialize parameters:
Rmax <- 37.9
hbite <- 0.011
sdbite <- 8.58

#Optimization: find the parameters from minimum negative log likelihood:
bite.par <- c(Rmax, hbite, sdbite)
bite.fit <- optim (bite.par, bite.nll, Sobs=Sobs, Iobs=Iobs)
bite.fit

#-------Fit the density model to the data

#Initialize parameters:
Vmax <- 1000
hdens <- 0.01
a <- 1
sddens <- 16.1

#Optimization: find the parameters from minimum log likelihood:
dens.par <- c(Vmax, hdens, sddens)
dens.fit <- optim (dens.par, dens.nll, a=a, Sobs=Sobs, Dobs=Dobs, Iobs=Iobs)
dens.fit

#-------Plot the data against predicted values

par(mfrow=c(1,2))

#Calculate fitted model dynamics for best parameter values in bite model:
fittedIbite <- bite.mod(Rmax=bite.fit$par[1], hbite=bite.fit$par[2], Sobs=Sobs)

#Plot fitted bite model over the data:
plot(Sobs, Iobs, xlab="Bite mass (S)", ylab="Intake rate (I)", main="Error fit of bite model")
points(sort(Sobs), sort(fittedIbite), col="red", type="l")
segments(Sobs, fittedIbite, Sobs, Iobs, col="green")

#Plot observed vs. fitted values:
plot(fittedIbite, Iobs, xlab="Predicted intake rate", ylab="Observed intake rate", main="Bite model")
abline(0,1,col="red")


#Calculate fitted model dynamics for best parameter values in density model:
fittedIdens <- dens.mod(Vmax=dens.fit$par[1], hdens=dens.fit$par[2], a=a, Dobs=Dobs, Sobs=Sobs)

#Plot fitted density model over the data:
plot(Sobs, Iobs, xlab="Bite mass (S)", ylab="Intake rate (I)", main="Error fit of density model for S values")
points(sort(Sobs), sort(fittedIdens), col="blue", type="l")
segments(Sobs, fittedIdens, Sobs, Iobs, col="green")

#Plot observed vs. fitted values:
plot(fittedIdens, Iobs, xlab="Predicted intake rate", ylab="Observed intake rate", main="Density model")
abline(0,1,col="red")

#----------CONFIDENCE INTERVALS-----------------------------------

source("C:/Users/Cheng/Downloads/confintprof.r") 

par(mfrow=c(1,3))

biteconfint <- confintprof( pars=bite.fit$par,
                            profpar=1,
                            nllfn=bite.nll,
                            plim=c(32,45),
                            conf_lev=0.95,
                            n=50,
                            pname="Rmax",
                            Sobs=Sobs,
                            Iobs=Iobs )

biteconfint <- confintprof( pars=bite.fit$par,
                            profpar=2,
                            nllfn=bite.nll,
                            plim=c(0.005,0.018),
                            conf_lev=0.95,
                            n=50,
                            pname="h",
                            Sobs=Sobs,
                            Iobs=Iobs )

biteconfint <- confintprof( pars=bite.fit$par,
                            profpar=3,
                            nllfn=bite.nll,
                            plim=c(7,11),
                            conf_lev=0.95,
                            n=50,
                            pname="sd",
                            Sobs=Sobs,
                            Iobs=Iobs )

# densconfint <- confintprof( pars=dens.fit$par,
#                             profpar=1,
#                             nllfn=dens.nll,
#                             plim=c(1,1000),
#                             conf_lev=0.95,
#                             n=50,
#                             pname="Vmax",
#                             a=a,
#                             Sobs=Sobs,
#                             Dobs=Dobs,
#                             Iobs=Iobs )



#--------------------AKAIKE-------------------------------------------

#---------AIC()

AIC <- function(nll, k) {
    aic <- (2 * nll) + 2 * k
    return(aic)
}

biteAIC <- AIC(nll=bite.fit$value, k=3)
biteAIC
densAIC <- AIC(nll=dens.fit$value, k=3)
densAIC

aic_values <- c(biteAIC, densAIC)

#---------deltaAIC

deltaAICbite <- biteAIC - biteAIC
deltaAICdens <- densAIC - biteAIC

delta_values <-c(deltaAICbite, deltaAICdens)

#---------AICc()

AICc <- function(aic, k, n) {
    aicc <- aic + ((2 * k * (k + 1))/(n - k - 1))
    return(aicc)
}

biteAICc <- AICc(aic=biteAIC, k=3, (n=length(Iobs)))
biteAICc

densAICc <- AICc(aic=densAIC, k=3, (n=length(Iobs)))
densAICc

AICc_values <- c(biteAICc, densAICc)

#---------Ak_weights()

Ak_weights <- function(aic_values) {
    best_aic <- min(aic_values)
    deltas <- aic_values - best_aic
    weights <- rep(NA, length(aic_values))
    for (i in (1:length(aic_values))) {
        weights[i] <- exp(-.5 * deltas[i])/(sum(exp(-.5 * deltas)))
    }
    return(weights)
}

weights <- Ak_weights(aic_values) 
weights

list(aic_values, delta_values, AICc_values, weights)
table(aic_values, delta_values, AICc_values, weights)


#-------------BITE MODEL DIAGNOSTICS----------------

par(mfrow=c(1,1))

#-------Histogram of residuals:

residuals <- fittedIbite - Iobs
hist(residuals, ylim= c(0,0.05), xlab = "Residuals", main = "Histogram of residuals",freq=FALSE)
rr <- seq(min(residuals), max(residuals), length.out=100)
lines(rr, dnorm(rr, mean=0, sd=sd(residuals)), col="red")
box()

#-------Residuals vs predicted (fitted):

plot(fittedIbite, residuals, ylab = "Residuals", xlab = "Fitted values", main = "Residuals vs fitted")
abline(h=0,col="red")

#-------Q-Q plot:

qqnorm(residuals)
qqline(residuals, col="red")


