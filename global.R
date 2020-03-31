library(shiny)
library(ggplot2)
library(RColorBrewer)
library(tidyr) # for "spread"
library(shinyjs) # for reset button
library(shinythemes) # for "spacelab" theme
library(shinycssloaders) # for "calculating" spinners
library(gridExtra) # for arranging sub-plots
library(Cairo) # for better graphics resolution
options(shiny.usecairo=T)
library(markdown) # for "details" and "about" pages
library(tidyverse) # for handling m 
library(viridis) # for contour plot

options(shiny.trace = FALSE)

### initialise 
# x-axis values used for plotting
tau.test <- seq(from=0, to=13, by = 0.1)

### global functions

# probability symptomatic
prob.symp <- function(tau, incper.meanlog, incper.sdlog) {
    plnorm(tau, meanlog = incper.meanlog, sdlog = incper.sdlog)
}

# probability asymptomatic
prob.asymp <- function(tau, incper.meanlog, incper.sdlog) {
    1 - plnorm(tau, meanlog = incper.meanlog, sdlog = incper.sdlog)
}

# serial interval
serint <- function(x, serint.shape, serint.scale) {
    dweibull(x = x, shape = serint.shape, scale = serint.scale)
}

# environmental contribution
p.el <- function(l, env.decay.rate, env.constant.duration, env.infectiousness.type) {
  if (env.infectiousness.type == "constant") {
    return(l < env.constant.duration)
  } else if (env.infectiousness.type == "exp.decay") {
    return(exp(-env.decay.rate * l))
  }
}

model.gen.f <- function(tau, incper.meanlog, incper.sdlog, P.a, xp) {
  s.of.tau <- prob.symp(tau = tau, 
                        incper.meanlog = incper.meanlog,
                        incper.sdlog = incper.sdlog)  
  (1 - P.a) * (s.of.tau + xp * (1 - s.of.tau))
}

model.gen.beta.s.div.by.RSorP <- function(tau, incper.meanlog, incper.sdlog,
                                          serint.scale, serint.shape,
                                          P.a, xp) {
  s.of.tau <- prob.symp(tau = tau, incper.meanlog = incper.meanlog,
                        incper.sdlog = incper.sdlog)
  
  serint(x = tau, serint.scale = serint.scale,
         serint.shape = serint.shape) /
    ((1 - P.a) * (s.of.tau + xp * (1 - s.of.tau)))
}

model.gen.beta.presym.tot <- function(tau, incper.meanlog, incper.sdlog, serint.scale, serint.shape, P.a, xp, RSorP) {
    xp * (1 - P.a) * RSorP * model.gen.beta.s.div.by.RSorP(tau = tau,
                                                           incper.meanlog = incper.meanlog,
                                                           incper.sdlog = incper.sdlog,
                                                           serint.scale = serint.scale,
                                                           serint.shape = serint.shape,
                                                           P.a = P.a,
                                                           xp = xp) *
    prob.asymp(tau = tau, 
               incper.meanlog = incper.meanlog, 
               incper.sdlog = incper.sdlog)
}

model.gen.beta.sym.tot <- function(tau, incper.meanlog, incper.sdlog, serint.scale, serint.shape, P.a, xp, RSorP) {
    (1 - P.a) * RSorP * model.gen.beta.s.div.by.RSorP(tau = tau,
                                                      incper.meanlog = incper.meanlog,
                                                      incper.sdlog = incper.sdlog,
                                                      serint.scale = serint.scale,
                                                      serint.shape = serint.shape,
                                                      P.a = P.a,
                                                      xp = xp) *
      prob.symp(tau = tau, 
                incper.meanlog = incper.meanlog, 
                incper.sdlog = incper.sdlog)
}


### Fixed values for the "interventions" tab for discretization of the integral
nmax<-100
ndiscr<-10
v<-c(1:nmax)/ndiscr
Y<-rep(1,nmax)/nmax
n<-9
Yerror<-0.1
miniter<-1
maxiter<-3