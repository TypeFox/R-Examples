################################################################################
# Working with dynamic models for agriculture
# Francois Brun (ACTA), Juliette Adrian (ACTA)
# version : 2013-06-09
############################### MAIN PROGRAM ###################################
# Sensitivity Analysis
library(ZeBook)
library(sensitivity)

list_individuals=carcass_data$list_individuals
energy=carcass_data$energy
init_condition=carcass_data$init_condition
observation_dynamic=carcass_data$observation_dynamic
observation_slaughter=carcass_data$observation_slaughter

################################################################################
# 1) definition of the distribution of parameter values
#show parameter value : nominal, minimum and maximum
param.default=carcass.define.param()

# We fix the borne of the other parameters considered as constant for the sensitivity analysis
param=param.default
param["binf",colnames(param.default)[(is.na(param.default["binf",]))]]=0.85*param.default["nominal",colnames(param.default)[(is.na(param.default["binf",]))]]
param["bsup",colnames(param.default)[(is.na(param.default["bsup",]))]]=1.15*param.default["nominal",colnames(param.default)[(is.na(param.default["bsup",]))]]

# simulation characteristics
id=list_individuals[2,]
energyI=subset(energy,Individu==id)
energyInd=energyI$energie/7*10
day=energyI$time*7-6
energy_interpol<-approx(day,energyInd,xout=1:max(day),method = "constant")
duration=length(energy_interpol$y)
#carcass.EMI.simule(param, energy_interpol,init_condition[init_condition$Individu==id,"Pvi"],duration)

################################################################################
# 2) Morris SA
# MORRIS's method, see help(morris)
# 2a) building the plan of simulation
nfac=dim(param)[2]
set.seed(123)
plan.morris <- morris(model=NULL , factors=row.names(t(param))[1:nfac],
r = 100, design = list(type = "oat", levels = 6 , grid.jump = 3), scale=T,
binf=as.vector(param["binf",]), bsup=as.vector(param["bsup",]))

# 2b)run the model for all the values of parameters
system.time(simX<-carcass.EMI.simule(plan.morris$X, energy_interpol,init_condition[init_condition$Individu==id,"Pvi"],duration))

# 2c) compute the indices for output of interest - here with MORRIS
choice_var="PV"
choice_time=100
simX_choice = simX[simX[,"time"]==choice_time,choice_var]
output.morris<-tell(plan.morris,simX_choice)

mu <- apply(output.morris$ee, 2, mean)
mu.star <- apply(output.morris$ee, 2, function(x) mean(abs(x)))
sigma <- apply(output.morris$ee, 2, sd)

table.morris =data.frame(output.morris$factors,mu,mu.star,sigma)
table.morris[order(table.morris$mu.star,decreasing=TRUE),]

plot(output.morris,xlim=c(0,max(mu.star)*1.05))
# end of file