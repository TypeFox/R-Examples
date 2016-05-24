################################################################################
# "Working with dynamic models for agriculture"
# Francois Brun (ACTA), Juliette Adrian (ACTA)
# version : 2013-05-22
############################### MAIN PROGRAM ###################################
# Sensitivity Analysis
library(ZeBook)
library(sensitivity)

################################################################################
# 1) definition of the distribution of parameter values
#show parameter value : nominal, minimum and maximum
param<-lactation.define.param()
param["bsup",]=1.15*param["nominal",]
param["binf",]=0.85*param["nominal",]

# simulation characteristics
duration = 10*7-1
dt=0.1
#lactation.calf.simule(param,duration, dt)
################################################################################
# 2) Morris SA
# MORRIS's method, see help(morris)
# 2a) building the plan of simulation
nfac=dim(param)[2]
set.seed(123)
plan.morris <- morris(model=NULL , factors=row.names(t(param))[1:nfac],
r = 50, design = list(type = "oat", levels = 6 , grid.jump = 3), scale=T,
binf=as.vector(param["binf",]), bsup=as.vector(param["bsup",]))

# 2b)run the model for all the values of parameters
system.time(simX<-lactation.calf.simule(plan.morris$X,duration, dt))

# 2c) compute the indices for output of interest - here with ANOVA
choice_var="RM"
choice_week=3 # correspond to around the max of lactation
simX_choice = simX[simX[,"week"]==choice_week,choice_var]
output.morris<-tell(plan.morris,simX_choice)

mu <- apply(output.morris$ee, 2, mean)
mu.star <- apply(output.morris$ee, 2, function(x) mean(abs(x)))
sigma <- apply(output.morris$ee, 2, sd)

table.morris =data.frame(output.morris$factors,mu,mu.star,sigma)
table.morris[order(table.morris$mu.star,decreasing=TRUE),]

plot(output.morris,xlim=c(0,max(mu.star)+0.2))


################################################################################
# 3) complete factorial plan with ANOVA
# 3a) building the plan of simulation
choice_param=c("km","cu","mum","ks","kdiv","mm")
n=3 # up to  minutes
levels.param = apply(param[,choice_param],2, function(v) {seq.int(v['binf'], v['bsup'], length.out=n)} )
list.levels.param=list()
for (p in colnames(levels.param)){list.levels.param=c(list.levels.param,list(levels.param[,p]))}
names(list.levels.param)= colnames(levels.param)
param.mat <-expand.grid(list.levels.param)

param.mat=as.matrix(cbind(param.mat,t(param["nominal", -match(colnames(param.mat),names(param["nominal",]))]))[,names(param["nominal",])])

# 3b)run the model for all the values of parameters
system.time(simX<-lactation.calf.simule((param.mat),duration, dt))

# 3c) compute the indices for output of interest - here with ANOVA
choice_var="RM"
choice_week=3
simX_choice = simX[simX[,"week"]==choice_week,choice_var]

#for (p in colnames(param.mat)){param.mat[,p]=as.factor(param.mat[,p])}
formula =paste("simX_choice~",choice_param[1])
for (p in choice_param[-1]) formula= paste(formula, "*",p)
formula=eval(parse(text=formula))
Fit <- summary(aov(formula, data=cbind(simX_choice, as.data.frame(param.mat[,choice_param]))))
print(Fit)
SumSq <- Fit[[1]][,2]
Total <- (n^dim(levels.param)[2]-1)*var(simX_choice)
Indices <- 100*SumSq/Total
TabIndices <- cbind(Fit[[1]],Indices)[order(Indices, decreasing=T),]

# print only influent factors and interactions, superior to 0.1%
subset(TabIndices, Indices>0.1)


# visualization of the curve
plot(simX[,c("week","RM")])
# end of file