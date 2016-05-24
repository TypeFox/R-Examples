################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-22
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
library(ZeBook)
library(sensitivity)
list_n_sy=unique(maize.data_EuropeEU$sy)

################################################################################
# 1) definition of the distribution of parameter values
#show parameter value : nominal, minimum and maximum
param<-maize.define.param()

################################################################################
# 2a) complete factorial plan on one particular site-year
n=3
levels.param = apply(param,2, function(v) {seq.int(v['binf'], v['bsup'], length.out=n)} )
list.levels.param=list()
for (p in colnames(levels.param)){list.levels.param=c(list.levels.param,list(levels.param[,p]))}
names(list.levels.param)= colnames(levels.param)
param.mat <-expand.grid(list.levels.param) 

# on a single site-year
sy="18-2006"
weather = maize.weather(working.year=strsplit(sy,"-")[[1]][2], working.site=strsplit(sy,"-")[[1]][1],weather_all=weather_EuropeEU)
    
# run the model for all the values of parameters
system.time(simX240<-maize.simule240(param.mat,  weather=weather, sdate=100, ldate=250))

for (p in colnames(param.mat)){param.mat[,p]=as.factor(param.mat[,p])}
Fit <- summary(aov(simX240~Tbase*RUE*K*alpha*LAImax*TTM*TTL, data=cbind(simX240, param.mat)))
print(Fit)
SumSq <- Fit[[1]][,2]
Total <- (n^dim(levels.param)[2]-1)*var(simX240)
Indices <- 100*SumSq/Total
TabIndices <- cbind(Fit[[1]],Indices)[order(Indices, decreasing=T),]

# print only influent factors and interactions, superior to 0.1%
subset(TabIndices, Indices>0.1)

# 2b) complete factorial plan mean of B240 for the n sites-years  (seem to be the more adequate for us)
n=2 # up to 20 minutes for n=2, and up to 6 hours for n=3 !
levels.param = apply(param,2, function(v) {seq.int(v['binf'], v['bsup'], length.out=n)} )
list.levels.param=list()
for (p in colnames(levels.param)){list.levels.param=c(list.levels.param,list(levels.param[,p]))}
names(list.levels.param)= colnames(levels.param)
param.mat <-expand.grid(list.levels.param)

# run the model for all the values of parameters
system.time(simX240<-maize.simule_multisy240(param.mat,list_n_sy , sdate=100, ldate=250, all=FALSE))

for (p in colnames(param.mat)){param.mat[,p]=as.factor(param.mat[,p])}
Fit <- summary(aov(simX240~Tbase*RUE*K*alpha*LAImax*TTM*TTL, data=cbind(simX240, param.mat)))
print(Fit)
SumSq <- Fit[[1]][,2]
Total <- (n^dim(levels.param)[2]-1)*var(simX240)
Indices <- 100*SumSq/Total
TabIndices <- cbind(Fit[[1]],Indices)[order(Indices, decreasing=T),]

# print only influent factors and interactions, superior to 0.1%
subset(TabIndices, Indices>0.1)

################################################################################
# 3a) Morris SA on one particular site-year
sy="18-2006"
weather = maize.weather(working.year=strsplit(sy,"-")[[1]][2], working.site=strsplit(sy,"-")[[1]][1],weather_all=weather_EuropeEU)

# MORRIS's method, see help(morris), with 7 parameters
nfac=7
set.seed(123)
system.time(output.morris <- morris(model=maize.simule240 , factors=row.names(t(param))[1:nfac],
r = 20, design = list(type = "oat", levels = 6 , grid.jump = 3), scale=T,
binf=as.vector(param["binf",]), bsup=as.vector(param["bsup",]), weather=weather,
sdate=100, ldate=250))

plot(output.morris, xlim=c(0,1))

table.morris = print(output.morris)
table.morris[order(table.morris$mu.star,decreasing=TRUE),]

# 3b) Morris SA on mean of B240 for the n sites-years (seem to be the more adequate for us)
# MORRIS's method, see help(morris), with 7 parameters
nfac=7

set.seed(123)
system.time(output.morris <- morris(model=maize.simule_multisy240 , factors=row.names(t(param))[1:nfac],
r = 20, design = list(type = "oat", levels = 6 , grid.jump = 3), scale=T,
binf=as.vector(param["binf",]), bsup=as.vector(param["bsup",]), liste_sy=list_n_sy,
sdate=100, ldate=250))

plot(output.morris, xlim=c(0,1.1))

table.morris = print(output.morris)
table.morris[order(table.morris$mu.star,decreasing=TRUE),]

# end of file