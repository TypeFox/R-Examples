################################################################################
# "Working with dynamic models for agriculture and the environment"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-05
# Be carefull, running all this script can be very long (up to 30 min)
library(sensitivity)
library(ZeBook)

# Definition of weed.factorseter nominal, min and max value
weed.factors = weed.define.param()
weed.factors

# Definition of decision variables for Soil, Crop et Herbicide on 10 years
# Case 1. Herbicid application each year.
weed.deci.1=data.frame(Soil=c(1,0,1,0,1,0,1,0,1,0), Crop=rep(1,10), Herb=c(1,1,1,1,1,1,1,1,1,1)) 
# Case 2. No herbicide application on year 3
weed.deci.2=data.frame(Soil=c(1,0,1,0,1,0,1,0,1,0), Crop=rep(1,10), Herb=c(1,1,0,1,1,1,1,1,1,1)) 

# Simulations of yield for case 1 and 2 with nominal values of weed.factorseters.
pred.1<-weed.model(weed.factors["nominal",],weed.deci.1)
pred.2<-weed.model(weed.factors["nominal",],weed.deci.2)

#png("weed.pred_1&2.png", width = 17, height = 8.5,units = "cm", res = 300,pointsize = 12)
par(mfrow=c(1,2),mar=c(4.1, 4.1, 1.1, 1.1))
plot(pred.1$year,pred.1$Yield, xlab="year", ylab="yield (t/ha)", type="b")
plot(pred.2$year,pred.2$Yield, xlab="year", ylab="yield (t/ha)", type="b")
#dev.off()

################################################################################
# Uncertainty Analysis on yield of year 3
N<-1000
duration = length(weed.deci.1$Soil)
Yield.mat.1<-matrix(NA, nrow=N,ncol=(duration+1),dimnames = list(NULL,paste("year",0:duration,sep="")))
Yield.mat.2<-Yield.mat.1

# Building the virtual experiment plan
weed.factors.mat <- matrix(ncol=length(weed.factors["nominal",]),nrow= N, dimnames = list(NULL,colnames(weed.factors)))
for (p in colnames(weed.factors.mat)) {weed.factors.mat[,p]<-runif(N,weed.factors["binf",p], weed.factors["bsup",p])}

# Plotting sample of parameter
#png("weed.sampling_Ymax_mh.png", width = 17, height = 8.5,units = "cm", res = 300,pointsize = 12)
par(mfrow=c(1,2),mar=c(4.1, 4.1, 1.1, 1.1))
for (k in c(100,1000)){
plot(weed.factors.mat[1:k,c("mh","Ymax")],xlim=weed.factors[c("binf","bsup"),"mh"],ylim=weed.factors[c("binf","bsup"),"Ymax"],sub=paste("sample of ",k) )
abline(v=weed.factors[c("binf","bsup"),"mh"],lty=2)
abline(h=weed.factors[c("binf","bsup"),"Ymax"],lty=2)
}
#dev.off()
# Simulate the plan
for (k in 1:N) {
  Yield.mat.1[k,]<-weed.model(weed.factors.mat[k,], weed.deci.1)$Yield
  Yield.mat.2[k,]<-weed.model(weed.factors.mat[k,], weed.deci.2)$Yield
}

#png("weed.uncertainty.png", width = 17, height = 17,units = "cm", res = 300,pointsize = 12)
par(mfrow=c(2,2),mar=c(4.1, 4.1, 1.1, 1.1))
hist(Yield.mat.1[,"year3"], main=" ", xlab="Yield with herbicide (t/ha)")
hist(Yield.mat.2[,"year3"], main=" ", xlab="Yield without herbicide (t/ha)")
hist(Yield.mat.1[,"year3"]-Yield.mat.2[,"year3"], main=" ", xlab="Yield loss (t/ha)")
hist(100*(Yield.mat.1[,"year3"]-Yield.mat.2[,"year3"])/Yield.mat.1[,"year3"], main=" ", xlab="relative Yield loss (%)")
#dev.off()
# Confidence intervals
quantile(Yield.mat.1[,"year3"],probs=c(0.01,0.05,0.95,0.99))
quantile(Yield.mat.2[,"year3"],probs=c(0.01,0.05,0.95,0.99))

#statistics on yield loss (with herb-without herb)
summary_diffYield1_2=summary(Yield.mat.1[,"year3"]-Yield.mat.2[,"year3"])
summary_diffYield1_2

quantile((Yield.mat.1[,"year3"]-Yield.mat.2[,"year3"]),probs=c(0.01,0.05,0.10,0.9,0.95,0.99) )

## Graphic : mean value and quantile estimation depending on size N of sample
N_mean_quant=data.frame()
for (k in c(1:100, seq(100,N,by=10))){
    V=(Yield.mat.1[,"year3"]-Yield.mat.2[,"year3"])[1:k]
    N_mean_quant=rbind(N_mean_quant,c("N"=k,"mean"=mean(V), quantile(V,probs=c(0.01,0.05,0.10,0.90,0.95,0.99))))
}
names(N_mean_quant)=c("N","mean","1%","5%","10%","90%","95%","99%")

#png("weed.uncertainty_SizeSample.png", width = 8.5, height = 8.5,units = "cm", res = 300,pointsize = 10)
par(mfrow=c(1,1),mar=c(4.1, 4.1, 1.1, 1.1))
plot(c(1,10000),c(min(V),max(V)), log="x", xlab="N", ylab="Yield loss (t/ha)", type="n", xaxp=c(1,10,1))
lines(N_mean_quant[,"N"],N_mean_quant[,"mean"],lwd=2)
lines(N_mean_quant[,"N"],N_mean_quant[,"1%"],lwd=1,lty=2)
lines(N_mean_quant[,"N"],N_mean_quant[,"99%"],lwd=1,lty=2)
lines(N_mean_quant[,"N"],N_mean_quant[,"5%"],lwd=1,lty=3)
lines(N_mean_quant[,"N"],N_mean_quant[,"95%"],lwd=1,lty=3)
#dev.off()

################################################################################
# Sensibility analysis on yield for year 3.
# MORRIS's method
# help(morris)
# 16 parameters
paraNames<-  colnames(weed.factors)
set.seed(2)

#Case 1. With herbicide
output.morris1 = morris(model =  weed.simule , factors=paraNames, 
r = 500, design = list(type = "oat", levels = 4 , grid.jump = 2), scale=T,
binf= weed.factors["binf",], bsup= weed.factors["bsup",], weed.deci=weed.deci.1)

#png("weed.morris_1.png", width = 8.5, height = 8.5,units = "cm", res = 300,pointsize = 10)
plot(output.morris1, xlim=c(-0.1,1.2),  main="with herbicide")
#dev.off()
#write.table(print(output.morris1),file="weed.morris1.csv",sep=";")
print(output.morris1)

#Case 2. Without herbicide on year 3
output.morris2 = morris(model =  weed.simule , factors=paraNames, 
r = 500, design = list(type = "oat", levels = 4 , grid.jump = 2), scale=T,
binf= weed.factors["binf",], bsup= weed.factors["bsup",], weed.deci=weed.deci.2)

#png("weed.morris_2.png", width = 8.5, height = 8.5,units = "cm", res = 300,pointsize = 10)
plot(output.morris2, xlim=c(-0.1,1.2), main = "without herbicide on year 3")
#dev.off()
#write.table(print(output.morris2),file="weed.morris2.csv",sep=";")
print(output.morris2)

################################################################################
# FAST's method
# help(fast99)
# 16 parameters
paraNames<-  colnames(weed.factors)
set.seed(2)

q.arg.fast = q.arg.fast.runif(weed.factors)

#Case 1. With herbicide
output.fast99_1 = fast99(model=weed.simule, factors=paraNames,
n = 1000, q = "qunif", q.arg =q.arg.fast , weed.deci=weed.deci.1)
#png("weed.fast99_1.png", width = 8.5, height = 8.5,units = "cm", res = 300,pointsize = 10)
plot(output.fast99_1)
#dev.off()
#write.table(print(output.fast99_1),file="weed.fast99_1.csv",sep=";")
print(output.fast99_1)

#Case 2. Without herbicide on year 3
output.fast99_2 = fast99(model=weed.simule, factors=paraNames,
n = 1000, q = "qunif", q.arg =q.arg.fast , weed.deci=weed.deci.2)
#png("weed.fast99_2.png", width = 8.5, height = 8.5,units = "cm", res = 300,pointsize = 10)
plot(output.fast99_2, cex.lab=0.75)
#dev.off()
#write.table(print(output.fast99_2),file="weed.fast99_2.csv",sep=";")
print(output.fast99_2)

################################################################################
# SOBOL2002's method
# help(sobol2002)
# 16 parameters
N=20000
set.seed(2)
X1=param.runif(weed.factors,N)
X2=param.runif(weed.factors,N)

# Running Sobol2002 : (N+2)*p simulation
#Case 1. With herbicide
system.time(output.sobol2002_1 <- sobol2002(model = weed.simule, X1, X2, nboot = 0, conf = 0.95,weed.deci=weed.deci.1))
#png("weed.sobol2002_1.png", width = 8.5, height = 8.5,units = "cm", res = 300, pointsize = 10)
par(mfrow=c(1,1))
plot(output.sobol2002_1)
#dev.off()
#write.table(print(output.sobol2002_1),file="weed.sobol2002_1.csv",sep=";")
print(output.sobol2002_1)

#Case 2. Without herbicide on year 3
system.time(output.sobol2002_2 <- sobol2002(model = weed.simule, X1, X2, nboot = 0, conf = 0.95,weed.deci=weed.deci.2))
#png("weed.sobol2002_2.png", width = 8.5, height = 8.5,units = "cm", res = 300, pointsize = 10)
par(mfrow=c(1,1))
plot(output.sobol2002_2)
#dev.off()
#write.table(print(output.sobol2002_2),file="weed.sobol2002_2.csv",sep=";")
print(output.sobol2002_2)

################################################################################
# ANOVA on the third year yield
# Function to run the weed model with modification of only the 5 parameters of interest returning the vector of the 3rd year yield for each parameter combination.
weed.simul5p<-function(param, para.mat, weed.deci){
	result <- c()
	for (i in 1:nrow(para.mat)){   
		param["mu"] <- para.mat$mu[i]
		param["beta.1"] <- para.mat$beta.1[i]
		param["beta.0"] <- para.mat$beta.0[i]
		param["mh"] <- para.mat$mh[i]
		param["Ymax"] <- para.mat$Ymax[i]	
		Y <- weed.model(param, weed.deci)	
		Y <- Y[Y$year ==3,"Yield"]
		result[i] <- Y
	}	
	return(result)
}

# Function
weed.annova=function(param,weed.factors,nlevel=2, npar=5, para.mat, weed.deci){
para.mat <- expand.grid(mu = seq.int(weed.factors['binf', 'mu'], weed.factors['bsup', 'mu'], length.out=nlevel), 
beta.1 = seq.int(weed.factors['binf', 'beta.1'], weed.factors['bsup', 'beta.1'], length.out=nlevel),
beta.0 = seq.int(weed.factors['binf', 'beta.0'], weed.factors['bsup', 'beta.0'], length.out=nlevel), 
mh = seq.int(weed.factors['binf', 'mh'], weed.factors['bsup', 'mh'], length.out=nlevel), 
Ymax = seq.int(weed.factors['binf', 'Ymax'], weed.factors['bsup', 'Ymax'], length.out=nlevel))

mu <- as.factor(para.mat$mu)
beta.1 <- as.factor(para.mat$beta.1)
beta.0 <- as.factor(para.mat$beta.0)
mh <- as.factor(para.mat$mh)
Ymax <- as.factor(para.mat$Ymax)
# for (p in colnames(para.mat)){para.mat[,p]=as.factor(para.mat[,p])}

Yield <- weed.simul5p(param, para.mat, weed.deci)
Fit <- summary(aov(Yield~mu*beta.1*beta.0*mh*Ymax))
print(Fit)
SumSq <- Fit[[1]][,2]
Total <- (nlevel^npar-1)*var(Yield)
Indices <- 100*SumSq/Total
print(Indices)
TabIndices <- cbind(Fit[[1]],Indices)[order(Indices, decreasing=T),]
print(TabIndices) 
return(TabIndices)    
}
weed.plot_anova=function(title){
# graph
#png(title, width = 17, height = 8.5,units = "cm", res = 300,pointsize = 12)
par(mfrow=c(1,2),mar=c(2.1, 10.1, 0.1, 0.1), adj=0,cex=0.5)
b=barplot(TabIndices.1[,"Indices"],  las=1, horiz=TRUE,xlim=c(0,105))
axis(side=2, at = b, line=9,labels = row.names(TabIndices.1), tick = FALSE,
lty = "solid",  lwd = 1,  hadj = 0, las=1)
b=barplot(TabIndices.2[,"Indices"],  las=1, horiz=TRUE,xlim=c(0,105))
axis(side=2, at = b, line=9,labels = row.names(TabIndices.2), tick = FALSE,
lty = "solid",  lwd = 1,  hadj = 0, las=1)
#dev.off()
}

######
nlevel=2
npar=5

# Scenario 1
TabIndices.1=weed.annova(weed.factors["nominal",],weed.factors,nlevel, npar, para.mat,weed.deci.1)
# Scenario 2
TabIndices.2=weed.annova(weed.factors["nominal",],weed.factors,nlevel, npar, para.mat, weed.deci.2)
weed.plot_anova("anova_5par_2Lev.png")
######
nlevel=3
npar=5
# Scenario 1
TabIndices.1=weed.annova(weed.factors["nominal",],weed.factors,nlevel, npar, para.mat, weed.deci.1)
# Scenario 2
TabIndices.2=weed.annova(weed.factors["nominal",],weed.factors,nlevel, npar, para.mat, weed.deci.2)
weed.plot_anova("anova_5par_3Lev.png")

######
nlevel=5
npar=5
# Scenario 1
TabIndices.1=weed.annova(weed.factors["nominal",],weed.factors,nlevel, npar, para.mat, weed.deci.1)
# Scenario 2
TabIndices.2=weed.annova(weed.factors["nominal",],weed.factors,nlevel, npar, para.mat, weed.deci.2)
weed.plot_anova("anova_5par_5Lev.png")

# End of file