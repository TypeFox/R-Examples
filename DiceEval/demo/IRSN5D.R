library(gam)
library(mda)
library(polspline)

#----------------------------------
# learning data
data(dataIRSN5D)
straussX <- dataIRSN5D[,1:5]
straussY <- dataIRSN5D[,6]
# test data
data(testIRSN5D)
data324 <- testIRSN5D
design324 <- as.matrix(data324[,1:5])
Y324  <- as.matrix(data324[,6])

# Simple auxilary function to transform the experimental design in [-1,1]^5
reduc <- function(V){
	return(2*(V-0.5))
}
Xtest <- apply(design324,2,reduc)
X <- apply(straussX,2,reduc)

# dimensionless response
mY <- mean(straussY)
sY <- sd(straussY)
Y <- as.matrix((straussY-mY)/sY)
Ytest <- (Y324-mY)/sY

# GRAPHIC: study of the response 
layout(matrix(c(1,2),1,2))
h<-hist(Y,freq=FALSE,20,xlab = "keff_DOE",xlim = c(-3,5),main = "")
lines(density(Y),col="red",lwd=2)
h<-hist(Ytest,freq=FALSE,40,xlab = "keff_TEST",xlim = c(-3,5),main = "")
lines(density(Ytest),col="red",lwd=2)
  
# BOXPLOT: training design (X) and test design (Xtest)
layout(matrix(c(1,2),1,2))
boxplot(as.data.frame(X),main = "Design of Experiments")
boxplot(as.data.frame(Xtest),main = "Test design")

# GRAPHIC: projection of the Strauss design on the principal axis
layout(matrix(1:6,3,2))
for (i in 1:5){
	hist(X[,i],freq=FALSE,20,xlab = paste("X",i,sep=""),main = "", cex.lab = 1.5)
	lines(density(X[,i]),col= "red", lwd=2)
}

# GRAPHIC: projection of the test points on the axis
layout(matrix(1:6,3,2))
for (i in 1:5){
	hist(Xtest[,i],freq=FALSE,20,xlab = paste("X",i,sep=""),main = "", cex.lab = 1.5)
}

# GRAPHIC: output vs inputs
layout(matrix(1:6,3,2))
for (i in 1:5){
	plot(Xtest[,i],Ytest,xlab = paste("X",i,sep=""), cex.lab = 1.5,ylab = "keff")
	points(X[,i],Y,col="red",pch = 19,bg = "red")
}

#------------------------------------------------------------------------
# Construction of the metamodel and quality criteria
#------------------------------------------------------------------------
# Fitting a model form training data
modLm <- modelFit(X,Y,type = "Linear",formula=Y~.)
summary(modLm$model)

# Prediction on test data
Ytest_lm <- modelPredict(modLm,Xtest)

# Quality criteria (R2)
cat("R2Lm =",  R2(Y,modLm$model$fitted.values),"\n",
    "Q2_50_Lm =", crossValidation(modLm,K=50)$Q2,"\n",
    "Q2_10_Lm =", crossValidation(modLm,K=10)$Q2,"\n",
    "R2Lmtest =",  R2(Ytest,Ytest_lm),"\n",sep="") 
	
# Other criteria of performance 
cat("RMA_DOE =",  RMA(Y,modelPredict(modLm,X))$max.value,"\n",
    "RMA_TEST =",  RMA(Ytest,modelPredict(modLm,Xtest))$max.value,"\n",
 	"MAE_DOE =",  MAE(Y,modelPredict(modLm,X)),"\n",
    "MAE_TEST =",  MAE(Ytest,modelPredict(modLm,Xtest)),"\n",
	"RMSE_DOE =",  RMSE(Y,modelPredict(modLm,X)),"\n",
    "RMSE_TEST =",  RMSE(Ytest,modelPredict(modLm,Xtest)),"\n")

# Focus on the Cross-Validation procedure
par(mfrow=c(1,1))
testCrossValidation(modLm,Kfold=c(2,5,10,20,30,40,dim(modLm$data$X)[1]))

#------------------------------------------------------------------------
# Study of the resiuals of an Additive Model
#------------------------------------------------------------------------
modAm <- modelFit(X,Y,type = "Additive",formula=formulaAm(X,Y))
summary(modAm)
residualsStudy(modAm)

#------------------------------------------------------------------------
# Focus on the PolyMARS procedure
#------------------------------------------------------------------------
Crit <- penaltyPolyMARS(X,Y,test=data.frame(Xtest,Ytest),graphic=TRUE)

#------------------------------------------------------------------------
# Focus on the stepwise model
#------------------------------------------------------------------------
out <- stepEvolution(X,Y,Y~.^2,P=c(1,2,5,10,20,30))

#------------------------------------------------------------------------
# Focus on the Kriging model (see DiceKriging package for more details)
#------------------------------------------------------------------------
mKm <- modelFit(X,Y, type="Kriging",covtype="powexp",control=list(trace=FALSE))
K <- 10
out   <- crossValidation(mKm, K)
par(mfrow=c(2,2))
plot(c(0,1:K),c(mKm$model@covariance@range.val[1],out$theta[,1]),xlab='',ylab='Theta1')
plot(c(0,1:K),c(mKm$model@covariance@range.val[2],out$theta[,2]),xlab='',ylab='Theta2')
plot(c(0,1:K),c(mKm$model@covariance@range.val[1],out$shape[,1]),xlab='',ylab='p1')
plot(c(0,1:K),c(mKm$model@covariance@range.val[2],out$shape[,2]),xlab='',ylab='p2')
par(mfrow=c(1,1))

#------------------------------------------------------------------------
# Comparison of metamodels
#------------------------------------------------------------------------
crit  <- modelComparison(X,Y,type=c("Linear","Additive","MARS","PolyMARS","Kriging"),
	test=data.frame(Xtest,Ytest),degree=2,gcv=4,
	formula=c(Y~.^2,formulaAm(X,Y),Y~.))
