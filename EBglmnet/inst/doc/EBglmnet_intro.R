## ---- eval=FALSE---------------------------------------------------------
#  install.packages("EBglmnet", repos = "http://cran.us.r-project.org")

## ------------------------------------------------------------------------
rm(list = ls())
library(EBglmnet)

## ------------------------------------------------------------------------
varNames = colnames(state.x77);
varNames
y= state.x77[,"Life Exp"]
xNames = c("Population","Income","Illiteracy", "Murder","HS Grad","Frost","Area")
x = state.x77[,xNames]

## ------------------------------------------------------------------------
set.seed(1)
output = EBglmnet(x,y,hyperparameters = c(0.1, 0.1))

## ------------------------------------------------------------------------
glmfit = output$fit
variables = xNames[glmfit[,1,drop=FALSE]]
cbind(variables,as.data.frame(round(glmfit[,3:6,drop=FALSE],4)))

## ------------------------------------------------------------------------
cvfit = cv.EBglmnet(x, y)

## ------------------------------------------------------------------------
cvfit$CrossValidation

## ------------------------------------------------------------------------
cvfit$hyperparameters
cvfit$fit

## ------------------------------------------------------------------------
output$Intercept
output$residual

## ------------------------------------------------------------------------
yy = y>mean(y);
output = EBglmnet(x,yy,family="binomial", hyperparameters = c(0.1, 0.1))

## ------------------------------------------------------------------------
output = EBglmnet(x,yy,family="binomial", prior = "elastic net", hyperparameters = c(0.1, 0.1))

## ------------------------------------------------------------------------
output = EBglmnet(x,yy,family="binomial", prior = "elastic net", 
                  hyperparameters = c(0.1, 0.1),Epis = TRUE)
output$fit

## ------------------------------------------------------------------------
data(BASIS)#this is the genotype of the the F2 population
N = nrow(BASIS)
p = ncol(BASIS)
j = sample((p-2),1)
cor(BASIS[,c(j,j+1,j+2)]) #Correlation structure among neighboring markers

## ------------------------------------------------------------------------
set.seed(1);
Mu = 100; #population mean;
nTrue = 10; # we assume 10 out of the 481 effects are true QTLs
trueLoc	= sort(sample(p,nTrue));
trueEff = runif(nTrue,2,3); #effect size from 2-3
xbeta = BASIS[,trueLoc]%*%trueEff;
s2 =  var(xbeta)*0.1/0.9 #residual variance with 10% noise
residual = rnorm(N,mean=0,sd=sqrt(s2))
y = Mu + xbeta + residual;  

## ------------------------------------------------------------------------
n = 300;
index = sample(N,n); 
CV = cv.EBglmnet(x=BASIS[index,],y=y[index],family="gaussian",prior= "lassoNEG",nfold= 5)

## ------------------------------------------------------------------------
CV$fit
trueLoc

## ------------------------------------------------------------------------
n = 300;
set.seed(1)
index = sample(nrow(BASIS),n)
p = ncol(BASIS);
m = p*(p+1)/2;
#1. simulate true QTL locations
nMain = 10;
nEpis  = 10;
mainLoc = sample(p,nMain);
episLoc = sample(seq((p+1),m,1),nEpis);
trueLoc = sort(c(mainLoc,episLoc)); #a vector in [1,m]
nTrue = length(trueLoc);
trueLocs = ijIndex(trueLoc, p); #two columns denoting the pair (i,j)
#2. obtain true QTL genotype
basis = matrix(0,n,nTrue);
for(i in 1:nTrue)
{
	if(trueLocs[i,1]==trueLocs[i,2])
	{
		basis[,i] = BASIS[index,trueLocs[i,1]]
	}else
	{
		basis[,i] = BASIS[index,trueLocs[i,1]]*BASIS[index,trueLocs[i,2]]
	}
}
#3. simulate true QTL effect size	
trueEff  = runif(nTrue,2,3);
#4. simulate phenotype
xbeta = basis%*%trueEff;
vary = var(xbeta);
Pr = 1/(1+ exp( -xbeta));
y = rbinom(n,1,Pr);

## ----eval=FALSE----------------------------------------------------------
#  CV = cv.EBglmnet(x=BASIS[index,],y=y,family="binomial",prior="lasso",nfold=5,Epis =TRUE)
#  ind = which(CV$fit[,6]<=0.1)#p-value cutoff
#  CV$fit[ind,]

## ------------------------------------------------------------------------
X = matrix(0,n,m);
X[,1:p] = BASIS[index,];
kk = p + 1;
for(i in 1:(p-1))
{
	for(j in (i+1):p)
	{
		X[,kk] = BASIS[index,i] * BASIS[index,j];
		kk = kk + 1;
	}
}

## ------------------------------------------------------------------------
library(glmnet);
alpha = 1
lambdaRatio = 1e-4; #same as in EBlasso
cv = cv.glmnet(X, y, alpha = alpha,family="binomial",nfolds = 5,lambda.min.ratio=lambdaRatio)
nLambda = length(cv$lambda)
nLambda
nbeta = rep(0,nLambda);
fit0 = cv$glmnet.fit;
for(i in 1:nLambda)
{
  nbeta[i] = length(which(fit0$beta[,i]!=0))
}
plot(nbeta,xlab=expression(paste(lambda, " in lasso selection path(n=300,p=115,921)")),
     ylab="No. of nonzero effects",xaxt="n")#
ticks = seq(1,nLambda,10)
axis(side=1, at= ticks,labels=round(cv$lambda[ticks],5), las=1,cex.axis = 0.5)
title("Number of nonzero effects in lasso selection path")

## ------------------------------------------------------------------------
lambda= cv$lambda.min
coefs = fit0$beta
ind = which(cv$lambda==cv$lambda.min)
beta = coefs[,ind]
betaij = which(beta!=0)
Xdata = X[,betaij];
colnames(Xdata) = betaij; 
refit = glm(y ~ Xdata, family = binomial(link = "logit"))#separation occurs

## ----eval=FALSE----------------------------------------------------------
#  #1. create the hyperparameters to be evaluated
#  familyPara = "gaussian";
#  priorPara = "elastic net";
#  epis =TRUE;
#  lambda_Max = lambdaMax(BASIS[index,],y,epis);
#  nStep = 9;#steps from lambda_max to 0.0001*lambda_min
#  lambda_Min = 0.0001*lambda_Max;
#  step = (log(lambda_Max) - log(lambda_Min))/nStep;
#  Lambda	= exp(seq(from = log(lambda_Max),to=log(lambda_Min),by= -step))
#  N_step	= length(Lambda);
#  Alpha = seq(from = 1, to = 0.1, by = -0.1)# values of alpha
#  nAlpha = length(Alpha);	
#  nPara = nAlpha * N_step;
#  HyperPara = matrix(0,nPara,2);
#  for(j in 1:nAlpha)
#  {
#  	strt = (j-1)*N_step + 1;
#  	ed 	= (j-1)*N_step + N_step;
#  	HyperPara[strt:ed,1] = rep(Alpha[j],N_step);
#  	HyperPara[strt:ed,2] = Lambda;
#  }
#  
#  #2. create a function for parallel computaiton
#  EBglmnet.CVonePar <-function(iHyper,X,y,nFolds,foldId,hyperPara,...)
#  {
#    parameters= hyperPara[iHyper,];
#    SSE = CVonePair(X,y,nFolds,foldId,parameters,...);
#    return(SSE);
#  }
#  
#  #3. setup parallel computation
#  library(snow)
#  library(EBglmnet)
#  ncl = 4; #use 4 CPUs
#  cl<-makeCluster(ncl,type="SOCK")
#  clusterEvalQ(cl,{library(EBglmnet)})
#  iPar = matrix(seq(1,nPara,1),nPara,1);
#  nFolds = 5;#5 fold CV
#  if(n%%nFolds!=0){
#  	foldId= sample(c(rep(1:nFolds,floor(n/nFolds)),1:(n%%nFolds)),n);
#  }else{
#  	foldId= sample(rep(1:nFolds,floor(n/nFolds)),n);
#  }
#  #call parRapply to perform parallel computation
#  SSE = parRapply(cl,iPar,EBglmnet.CVonePar,BASIS[index,],y,nFolds,foldId,HyperPara,
#                  Epis = epis, prior = priorPara, family= familyPara);	
#  #collect the result in CVresult
#  CVresult =matrix(SSE,nPara,4,byrow=TRUE)#4 columns of cv$
#  stopCluster(cl)

