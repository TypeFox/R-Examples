"ss.aipe.reg.coef.sensitivity" <-function(True.Var.Y=NULL,True.Cov.YX=NULL, True.Cov.XX=NULL, 
Estimated.Var.Y=NULL, Estimated.Cov.YX=NULL, Estimated.Cov.XX=NULL, Specified.N=NULL, 
which.predictor=1, w=NULL, Noncentral=FALSE, Standardize=FALSE, conf.level=.95, 
degree.of.certainty=NULL, assurance=NULL, certainty=NULL, G=1000, print.iter=TRUE)
{
if(!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")
  
if(!is.null(certainty)& is.null(degree.of.certainty)&is.null(assurance)) degree.of.certainty<-certainty
if (is.null(assurance) && !is.null (degree.of.certainty)& is.null(certainty)) assurance <-degree.of.certainty
if (!is.null(assurance) && is.null (degree.of.certainty)& is.null(certainty)) assurance -> degree.of.certainty

if(!is.null(assurance) && !is.null (degree.of.certainty) && assurance!=degree.of.certainty) 
stop("The arguments 'assurance' and 'degree.of.certainty' must have the same value.")

if(!is.null(assurance) && !is.null (certainty) && assurance!=certainty) 
stop("The arguments 'assurance' and 'certainty' must have the same value.")

if(!is.null(degree.of.certainty) && !is.null (certainty) && degree.of.certainty!=certainty) 
stop("The arguments 'degree.of.certainty' and 'certainty' must have the same value.")

if(Noncentral==TRUE & is.null(True.Var.Y)) True.Var.Y <- 1

if(is.null(w)) stop("You must specify \'w\' (i.e., a confidence interval width).")
width <- w
if(is.null(conf.level)) stop("You must specify a confidence level (i.e., 1 - Type I error rate).")
if(is.null(G)) stop("You must specify 'G/' (i.e., the number of generations of the simulation).")

if(is.null(True.Cov.XX)) stop("You must specify 'True.Cov.XX' (i.e., the covariance matrix of the predictors).")
if(is.null(True.Cov.YX)) stop("You must specify 'True.Cov.YX' (i.e., the covariance vector of the predictors with the dependent variable).")

if((sum(round(True.Cov.XX,5)==round(t(True.Cov.XX), 5)))!=(dim(True.Cov.XX)[1]*dim(True.Cov.XX)[2])) stop("The correlation matrix, \'True.Cov.XX\' should be symmetric.")


p <- dim(True.Cov.XX)[1]

if(!is.null(Estimated.Cov.XX))
{
if((sum(Estimated.Cov.XX==t(Estimated.Cov.XX)))==(dim(Estimated.Cov.XX)[1]*dim(Estimated.Cov.XX)[2])) stop("The covariance matrix, \'True.Cov.XX\' should be symmetrix.")
}

if(is.null(Estimated.Var.Y)) Estimated.Var.Y <- True.Var.Y
if(is.null(Estimated.Cov.XX)) Estimated.Cov.XX <- True.Cov.XX
if(is.null(Estimated.Cov.YX)) Estimated.Cov.YX <- True.Cov.YX

Estimated.Sigma <- cbind(c(Estimated.Var.Y, Estimated.Cov.YX), rbind(Estimated.Cov.YX, Estimated.Cov.XX))
sigma.Y<- sqrt(Estimated.Sigma[1,1])
sigma.X<- sqrt(Estimated.Sigma[(1+which.predictor),(1+which.predictor)])
Estimated.Rho2.Y_X <- (Estimated.Cov.YX%*%solve(Estimated.Cov.XX)%*%Estimated.Cov.YX)/(sigma.Y^2)
Estimated.Rho2.j_X.without.j <- 1 - ((solve(Estimated.Cov.XX)[which.predictor,which.predictor]*Estimated.Cov.XX)[which.predictor,which.predictor])^(-1)
Estimated.b.j <- (solve(Estimated.Cov.XX)%*%Estimated.Cov.YX)[which.predictor]

# Covariance structure.
True.Sigma <- cbind(c(True.Var.Y, True.Cov.YX), rbind(True.Cov.YX, True.Cov.XX))

True.Rho2.Y_X <- (True.Sigma[1,-1]%*%solve(True.Sigma[-1,-1])%*%True.Sigma[-1,1])/(True.Sigma[1,1])
True.Rho2.j_X.without.j <- 1 - ((solve(True.Cov.XX)[which.predictor,which.predictor]*True.Cov.XX)[which.predictor,which.predictor])^(-1)
True.b.j <- (solve(True.Sigma[-1,-1])%*%True.Sigma[-1,1])[which.predictor]

if(True.Rho2.Y_X>1) stop("You have specified an impossible correlational structure of \'True.Cov.XX\' and/or \'True.Cov.YX\' (the multiple R square is above 1).")
if(Estimated.Rho2.Y_X>1) stop("You have specified an impossible correlational structure of \'Estimated.Cov.XX\' and/or \'Estimated.Cov.YX\' (the multiple R square is above 1).")

# See if this needs to be modified
if(is.null(Specified.N))
{
Estimated.Sigma.as.Cor <- cov2cor(Estimated.Sigma)
N <- ss.aipe.reg.coef(width=width, RHO.XX=Estimated.Sigma.as.Cor[2:(p+1),2:(p+1)], Rho.YX=Estimated.Sigma.as.Cor[1,2:(p+1)], which.predictor=which.predictor, 
conf.level=conf.level, Noncentral=Noncentral, degree.of.certainty=degree.of.certainty, sigma.Y=Estimated.Sigma[1,1]^.5, sigma.X=(Estimated.Sigma[(1+which.predictor),(1+which.predictor)])^.5)
} else N <- Specified.N



# Means (arbitrary)
MU <- rep(0, p+1) 


# Begin simulation.
Results <- matrix(NA, G, 6)
colnames(Results) <- c("b.j", "LL.CI.beta.j", "UL.CI.beta.j", "R.Square", "SE.b.j", "t.for.b.j")
for(i in 1:G)
{
if(print.iter==TRUE) cat(c(i),"\n")
DATA <- MASS::mvrnorm(N, mu=MU, Sigma=True.Sigma)

if(Standardize==TRUE) DATA <- scale(DATA)

Regression.Results <- lm(DATA[,1] ~ DATA[,-1])
Summary.Results <- summary(Regression.Results)

b.j <- coef(Summary.Results)[(which.predictor+1),1]
SE.b.j <- coef(Summary.Results)[(which.predictor+1),2]

if(Noncentral==FALSE) CI.Lims <- ci.reg.coef(b.j=b.j, SE.b.j=SE.b.j, s.Y=(var(DATA[,1]))^.5, s.X=(var(DATA[,1+which.predictor]))^.5, N=dim(DATA)[1], p=(dim(DATA)[2]-1), R2.Y_X=NULL, R2.j_X.without.j=NULL, conf.level=conf.level, R2.Y_X.without.j=NULL, t.value=NULL, alpha.lower=NULL, alpha.upper=NULL, Noncentral=FALSE, Suppress.Statement=TRUE)
if(Noncentral==TRUE) CI.Lims <- ci.reg.coef(b.j=b.j, SE.b.j=SE.b.j, s.Y=(var(DATA[,1]))^.5, s.X=(var(DATA[,1+which.predictor]))^.5, N=dim(DATA)[1], p=(dim(DATA)[2]-1), R2.Y_X=NULL, R2.j_X.without.j=NULL, conf.level=conf.level, R2.Y_X.without.j=NULL, t.value=NULL, alpha.lower=NULL, alpha.upper=NULL, Noncentral=TRUE, Suppress.Statement=TRUE)


Results[i,1] <- b.j
Results[i,2] <- CI.Lims$Lower
Results[i,3] <- CI.Lims$Upper
Results[i,4] <- Summary.Results$r.squared
Results[i,5] <- SE.b.j
Results[i,6] <- b.j/SE.b.j
}
Results <- as.data.frame(Results)
# End Simulation.

Summary.of.Results <- list(Mean.b.j=mean(Results[,1]), Median.b.j=median(Results[,1]), SD.b.j=(var(Results[,1]))^.5, 
Mean.CI.width=mean(Results[,3]-Results[,2]), Median.CI.width=median(Results[,3]-Results[,2]), SD.CI.width=(var(Results[,3]-Results[,2]))^.5, 
Pct.CI.Less.w=mean((Results[,3]-Results[,2])<=w)*100,Pct.CI.Miss.Low=mean(True.b.j < Results[,2])*100, Pct.CI.Miss.High=mean(True.b.j > Results[,3])*100, Total.Type.I.Error=(mean((True.b.j < Results[,2]) | (True.b.j > Results[,3])))*100,
Mean.R2=mean(Results[,4]), Median.R2=median(Results[,4]), sd.R2=(var(Results[,4]))^.5)

###################################################################################################
# Vector of specification values.
if(is.null(degree.of.certainty)) degree.of.certainty <- 0
Specifications <- list(Sample.Size=round(N), True.Rho2.Y_X=round(True.Rho2.Y_X, 5), Estimated.Rho2.Y_X=round(Estimated.Rho2.Y_X, 5),
True.Rho2.j_X.without.j=round(True.Rho2.j_X.without.j, 5), Estimated.Rho2.j_X.without.j=round(Estimated.Rho2.j_X.without.j, 5),
True.b.j=round(True.b.j, 5), Estimated.b.j=round(Estimated.b.j, 5), width.specified=round(width, 5), sigma.Y=round(sigma.Y, 5),  sigma.X=round(sigma.X, 5), 
Noncentral=Noncentral, Standardize=Standardize, conf.level=round(conf.level), degree.of.certainty=round(degree.of.certainty), G=round(G))

return(list(Data.from.Simulation=Results, Specifications=Specifications, Summary.of.Results=Summary.of.Results))
}
