ss.aipe.R2.sensitivity <- function(True.R2=NULL, Estimated.R2=NULL, w=NULL, p=NULL, 
Random.Predictors=TRUE, Selected.N=NULL, degree.of.certainty=NULL, assurance=NULL, certainty=NULL, conf.level=.95, 
Generate.Random.Predictors=TRUE, rho.yx=.3, rho.xx=.3, G=10000, print.iter=TRUE, ...)
{
  requireNamespace("MASS", quietly = TRUE)
  
if(!is.null(certainty)& is.null(degree.of.certainty)&is.null(assurance)) degree.of.certainty<-certainty
if (is.null(assurance) && !is.null (degree.of.certainty)& is.null(certainty)) assurance <-degree.of.certainty
if (!is.null(assurance) && is.null (degree.of.certainty)& is.null(certainty)) assurance -> degree.of.certainty

if(!is.null(assurance) && !is.null (degree.of.certainty) && assurance!=degree.of.certainty) 
stop("The arguments 'assurance' and 'degree.of.certainty' must have the same value.")

if(!is.null(assurance) && !is.null (certainty) && assurance!=certainty) 
stop("The arguments 'assurance' and 'certainty' must have the same value.")

if(!is.null(degree.of.certainty) && !is.null (certainty) && degree.of.certainty!=certainty) 
stop("The arguments 'degree.of.certainty' and 'certainty' must have the same value.")

if(True.R2>=1 | True.R2<=0) stop("The values of \'True.R2\' (i.e., the squared multiple correlation coefficient (R^2)) must be between zero and one.")
if(w==0 | w>=1) stop("The width is not specified correctly.")
if(w==0 | w>=1) stop("The width is not specified correctly.")

if(is.null(Estimated.R2) & is.null(Selected.N)) stop("You must specify either \'Estimated.R2\' or \'Selected.N\'.", call.=FALSE)

if(!is.null(degree.of.certainty))
{
if(degree.of.certainty<0 | degree.of.certainty>1) stop("You must specify either \'degree.of.certainty\' to be a value between 0 and 1.", call.=FALSE)
}

options(warn=-1)

if(!is.null(Estimated.R2))
{
if(Estimated.R2>=1 | Estimated.R2<=0) stop("The values of \'Estimated.R2\' (i.e., the squared multiple correlation coefficient (R^2)) must be between zero and one.")
N <- ss.aipe.R2(Population.R2=Estimated.R2, conf.level=conf.level, width=w, which.width="Full", p=p, degree.of.certainty=degree.of.certainty, Random.Predictors=Random.Predictors)$Required.Sample.Size
}
else
{
N <- Selected.N
}
#############################################################################################################################

# Means (arbitrary)
MU <- rep(0, p+1) 

# Correlation between Y and the X variables (arbitrary for plausible scenarios)
sigma.YX <- rbind(rep(rho.yx, p)); sigma.XY <- t(sigma.YX)

# Correlation among the predictors (arbitrary for plausible scenarios).
Sigma.XX <- matrix(rep(rho.xx, p^2), nrow=p, ncol=p); diag(Sigma.XX) <- 1

# Defines the numerator so that the desired P^2 (Rho Squared; Population multiple correlation coefficient) can be obtained.
Numerator.P.Square <- (sigma.YX%*%solve(Sigma.XX)%*%sigma.XY)

# Define the variance of Y so that the pop. mult. cor. coef. is as specified.
sigma.Y <- Numerator.P.Square/True.R2

Sigma <- rbind(c(sigma.Y, sigma.YX), cbind(sigma.XY, Sigma.XX))
#############################################################################################################################


R.Square.Results <- matrix(NA, G, 3)
colnames(R.Square.Results) <- c("Lower.CI.Limit.R2", "Observed.R2", "Upper.CI.Limit.R2")

if(Generate.Random.Predictors==TRUE)
{
for(i in 1:G)
{
if(print.iter==TRUE) cat(c(i),"\n")
DATA <- MASS::mvrnorm(N, mu=MU, Sigma=Sigma)

Regression.Results <- lm(DATA[,1] ~ DATA[,-1])
Summary.Regression.Results <- summary(Regression.Results)

R.Square.Results[i,2] <- Summary.Regression.Results$r.squared

CI.Limits.R2 <- try(ci.R2(R2 = R.Square.Results[i,2], conf.level = conf.level, N = N, p = p, Random.Predictors=Random.Predictors))

R.Square.Results[i,1] <- CI.Limits.R2$Lower.Conf.Limit.R2
R.Square.Results[i,3] <- CI.Limits.R2$Upper.Conf.Limit.R2
}
}

if(Generate.Random.Predictors==FALSE)
{
DATA.Pop.Cov.Structure <- MASS::mvrnorm(N, mu=MU, Sigma=Sigma, empirical = TRUE)[,-1]
BETA <- cbind(c(sigma.YX%*%solve(Sigma.XX)))
True.Y <- DATA.Pop.Cov.Structure%*%BETA

for(i in 1:G)
{
if(print.iter==TRUE) cat(c(i),"\n")

# So, only Y is random from sample to sample.
Obs.Y <- True.Y + rnorm(N, 0, sqrt(sigma.Y*(1-True.R2)))

Regression.Results <- lm(Obs.Y ~ DATA.Pop.Cov.Structure)

Summary.Regression.Results <- summary(Regression.Results)

R.Square.Results[i,2] <- Summary.Regression.Results$r.squared
print(Summary.Regression.Results$r.squared)
CI.Limits.R2 <- try(ci.R2(R2 = R.Square.Results[i,2], conf.level = conf.level, N = N, p = p, Random.Predictors=Random.Predictors))

R.Square.Results[i,1] <- CI.Limits.R2$Lower.Conf.Limit.R2
R.Square.Results[i,3] <- CI.Limits.R2$Upper.Conf.Limit.R2
}
}


#Summary Section
Lower.Type.I.Error <- mean(True.R2 <= R.Square.Results[,1])
Upper.Type.I.Error <- mean(True.R2 >= R.Square.Results[,3])
Type.I.Error <- Lower.Type.I.Error + Upper.Type.I.Error

Lower.Width.CI <- R.Square.Results[,2] - R.Square.Results[,1]
Upper.Width.CI <- R.Square.Results[,3] - R.Square.Results[,2]
Width.CI <- Lower.Width.CI + Upper.Width.CI
#############################################################

Results <- list(Lower.Limit.R2=R.Square.Results[,1], R2=R.Square.Results[,2], Upper.Limit.R2=R.Square.Results[,3], 
Lower.Width.CI=Lower.Width.CI, Upper.Width.CI=Upper.Width.CI, Width.CI=Width.CI)

Num.Probs.with.CIs <- G-length(na.omit(Results$Width))

Specifications <- list(Desired.width=w, True.R2=True.R2, Estimated.R2=Estimated.R2, Num.of.Predictors=p, N=N, degree.of.certainty=degree.of.certainty, Num.of.Replications=G, Conf.Level=conf.level, Random.Predictors=Random.Predictors, Generate.Random.Predictors=Generate.Random.Predictors)

Summary <- list(mean.low.lim.R2=mean(Results$Lower.Limit.R2, na.rm=TRUE), median.low.lim.R2=median(Results$Lower.Limit.R2, na.rm=TRUE), sd.low.lim.R2=sqrt(var(Results$Lower.Limit.R2, na.rm=TRUE)),
mean.up.lim.R2=mean(Results$Upper.Limit.R2, na.rm=TRUE), median.up.lim.R2=median(Results$Upper.Limit.R2, na.rm=TRUE), sd.up.lim.R2=sqrt(var(Results$Upper.Limit.R2, na.rm=TRUE)),
mean.R2=mean(Results$R2, na.rm=TRUE), median.R2=median(Results$R2, na.rm=TRUE), sd.R2=sqrt(var(Results$R2, na.rm=TRUE)),mean.lower.CI.width.R2=mean(Results$Lower.Width.CI, na.rm=TRUE), median.lower.CI.width.R2=median(Results$Lower.Width.CI, na.rm=TRUE), sd.lower.CI.width.R2=sqrt(var(Results$Lower.Width.CI, na.rm=TRUE)),
mean.upper.CI.width.R2=mean(Results$Upper.Width.CI, na.rm=TRUE), median.upper.CI.width.R2=median(Results$Upper.Width.CI, na.rm=TRUE), sd.upper.CI.width.R2=sqrt(var(Results$Upper.Width.CI, na.rm=TRUE)),
mean.CI.width.R2=mean(Results$Width.CI, na.rm=TRUE), median.CI.width.R2=median(Results$Width.CI, na.rm=TRUE), sd.CI.width.R2=sqrt(var(Results$Width.CI, na.rm=TRUE)), Pct.Less.w=mean(Width.CI<=w, na.rm=TRUE), Num.Probs.with.CIs=Num.Probs.with.CIs)

return(list(Results=Results, Specifications=Specifications, Summary=Summary))
options(warn=1)
}
