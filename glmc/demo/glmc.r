#
pause <- function(){readline(prompt="Pause. Press <Enter> to continue...");invisible()}
# In this example we apply the methodology to combine survey data
# from the British Household Panel Survey (BHPS)(Taylor et al.,
# 1995) with population level information from the birth
# registration system on the General Fertility Rate (GFR) to
# estimate annual birth probabilities by parity. This situation was
# considered by Handcock, Huovilainen, and Rendall (2000).

#Specify the data. 

n <- rbind(c(5903,230),c(5157,350))
dimnames(n) <- list(c("no child","child"), c("no birth","birth"))
#
n
pause()
#
# Create indicator covariates
#
mat <- matrix(0,nrow=sum(n),ncol=2)
mat <- rbind(matrix(1,nrow=n[1,1],ncol=1)%*%c(0,0),
             matrix(1,nrow=n[1,2],ncol=1)%*%c(0,1),
             matrix(1,nrow=n[2,1],ncol=1)%*%c(1,0),
             matrix(1,nrow=n[2,2],ncol=1)%*%c(1,1))

#Specifying the population constraints.
#
# From the combination of birth registration numerator and
# population estimate denominator we can determine the general
# fertility rate (GFR) for the years 1992 to 1996, of England and
# Wales, which we assume to be measured without error (Office for
# National Statistics, 1998). The population level value of the GFR
# was found to be 0.06179.

gfr <- .06179*matrix(1,nrow=nrow(mat),ncol=1)
g <- matrix(1,nrow=nrow(mat),ncol=1)
amat <- matrix(mat[,2]*g-gfr,ncol=1)

# Method 1. Defining constraints in the data frame.

hrh <- data.frame(birth=mat[,2], child=mat[,1], constraints=amat)

gfit <- glmc(birth~child, data=hrh, family="binomial",emplik.method="Owen",
             control=glmc.control(maxit.glm=10,maxit.weights=200,
             itertrace.weights=TRUE,gradtol.weights=10^(-6)))

summary(gfit)
pause()

# Method 2. Defining constraints through a matrix.

gfit<- glmc(birth~child,data=hrh,family=binomial(link=logit),
            emplik.method="Owen",control=glmc.control(maxit.glm=10,
            maxit.weights=200,itertrace.weights=TRUE,gradtol.weights=10^(-10)),
            Amat=amat,confn=NULL)

summary(gfit)
pause()

# Method 3. Defining constraints through a function.

fn <- function(t,Y,X){
 grf <- .06179*matrix(1,nrow=nrow(as.matrix(X)),ncol=1)
 g <- matrix(1,nrow=nrow(X),ncol=1)
 amat <- matrix(Y*g-grf,ncol=1)
 return(amat)
}

gfit <- glmc(birth~child,data=hrh,family=binomial(link=logit),
             optim.method="BFGS",emplik.method="Owen",
             control=glmc.control(maxit.glm=10,maxit.optim=10^(8),
             reltol.optim=10^(-10),trace.optim=9,REPORT.optim=1,
             maxit.weights=200,gradtol.weights=10^(-6),itertrace.weights=FALSE),
             optim.hessian=TRUE,Amat=NULL,confn=fn)

summary(gfit)
