## This is testing whether thresholds are increasing when the coverage
## probability requested during calibration is increasing
context("Monotonicity of threshold in Coverage Proability")

chart <- new("SPCCUSUM",model=SPCModelNormal(Delta=1))
chartnp <- new("SPCCUSUM",model=SPCModelNonparCenterScale(Delta=1))

rX <- function(x,truecoeff) {
  X1 <- rbinom(x,1,0.4)
  X2 <- runif(x,0,1)
  X3 <- rnorm(x)
  data.frame(y=truecoeff[2]*X1+truecoeff[3]*X2+truecoeff[4]*X3+rnorm(x,mean=truecoeff[1]),x1=X1,x2=X2,x3=X3)
}
Xlinreg<- rX(1000, truecoeff=c(0,1,1,1))
chartlm <- new("SPCCUSUM",model=SPCModellm(formula="y~x1+x2+x3",Delta=1))


rX <- function(x,truecoeff) {
  X1 <- rbinom(x,1,0.4)
  X2 <- runif(x,0,1)
  X3 <- rnorm(x)
  xbeta <- truecoeff[1]+truecoeff[2]*X1+truecoeff[3]*X2+truecoeff[4]*X3
  p <- exp(xbeta)/(1+exp(xbeta))
  y <- rbinom(x,1,p)
  data.frame(y=y,x1=X1,x2=X2,x3=X3)
}
Xlogreg<- rX(1000, truecoeff=c(-1,100,1,1))
chartlogreg <- new("SPCCUSUM",model=SPCModellogregLikRatio(formula="y~x1+x2+x3",Delta=1))


model <- SPCModelNormal(Delta=1)
model$Pofdata <- function(data){
    list(mu= median(data), sd= mad(data), m=length(data))
}
chartrobust <- new("SPCCUSUM",model=model)
SPCModelExponential=function(Delta=1){
    structure(
        list(
            Pofdata=function(data){
                list(lambda=1/mean(data), n=length(data))
            },
            xiofP=function(P) P$lambda,
            resample=function(P) rexp(P$n,rate=P$lambda),
            getcdfupdates=function(P, xi) {
                function(x){ if(Delta<1)
                                 pmax(0,1-exp(-P$lambda*(x-log(Delta))/(xi*(1-Delta))))
                else
                    pmin(1,exp(-P$lambda*(log(Delta)-x)/(xi*(Delta-1))))
                         }
            },
            updates=function(xi,data) log(Delta)-xi*(Delta-1)*data),
        class="SPCDataModel")
}

ExpCUSUMchart=new("SPCCUSUM",model=SPCModelExponential(Delta=1.25))



testmonotonicity <- function(covprobs=c(0.1,0.5,0.95),property="ARL",
                               target=100,nofreps=50,testchart,X,...){
    ## check ordering of thresholds.
    thresholds <- sapply(covprobs, function(k){
        SPCproperty(data=X,
                    nrep=nofreps,
                    property=paste("cal",property,sep=""),
                    chart=testchart,
                    params=list(target=target,...),covprob=k,quiet=TRUE)@res}
                         )
    for (i in 2:length(covprobs))
        expect_less_than(thresholds[i-1],thresholds[i])
}

test_that("check monotonicity thresholds  Exp",{  
    skip_on_cran()
    X <- rexp(1000)
    testmonotonicity(target=1000,testchart=ExpCUSUMchart,X=X)
    testmonotonicity(target=0.1,testchart=ExpCUSUMchart,X=X,property="hitprob",nsteps=50)
})

test_that("check monotonicity thresholds standard chart",{
    skip_on_cran()
    X <-  rnorm(1000)
    
    testmonotonicity(target=100,testchart=chart,X=X)
    testmonotonicity(target=100,testchart=chartnp,X=X)
    testmonotonicity(target=100,testchart=chartrobust,X=X)
})
