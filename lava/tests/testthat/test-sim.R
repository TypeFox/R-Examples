context("Simulation")

test_that("Missing", {
    m <- lvm(y~1)
    m <- Missing(m,y~1,r~x)
    set.seed(1)
    d <- simulate(m,1e3,seed=1)
    expect_equal(sum(d$r),sum(!is.na(d$y0)))

    g <- glm(r~x,data=d,family=binomial)
    expect_true(all.equal(coef(g),c(0,1),tolerance=0.2,check.attributes=FALSE))
})


test_that("sim.default I", {
    m <- lvm(y~x+e)
    distribution(m,~y) <- 0
    distribution(m,~x) <- uniform.lvm(a=-1.1,b=1.1)
    transform(m,e~x) <- function(x) (1*x^4)*rnorm(length(x),sd=1)

    onerun <- function(iter=NULL,...,n=2e3,b0=1,idx=2) {
        d <- sim(m,n,p=c("y~x"=b0))
        l <- lm(y~x,d)
        res <- c(coef(summary(l))[idx,1:2],
                 confint(l)[idx,],
                 estimate(l,only.coef=TRUE)[idx,2:4])
        names(res) <- c("Estimate","Model.se","Model.lo","Model.hi",
                        "Sandwich.se","Sandwich.lo","Sandwich.hi")
        res
    }

    val <- sim(onerun,R=2,b0=1,n=10,messages=0,mc.cores=1)
    expect_true(nrow(val)==2)
    val <- sim(val,R=2,b0=1,n=10,type=0) ## append results
    expect_true(nrow(val)==4)

    s1 <- summary(val,estimate=c(1,1),confint=c(3,4,6,7),true=c(1,1),names=c("Model","Sandwich"))
    expect_true(length(grep("Coverage",rownames(s1)))>0)
    expect_equivalent(colnames(s1),c("Model","Sandwich"))
    
    val <- sim(onerun,R=2,cl=TRUE,seed=1,messages=0,mc.cores=2)
    expect_true(val[1,1]!=val[1,2])
        
    onerun2 <- function(a,b,...) {
        return(cbind(a=a,b=b,c=a-1,d=a+1))
    }
    R <- data.frame(a=1:2,b=3:4)
    dm <- capture.output(val2 <- sim(onerun2,R=R,messages=1,mc.cores=2))
    expect_true(all(R-val2[,1:2]==0))
    res <- summary(val2)
    expect_equivalent(res["Mean",],c(1.5,3.5,0.5,2.5))

    expect_output(print(val2[1,]),"a b c d")
    expect_output(print(val2[1,]),"1 3 0 2")
       
    res <- summary(val2,estimate="a",se="b",true=1.5,confint=c("c","d"))
    expect_true(res["Coverage",]==1)
    expect_true(res["SE/SD",]==mean(val2[,"b"])/sd(val2[,"a"]))
      
})


test_that("distributions", {
    m <- lvm(y1~x)
    distribution(m,~y1) <- binomial.lvm("probit")
    distribution(m,~y2) <- poisson.lvm()
    distribution(m,~y3) <- normal.lvm(mean=1,sd=2)
    distribution(m,~y3) <- lognormal.lvm()
    distribution(m,~y3) <- pareto.lvm()
    distribution(m,~y3) <- loggamma.lvm()
    distribution(m,~y3) <- weibull.lvm()
    distribution(m,~y3) <- chisq.lvm()
    distribution(m,~y3) <- student.lvm(mu=1,sigma=1)    

    expect_output(print(distribution(m)$y2),"Family: poisson")
    expect_output(print(distribution(m)$y1),"Family: binomial")
    latent(m) <- ~u    
    expect_output(print(m),"binomial\\(probit\\)")
    expect_output(print(m),"poisson\\(log\\)")

    ## Generator:
    m <- lvm()
    distribution(m,~y,TRUE) <- function(n,...) {
        res <- exp(rnorm(n)); res[seq(min(n,5))] <- 0
        return(res)
    }
    d <- sim(m,10)
    expect_true(all(d[1:5,1]==0))
    expect_true(all(d[6:10,1]!=0))

    m <- lvm()
    distribution(m,~y,"a",init.par=2) <- function(n,a,...) {
        rep(1,n)*a
    }
    expect_true(all(sim(m,2)==2))
    expect_true(all(sim(m,2,p=c(a=10))==10))
    expect_equivalent(sim(m,2,p=c(a=10)),sim(m,2,a=10))

    ## Multivariate distribution
    m <- lvm()
    rmr <- function(n,rho,...) rmvn(n,sigma=diag(2)*(1-rho)+rho)
    distribution(m,~y1+y2,rho=0.9) <- rmr
    expect_equivalent(c("y1","y2"),colnames(d <- sim(m,5)))

    ## Special 'distributions'
    m <- lvm()
    distribution(m,~x1) <- sequence.lvm(int=TRUE)
    distribution(m,~x2) <- sequence.lvm(a=1,b=2)
    distribution(m,~x3) <- sequence.lvm(a=NULL,b=2)
    distribution(m,~x4) <- sequence.lvm(a=2,b=NULL)
    ex <- sim(m,5)
    expect_equivalent(ex$x1,1:5)
    expect_equivalent(ex$x2,seq(1,2,length.out=5))
    expect_equivalent(ex$x3,seq(-2,2))
    expect_equivalent(ex$x4,seq(2,6))

    m <- lvm()
    distribution(m,~x1) <- ones.lvm()
    distribution(m,~x2) <- ones.lvm(p=0.5)
    distribution(m,~x3) <- ones.lvm(interval=c(0.4,0.6))
    ex <- sim(m,10)
    expect_equivalent(ex$x1,rep(1,10))
    expect_equivalent(ex$x2,c(rep(0,5),rep(1,5)))
    expect_equivalent(ex$x3,c(0,0,0,1,1,1,0,0,0,0))

    m <- lvm()
    expect_error(distribution(m,~y) <- threshold.lvm(p=c(0.5,.75)))
    distribution(m,~y) <- threshold.lvm(p=c(0.25,.25))
    set.seed(1)
    expect_equivalent(1:3,sort(unique(sim(m,200))[,1]))

    ## distribution(m,~y) <- threshold.lvm(p=c(0.25,.25),labels=letters[1:3])
    ## expect_equivalent(c("a","b","c"),sort(unique(sim(m,200))[,1]))
        
})


test_that("eventTime", {
    m <- lvm(eventtime~x)
    distribution(m,~eventtime) <- coxExponential.lvm(1/100)
    distribution(m,~censtime) <- coxWeibull.lvm(1/500)
    eventTime(m) <- time~min(eventtime=1,censtime=0)

    set.seed(1)
    d <- sim(m,100)
    expect_equivalent((d$time<d$cens)*1L,d$status)

    ## TODO
    plot(m)
    expect_output(print(m),"Event History Model")

    ## Time varying effect
    m <- lvm(y~1)
    distribution(m,~z1) <- ones.lvm(0.5)
    R <- log(cbind(c(0.2,0.7,0.9),c(0.5,0.3,0.3)))
    m <- timedep(m,y~z1,timecut=c(0,3,5),rate=R)
   
    ## sim(m,100)
    ## d <- sim(m,1e4); d$status <- TRUE
    ## dd <- mets::lifetable(survival::Surv(y,status)~z1,data=d,breaks=c(0,3,5,Inf));
    ## exp(coef(glm(events ~ offset(log(atrisk)) + -1 + interval+z1:interval, dd, family=poisson)))

})


