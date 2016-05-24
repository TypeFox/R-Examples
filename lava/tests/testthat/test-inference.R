context("Inference")

test_that("Effects",{
    m <- lvm()
    regression(m) <- c(y1,y2,y3)~u; latent(m) <- ~u
    regression(m) <- c(z1,z2,z3)~v; latent(m) <- ~v
    regression(m) <- u~v
    regression(m) <- c(u,v,z3,y1)~x
    d <- sim(m,100)
    start <- c(rep(0,6),rep(1,17))
    suppressWarnings(e <- estimate(m,d,control=list(iter.max=0,start=start)))
    f <- coef(ef <- effects(e,y1~x))
    expect_true(all(f[,2]>0)) ## Std.err
    expect_equal(f["Total",1],3) 
    expect_equal(f["Direct",1],1)
    f2 <- coef(effects(e,u~v))
    expect_equal(f2["Total",1],1)
    expect_equal(f2["Direct",1],1)
    expect_equal(f2["Indirect",1],0)

    expect_output(print(ef),"Indirect effects")
    expect_equivalent(confint(ef)["Direct",],
                      confint(e)["y1~x",])

    expect_equivalent(totaleffects(e,y1~x)[,1:4],f["Total",])
    
    g <- graph::updateGraph(plot(m,noplot=TRUE))
    expect_equivalent(path(g,y1~x),path(m,y1~x))    
})

test_that("Profile confidence limits", {
    m <- lvm(y~b*x)
    constrain(m,b~psi) <- identity
    set.seed(1)
    d <- sim(m,100)
    e <- estimate(m, d)
    ci0 <- confint(e,3)
    ci <- confint(e,3,profile=TRUE)
    expect_true(mean((ci0-ci)^2)<0.1)
})

test_that("IV-estimator", {
    m <- lvm(c(y1,y2,y3)~u); latent(m) <- ~u    
    set.seed(1)
    d <- sim(m,100)
    e0 <- estimate(m,d)
    e <- estimate(m,d,estimator="iv") ## := MLE
    expect_true(mean((coef(e)-coef(e0))^2)<1e-9)
})

test_that("glm-estimator", {         
    m <- lvm(y~x+z)
    regression(m) <- x~z
    distribution(m,~y+z) <- binomial.lvm("logit")
    set.seed(1)
    d <- sim(m,1e3)
    head(d)
    e <- estimate(m,d,estimator="glm")
    c1 <- coef(e,2)[c("y","y~x","y~z"),1:2]
    c2 <- estimate(glm(y~x+z,d,family=binomial))$coefmat[,1:2]  
    expect_equivalent(c1,c2)
})


test_that("gaussian", {
    m <- lvm(y~x)
    d <- simulate(m,100,seed=1)
    S <- cov(d[,vars(m),drop=FALSE])
    mu <- colMeans(d[,vars(m),drop=FALSE])
    f <- function(p) lava:::gaussian_objective.lvm(p,x=m,S=S,mu=mu,n=nrow(d))
    g <- function(p) lava:::gaussian_score.lvm(p,x=m,n=nrow(d),data=d,indiv=TRUE)
    s1 <- numDeriv::grad(f,c(0,1,1))
    s2 <- g(c(0,1,1))
    expect_equal(s1,-colSums(s2),tolerance=0.1)
})


test_that("Association measures", {
    P <- matrix(c(0.25,0.25,0.25,0.25),2)
    a1 <- lava:::assoc(P)
    expect_equivalent(-log(0.25),a1$H)
    expect_true(with(a1, all(c(kappa,gamma,MI,U.sym)==0)))
    
    p <- lava:::prob.normal(sigma=diag(nrow=2),breaks=c(-Inf,0),breaks2=c(-Inf,0))[1]
    expect_equal(p[1],0.25)
    ## q <- qnorm(0.75)
    ## m <- ordinal(lvm(y~x),~y, K=3)#, breaks=c(-q,q))
    ## normal.threshold(m,p=c(0,1,2))
})


test_that("equivalence", {
    m <- lvm(c(y1,y2,y3)~u,u~x,y1~x)
    latent(m) <- ~u
    d <- sim(m,100)
    cancel(m) <- y1~x
    regression(m) <- y2~x
    e <- estimate(m,d)
    ##eq <- equivalence(e,y1~x,k=1)
    dm <- capture.output(eq <- equivalence(e,y2~x,k=1))
    expect_output(print(eq),"y1,y2")
    expect_true(all(c("y1","y3")%in%eq$equiv[[1]][1,]))    
})

test_that("multiple testing", {
    expect_equivalent(lava:::holm(c(0.05/3,0.025,0.05)),rep(0.05,3))
    
    ci1 <- scheffe(l <- lm(1:5~c(0.5,0.7,1,1.3,1.5)))
    ci2 <- predict(l,interval="confidence")
    expect_equivalent(ci1[,1],ci2[,1])
    expect_true(all(ci1[,2]<ci2[,2]))
    expect_true(all(ci1[,3]>ci2[,3]))
})


test_that("modelsearch and GoF", {
    m <- lvm(c(y1,y2)~x)
    d <- sim(m,100)
    e <- estimate(lvm(c(y1,y2)~1,y1~x),d)
    e0 <- estimate(lvm(c(y1,y2)~x,y1~~y2),d)

    s1 <- modelsearch(e,silent=TRUE)
    expect_true(nrow(s1$res)==2)
    s2 <- modelsearch(e0,silent=TRUE,dir="backward")
    expect_true(nrow(s2$res)==3)
    e00 <- estimate(e0,vcov=vcov(e0))$coefmat
    ii <- match(s2$res[,"Index"],rownames(e00))
    expect_equivalent(e00[ii,5],s2$test[,2])

    s3 <- modelsearch(e0,dir="backward",k=3)
    expect_true(nrow(s3$res)==1)
    
    ee <- modelsearch(e0,dir="backstep",messages=FALSE)
    expect_true(inherits(ee,"lvm"))

    ## TODO
    gof(e,all=TRUE)    
    r <- rsq(e)[[1]]
    expect_true(abs(summary(lm(y1~x,d))$r.square-r["y1"])<1e-5)
    
})

test_that("Bootstrap", {
    y <- rep(c(0,1),each=5)
    x <- 1:10
    e <- estimate(y~x)
    B1 <- bootstrap(e,R=2,silent=TRUE,mc.cores=1,sd=TRUE)    
    B2 <- bootstrap(e,R=2,silent=TRUE,bollenstine=TRUE,mc.cores=1)
    expect_false(B1$bollenstine)
    expect_true(B2$bollenstine)
    expect_true(nrow(B1$coef)==2)
    expect_output(print(B1),"Standard errors:")
    dm <- capture.output(B3 <- bootstrap(e,R=2,fun=function(x) coef(x)[2]^2+10))
    expect_true(all(mean(B3$coef)>10))

    y <- rep(c(0,1),each=5)
    x <- 1:10
    m <- lvm(y~b*x)
    constrain(m,alpha~b) <- function(x) x^2
    e <- estimate(m,data.frame(y=y,x=x))
    b <- bootstrap(e,R=1,silent=TRUE)
    expect_output(print(b),"alpha")
})


test_that("Survreg", {
    m <- lvm(y0~x)   
    transform(m,y~y0) <- function(x) pmin(x,2)
    transform(m,status~y0) <- function(x) x<2
    d <- simulate(m,100,seed=1)
    require('survival')
    m <- survreg(Surv(y,status)~x,data=d,dist='gaussian')
    s <- score(m)
    expect_true(length(pars(m))==length(coef(m))+1)
    expect_true(abs(attr(score(m,pars(m)),'logLik')-logLik(m))<1e-9)
    expect_true(mean(colSums(s)^2)<1e-6)
    expect_equivalent(vcov(m),attr(s,'bread'))
})



test_that("Combine", { ## Combine model output
    data(serotonin)
    m1 <- lm(cau ~ age*gene1 + age*gene2,data=serotonin)
    m2 <- lm(cau ~ age + gene1,data=serotonin)

    cc <- Combine(list('model A'=m1,'model B'=m2),fun=function(x) c(R2=format(summary(x)$r.squared,digits=2)))
    expect_true(nrow(cc)==length(coef(m1))+1)
    expect_equivalent(colnames(cc),c('model A','model B'))
    expect_equivalent(cc['R2',2],format(summary(m2)$r.squared,digits=2))
})



test_that("zero-inflated binomial regression (zib)", {
    set.seed(1)
    n <- 1e3
    x <- runif(n,0,20)
    age <- runif(n,10,30)
    z0 <- rnorm(n,mean=-1+0.05*age)
    z <- cut(z0,breaks=c(-Inf,-1,0,1,Inf))
    p0 <- lava::expit(model.matrix(~z+age) %*% c(-.4, -.4, 0.2, 2, -0.05))
    y <- (runif(n)<lava:::tigol(-1+0.25*x-0*age))*1
    u <- runif(n)<p0
    y[u==0] <- 0
    d <- data.frame(y=y,x=x,u=u*1,z=z,age=age)

    ## Estimation
    e0 <- zibreg(y~x*z,~1+z+age,data=d)    
    e <- zibreg(y~x,~1+z+age,data=d)
    compare(e,e0)
    expect_equivalent(score(e,data=d),
                      colSums(score(e,indiv=TRUE)))
    expect_equivalent(logLik(e,data=d),
                      logLik(e)) 
    expect_equivalent(vcov(e), Inverse(information(e,type="obs",data=d)))

    expect_output(print(e), "Prevalence probabilities:")
    
    PD(e0,intercept=c(1,3),slope=c(2,6))

    B <- rbind(c(1,0,0,0,20),
               c(1,1,0,0,20),
               c(1,0,1,0,20),
               c(1,0,0,1,20))
    prev <- summary(e,pr.contrast=B)$prevalence

    x <- seq(0,100,length.out=100)
    newdata <- expand.grid(x=x,age=20,z=levels(d$z))
    fit <- predict(e,newdata=newdata)
})

test_that("Optimization", {
    m <- lvm(y~x+z)
    d <- simulate(m,20,seed=1)
    e1 <- estimate(m,d,control=list(method="nlminb0"))
    e2 <- estimate(m,d,control=list(method="NR"))
    expect_equivalent(round(coef(e1),3),round(coef(e2),3))

    f <- function(x) x^2*log(x) # x>0
    df <- function(x) 2*x*log(x) + x
    df2 <- function(x) 2*log(x) + 3
    op <- NR(5,f,df,df2,control=list(tol=1e-40)) ## Find root
    expect_equivalent(round(op$par,digits=7),.6065307)
    op2 <- estfun0(5,gradient=df)
    op3 <- estfun(5,gradient=df,hessian=df2,control=list(tol=1e-40))
    expect_equivalent(op$par,op2$par)
    expect_equivalent(op$par,op3$par)
})


test_that("Prediction with missing data, random intercept", {
    ## Random intercept model
    m <- lvm(c(y1,y2,y3)~u[0])
    latent(m) <- ~u
    regression(m) <- y1~x1
    regression(m) <- y2~x2
    regression(m) <- y3~x3

    d <- simulate(m,1e3,seed=1)
    dd <- reshape(d,varying=list(c('y1','y2','y3'),c('x1','x2','x3')),direction='long',v.names=c('y','x'))
    
    ##system.time(l <- lme4::lmer(y~x+(1|id), data=dd, REML=FALSE))
    system.time(l <- nlme::lme(y~x,random=~1|id, data=dd, method="ML"))
    m0 <- lvm(c(y1[m:v],y2[m:v],y3[m:v])~1*u[0])
    latent(m0) <- ~u
    regression(m0,y=c('y1','y2','y3'),x=c('x1','x2','x3')) <- rep('b',3)
    system.time(e <- estimate(m0,d))

    mytol <- 1e-6
    mse <- function(x,y=0) mean(na.omit(x-y)^2)
    expect_true(mse(logLik(e),logLik(l))<mytol)
    expect_true(mse(nlme::fixef(l),coef(e)[1:2])<mytol)
    u1 <- nlme::ranef(l)##[[1]][,1]
    u2 <- predict(e,endogenous(e))
    expect_true(mse(u1,u2)<1e-9)

    ## Missing data
    idx <- sample(seq(nrow(dd)),nrow(dd)*0.5)
    dd0 <- dd[idx,,drop=FALSE]
    d0 <- mets::fast.reshape(subset(dd0,select=-u),id='id',num='time')
    
    system.time(e0 <- estimate(m0,d0,missing=TRUE))
    ##system.time(l0 <- lme4::lmer(y~x+(1|id), data=dd0, REML=FALSE))
    system.time(l0 <- nlme::lme(y~x,random=~1|id, data=dd0, method="ML"))
    expect_true(mse(logLik(e0),logLik(l0))<mytol)
    expect_true(mse(nlme::fixef(l0),coef(e0)[1:2])<mytol)
    u01 <- nlme::ranef(l0)##[[1]][,1]
    u02 <- predict(e0,endogenous(e0))
    expect_true(mse(u01,u02)<1e-9)

    s <- summary(e0)
    expect_output(print(s),paste0("data=",nrow(d0)))
    expect_true(inherits(e0$estimate,"multigroupfit"))
    expect_output(print(e0$estimate),"Group 1")
    expect_output(print(summary(e0$estimate)),paste0("observations = ",nrow(d0)))

})



test_that("Random slope model", {

    set.seed(1)
    m <- lvm()
    regression(m) <- y1 ~ 1*u+1*s
    regression(m) <- y2 ~ 1*u+2*s
    regression(m) <- y3 ~ 1*u+3*s
    latent(m) <- ~u+s
    dw <- sim(m,20)

    dd <- mets::fast.reshape(dw)
    dd0 <- dd[-c(1:2*3),]
    library(lme4)
    l <- lmer(y~ 1+num +(1+num|id),dd,REML=FALSE)
    sl <- lava:::varcomp(l,profile=FALSE)

    d <- mets::fast.reshape(dd,id="id")
    d0 <- mets::fast.reshape(dd0,id="id")
  
    m0 <- lvm(c(y1[0:v],y2[0:v],y3[0:v])~1*u)
    addvar(m0) <- ~num1+num2+num3
    covariance(m0) <- u~s
    latent(m0) <- ~s+u
    regression(m0) <- y1 ~ num1*s
    regression(m0) <- y2 ~ num2*s
    regression(m0) <- y3 ~ num3*s
    system.time(e <- estimate(m0,d,param="none",control=list(trace=0,constrain=FALSE)))
    
    expect_true(mean(sl$coef-coef(e)[c("u","s")])^2<1e-5)
    expect_true((logLik(l)-logLik(e))^2<1e-5)    
    expect_true(mean(diag(sl$varcomp)-coef(e)[c("u,u","s,s")])^2<1e-5)

    ## missing
    expect_output(e0 <- estimate(m0,d0,missing=TRUE,param="none",control=list(method="NR",constrain=FALSE,start=coef(e),trace=1)),
                  "Iter=")
    l0 <- lmer(y~ 1+num +(1+num|id),dd0,REML=FALSE)    
    expect_true((logLik(l0)-logLik(e0))^2<1e-5)    
    
    m1 <- lvm(c(y1[0:v],y2[0:v],y3[0:v])~1*u)
    addvar(m1) <- ~num1+num2+num3
    covariance(m1) <- u~s
    latent(m1) <- ~s+u
    regression(m1) <- y1 ~ b1*s
    regression(m1) <- y2 ~ b2*s
    regression(m1) <- y3 ~ b3*s
    constrain(m1,b1~num1) <- function(x) x
    constrain(m1,b2~num2) <- function(x) x
    constrain(m1,b3~num3) <- function(x) x
    system.time(e1 <- estimate(m1,d,param="none",p=coef(e)))
    expect_true((logLik(e1)-logLik(e))^2<1e-5)

    ## TODO    
    ## missing    
    ## system.time(e10 <- estimate(m1,d0,missing=TRUE,param="none",
    ##                             control=list(trace=coef(e0))))
    
})


test_that("multinomial", {
    set.seed(1)
    breaks <- c(-Inf,-1,0,Inf)
    m <- lvm(); covariance(m,pairwise=TRUE) <- ~y1+y2+y3+y4
    d <- transform(sim(m,5e2),
                   z1=cut(y1,breaks=breaks),
                   z2=cut(y2,breaks=breaks),
                   z3=cut(y3,breaks=breaks),
                   z4=cut(y4,breaks=breaks))

    m <- multinomial(d[,5])
    lev <- levels(d[,5])    
    e <- estimate(l <- lm(d[,5]==lev[1]~1))
    expect_true(abs(vcov(e)[1]-vcov(m)[1])<1e-9)
        
    (a1 <- multinomial(d[,5:6],marginal=TRUE))
    K1 <- kappa(a1) ## Cohen's kappa
    P <- a1$P
    marg1 <- rowSums(P)
    marg2 <- colSums(P)
    expect_equivalent(K1$coef,sum(diag(P)-marg1*marg2)/(1-sum(marg1*marg2)))
    
    K2 <- kappa(d[,7:8])
    ## Testing difference K1-K2:
    e1 <- estimate(merge(K1,K2,id=TRUE),diff)
    e2 <- estimate(merge(K1,K2,id=NULL),diff)
    expect_true(vcov(e1)!=vcov(e2))
    expect_equivalent(vcov(e2),(vcov(K1)+vcov(K2)))

    g1 <- gkgamma(d[,5:6])
    g2 <- gkgamma(table(d[,5:6]))
    g3 <- gkgamma(multinomial(d[,5:6]))
    expect_equivalent(g1$coefmat,g2$coefmat)
    expect_equivalent(g3$coefmat,g2$coefmat)

    
    ## TODO
    lava:::independence(d[,5:6])
    information(d[,5:6])
    ## pcor

    if (requireNamespace("polycor",quietly=TRUE)) {
        require('mvtnorm')
        system.time(rho <- pcor(d[,5],d[,6]))
        rho2 <- polycor::polychor(d[,5],d[,6],ML=TRUE,std.err=TRUE)
        expect_true(abs(rho$coef[1]-rho2$rho)^2<1e-5)
        expect_true(abs(rho$vcov[1]-rho2$var[1])^2<1e-5)
        expect_true(mean(score(rho))^2<1e-5)
    }
})


test_that("predict,residuals", {
    m <- lvm(c(y1,y2,y3)~u,u~x)
    latent(m) <- ~u
    set.seed(1)
    d <- sim(m,100,'y1~u'=1,'y3~u'=3)
    e <- estimate(m,d)
    
    l <- lm(y3~x,data=d)
    e3 <- estimate(y3~x,d)
    expect_true(mean((residuals(l)-residuals(e3))^2)<1e-12)
    expect_true(mean(var(residuals(e3,std=TRUE))[1]*99/100-1)<1e-3)
   
    r <- residuals(e)
    expect_true(ncol(r)==3)         
})


test_that("partialcor", {
    m <- lvm(c(y1,y2,y3)~x1+x2)
    covariance(m) <- c(y1,y2,y3)~y1+y2+y3
    set.seed(1)
    d <- sim(m,500)
    c1 <- partialcor(~x1+x2,d)
    e <- estimate(m,d)
    c2 <- correlation(e)    
    expect_true(mean(c1[,1]-c2[,1])^2<1e-9)
    ## CI, note difference var(z)=1/(n-k-3) vs var(z)=1/(n-3)
    expect_true(mean(c1[,4]-c2[,3])^2<1e-3)
    expect_true(mean(c1[,5]-c2[,4])^2<1e-3)

})

test_that("partialgamma", {
      
})

test_that("multipletesting", {
    ## TODO
})

