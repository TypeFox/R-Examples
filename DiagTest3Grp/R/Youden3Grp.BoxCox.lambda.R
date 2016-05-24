Youden3Grp.BoxCox.lambda <-
function(xx,yy,zz,nboot=30,lambda0=0.00001)
  {

    ################This function confirms the box-cox transformation parameter lambda being 0 or not:
    ################For xx,  yy, zz from 3 ordinal diagnostic groups, first estimate lambda (lambda.hat) from TN maximizing profile likelihood
    ################ if sample size>=2000, do asymptotic Chisq(df=1) test on H0:lambda.hat=lambda0 (a close to 0 value, eg., 0.0001)
    ################ else do nboot (default=50) bootstrap samples, repeatedly estimate lambda.hat and obtain 0.025~0.975 percentile CI,
    #####                  if the CI contains 0, then lambda.hat is set to be 0 and log transformation will be used
    
    ###For details on the Chisq(df=1) test, See the slides on web---" Box-Cox Transformations: An Overview"
    ###by PENGFEI LI, Depart. of Statistics, Univ. of Connecticut, 04/11/2005
    
    ###xx,yy,zz: marker values from the three ordinal groups
    ###nboot:number of bootstrap samples
    ###lambda0 : the close-to-0 lambda used in the hypothesis testing

    
    ##these steps have been delt in Youden3Grp.boxcox.fit() function
    if(min(xx)<=0) new.xx <- xx+abs(min(xx))+1
    else new.xx <- xx

    if(min(yy)<=0) new.yy <- yy+abs(min(yy))+1
    else new.yy <- yy

    if(min(zz)<=0) new.zz <- zz+abs(min(zz))+1
    else new.zz <- zz


    fit0 <- Youden3Grp.boxcox.fit(object1=new.xx,object2=new.yy,object3=new.zz)
    
    lambda.est <- fit0$lambda[1]
    
    #if(!is.na(fit0$loglik) & !is.infinite(fit0$loglik))    lambda.est <- fit0$lambda[1]
    #else return(NA)
    
    #cat("lambda.est=",lambda.est,"\n")
    
    ########Determine whether use lambda.est as per se or view it equivalent as 0
    #######by bootstrap for small samples or by chisq test if large samples >=2000
    n1 <- length(new.xx)
    n2 <- length(new.yy)
    n3 <- length(new.zz)
    
    if(min(n1,n2,n3)>=2000)##asymptotic chisq(df=1) test
      {
        #loglik.h1 <- Youden3Grp.boxcox.Loglike.GivenLambda(dat1=new.xx,dat2=new.yy,dat3=new.zz,lambda.est=lambda.est)##loglike under alternative hypotheiss:lambda=lambda.est
        #loglik.h0 <- Youden3Grp.boxcox.Loglike.GivenLambda(dat1=new.xx,dat2=new.yy,dat3=new.zz,lambda.est=lambda0)##loglike under alternative hypotheiss:lambda=lambda0~0
        loglik.h1 <- negloglik.boxcox.3Grp(lambda.val=lambda.est, data1=new.xx, xmat1=as.matrix(rep(1,n1)),data2=new.yy,xmat2=as.matrix(rep(1,n2)),data3=new.zz,xmat3=as.matrix(rep(1,n3)), lik.method = "ML")
        loglik.h1 <- -loglik.h1
        loglik.h0 <- negloglik.boxcox.3Grp(lambda.val=lambda0, data1=new.xx, xmat1=as.matrix(rep(1,n1)),data2=new.yy,xmat2=as.matrix(rep(1,n2)),data3=new.zz,xmat3=as.matrix(rep(1,n3)), lik.method = "ML")
        loglik.h0 <- -loglik.h0
        
        #print(c(loglik.h1,loglik.h0))
        W <- 2*(loglik.h1-loglik.h0)
        
        p.value <- 1-pchisq(W,df=1,lower.tail=T)##higher W indicate larger likelihood ratio of lambda vs. lambda0, ->prefer lambda
        if(p.value>0.05) lambda.est <- 0 
      }
    else##do bootstrap sampling , if 95% CI cross zero,take zero
      {

        all.lambda <- sapply(1:nboot,function(jj)
               {
                 new.xx2 <- bootSample(new.xx,seed0=jj*10)
                 new.yy2 <- bootSample(new.yy,seed0=jj*10+1)
                 new.zz2 <- bootSample(new.zz,seed0=jj*10+2)
                 Youden3Grp.boxcox.fit(object1=new.xx2,object2=new.yy2,object3=new.zz2)$lambda
               })
        all.lambda <- all.lambda[1,]
        lower <- quantile(all.lambda,prob=0.025)
        upper <- quantile(all.lambda,prob=0.975)

        
        if(0>=lower && 0<=upper) lambda.est <- 0

      }  
    return(lambda.est)
          
  }

