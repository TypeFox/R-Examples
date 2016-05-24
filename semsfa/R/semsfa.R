semsfa<-function(formula, data=list(), sem.method="gam", var.method="fan", ineffDecrease=TRUE, tol = 1e-05, n.boot=0,  ...){     
    returnObj<-list()
    dati<-data
    if (!is.data.frame(data)) {
       stop("Data shoud be a data frame object")
    }
    
    if (sem.method!="gam" & sem.method!="kernel" & sem.method!="loess") {
       stop("admitted values for sem.method are 'gam' 'kernel' or 'loess' ")
    }
    if (var.method!="fan" & var.method!="mm") {
       stop("admitted values for var.method are 'fan' or 'mm'")
    }
    
    if(sem.method=="gam"){#generalized additive model
      reg <- gam(formula=formula,data=data,...)
      df <- sum(reg$edf)
      resp<-reg$model[,1]
    }

    if(sem.method=="kernel"){#kernel regression  
      dati<-data.frame(data)
      #reg<-npregbw(formula=formula,data=dati,...)
      reg<-npregbw(formula=formula,data=dati,regtype="ll",bwmethod="cv.aic",...) 
      reg<-npreg(reg)  
      df <-reg$ndim
      resp<-residuals(reg)+fitted(reg)
    }

    if(sem.method=="loess"){#local polynomial regression
      reg <- loess(formula=formula,data=data,...)
      df <- reg$enp
      resp<-reg$y
    }
             
    if (ineffDecrease & skewness(residuals(reg))>0) {
         warning("the residuals of the FirstStep estimates are left skewed \n this might indicate that there is no inefficiency or that the model is misspecified")
    }    
    if (!ineffDecrease & skewness(residuals(reg))<0) {
         warning("the residuals of the FirstStep estimates are right skewed \n this might indicate that there is no inefficiency or that the model is misspecified")
    }    

    reg.fitted<-fitted(reg)
    if(var.method=="fan"){
        sec.step<-optimize(f=fan,ineffD=ineffDecrease,resp=resp,Ey=reg.fitted,interval=c(0, 10), tol = tol,maximum=TRUE)
        lambda<-sec.step$maximum
        sigma<-sqrt(mean((resp-reg.fitted)^2)/(1-(2*lambda^2)/(pi*(1+lambda^2))))
        mu=sqrt(2)*sigma*lambda/sqrt(pi*(1+lambda^2))
        logL<-sec.step$objective
        }
    if(var.method=="mm"){
        if(ineffDecrease){sigma_sq_u=(moment(residuals(reg),order=3)/(sqrt(2/pi)*(1-4/pi)))^(2/3)}else{
        sigma_sq_u=(moment(residuals(reg),order=3)/(sqrt(2/pi)*(4/pi-1)))^(2/3)}
        if(sigma_sq_u<0){stop("efficiency variance estimate is negative")}
                else{
                sigma_sq_v=moment(residuals(reg),order=2)- (1-2/pi)*sigma_sq_u
                lambda<- sqrt(sigma_sq_u/sigma_sq_v)
                sigma<-sqrt(sigma_sq_u + sigma_sq_v)
                mu=sqrt(sigma_sq_u*2/pi)
                }
        if(ineffDecrease){
        logL<-+(length(resp)/2)*log(2/pi)-length(resp)*log(sigma) + sum(pnorm(-(resp-reg.fitted-mu)*lambda/sigma,log.p=TRUE))
          -(1/(2*sigma^2))*sum((resp-reg.fitted-mu)^2)
        }else{
        logL<- +(length(resp)/2)*log(2/pi)-length(resp)*log(sigma) + sum(pnorm((resp-reg.fitted+mu)*lambda/sigma,log.p=TRUE))
          -(1/(2*sigma^2))*sum((resp-reg.fitted+mu)^2)}
        
    }
    
    if(ineffDecrease){fitted=reg.fitted+mu}else{fitted=reg.fitted-mu}
        
    returnObj$formula<-formula
    returnObj$y<-resp
    returnObj$data<-data 
    returnObj$call <- match.call() 
    returnObj$sem.method<-sem.method
    returnObj$var.method<-var.method
    returnObj$ineffDecrease<-ineffDecrease
    returnObj$reg<-reg
    returnObj$reg.fitted<-reg.fitted
    returnObj$regkewness <- skewness(residuals(reg))
    returnObj$lambda<-lambda
    returnObj$sigma <- sigma
    returnObj$fitted<-fitted
    returnObj$tol <-tol
    returnObj$residual.df<- length(resp)-df -2   
    returnObj$bic <-  -2*logL +log(length(resp))*(df+2)    

    returnObj$n.boot<-n.boot    
    if(n.boot>0){
        sup <- function(w) if( any( grepl( "the residuals of the FirstStep", w) ) ) invokeRestart( "muffleWarning" ) 
        r <- foreach(icount(n.boot), .combine=rbind) %do% {
        ind <- sample(nrow(data),size=nrow(data),replace=TRUE)
        boot<- withCallingHandlers(semsfa(formula=formula,data=data[ind,],
               sem.method=sem.method,var.method=var.method,ineffDecrease=ineffDecrease,tot=tol, ...), 
               warning = sup)
        c(boot$lambda,boot$sigma)
        }
        returnObj$boot.mat<-r
        returnObj$b.se<-apply(r,2,sd)
    }
    class(returnObj) <- "semsfa"
    return(returnObj)
}
