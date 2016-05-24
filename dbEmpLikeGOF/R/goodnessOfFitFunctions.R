
#
# returns  cutoff given target alpha 
#

returnCutoff <-function(sample.size,testcall=c("uniform", "normal", "distribution.equality"), targetalpha=.05, num.mc=200, delta=0.5, delta.equality=0.10, pvl.Table=FALSE, random.seed.flag=TRUE){

  if( (testcall != "uniform")  & (testcall != "normal") & (testcall != "distribution.equality")){
    stop("testcall not appropriately defined in function call \n  testcall should be either 'uniform', 'normal', or 'distribution.equality' \n")
  }else{
    # until table is produced calucluate then add this option
    if(pvl.Table){
      cut.out <- getCutoff(testcall, sample.size, targetalpha)
    }else{ # ends if estimated by table
      
      #myalpha<-c()
      #cut.out<-c()    

      vec <- 1:num.mc
      z <- mapply(FUN="testfun2",vec, MoreArgs=list(sample.size=sample.size, testcall=testcall,delta=delta,delta.equality=delta.equality, random.seed.flag=random.seed.flag))
      cut.out<-quantile(z,(1-targetalpha), na.rm=TRUE)
    }# ends else calculate
    return(cut.out)
  }# end if testcall appropriate   
}# end function



#
# Goodness of Fit function
#

dbEmpLikeGOF <- function(x, y=NA, testcall=c("uniform", "normal"), delta=0.50, delta.equality=0.10,num.mc=1000, pvl.Table=TRUE, vrb=TRUE, random.seed.flag=TRUE){


  # if comparing to normal or uniform 
  if(is.na(y[1])){
  
    if( (testcall != "uniform")  & (testcall != "normal") ){
      
      stop("testcall not appropriately defined in function call \n  testcall should be either 'uniform' or 'normal' \n")
      

    }else{
      
      nx<-length(x)

      if(testcall=="uniform"){
    
      #
      # get test stat
      #
        if(vrb) cat("\n...Working on teststat \n") 
        n<-length(x)
        
        #tomin<-c()

        x.o<-sort(x)
        m <- 1:n^(1-delta)
        tomin <- mapply(FUN="testfun", m=m, MoreArgs=list(x.o=x.o, n=n))
        
        #tomin =  apply(FUN="prod", MARGIN=2,tomin)
        #teststat = log(min(tomin))

        logTomin <- log(tomin)
        sumTomin <- apply(FUN="sum", MARGIN=2, logTomin)
        teststat <- min(sumTomin)


        
        
      #
      # get or estimate p value
      #
        if(pvl.Table){
          if(delta != .5) cat("Table calculated on a delta = 0.5 \n    To use delta=",delta," please change pvl.Table to FALSE \n")
          pvl <- getPval(testcall, teststat, x, y, vrb)          
        }else{        
          if(vrb) cat("...Working on p-value \n")
          
          mc.est<-c()
          #tomin<-c()

          for(j in 1:num.mc){
            
            if(random.seed.flag) {
              set.seed(j)
            }    
            x<-runif(nx)
            x.o<-sort(x)
            m <- 1:nx^(1-delta)
            tomin <- mapply(FUN="testfun", m=m, MoreArgs=list(x.o=x.o, n=nx))
            #tomin =  apply(FUN="prod", MARGIN=2,tomin)
            #mc.est[j] = log(min(tomin))
            logTomin <- log(tomin)
            sumTomin <- apply(FUN="sum", MARGIN=2, logTomin)
            mc.est[j] <- min(sumTomin)
          } #end j loop
  
        # count the number of observations on either end of spectrum
          estp.rh<-sum(mc.est>teststat)/num.mc
          #estp.lh=sum(mc.est<teststat)/num.mc
          #pvl = min(estp.rh,estp.lh)
          pvl <- estp.rh 
        }
        output <- list(teststat=teststat, pvalue=pvl)
        

      }else{ # if not uniform than normal


      #
      # get test stats
      #
        if(vrb) cat("\n...Working on teststat \n") 
        n<-length(x)
        m1 <- 1:round(n/2)
        m2 <- 1:n^(1-delta)
        x.o<-sort(x)
        VN2 <- mapply(FUN="testfun", m=m2, MoreArgs=list(x.o=x.o, n=n))
        #VN2 =  apply(FUN="prod", MARGIN=2,VN2)
        #tominVN2 = (2*pi*exp(1)*var(x))^(n/2)*VN2
        #teststat = log(min(tominVN2))
        logVN2 <- log(VN2)
        sumVN2 <- apply(FUN="sum", MARGIN=2, logVN2)
        teststat <- min((((n/2)*log(2*pi*exp(1)*var(x))) + sumVN2))

        
      #
      # get or estimate p value
      #
        if(pvl.Table){
          if(delta != .5) cat("Table calculated on a delta = 0.5 \n    To use delta=",delta," please change pvl.Table to FALSE \n")
          pvl <- getPval(testcall, teststat, x, y, vrb)          
        }else{        
          if(vrb) cat("...Working on p-value \n")

          #mc.estVN1<-mc.estVN2<-c()
          #tominVN1<-tominVN2<-c()
          mc.estVN2<-c()
          
          for(j in 1:num.mc){
            if(random.seed.flag) {
              set.seed(j)
            }    
            x<-rnorm(nx)
            
            m2 <- 1:nx^(1-delta)
            x.o<-sort(x)
            VN2 <- mapply(FUN="testfun", m=m2, MoreArgs=list(x.o=x.o, n=nx))
            
            #VN2 =  apply(FUN="prod", MARGIN=2,VN2)
            #tominVN2 = (2*pi*exp(1)*var(x))^(nx/2)*VN2
            #  mc.estVN2[j] = min(tominVN2)

            logVN2 <- log(VN2)
            sumVN2 <- apply(FUN="sum", MARGIN=2, logVN2)
            mc.estVN2[j] <- min((((n/2)*log(2*pi*exp(1)*var(x))) + sumVN2))
            #mc.estVN2[j] = log(min(tominVN2))
            
          }
          estp.rh<-sum(mc.estVN2>teststat)/num.mc
          ###estp.lh=sum(mc.estVN2<teststat)/num.mc
          ##pvl = min(estp.rh,estp.lh)
          pvl <- estp.rh
        }
        output <- list(teststat=teststat, pvalue=pvl)
        
      }
      
     } # if  appropriate call type  

  }else{ # if y is compared to uniform or normal else compare to another vector 
    
    if( (delta.equality>0.25) & (delta.equality<0) ) stop("delta.equality must be between 0.0 and 0.25")

    testcall<-"distribution.equality"
    
    #
    # get test stats
    #
    if(vrb) cat("\n...Working on teststat \n") 
    n1 <- length(x)
    n2 <- length(y)
  
    l1 <- floor(n1^(0.5+delta.equality))
    u1 <- ceiling(min((n1^(1-delta.equality)),(n1/2)))

    mvc <- seq(l1,u1,by=1)
    j1 <- 1:n1
    i<-1
    #tomin1 <- c()
    tomin1 <- mapply(FUN="testfun5", m=mvc, MoreArgs=list(x=x,y=y,i=i, j=j1))
    ans.pt1 <- min(tomin1)
        
    l2 <- floor(n2^(0.5+delta.equality))
    u2 <- ceiling(min((n2^(1-delta.equality)),(n2/2)))
  
    vvc <- seq(l2,u2, by=1)
    j2 <- 1:n2
    i<-2
    #tomin2 <- c()
    tomin2 <- mapply(FUN="testfun5", m=vvc, MoreArgs=list(x=x,y=y,i=i, j=j2))
    ans.pt2 <- min(tomin2)

    teststat <- log((ans.pt1*ans.pt2))

 
    #
    # get or estimate p value
    #
    if(pvl.Table){
      if(delta.equality != .1) cat("Table calculated on a delta = 0.1 \n    To use delta.equality=",delta.equality," please change pvl.Table to FALSE \n")
       pvl <- getPval(testcall, teststat, x, y, vrb)          
    }else{        
      if(vrb) cat("...Working on p-value \n")
      toMinpvec <- c()
      for(j in 1:num.mc){
        if(random.seed.flag) {
          set.seed(j)
        }    
        v1<-rnorm(n1,mean=0,sd=1)
        if(random.seed.flag) {
          set.seed(j+1)
        }    
        v2<-rnorm(n2,mean=0,sd=1)
        l1 <- floor(n1^(0.5+delta.equality))
        u1 <- ceiling(min((n1^(1-delta.equality)),(n1/2)))
        mvc <- seq(l1,u1,by=1)
        j1 <- 1:n1
        i<-1
        #tomin1 <- c()
        tomin1 <- mapply(FUN="testfun5", m=mvc, MoreArgs=list(x=v1,y=v2,i=i, j=j1))
        ans.pt1 <- min(tomin1)
        l2 <- floor(n2^(0.5+delta.equality))
        u2 <- ceiling(min((n2^(1-delta.equality)),(n2/2)))
        vvc <- seq(l2,u2, by=1)
        j2 <- 1:n2
        i<-2
        #tomin2 <- c()
        tomin2 <- mapply(FUN="testfun5", m=vvc, MoreArgs=list(x=v1,y=v2,i=i, j=j2))
        ans.pt2 <- min(tomin2)
        toMinpvec[j] <- log((ans.pt1*ans.pt2))     
      }
      estp.rh<-sum(toMinpvec>teststat)/num.mc
      #estp.lh=sum(toMinpvec<teststat)/num.mc
      #pvl = min(estp.rh,estp.lh)
      pvl <- estp.rh
    }
    
    output <- list(teststat=teststat, pvalue=pvl)
   }
  
  return(output)
 
}
  


