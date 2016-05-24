
#
#  mapply internal function 
#
testfun <- function(m, n, x.o){
   
  #tempx.o=c()
  nvl <- 1:n
  lep <- rep(NA, n)
  rep <- rep(NA, n)
  
  i.th <- 1+m
  dx1 <- which(nvl <= i.th)
  dx2 <- setdiff(nvl, dx1)
  lep[dx1] <- rep(x.o[1], length(dx1))
  lep[dx2] <- x.o[dx2-m]  
  
  i.th2 <- n-m
  dx1 <- which(nvl >= i.th2)
  dx2 <- setdiff(nvl, dx1)
  rep[dx1] <-  rep(x.o[n], length(dx1))
    rep[dx2] <-  x.o[dx2+m]
  
  tempx.o <- 2*m/(n*(rep-lep))
  return(tempx.o)
}


#
# mapply function for equation2.16
#
testfun5 <- function(x,y,j,m,i){

  
  eqn <- mapply("equation210", j=j,  MoreArgs=list(x=x,y=y,m=m,i=i))
  if(i==1) n1 <- length(x)
  if(i==2) n1 <- length(y)  
  eqn <- (2*m)/(eqn*n1)
  return(prod(eqn))
  
}




#
# equation 2.10 in paper to be used in testfun5 for equation 2.16 
#  Vexler A, Yu J, Tian L, Liu S (2010). “Two-sample nonparametric likelihood inference based
#  on incomplete data with an application to a pneumonia study.” Biometrical Journal, 52(3),
#  348–361.


equation210 <- function(x,y,m,i,j){

  if( (i!=1) & (i!=2) ) stop(" i must be of value 1 or 2 \n")
  
  n1 <- length(x)
  n2 <- length(y)
  x.o1 <- sort(x)
  x.o2 <- sort(y)
  
  # first vector
  dx <- j+m
  if((j+m)<=0) dx <- 1
  if(i == 1){
    if((j+m)> n1) dx <- n1
    cut1p <- x.o1[dx]
  }else{
    if((j+m)> n2) dx <- n2
    cut1p <- x.o2[dx]
  }
  dx <- j-m
  if((j-m)<=0) dx <- 1
  if(i == 1){
    if((j-m)> n1) dx <- n1
    cut1n <- x.o1[dx]
  }else{
    if((j-m)> n2) dx <- n2
    cut1n <- x.o2[dx]
  }


  tosum <- c()
  for(k in 1:2){
         
    if(k == 1){

      if(i == 1){ 

        dx1 <- length(which(x.o1 <= cut1p))
        if(dx1!=n1){pt1 <- c(rep(1,dx1), rep(0,length((dx1+1):n1)))}else{pt1 <- c(rep(1,dx1))}
        dx2 = length(which(x.o1 <= cut1n))
        if(dx2!=n1){pt2 <- c(rep(1,dx2), rep(0,length((dx2+1):n1)))}else{pt2 <- c(rep(1,dx2))}
        vl1 <- sum((pt1-pt2))
              
      }else{# i is 2
        
        dx3 <- length(which(x.o1 <= cut1p))
        if(dx3!=n1){pt3 <- c(rep(1,dx3), rep(0,length((dx3+1):n1)))}else{pt3 <- c(rep(1,dx3))}
        dx4 <- length(which(x.o1 <= cut1n))
        if(dx4!=n1){pt4 <- c(rep(1,dx4), rep(0,length((dx4+1):n1)))}else{pt4 <- c(rep(1,dx4))}
        vl1 <- sum((pt3-pt4))
        
      }

      tosum[k] <- vl1
      
    }else{ # k is 2


      if(i == 1){ 

        dx1 <- length(which(x.o2 <= cut1p))
        if(dx1!=n2){pt1 <- c(rep(1,dx1), rep(0,length((dx1+1):n2)))}else{pt1 <- c(rep(1,dx1))}
        dx2 <- length(which(x.o2 <= cut1n))
        if(dx2!=n2){pt2 <- c(rep(1,dx2), rep(0,length((dx2+1):n2)))}else{pt2 <- c(rep(1,dx2))}
        vl2 <- sum((pt1-pt2))
              
      }else{# i is 2
        
        dx3 <- length(which(x.o2 <= cut1p))
        if(dx3!=n2){pt3 <- c(rep(1,dx3), rep(0,length((dx3+1):n2)))}else{pt3 <- c(rep(1,dx3))}
        dx4 <- length(which(x.o2 <= cut1n))
        if(dx4!=n2){pt4 <- c(rep(1,dx4), rep(0,length((dx4+1):n2)))}else{pt4 <- c(rep(1,dx4))}
        vl2 <- sum((pt3-pt4))        
      }

      tosum[k] <- vl2
    }
  
  }

  ans <- sum(tosum)/(n1+n2) 
  return(ans)
}



# used in getPval 
helperPval <- function(rowdx, colvls, teststat, testmat){


  if(rowdx == 1){
    if(teststat >= colvls[1]){
      pvl <- as.numeric(rownames(testmat))[1]
    }else{
      tempx <- c(colvls[1], colvls[2])
      tempy <- c(as.numeric(rownames(testmat))[1], as.numeric(rownames(testmat))[2])
      pvl <- NA
    }
  }
  if(rowdx == length(colvls)){
    if(teststat <= colvls[length(colvls)]){
      pvl <- as.numeric(rownames(testmat))[length(colvls)]
    }else{
      tempx <- c(colvls[(length(colvls)-1)], colvls[length(colvls)])
      tempy <- c(as.numeric(rownames(testmat))[(length(colvls)-1)], as.numeric(rownames(testmat))[length(colvls)])
      pvl <- NA
    }
  } 
  if(!exists("pvl")){
    
    if( (teststat > colvls[rowdx]) & (teststat < colvls[rowdx-1]) ){
      tempx <- c(colvls[rowdx-1], colvls[rowdx])
      tempy <- c(as.numeric(rownames(testmat))[rowdx-1], as.numeric(rownames(testmat))[rowdx])
      pvl <- NA
    }else{
      tempx <- c(colvls[rowdx], colvls[rowdx+1])
      tempy <- c(as.numeric(rownames(testmat))[rowdx], as.numeric(rownames(testmat))[rowdx+1])
      pvl <- NA
    }
  }
  if(is.na(pvl)){
    modelest <- lm(tempy~tempx, data=as.data.frame(list(tempy=tempy, tempx=tempx)))
    pvl <-   as.numeric(modelest$coefficients[1]) + (teststat*as.numeric(modelest$coefficients[2]))      
  }
  return(pvl)
}




#
# estimate pval by premade table
#
getPval <-function(testcall, teststat, x, y, vrb=TRUE){
  
  if(vrb) cat("estimating pvalue based on table \n")

  # slightly different if normal/uniform to equality so split 
  if( (testcall=="normal")  | (testcall=="uniform")  ){
    
    # load preset data matrix 
    if(testcall == "uniform"){
      #data(uniformData, envir=environment())
      testmat <- unifCut
    }
    if(testcall == "normal"){
      #data(normalData, envir=environment())
      testmat <- normCut
    }
    # get sample size and check if that was and estimated value 
    smp.size <- length(x)
    coldx <- which(colnames(testmat) == smp.size)

    # if it was an estimated value only working with one column and do not
    # need to interpolate between columns 
    if(length(coldx)>0){
      colvls <- testmat[,coldx]
      difvls <- colvls-teststat
      absvls <- abs(difvls)
      # is there an exact match return that value
      matchVl <- which(absvls == 0)
      if(length(matchVl >0)){
        pvl <- as.numeric(rownames(testmat)[matchVl])
      # estimate between rows
      }else{        
        rowdx <- which(absvls == min(absvls))
        pvl <- helperPval(rowdx, colvls, teststat, testmat)
      }
    # between columns as sample size was not exactly in data matrix      
    }else{
      
      if((smp.size < 10) | (smp.size>10000)){

        if(smp.size < 10){
          cat("Table generated on sample size between 10 and 10000. \n Estimating based on sample size of 10 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
          coldx <- 1
          colvls <- testmat[,coldx]
          difvls <- colvls-teststat
          absvls <- abs(difvls)
          # is there an exact match return that value
          matchVl <- which(absvls == 0)
          if(length(matchVl >0)){
            pvl <- as.numeric(rownames(testmat)[matchVl])
           # estimate between rows
          }else{        
            rowdx <- which(absvls == min(absvls))
            pvl <- helperPval(rowdx, colvls, teststat, testmat)
          }
        }else{
          cat("Table generated on sample size between 10 and 10000. \n Estimating based on sample size of 10000 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
          coldx <- 13
          colvls <- testmat[,coldx]
          difvls <- colvls-teststat
          absvls <- abs(difvls)
          # is there an exact match return that value
          matchVl <- which(absvls == 0)
          if(length(matchVl >0)){
            pvl <- as.numeric(rownames(testmat)[matchVl])
           # estimate between rows
          }else{        
            rowdx <- which(absvls == min(absvls))
            pvl <- helperPval(rowdx, colvls, teststat, testmat)
          }

        }
        
      }else{
        coldx1 <- max(which(as.numeric(colnames(testmat))<smp.size))
        coldx2 <- min(which(as.numeric(colnames(testmat))>smp.size))
        # first vls
        colvls <- testmat[,coldx1]
        difvls <- colvls-teststat
        absvls <- abs(difvls)
        # is there an exact match return that value
        matchVl <- which(absvls == 0)
        if(length(matchVl >0)){
          pvl1 <- as.numeric(rownames(testmat)[matchVl])
        # estimate between rows
        }else{
          rowdx <- which(absvls == min(absvls))
          pvl1 <- helperPval(rowdx, colvls, teststat, testmat)       
        }              
          
        # second value    
        colvls <- testmat[,coldx2]
        difvls <- colvls-teststat
        absvls <- abs(difvls)
        # is there an exact match return that value
        matchVl <- which(absvls == 0)
        if(length(matchVl >0)){
          pvl2 <- as.numeric(rownames(testmat)[matchVl])
        # estimate between rows
        }else{        
          rowdx <- which(absvls == min(absvls))
          pvl2 <- helperPval(rowdx, colvls, teststat, testmat)       
        }              
        tempx <- c(as.numeric(colnames(testmat)[coldx1]),as.numeric(colnames(testmat)[coldx2]))
        tempy <- c(pvl1, pvl2)
        datamat <- as.data.frame(cbind(tempx,tempy))
        modelest <- lm("tempy~tempx", data=datamat) 
        pvl <-  as.numeric(modelest$coefficients[1]) + (smp.size*as.numeric(modelest$coefficients[2]))
      }

      
    }
  }else{  # if not normal or uniform then distribution equality  
   
    #data(distributionEqualityData, envir=environment())
    testmat <- distEqCut
    # start with x length and estimate 
    # get sample size and check if that was and estimated value 
    smp.size <- length(x)
    coldx <- which(colnames(testmat) == smp.size)
    # if it was an estimated value only working with one column and do not
    # need to interpolate between columns 
    if(length(coldx)>0){
      colvls <- testmat[,coldx]
      difvls <- colvls-teststat
      absvls <- abs(difvls)
      # is there an exact match return that value
      matchVl <- which(absvls == 0)
      if(length(matchVl >0)){
        pvlx <- as.numeric(rownames(testmat)[matchVl])
      # estimate between rows
      }else{        
        rowdx <- which(absvls == min(absvls))
        pvlx <- helperPval(rowdx, colvls, teststat, testmat)
      }
    # between columns as sample size was not exactly in data matrix      
    }else{
      if((smp.size < 10) | (smp.size>500)){
        

        if(smp.size < 10){
          cat("Table generated on sample size between 10 and 500. \n Estimating based on sample size of 10 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
          coldx <- 1
          colvls <- testmat[,coldx]
          difvls <- colvls-teststat
          absvls <- abs(difvls)
          # is there an exact match return that value
          matchVl <- which(absvls == 0)
          if(length(matchVl >0)){
            pvl <- as.numeric(rownames(testmat)[matchVl])
           # estimate between rows
          }else{        
            rowdx <- which(absvls == min(absvls))
            pvl <- helperPval(rowdx, colvls, teststat, testmat)
          }
        }else{
          cat("Table generated on sample size between 10 and 500. \n Estimating based on sample size of 500 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
          coldx <- 8
          colvls <- testmat[,coldx]
          difvls <- colvls-teststat
          absvls <- abs(difvls)
          # is there an exact match return that value
          matchVl <- which(absvls == 0)
          if(length(matchVl >0)){
            pvl <- as.numeric(rownames(testmat)[matchVl])
           # estimate between rows
          }else{        
            rowdx <- which(absvls == min(absvls))
            pvl <- helperPval(rowdx, colvls, teststat, testmat)
          }

        }

        
      }else{
        coldx1 <- max(which(as.numeric(colnames(testmat))<smp.size))
        coldx2 <- min(which(as.numeric(colnames(testmat))>smp.size))
        # first vls
        colvls <- testmat[,coldx1]
        difvls <- colvls-teststat
        absvls <- abs(difvls)
        # is there an exact match return that value
        matchVl <- which(absvls == 0)
        if(length(matchVl >0)){
          pvl1 <- as.numeric(rownames(testmat)[matchVl])
        # estimate between rows
        }else{
          rowdx <- which(absvls == min(absvls))
          pvl1 <- helperPval(rowdx, colvls, teststat, testmat) 
        }
        # second value    
        colvls <- testmat[,coldx2]
        difvls <- colvls-teststat
        absvls <- abs(difvls)
        # is there an exact match return that value
        matchVl <- which(absvls == 0)
        if(length(matchVl >0)){
          pvl2 <- as.numeric(rownames(testmat)[matchVl])
        # estimate between rows
        }else{
          rowdx <- which(absvls == min(absvls))
          pvl2 <- helperPval(rowdx, colvls, teststat, testmat) 
        }
        tempx <- c(as.numeric(colnames(testmat)[coldx1]),as.numeric(colnames(testmat)[coldx2]))
        tempy <- c(pvl1, pvl2)
        datamat <- as.data.frame(cbind(tempx,tempy))
        modelest <- lm("tempy~tempx", data=datamat) 
        pvlx <-  as.numeric(modelest$coefficients[1]) + (smp.size*as.numeric(modelest$coefficients[2]))
      }
    }
    # Now if x and y are same length return pvlx as value
    # Else calculate pvl for y and average 
    if(length(y) == length(x)){ pvl = pvlx }else{
      smp.size <- length(y)
      coldx <- which(colnames(testmat) == smp.size)
      # if it was an estimated value only working with one column and do not
      # need to interpolate between columns 
      if(length(coldx)>0){
        colvls <- testmat[,coldx]
        difvls <- colvls-teststat
        absvls <- abs(difvls)
        # is there an exact match return that value
        matchVl <- which(absvls == 0)
        if(length(matchVl >0)){
          pvly <- as.numeric(rownames(testmat)[matchVl])
        # estimate between rows
        }else{        
          rowdx <- which(absvls == min(absvls))
          pvly <- helperPval(rowdx, colvls, teststat, testmat)                  
        }
      # between columns as sample size was not exactly in data matrix      
      }else{
        if((smp.size < 10) | (smp.size>500)){
      

          if(smp.size < 10){
            cat("Table generated on sample size between 10 and 500. \n Estimating based on sample size of 10 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
            coldx <- 1
            colvls <- testmat[,coldx]
            difvls <- colvls-teststat
            absvls <- abs(difvls)
                                        # is there an exact match return that value
            matchVl <- which(absvls == 0)
            if(length(matchVl >0)){
              pvl <- as.numeric(rownames(testmat)[matchVl])
                                        # estimate between rows
            }else{        
              rowdx <- which(absvls == min(absvls))
              pvl <- helperPval(rowdx, colvls, teststat, testmat)
            }
          }else{
            cat("Table generated on sample size between 10 and 500. \n Estimating based on sample size of 500 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
            coldx <- 8
            colvls <- testmat[,coldx]
            difvls <- colvls-teststat
            absvls <- abs(difvls)
                                        # is there an exact match return that value
            matchVl <- which(absvls == 0)
            if(length(matchVl >0)){
              pvl <- as.numeric(rownames(testmat)[matchVl])
                                        # estimate between rows
            }else{        
              rowdx <- which(absvls == min(absvls))
              pvl <- helperPval(rowdx, colvls, teststat, testmat)
            }

          }



          
        }else{
          coldx1 <- max(which(as.numeric(colnames(testmat))<smp.size))
          coldx2 <- min(which(as.numeric(colnames(testmat))>smp.size))
          # first vls
          colvls <- testmat[,coldx1]
          difvls <- colvls-teststat
          absvls <- abs(difvls)
          # is there an exact match return that value
          matchVl <- which(absvls == 0)
          if(length(matchVl >0)){
            pvl1 <- as.numeric(rownames(testmat)[matchVl])
          # estimate between rows
          }else{        
            rowdx <- which(absvls == min(absvls))
            pvl1 <- helperPval(rowdx, colvls, teststat, testmat)
          } 
          # second value    
          colvls <- testmat[,coldx2]
          difvls <- colvls-teststat
          absvls <- abs(difvls)
          # is there an exact match return that value
          matchVl <- which(absvls == 0)
          if(length(matchVl >0)){
            pvl2 <- as.numeric(rownames(testmat)[matchVl])
          # estimate between rows
          }else{
            rowdx <- which(absvls == min(absvls))
            pvl2 <-  helperPval(rowdx, colvls, teststat, testmat)       
          }
          tempx <- c(as.numeric(colnames(testmat)[coldx1]),as.numeric(colnames(testmat)[coldx2]))
          tempy <- c(pvl1, pvl2)
          datamat <- as.data.frame(cbind(tempx,tempy))
          modelest <- lm("tempy~tempx", data=datamat) 
          pvly <-  as.numeric(modelest$coefficients[1]) + (smp.size*as.numeric(modelest$coefficients[2]))
        }
      }
      pvl <- mean(c(pvlx, pvly))
    }
  }# else for normal/uniform vs distribution equality
  return(pvl)
}# end function 


#
# function to just return value not pvalue 
#


emplikeGOFinternal <- function(x, testcall=c("uniform", "normal"), delta=0.50){

  if( (testcall != "uniform")  & (testcall != "normal") ){

    stop("testcall not appropriately defined in function call \n  testcall should be either 'uniform' or 'normal' \n")
    

  }else{
    nx <- length(x)
    if(testcall=="uniform"){
    
      #
      # get test stat
      #
      n <- length(x)
      #tomin=c()
      x.o <- sort(x)
      m <- 1:n^(1-delta)
      tomin <- mapply("testfun", m=m, MoreArgs=list(x.o=x.o, n=n))
      #tomin =  apply(FUN="prod", MARGIN=2,tomin)
      #unif.stat = min(tomin)
      #output = unif.stat
      logTomin <- log(tomin)
      sumTomin <- apply(FUN="sum", MARGIN=2, logTomin)
      output <- min(sumTomin)


    }else{ # if not uniform than normal


      #
      # get test stats
      #
      n <- length(x)
      m1 <- 1:round(n/2)
      m2 <- 1:n^(1-delta)
      x.o <- sort(x)
      VN2 <- mapply("testfun", m=m2, MoreArgs=list(x.o=x.o, n=n))
      #VN2 =  apply(FUN="prod", MARGIN=2,VN2)
      #tominVN2 = (2*pi*exp(1)*var(x))^(n/2)*VN2
      #Vn2 = min(tominVN2)
      #output = Vn2
      logVN2 <- log(VN2)
      sumVN2 <- apply(FUN="sum", MARGIN=2, logVN2)
      output <- min((((n/2)*log(2*pi*exp(1)*var(x))) + sumVN2))

      
    }

    return(output)

  } # if  appropriate call type  
  
}

distEqInternal <-function(x,y,delta.equality=0.10){

  n1 <- length(x)
  n2 <- length(y)
  
  l1 <- floor(n1^(0.5+delta.equality))
  u1 <- ceiling(min((n1^(1-delta.equality)),(n1/2)))
  
  mvc <- seq(l1,u1,by=1)
  j1 <- 1:n1
  i <- 1
  #tomin1 = c()
  tomin1 <- mapply("testfun5", m=mvc, MoreArgs=list(x=x,y=y,i=i, j=j1))
  ans.pt1 <- min(tomin1)
  
  l2 <- floor(n2^(0.5+delta.equality))
  u2 <- ceiling(min((n2^(1-delta.equality)),(n2/2)))
  
  vvc <- seq(l2,u2, by=1)
  j2 <- 1:n2
  i <- 2
  #tomin2 = c()
  tomin2 <- mapply("testfun5", m=vvc, MoreArgs=list(x=x,y=y,i=i, j=j2))
  ans.pt2 <- min(tomin2)
  
  teststat <- (ans.pt1*ans.pt2)

  return(teststat)

}


#
#  mapply internal function 
#
testfun2 <- function(iter, sample.size, testcall, delta=0.5, delta.equality=0.1, random.seed.flag=FALSE){

  if(testcall=="uniform"){
    if(random.seed.flag) {
      set.seed(iter)
    }
    x <- runif(sample.size[1])
    vl <- emplikeGOFinternal(x, testcall, delta=delta)
    #z=log(vl)
    z <- vl
  }else{
    if(testcall=="normal"){
      if(random.seed.flag) {
        set.seed(iter)
      }    
      x <- rnorm(sample.size[1])
      vl <- emplikeGOFinternal(x, testcall, delta=delta)
      #z=log(vl)
      z <- vl
    }else{
      if(length(sample.size)==1) sample.size = c(sample.size,sample.size)
      if(random.seed.flag) {
        set.seed(iter)
      }
      x <- rnorm(sample.size[1])
      if(random.seed.flag) {
        set.seed(iter+1)
      }
      y <- rnorm(sample.size[2])
      vl <- distEqInternal(x,y,delta.equality=delta.equality)
      z <- log(vl)
    }
  }        
  return(z)
}


#
# return teststatistic/cutoff 
#

getCutoff <- function(testcall, smp.size, targetalpha){
   
  cat("estimating cutoff based on table \n")
  
  # slightly different if normal/uniform to equality so split 
  if( (testcall=="normal")  | (testcall=="uniform")  ){
    
    # load preset data matrix 
    if(testcall == "uniform"){
      #data(uniformData, envir=environment())
      testmat <- unifCut
    }
    if(testcall == "normal"){
      #data(normalData, envir=environment())
      testmat <- normCut
    }

    coldx <- which(colnames(testmat) == smp.size)
    rowdx <- which(rownames(testmat) == targetalpha)
    if(length(rowdx)==0){
      cat("Rounding target alpha: ", targetalpha, " to ")
      targetalpha <- round(targetalpha, digits=3)
      cat(targetalpha, ".\n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
    }
    rowdx <- which(rownames(testmat) == targetalpha)
    # if it was an estimated value only working with one column and do not
    # need to interpolate between columns 
    if(length(coldx)>0){
      cutoff <- testmat[rowdx,coldx]
    # between columns as sample size was not exactly in data matrix      
    }else{
      # if range issue
      if((smp.size < 10) | (smp.size>10000)){
        if(smp.size < 10){
          cat("Table generated on sample size between 10 and 10000. \n Estimating based on sample size of 10 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
          coldx <- 1
          cutoff <- testmat[rowdx,coldx]
        }else{
          cat("Table generated on sample size between 10 and 10000. \n Estimating based on sample size of 10000 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
          coldx <- 13
          cutoff <- testmat[rowdx,coldx]
        }        
      }else{# inbetween columns
        coldx1 <- max(which(as.numeric(colnames(testmat))<smp.size))
        coldx2 <- min(which(as.numeric(colnames(testmat))>smp.size))
        tempx <- as.numeric(colnames(testmat)[c(coldx1,coldx2)])
        tempy <- c(testmat[rowdx,coldx1], testmat[rowdx,coldx2])
        modelest <- lm(tempy~tempx, data=as.data.frame(list(tempy=tempy, tempx=tempx)))
        cutoff <- as.numeric(modelest$coefficients[1]) + (smp.size*as.numeric(modelest$coefficients[2])) 
      }    
    }    
  }else{# ends if normal or uniform 

    #data(distributionEqualityData, envir=environment())
    testmat <- distEqCut
    
    coldx <- which(colnames(testmat) == smp.size[1])
    rowdx <- which(rownames(testmat) == targetalpha)
    if(length(rowdx)==0){
      cat("Rounding target alpha: ", targetalpha, " to ")
      targetalpha <- round(targetalpha, digits=3)
      cat(targetalpha, ".\n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
    }
    rowdx <- which(rownames(testmat) == targetalpha)
    # if it was an estimated value only working with one column and do not
    # need to interpolate between columns 
    if(length(coldx)>0){
      cutoff1 <- testmat[rowdx,coldx]
    # between columns as sample size was not exactly in data matrix      
    }else{
      # if range issue
      if((smp.size[1] < 10) | (smp.size[1] > 500)){
        if(smp.size[1] < 10){
          cat("Table generated on sample size between 10 and 500. \n Estimating based on sample size of 10 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
          coldx <- 1
          cutoff1 <- testmat[rowdx,coldx]
        }else{
          cat("Table generated on sample size between 10 and 500. \n Estimating based on sample size of 500 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
          coldx <- 8
          cutoff1 <- testmat[rowdx,coldx]
        }        
      }else{# inbetween columns
        coldx1 <- max(which(as.numeric(colnames(testmat))<smp.size[1]))
        coldx2 <- min(which(as.numeric(colnames(testmat))>smp.size[1]))
        tempx <- as.numeric(colnames(testmat)[c(coldx1,coldx2)])
        tempy <- c(testmat[rowdx,coldx1], testmat[rowdx,coldx2])
        modelest <- lm(tempy~tempx, data=as.data.frame(list(tempy=tempy, tempx=tempx)))
        cutoff1 <- as.numeric(modelest$coefficients[1]) + (smp.size[1]*as.numeric(modelest$coefficients[2])) 
      }    
    }    
    #
    # if smp.size for both sets are the same return cutoff1
    # else
    # calculate second cutoff for varying smp.size and average
    if((length(smp.size) == 1) | (smp.size[1] == smp.size[2])){
      cutoff <- cutoff1
    }else{

      coldx <- which(colnames(testmat) == smp.size[2])
      # if it was an estimated value only working with one column and do not
      # need to interpolate between columns 
      if(length(coldx)>0){
        cutoff2 <- testmat[rowdx,coldx]
      # between columns as sample size was not exactly in data matrix      
      }else{
        # if range issue
        if((smp.size[2] < 10) | (smp.size[2] > 500)){
          if(smp.size[2] < 10){
            cat("Table generated on sample size between 10 and 500. \n Estimating based on sample size of 10 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
            coldx <- 1
            cutoff2 <- testmat[rowdx,coldx]
          }else{
            cat("Table generated on sample size between 10 and 500. \n Estimating based on sample size of 500 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
            coldx <- 8
            cutoff2 <- testmat[rowdx,coldx]
          }        
        }else{# inbetween columns
          coldx1 <- max(which(as.numeric(colnames(testmat))<smp.size[2]))
          coldx2 <- min(which(as.numeric(colnames(testmat))>smp.size[2]))
          tempx <- as.numeric(colnames(testmat)[c(coldx1,coldx2)])
          tempy <- c(testmat[rowdx,coldx1], testmat[rowdx,coldx2])
          modelest <- lm(tempy~tempx, data=as.data.frame(list(tempy=tempy, tempx=tempx)))
          cutoff2 <- as.numeric(modelest$coefficients[1]) + (smp.size[2]*as.numeric(modelest$coefficients[2])) 
        }    
      }
      cutoff <- mean(c(cutoff1, cutoff2))
      
    }# ends if same length 
  }# ends if distribution equality 
  return(cutoff)
 
}# ends getCutoff
