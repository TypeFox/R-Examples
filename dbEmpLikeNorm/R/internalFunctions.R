# delta must be between 0 and 1
computeValue <- function(Lst, delta){


###########################
#
# Term A for all K
#
###########################


  # mean of each k 
  avg = unlist(lapply(Lst, FUN="mean"))
  # length of each k
  lenVec = unlist(lapply(Lst, FUN="length"))
  # initialize vector
  sigmaSq=rep(NA,length(avg))
  # for each k compute estimated sigma squared value
  for(le in 1:length(Lst)){
    sigmaSq[le] = sum(((Lst[[le]] - avg[le])^2)/length(Lst[[le]]))
  }
  # computer term A for each k
  #termA = (2*pi*exp(1)*sigmaSq)^(lenVec/2)
  termLogA = (lenVec/2)*log((2*pi*exp(1)*sigmaSq))

  # sort each k in  list 
  Lst = lapply(Lst, FUN="sort")


###########################
#
# Term B for all K
#
###########################

  termB = rep(NA,length(avg))

  for(k in 1:length(Lst)){

    # first determine which values of m
    # m must be integer
    vecM = 1:floor((lenVec[k]^(1.0-delta)))

    matN = matrix(NA,nrow=length(vecM), ncol=lenVec[k])
    # loop of m values 
    for(mdx in 1:length(vecM)){
      m = vecM[mdx]
      j = 1:lenVec[k]
      # find j+m index 
      ldx = j+m
      ldx[which(ldx>lenVec[k])] = lenVec[k]
      # find j-m index 
      udx = j-m
      udx[which(udx<1)] = 1

      # compute the value for each n in k for given m 
      #matN[mdx,] = (2*m)/((Lst[[k]][ldx]-Lst[[k]][udx])*lenVec[k])

      matN[mdx,] = exp(log((2*m)) - log(((Lst[[k]][ldx]-Lst[[k]][udx])*lenVec[k])))
      
    }
    # take product across across rows (for each m)
    # find minumum 
    termB[k] = min(apply(matN, MARGIN=1,FUN="prod"))
  }

###########################
#
# finish statistic 
#
###########################

  #valK = termA+termB
  valK = termLogA+log(termB)
  #T = prod(valK)
  T = sum(valK)
  #return(log(T))
  return(T)
  
}# end computeValue




computePval <- function(teststat, Lst, vrb=TRUE){

   if(vrb) cat("estimating pvalue based on table \n")
   
   len = length(Lst)
   if(len > 3){
     if(vrb) cat("Exiting: tables based on two or three group sets of equal sample size \n    Please run using pvl.Table=FALSE \n")
     pvl = NA
   }else{ 
     if(len == 2){
       #data(twoGroup)
       testmat = twoMat
     }else{
       if(len == 3){
         #data(threeGroup)
         testmat = threeMat
       }
     }
     if(vrb) cat("Tables based on sets of equal sample size \n    Continuing with sample size of first element in list \n") 
     smp.size = length(Lst[[1]])
     coldx = which(colnames(testmat) == smp.size)

     # if it was an estimated value only working with one column and do not
     # need to interpolate between columns 
     if(length(coldx)>0){
       colvls = testmat[,coldx]
       difvls = colvls-teststat
       absvls = abs(difvls)
       # is there an exact match return that value
       matchVl = which(absvls == 0)
       if(length(matchVl) >0){
         pvl = as.numeric(rownames(testmat)[matchVl])
       # estimate between rows
       }else{        
         rowdx = which(absvls == min(absvls))
         pvl = helperPval(rowdx, colvls, teststat, testmat)
       }
     # between columns as sample size was not exactly in data matrix      
     }else{

       if((smp.size < 10) | (smp.size>300)){
         if(smp.size < 10){
           cat("Table generated on sample size between 10 and 300. \n Estimating based on sample size of 10 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
           coldx = 1
           colvls = testmat[,coldx]
           difvls = colvls-teststat
           absvls = abs(difvls)
           # is there an exact match return that value
           matchVl = which(absvls == 0)
           if(length(matchVl) >0){
             pvl = as.numeric(rownames(testmat)[matchVl])
           # estimate between rows
           }else{        
             rowdx = which(absvls == min(absvls))
             pvl = helperPval(rowdx, colvls, teststat, testmat)
           }
         }else{
           cat("Table generated on sample size between 10 and 300. \n Estimating based on sample size of 300 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
           coldx = 13
           colvls = testmat[,coldx]
           difvls = colvls-teststat
           absvls = abs(difvls)
           # is there an exact match return that value
           matchVl = which(absvls == 0)
           if(length(matchVl) >0){
             pvl = as.numeric(rownames(testmat)[matchVl])
            # estimate between rows
           }else{        
             rowdx = which(absvls == min(absvls))
             pvl = helperPval(rowdx, colvls, teststat, testmat)
           }
         }
                   
       }else{
         coldx1 = max(which(as.numeric(colnames(testmat))<smp.size))
         coldx2 = min(which(as.numeric(colnames(testmat))>smp.size))
         # first vls
         colvls = testmat[,coldx1]
         difvls = colvls-teststat
         absvls = abs(difvls)
         # is there an exact match return that value
         matchVl = which(absvls == 0)
         if(length(matchVl) >0){
           pvl1 = as.numeric(rownames(testmat)[matchVl])
         # estimate between rows
         }else{
           rowdx = which(absvls == min(absvls))
           pvl1 = helperPval(rowdx, colvls, teststat, testmat)       
         }              
          
         # second value    
         colvls = testmat[,coldx2]
         difvls = colvls-teststat
         absvls = abs(difvls)
         # is there an exact match return that value
         matchVl = which(absvls == 0)
         if(length(matchVl) >0){
           pvl2 = as.numeric(rownames(testmat)[matchVl])
         # estimate between rows
         }else{        
           rowdx = which(absvls == min(absvls))
           pvl2 = helperPval(rowdx, colvls, teststat, testmat)       
         }              
         tempx = c(as.numeric(colnames(testmat)[coldx1]),as.numeric(colnames(testmat)[coldx2]))
         tempy = c(pvl1, pvl2)
         datamat = as.data.frame(cbind(tempx,tempy))
         modelest = lm("tempy~tempx", data=datamat) 
         pvl =  as.numeric(modelest$coefficients[1]) + (smp.size*as.numeric(modelest$coefficients[2]))
       }
     }      
   }# end if 2 or 3 group
   
   return(pvl)
   
 }




getCutoffValue <- function(numberOfgroups, sample.size, targetalpha){

  if(numberOfgroups == 2){
    #data(twoGroup)
    testmat = twoMat  
  }else{
    #data(threeGroup)
    testmat = threeMat  
  }

  if(length(sample.size) >1){
    cat("Tables are based on equal length groups. Continuing using the first sample size \n   For a more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
  }
  smp.size = sample.size[1]
  
  coldx = which(colnames(testmat) == as.character(smp.size))
  rowdx = which(rownames(testmat) == as.character(targetalpha))
  if(length(rowdx)==0){
    cat("Rounding target alpha: ", targetalpha, " to ")
    targetalpha = round(targetalpha, digits=3)
    if(targetalpha == 0) targetalpha = .001
    if(targetalpha == 1) targetalpha = .999
    cat(targetalpha, ".\n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
  }
  rowdx = which(rownames(testmat) == as.character(targetalpha))
    # if it was an estimated value only working with one column and do not
    # need to interpolate between columns 
  if(length(coldx)>0){
    cutoff = testmat[rowdx,coldx]
    # between columns as sample size was not exactly in data matrix      
  }else{
                                        # if range issue
    if((smp.size < 10) | (smp.size>300)){
      if(smp.size < 10){
        cat("Table generated on sample size between 10 and 300. \n Estimating based on sample size of 10 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
        coldx = 1
        cutoff = testmat[rowdx,coldx]
      }else{
        cat("Table generated on sample size between 10 and 300. \n Estimating based on sample size of 300 \n For more accurate calculation use monte carlo method [pvl.Table=FALSE] \n")
        coldx = 13
        cutoff = testmat[rowdx,coldx]
      }        
    }else{# inbetween columns
      coldx1 = max(which(as.numeric(colnames(testmat))<smp.size))
      coldx2 = min(which(as.numeric(colnames(testmat))>smp.size))
      tempx = as.numeric(colnames(testmat)[c(coldx1,coldx2)])
      tempy = c(testmat[rowdx,coldx1], testmat[rowdx,coldx2])
      modelest = lm(tempy~tempx, data=as.data.frame(list(tempy=tempy, tempx=tempx)))
      cutoff = as.numeric(modelest$coefficients[1]) + (smp.size*as.numeric(modelest$coefficients[2])) 
    }    
  }   
 
  return(cutoff)
}






####################################
####################################
####################################
####################################





BayesWays <- function(numberOfgroups, sample.size, alpha, nsims=200, v.threshold=NA){

  # this already should have been checked in returnCutoffValue
  #if(numberOfgroups > 3){
  #  cat("Tables are based on 2 or 3 groups. To get a cutoff for greater than 3 groups change [pvl.Table=FALSE]  \n")
  #  numberOfgroups = 3
  #}
  # fix sample.size list if not equal to number of groups
  if(length(sample.size)<numberOfgroups) {
    dif <- numberOfgroups-length(sample.size)
    sample.size <- c(sample.size,rep(sample.size[length(sample.size)], dif))
  }
  
  # Part 1 :  find value in table and estimate
  if(numberOfgroups == 2){
    part1 = getCutoffTableTwo(sample.size, alpha)
  }
  if(numberOfgroups == 3){
    part1 = getCutoffTableThree(sample.size, alpha)
  }
  
  # set threshold for varience 
  if(is.na(v.threshold)) v.threshold = part1[2]

  # keep original as nsims will change in loop
  o.nsims = nsims

  # Step 2: compute qhat
  vec.teststats = helperGetTstat(numberOfgroups, sample.size, nsims)
  qhat = computeQhat(vec.teststats, alpha, nsims, part1)
  
  # Step 3 compute variance estimator
  veTemp = varianceEstimator(numberOfgroups, sample.size, nsims, qhat, o.nsims)
  new.teststats = veTemp[[1]]
  var.est = veTemp[[2]]
  
  while((nsims<=35000) & (var.est>v.threshold)){
   
    # Step 4 compute new qhat 
    nsims = length(vec.teststats)+length(new.teststats)
    vec.teststats = c(vec.teststats, new.teststats)
    new.qhat = computeQhat(vec.teststats, alpha, nsims, part1)

    veTemp = varianceEstimator(numberOfgroups, sample.size, nsims, new.qhat, o.nsims)
    new.teststats = veTemp[[1]]
    var.est = veTemp[[2]]
   
  }
  
  if((nsims>35000) & (var.est>v.threshold)){
    cat("Exceded simulation limits without reaching acceptable variance \n")
    qhatC = NA
  }else{  
    # compute final qhat
    nsims = length(vec.teststats)+length(new.teststats)
    vec.teststats = c(vec.teststats, new.teststats)
    qhatC = computeQhat(vec.teststats, alpha, nsims, part1)
  }

  return(c(qhatC, var.est))
 
}





#
# Set for 2 groups
#
getCutoffTableTwo <- function(sample.size, alpha){
  
  #load("dataTable/twoGroup.RData")
  #data(twoGroup)
  mat <- twoMat
  cnames <- colnames(mat)
  rnames <- rownames(mat)

  n1 = sample.size[1]
  n2 = sample.size[2]
  
  NMequal = FALSE
  # determin if using same n and m index
  if(n1 == n2){
    NMequal = TRUE
  }else{
    # may not be equal but may use same index
    temp1 = which(as.numeric(cnames)<=n1)
    temp2 = which(as.numeric(cnames)>n1)
    coldx =  c(temp1[length(temp1)], temp2[1])

    temp3 = which(as.numeric(cnames)<=n2)
    temp4 = which(as.numeric(cnames)>n2)
    coldx2 =  c(temp3[length(temp3)], temp4[1])

    comp = (coldx==coldx2)
    if(length(comp[comp==FALSE])==0) NMequal = TRUE
  }
  
  q=n=m=alpha2=c()

  # if n and m are equal 
  if(NMequal){

    temp1 = which(as.numeric(rnames)<=alpha)
    temp2 = which(as.numeric(rnames)>alpha)
    rowdx =  c(temp1[length(temp1)], temp2[1])
    # adjust index if below or above minimum or maximum preset alpha values (ie. <.001 or >.999)
    if((length(rowdx)==1) & (rowdx[1]==1)) rowdx = c(1,2)
    if((length(rowdx)==2) & (is.na(rowdx[2]))) rowdx = c((length(rnames)-1),length(rnames))

    temp3 = which(as.numeric(cnames)<=n1)
    temp4 = which(as.numeric(cnames)>n1)
    coldx =  c(temp3[length(temp3)], temp4[1])
    # adjust index if below or above minimum or maximum preset column sample size (ie. <10 or >300)
    if((length(coldx)==1) & (coldx[1]==1)) coldx = c(1,2)
    if((length(coldx)==2) & (is.na(coldx[2]))) coldx = c((length(cnames)-1),length(cnames))

    # adjust index outward to collect 16 data points 
    # if on a boundary adjust to other side
    idx = 1
    maxD = rowdx[1]-idx
    maxU = rowdx[2]+idx
    if(maxD<=0){
      if(maxD==0){
        maxU=maxU+1
      }else{
        maxU = maxU+(abs(maxD)+1)
      }
      maxD = 1           
    }
    if(maxU>length(rnames)){
      extraDif  = maxU-length(rnames)
      maxD = maxD-extraDif
      maxU = length(rnames)
    }
    rowdx = maxD:maxU
    maxL = coldx[1]-idx
    maxR = coldx[2]+idx
    if(maxL<=0){
      if(maxL==0){
        maxR=maxR+1
      }else{
        maxR = maxR+(abs(maxL)+1)
      }
      maxL = 1           
    }
    if(maxR>length(cnames)){
      extraDif  = maxR-length(cnames)
      maxL = maxL-extraDif
      maxR = length(cnames)
      }
    coldx = maxL:maxR

    tempMat = mat[rowdx,coldx]
    for(i in 1:dim(tempMat)[2]){
      q = c(q,as.numeric(tempMat[,i]))
      n = c(n,rep(as.numeric(colnames(tempMat)[i]),dim(tempMat)[1]))
      alpha2 = c(alpha2,as.numeric(rownames(tempMat)))
    }


  #
  # ends if n=m
  # now if n!=m
  #
  }else{

    # find rows = alpha
    temp1 = which(as.numeric(rnames)<=alpha)
    temp2 = which(as.numeric(rnames)>alpha)
    rowdx =  c(temp1[length(temp1)], temp2[1])
    # adjust index if below or above minimum or maximum preset alpha values (ie. <.001 or >.999)
    if((length(rowdx)==1) & (rowdx[1]==1)) rowdx = c(1,2)
    if((length(rowdx)==2) & (is.na(rowdx[2]))) rowdx = c((length(rnames)-1),length(rnames))
    # adjust index outward  
    # if on a boundary adjust to other side
    idx = 1
    maxD = rowdx[1]-idx
    maxU = rowdx[2]+idx
    if(maxD<=0){
      if(maxD==0){
        maxU=maxU+1
      }else{
        maxU = maxU+(abs(maxD)+1)
      }
      maxD = 1           
    }
    if(maxU>length(rnames)){
      extraDif  = maxU-length(rnames)
      maxD = maxD-extraDif
      maxU = length(rnames)
    }
    rowdx = maxD:maxU
      
   
    # work on n1 = n
    temp3 = which(as.numeric(cnames)<=n1)
    temp4 = which(as.numeric(cnames)>n1)
    coldx =  c(temp3[length(temp3)], temp4[1])
    # adjust index if below or above minimum or maximum preset column sample size (ie. <10 or >300)
    if((length(coldx)==1) & (coldx[1]==1)) coldx = c(1,2)
    if((length(coldx)==2) & (is.na(coldx[2]))) coldx = c((length(cnames)-1),length(cnames))

    # adjust index outward to collect 16 data points 
    # if on a boundary adjust to other side
    maxL = coldx[1]-idx
    maxR = coldx[2]+idx
    if(maxL<=0){
      if(maxL==0){
        maxR=maxR+1
      }else{
        maxR = maxR+(abs(maxL)+1)
      }
      maxL = 1           
    }
    if(maxR>length(cnames)){
      extraDif  = maxR-length(cnames)
      maxL = maxL-extraDif
      maxR = length(cnames)
      }
    coldx = maxL:maxR
    
    tempMat = mat[rowdx,coldx]
    for(i in 1:dim(tempMat)[2]){
      q = c(q,as.numeric(tempMat[,i]))
      n = c(n,rep(as.numeric(colnames(tempMat)[i]),dim(tempMat)[1]))
      alpha2 = c(alpha2,as.numeric(rownames(tempMat)))
    }
 
    # work on n2 = m
    temp3 = which(as.numeric(cnames)<=n2)
    temp4 = which(as.numeric(cnames)>n2)
    coldx =  c(temp3[length(temp3)], temp4[1])
    # adjust index if below or above minimum or maximum preset column sample size (ie. <10 or >300)
    if((length(coldx)==1) & (coldx[1]==1)) coldx = c(1,2)
    if((length(coldx)==2) & (is.na(coldx[2]))) coldx = c((length(cnames)-1),length(cnames))

    # adjust index outward to collect 16 data points 
    # if on a boundary adjust to other side
    maxL = coldx[1]-idx
    maxR = coldx[2]+idx
    if(maxL<=0){
      if(maxL==0){
        maxR=maxR+1
      }else{
        maxR = maxR+(abs(maxL)+1)
      }
      maxL = 1           
    }
    if(maxR>length(cnames)){
      extraDif  = maxR-length(cnames)
      maxL = maxL-extraDif
      maxR = length(cnames)
      }
    coldx = maxL:maxR
    
    tempMat = mat[rowdx,coldx]
    for(i in 1:dim(tempMat)[2]){
      q = c(q,as.numeric(tempMat[,i]))
      n = c(n,rep(as.numeric(colnames(tempMat)[i]),dim(tempMat)[1]))
      alpha2 = c(alpha2,as.numeric(rownames(tempMat)))
    }
  # ends if not n=m  
  }

  nsq = n^2
 
  # fit model 
  eqn = summary(lm(q~n+nsq+alpha2))
  coeff = as.numeric(eqn$coefficients[,1])
  std.e = eqn$sigma
  vl = coeff[1] + coeff[2]*n1 + coeff[3]*(n1*n2) + coeff[4]*alpha + std.e

  return(c(vl, std.e^2))

  
} # end getCutoffTable







                           
helperGetTstat <- function(numberOfgroups, sample.size, num.mc=1000, delta=0.05){ 
  listvls <- c()
  for(i in 1:num.mc){
    testList<-c()
    #cat(i, "\n")
    for(j in 1:numberOfgroups){
      testList[[length(testList)+1]] <- rnorm(sample.size[j])
    }
    listvls <- c(listvls,computeValue(testList, delta))
   }
   return(listvls)
}






computeQhat <- function(vec.teststats, alpha, nsims, part1){

   # Part 2 B :  equation   
  Jmin = floor(max(c(1,((1-alpha)*nsims-sqrt(nsims)*log(nsims)))))
  Jmax = ceiling(min(c(nsims,((1-alpha)*nsims+sqrt(nsims)*log(nsims)))))
  Jint = as.numeric(Jmin:Jmax)
  num = exp((-nsims/(2*alpha*(1-alpha)))*((1-alpha - Jint/nsims)^2))*(sqrt(part1[2]/(2*pi))*(exp(-((vec.teststats[Jint-1]-part1[1])^2)/(2*part1[2]))-exp(-((vec.teststats[Jint]-part1[1])^2)/(2*part1[2])))+part1[1]*(pnorm(vec.teststats[Jint], part1[1], sqrt(part1[2]))-pnorm(vec.teststats[Jint-1], part1[1], sqrt(part1[2]))))
  den = exp((-nsims/(2*alpha*(1-alpha)))*((1-alpha - Jint/nsims)^2))*((pnorm(vec.teststats[Jint], part1[1], sqrt(part1[2]))-pnorm(vec.teststats[Jint-1], part1[1], sqrt(part1[2]))))
  value.cal = sum(num)/sum(den)

  return(value.cal)
  
}



varianceEstimator <- function(numberOfgroups, sample.size, nsims, qhat, o.nsims){

  # Part 3 A: get another 200 MC at n yielding 200 teststats
  part3 = helperGetTstat(numberOfgroups, sample.size, o.nsims)
  
  den2 = ((density(part3, from=qhat, to=qhat)$y[1])^2)*nsims
  # Part 3 B: get variance estimator
  F.fun =  ecdf(part3)
  num2 = F.fun(qhat)*(1-F.fun(qhat))
  var.est = num2/den2

  returnValue = list(part3, var.est)
  
  return(returnValue)
    
}

  

#
# Set for 3 groups
#
getCutoffTableThree <- function(sample.size, alpha){
  
  #load("dataTable/threeGroup.RData")
  #data(threeGroup)
  mat <- threeMat
  cnames <- colnames(mat)
  rnames <- rownames(mat)

  n1 = sample.size[1]
  n2 = sample.size[2]
  n3 = sample.size[3]
  
  NMTequal = FALSE
  # determin if using same n and m and t index
  if((n1 == n2) & (n2 == n3)){
    NMTequal = TRUE
  }else{
    # may not be equal but may use same index
    temp1 = which(as.numeric(cnames)<=n1)
    temp2 = which(as.numeric(cnames)>n1)
    coldx =  c(temp1[length(temp1)], temp2[1])

    temp3 = which(as.numeric(cnames)<=n2)
    temp4 = which(as.numeric(cnames)>n2)
    coldx2 =  c(temp3[length(temp3)], temp4[1])

    temp5 = which(as.numeric(cnames)<=n3)
    temp6 = which(as.numeric(cnames)>n3)
    coldx3 =  c(temp5[length(temp5)], temp6[1])
      
    comp = (coldx==coldx2)
    if(length(comp[comp==FALSE])==0){
      comp2 = (coldx2==coldx3)
      if(length(comp2[comp2==FALSE])==0){
        NMTequal = TRUE
      }
    }    
  }


  q=n=m=t=alpha2=c()

  # if n and m and t are equal 
  if(NMTequal){
    
    temp1 = which(as.numeric(rnames)<=alpha)
    temp2 = which(as.numeric(rnames)>alpha)
    rowdx =  c(temp1[length(temp1)], temp2[1])
    # adjust index if below or above minimum or maximum preset alpha values (ie. <.001 or >.999)
    if((length(rowdx)==1) & (rowdx[1]==1)) rowdx = c(1,2)
    if((length(rowdx)==2) & (is.na(rowdx[2]))) rowdx = c((length(rnames)-1),length(rnames))

    temp3 = which(as.numeric(cnames)<=n1)
    temp4 = which(as.numeric(cnames)>n1)
    coldx =  c(temp3[length(temp3)], temp4[1])
    # adjust index if below or above minimum or maximum preset column sample size (ie. <10 or >300)
    if((length(coldx)==1) & (coldx[1]==1)) coldx = c(1,2)
    if((length(coldx)==2) & (is.na(coldx[2]))) coldx = c((length(cnames)-1),length(cnames))

    # adjust index outward to collect 16 data points 
    # if on a boundary adjust to other side
    idx = 1
    maxD = rowdx[1]-idx
    maxU = rowdx[2]+idx
    if(maxD<=0){
      if(maxD==0){
        maxU=maxU+1
      }else{
        maxU = maxU+(abs(maxD)+1)
      }
      maxD = 1           
    }
    if(maxU>length(rnames)){
      extraDif  = maxU-length(rnames)
      maxD = maxD-extraDif
      maxU = length(rnames)
    }
    rowdx = maxD:maxU
    maxL = coldx[1]-idx
    maxR = coldx[2]+idx
    if(maxL<=0){
      if(maxL==0){
        maxR=maxR+1
      }else{
        maxR = maxR+(abs(maxL)+1)
      }
      maxL = 1           
    }
    if(maxR>length(cnames)){
      extraDif  = maxR-length(cnames)
      maxL = maxL-extraDif
      maxR = length(cnames)
      }
    coldx = maxL:maxR

    tempMat = mat[rowdx,coldx]
    for(i in 1:dim(tempMat)[2]){
      q = c(q,as.numeric(tempMat[,i]))
      n = c(n,rep(as.numeric(colnames(tempMat)[i]),dim(tempMat)[1]))
      alpha2 = c(alpha2,as.numeric(rownames(tempMat)))
    }
 
  #
  # ends if n=m
  # now if n, m or t differ
  #
  }else{

    # find rows = alpha
    temp1 = which(as.numeric(rnames)<=alpha)
    temp2 = which(as.numeric(rnames)>alpha)
    rowdx =  c(temp1[length(temp1)], temp2[1])
    # adjust index if below or above minimum or maximum preset alpha values (ie. <.001 or >.999)
    if((length(rowdx)==1) & (rowdx[1]==1)) rowdx = c(1,2)
    if((length(rowdx)==2) & (is.na(rowdx[2]))) rowdx = c((length(rnames)-1),length(rnames))
    # adjust index outward  
    # if on a boundary adjust to other side
    idx = 1
    maxD = rowdx[1]-idx
    maxU = rowdx[2]+idx
    if(maxD<=0){
      if(maxD==0){
        maxU=maxU+1
      }else{
        maxU = maxU+(abs(maxD)+1)
      }
      maxD = 1           
    }
    if(maxU>length(rnames)){
      extraDif  = maxU-length(rnames)
      maxD = maxD-extraDif
      maxU = length(rnames)
    }
    rowdx = maxD:maxU
  
    # work on n1 = n
    temp3 = which(as.numeric(cnames)<=n1)
    temp4 = which(as.numeric(cnames)>n1)
    coldx =  c(temp3[length(temp3)], temp4[1])
    # adjust index if below or above minimum or maximum preset column sample size (ie. <10 or >300)
    if((length(coldx)==1) & (coldx[1]==1)) coldx = c(1,2)
    if((length(coldx)==2) & (is.na(coldx[2]))) coldx = c((length(cnames)-1),length(cnames))

    # adjust index outward to collect 16 data points 
    # if on a boundary adjust to other side
    maxL = coldx[1]-idx
    maxR = coldx[2]+idx
    if(maxL<=0){
      if(maxL==0){
        maxR=maxR+1
      }else{
        maxR = maxR+(abs(maxL)+1)
      }
      maxL = 1           
    }
    if(maxR>length(cnames)){
      extraDif  = maxR-length(cnames)
      maxL = maxL-extraDif
      maxR = length(cnames)
      }
    coldx = maxL:maxR
    
    tempMat = mat[rowdx,coldx]
    for(i in 1:dim(tempMat)[2]){
      q = c(q,as.numeric(tempMat[,i]))
      n = c(n,rep(as.numeric(colnames(tempMat)[i]),dim(tempMat)[1]))
      alpha2 = c(alpha2,as.numeric(rownames(tempMat)))
    }
   
 
    
    # work on n2 = m
    temp3 = which(as.numeric(cnames)<=n2)
    temp4 = which(as.numeric(cnames)>n2)
    coldx =  c(temp3[length(temp3)], temp4[1])
    # adjust index if below or above minimum or maximum preset column sample size (ie. <10 or >300)
    if((length(coldx)==1) & (coldx[1]==1)) coldx = c(1,2)
    if((length(coldx)==2) & (is.na(coldx[2]))) coldx = c((length(cnames)-1),length(cnames))

    # adjust index outward to collect 16 data points 
    # if on a boundary adjust to other side
    maxL = coldx[1]-idx
    maxR = coldx[2]+idx
    if(maxL<=0){
      if(maxL==0){
        maxR=maxR+1
      }else{
        maxR = maxR+(abs(maxL)+1)
      }
      maxL = 1           
    }
    if(maxR>length(cnames)){
      extraDif  = maxR-length(cnames)
      maxL = maxL-extraDif
      maxR = length(cnames)
      }
    coldx = maxL:maxR
    
    tempMat = mat[rowdx,coldx]
    for(i in 1:dim(tempMat)[2]){
      q = c(q,as.numeric(tempMat[,i]))
      n = c(n,rep(as.numeric(colnames(tempMat)[i]),dim(tempMat)[1]))
      alpha2 = c(alpha2,as.numeric(rownames(tempMat)))
    }

    # work on n3 = t
    temp3 = which(as.numeric(cnames)<=n3)
    temp4 = which(as.numeric(cnames)>n3)
    coldx =  c(temp3[length(temp3)], temp4[1])
    # adjust index if below or above minimum or maximum preset column sample size (ie. <10 or >300)
    if((length(coldx)==1) & (coldx[1]==1)) coldx = c(1,2)
    if((length(coldx)==2) & (is.na(coldx[2]))) coldx = c((length(cnames)-1),length(cnames))

    # adjust index outward to collect 16 data points 
    # if on a boundary adjust to other side
    maxL = coldx[1]-idx
    maxR = coldx[2]+idx
    if(maxL<=0){
      if(maxL==0){
        maxR=maxR+1
      }else{
        maxR = maxR+(abs(maxL)+1)
      }
      maxL = 1           
    }
    if(maxR>length(cnames)){
      extraDif  = maxR-length(cnames)
      maxL = maxL-extraDif
      maxR = length(cnames)
      }
    coldx = maxL:maxR
    
    tempMat = mat[rowdx,coldx]
    for(i in 1:dim(tempMat)[2]){
      q = c(q,as.numeric(tempMat[,i]))
      n = c(n,rep(as.numeric(colnames(tempMat)[i]),dim(tempMat)[1]))
      alpha2 = c(alpha2,as.numeric(rownames(tempMat)))
    }
    
  # ends if not n=m  
  }
  
  ncub = n^3
 
  # fit model 
  eqn = summary(lm(q~n+ncub+alpha2))
  coeff = as.numeric(eqn$coefficients[,1])
  std.e = eqn$sigma
  vl = coeff[1] + coeff[2]*n1 + coeff[3]*(n1*n2*n3) + coeff[4]*alpha + std.e

  return(c(vl, std.e^2))

} # end getCutoffTable


