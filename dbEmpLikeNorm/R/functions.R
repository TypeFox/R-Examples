
#
# x is a list
# delta must be between 0 and 1 
#

dbELnorm <- function(x, delta=0.05, num.mc=1000, pvl.Table=TRUE, vrb=TRUE){

  if(class(x) == "list"){
  
    n<-length(x)
    if(n==1){
      if(vrb) cat("\n Only 1 group defined. Utilizing package dbEmpLikeGOF \n") 
      data <- x[[1]]
      dbEmpLikeGOF(x=data, testcall="normal", delta=delta, num.mc=num.mc, pvl.Table=pvl.Table, vrb-vrb)
    }else{
      
      #
      # get test statistic
      #
      if(vrb) cat("\n...Working on teststat \n")
      if(delta>=1 && delta<=0){
        cat("\n Delta value must be between 0 and 1. \n Changing to default delta=0.05 \n")
        delta<-0.05
      }
      teststat <- computeValue(Lst=x, delta=delta)
      
      #
      # get or estimate p value
      #
      if(pvl.Table){
        if(delta != .05) cat("Table calculated on a delta = 0.05 \n    To use delta=",delta," please change pvl.Table to FALSE \n")
        pvl <- computePval(teststat, x, vrb)          
      }else{        

        if(vrb) cat("...Working on p-value \n    This may take a few minutes \n")
        mc.est<-c()
        for(j in 1:num.mc){
          testList <- list()
          for(i in 1:n){
            testList[[length(testList)+1]] <- rnorm(length(x[[i]]))
          }
          mc.est[j] <- computeValue(testList, delta)
        }
        pvl <- sum((mc.est>teststat))/num.mc
       
      
      }# end  calculate pval


    
    }# end if list is greater than 2 groups 
  
  }else{ # end if not a list
    if(vrb) cat("\n argument x  must be a list. \n")
    teststat <- NA
    pvl <- NA
  }
  
  output <- list(teststat=teststat, pvalue=pvl)
  return(output)
 
}



returnCutoffValue <- function(numberOfgroups, sample.size, targetalpha=0.05, MC.Method=TRUE, Table.Method=FALSE, Bayes.Method=FALSE, num.mc=1000, delta=0.05,  nsims=200, v.threshold=NA){

  if(!MC.Method & !Table.Method & !Bayes.Method){
    cat("No Method Selected. Utilizing Monte Carlo Techniques \n")
    MC.Method = TRUE
  }
  
  cut.out1 = NA
  cut.out2 = NA
  cut.out3 = NA
  
  if(MC.Method){ 

    listvls <- c()
    for(i in 1:num.mc){
      testList <- list()
      if(length(sample.size)<numberOfgroups) {
        dif <- numberOfgroups-length(sample.size)
        sample.size <- c(sample.size,rep(sample.size[length(sample.size)], dif))
      }
      for(j in 1:numberOfgroups){
        testList[[length(testList)+1]] <- rnorm(sample.size[j])
      }
      listvls <- c(listvls,computeValue(testList, delta))
      
    }
    cut.out1 <- quantile(listvls, (1-targetalpha), na.rm=TRUE)
  }
    
  if(Table.Method){
    
    cat("estimating cutoff based on table \n")
    if(numberOfgroups == 1){
      cat("only 1 group defined. utilizing dbEmpLikeGOF getCutoff  \n")
      cut.out2 <- getCutoff(testcall="normal", smp.size=sample.size, targetalpha=targetalpha)
    }else{
      if(numberOfgroups > 3){
        cat("Tables are based on 2 or 3 groups. To get a cutoff for greater than 3 groups use [MC.Method=TRUE]  \n")
        numberOfgroups = 3
      }
      cut.out2 <- getCutoffValue(numberOfgroups, sample.size, targetalpha)
    }
  }

  if(Bayes.Method){

    cat("estimating cutoff based on bayesian method \n")
    if(numberOfgroups == 1){
      if(Table.Method){
        cut.out3 = cut.out2
      }else{
        cat("only 1 group defined. utilizing dbEmpLikeGOF getCutoff [non bayesian]  \n")
        cut.out3 <- getCutoff(testcall="normal", smp.size=sample.size, targetalpha=targetalpha)
      }
    }else{
      if(numberOfgroups > 3){
        cat("Tables are based on 2 or 3 groups. To get a cutoff for greater than 3 groups use [MC.Method=TRUE]  \n")
        numberOfgroups = 3
      }
      #  CHECK IF ROW AND COLUMN MATCH ... THEN JUST RETURN VALUE IN TABLE. NO BAYES
      if(length(sample.size)<numberOfgroups) {
        dif <- numberOfgroups-length(sample.size)
        sample.size <- c(sample.size,rep(sample.size[length(sample.size)], dif))
      }
      cnames = c(10, seq(25,300,25))
      rnames = seq(.001, .999, .001)
      if((length(rle(sample.size)$length) == numberOfgroups) & (length(which(cnames == sample.size[1]))>0) & (length(which(rnames == targetalpha))>0)){
        if(numberOfgroups == 2){
          #data(twoGroup)
          cut.out3 <- twoMat[which(rnames == targetalpha),which(cnames == sample.size[1])]
        }
        if(numberOfgroups == 3){
          #data(threeGroup)
          cut.out3 <- threeMat[which(rnames == targetalpha),which(cnames == sample.size[1])]
        }        
      }else{
        cut.out3 <- BayesWays(numberOfgroups, sample.size, targetalpha, nsims, v.threshold)[1]
      }
    }   
  }

  cut.out = list(MC.Method.rslt = as.numeric(cut.out1), Table.Method.rslt = as.numeric(cut.out2), Bayes.Method.rslt = as.numeric(cut.out3))
  if(length(which(is.na(cut.out)))>0){
    cut.rslt = cut.out[-which(is.na(cut.out))]
  }else{
    cut.rslt = cut.out
  }
  if(length(cut.rslt)==1) cut.rslt = unlist(cut.rslt)
  
  return(cut.rslt)
}
