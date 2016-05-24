############################
# permutation Test of symmetry (i.e. one sample test) 
# Y is the Nxp matrix of respoonces
# perms: number of permutations
# tail : vector of tails 1, -1 or 0
# permP.return, permT.return, permSpace.return : logical: shoul space of p-values, of statistic and of signs be returned?
############################
.symmetry.nptest <- function(data, perms=5000, statTest="t",  tail = NULL, testType="permutation",...){
	if(is.function(statTest)) {
		test<-statTest
	} else if(statTest%in%c("t","sum")){
		if (testType=="rotation") {
      test <- .t.rotation.nptest.1sample
      } else ## permutation test
        test <- .t.symmetry.nptest
	} else if(statTest%in%c("Wilcoxon","rank","Sign")){
		if (testType=="rotation") warning("Rotations are not allowed for Wilcoxon (i.e. ranks) test, permutations will be used instead.")
		 ## permutation test
		test <- .rank.symmetry.nptest
	} else if(statTest%in%c("McNemar")){
		if (testType=="rotation") warning("Rotations are not allowed for Mc Nemar test, permutations will be used instead.")
		 ## permutation test
		test <- .mcnemar.symmetry.nptest
	} else  {stop("This test statistic is not valid, nothing done."); return()}
  environment(test) <- sys.frame(sys.nframe())	
  out <- sys.frame(sys.nframe())
}

#####################################

.t.symmetry.nptest <- function(){
  
  if(is.null(data$W)){
    Ns=apply(!is.na(data$Y),2,sum)
    data$Y[is.na(data$Y)]=0
  } else { ### W can be a vector of length nrow(Y) or a matrix of same dim of Y
    Ns=apply((!is.na(data$Y) & !is.na(data$W)),2,sum)
    data$Y[(is.na(data$Y) | is.na(data$W))]=0
    data$Y=data$W * data$Y
#     perms <- make.signSpace(nrow(data$Y),perms)
#     permT <- rbind(rep(1,perms$n),perms$permID) %*% (data$W * data$Y)
  }
  perms <- make.signSpace(nrow(data$Y),perms)
  permT <- rbind(rep(1,perms$n),perms$permID) %*% data$Y
  colnames(permT) = .getTNames(data$Y)
  permT = rbind(permT,-permT[nrow(permT):1,,drop=FALSE])
  rownames(permT)=.getTRowNames(permT)
  
  if(statTest=="t") {
      M2s=apply(data$Y^2,2,sum)
    permT=permT/t(sqrt((M2s-t((permT)^2)/Ns)*((Ns)/(Ns-1))))
  }
  colnames(permT) = .getTNames(data$Y,permT=permT)
  return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test=statTest)))
}

######signed ranks test
.rank.symmetry.nptest <- function(){
  if(statTest=="Sign"){
    data$Y=sign(data$Y)  
    Test="Sign"
  } else {
    data$Y=apply(data$Y,2,function(y) 
    {values= !(is.na(y) | (y==0))
     values[values]=rank(abs(y[values]))*sign(y[values])
     values
    })
    Test="Wilcoxon Sign"
  }
	
	perms <- make.signSpace(nrow(data$Y),perms)
	
	if(is.null(data$W)){
		permT <- rbind(rep(1,perms$n),perms$permID) %*% data$Y
	} else { ### W can be a vector of length nrow(Y) or a matrix of same dim of Y
		permT <- rbind(rep(1,perms$n),perms$permID) %*% (data$W * data$Y)
	}
	
	permT = rbind(permT,-permT[nrow(permT):1,,drop=FALSE])
	rownames(permT)=.getTRowNames(permT)
	colnames(permT) = .getTNames(data$Y,permT=permT)
	return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test=Test)))
}

######McNemar test
.mcnemar.symmetry.nptest <- function(){
	data$Y=sign(data$Y)
	data$Y[is.na(data$Y)]=0

	perms <- make.signSpace(nrow(data$Y),perms)
	
	if(is.null(data$W)){
		permT <- rbind(rep(1,perms$n),perms$permID) %*% data$Y
	} else { ### W can be a vector of length nrow(Y) or a matrix of same dim of Y
		permT <- rbind(rep(1,perms$n),perms$permID) %*% (data$W * data$Y)
	}
	
	permT=scale(permT,scale=sqrt(colSums(abs(data$Y))))
	permT = rbind(permT,-permT[nrow(permT):1,,drop=FALSE])
	rownames(permT)=.getTRowNames(permT)
	colnames(permT) = .getTNames(data$Y,permT=permT)
	return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="McNemar")))
}

#############ROTATION

.t.rotation.nptest.1sample <- function(){
  data <- .orthoZ(data)
  
  naRows=apply(is.na(data$Y),1,sum)
  if(any(naRows>0)) {
    warning("Some NA on Y. In rotationTest observations are excluded row-wise")
    Y=Y[naRows==0,,drop=FALSE]
  }
  Ns=nrow(data$Y)
  M2s=apply(data$Y^2,2,sum)
  
  perms <- make.permSpace(1:Ns,perms,testType=testType)
  digitsK=trunc(log10(perms$B))+1
  
  permT=matrix(,perms$B, ncol(data$Y))
  permT[1,]=colSums(data$Y)
  for(i in 2:(perms$B)) { 
      if (i%%1000==0) {
        cat(rep("\b", 2*digitsK+3), i, " / ", perms$B, sep="")
        flush.console()
      }
      # R is random matrix of independent standard-normal entries 
      # Z shall be a random matrix with the same mean and covariance structure as Y 
      permT[i,] = colSums(perms$rotFunct())
  }
  cat(rep("\b", 2*digitsK+3));  flush.console()
  
  colnames(permT) = .getTNames(data$Y)
  rownames(permT)=.getTRowNames(permT)
  permT=permT/t(sqrt((M2s-t((permT)^2)/Ns)*((Ns)/(Ns-1))))
  return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="t")))
}