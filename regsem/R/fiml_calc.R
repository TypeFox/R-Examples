
fiml_calc = function(ImpCov,data,Areg,lambda,alpha,type,pen_vec,nvar){

# look at lav_objective.R from lavaan
  m = dim(ImpCov)[1]
  IntCol = which(colnames(Areg) == "1")
  IntCol2 = colnames(Areg[,1:nvar])
  fit = 0

  if(type=="none"){

    for(i in 1:nrow(data)){
      person1 = data[i,IntCol2]

      ind1 = which(is.na(person1)==T)
      misVar = colnames(data)[ind1]
      K = (nvar - sum(unique(ind1))) * log(2*pi)


      if(sum(is.na(person1)==T) >0){
        meanvec = Areg[(rownames(Areg)!="1"),IntCol]
        meanvec2 = meanvec[names(meanvec) != misVar]
        sub1 = c(person1[is.na(person1)==FALSE],0) - meanvec2
        indFit = K - log(det(ImpCov[-ind1,-ind1])) + t(sub1) %*% solve(ImpCov[-ind1,-ind1]) %*% sub1
      }else{
        meanvec = Areg[(rownames(Areg)!="1"),IntCol]
        sub1 = as.numeric(cbind(person1,0) - meanvec)
        indFit = K - log(det(ImpCov))  + t(sub1) %*% solve(ImpCov) %*% sub1
      }
      fit = fit + log(indFit)
    }


  }else{
    stop("Only type==none is currently supported")
  }


  fit

}
