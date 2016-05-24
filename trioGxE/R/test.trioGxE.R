test.trioGxE <- function(object,data=NULL,nreps,level=0.05,early.stop=FALSE,fix.sp=FALSE,
                         output=NULL,return.data=FALSE,return.object=FALSE,...){
  res=list()
  if(fix.sp)
    sp<-object$sp
  else
    sp<-NULL
  ## is object NULL?
  if(is.null(object)){
    if(!is.null(data)){
      # fit 'data' with trioGxE() ## other arguments must be passed through "...".
      object <- trioGxE(data=data,...)
    }
    else
      stop("Either \'trioGxE\' object or data must be provided.")
  }
  
  # 'object' is passed by the user
  else {
    ## check if the original data is contained the object
    data <- object$data
    if(is.null(data)) {
      stop("The original data set is necessary for testing GxE: 
           Please apply 'trioGxE()' to the data again 
           by setting 'return.data=TRUE'.")
    }
  }

  ## extract necessary information from the fitting object
  control = object$control
  pgenos = object$terms$pgenos ## is necessary?? ## for obtaining the maing type
  cgeno = object$terms$cgeno ## is necessary??
  cenv = object$terms$cenv ## need this to choose the column to permute
  penmod = object$penmod ## need this
  tem.k = object$smooth$bs.dim ## need this so that we fit the generated data using the same basis dimension
  coef = object$coef
  Vp = object$Vp
  Gp = object$Gp
  
  if(penmod == "codominant"){
    k = c(tem.k[1],tem.k[2])
  }
  
  else{
    k = tem.k[tem.k!=0][1]
  }

  GxE.stat <- test.stat(coef=coef, var.cov=Vp, k=k)
  if(!is.null(output)){
    write.table(GxE.stat,file=output,quote=FALSE,
                row.names=FALSE,col.names=FALSE,append=TRUE)
  }
  
  ## Permuting data begins here
  perm.dat = data
  ## is early-termination deployed? 
  if(early.stop){
    extreme.stat.num = ceiling(level*(nreps+1))
        
    cat("Early termination is being deployed:\n Sampling will be continued until ", 
        extreme.stat.num, " replicates (including itself) with test statistics greater", 
        "than or equal to the observed value are obtained. \n",sep="")
    
    extreme.stat.GxE=1; #including the observed test statistic itself    
    i=2 ## beginning for the second row
    while( ( ( extreme.stat.GxE < extreme.stat.num )  & ( ( i-1 ) <= nreps ) ) ) {
      ## permute data, calculate test statistic, write output
      perm.stat <- perm.test.stat(data=data,Gp=Gp,cenv=cenv,pgenos=pgenos,
                                  cgeno=cgeno,penmod=penmod,k=k,sp=sp,control=control)
      if(!is.null(output)){
        write.table(perm.stat,file=output,quote=FALSE,
                    row.names=FALSE,col.names=FALSE,append=TRUE)
      }
      GxE.stat = rbind(GxE.stat,perm.stat)

      if(GxE.stat[i,1]>=GxE.stat[1,1]) #contribute to the tail probability P(T>=Tobs)
        extreme.stat.GxE=extreme.stat.GxE+1
      else
        extreme.stat.GxE=extreme.stat.GxE
      i = i+1 ##
    }#while ends
  }#if(early.stop) ends
  
  else{ #not early stop
    for(i in 1:nreps){
      perm.stat <- perm.test.stat(data=data,Gp=Gp,cenv=cenv,pgenos=pgenos,
                                  cgeno=cgeno,penmod=penmod,k=k,sp=sp,control=control)
      if(!is.null(output)){
        write.table(perm.stat,file=output,quote=FALSE,
                    row.names=FALSE,col.names=FALSE,append=TRUE)
      }
      GxE.stat = rbind(GxE.stat,perm.stat)
    }
    
  }# not early stop
  res$GxE.stat = GxE.stat
  if(penmod=="codominant"){
  res$p.value = (1/nrow(GxE.stat))*c(sum(GxE.stat[,1]>=GxE.stat[1,1]),
                                      sum(GxE.stat[,2]>=GxE.stat[1,2]),
                                      sum(GxE.stat[,3]>=GxE.stat[1,3]))
  }
  else res$p.value = (1/nrow(GxE.stat))*sum(GxE.stat[,1]>=GxE.stat[1,1])
  res  
}

test.stat <- function(coef,var.cov,k){
  if(length(k)==2){
    k1 = k[1]
    k2 = k[2]
    sm.coef = coef[-c(1,(k1+1))] ## removing the intercept parameters
    sm1.coef = sm.coef[1:(k1-1)]
    sm2.coef = sm.coef[k1:(k1+k2-2)]
    
    V.sm.coef = var.cov[-c(1,(k1+1)),-c(1,(k1+1))] 
    V.sm1.coef=V.sm.coef[c(1:(k1-1)),c(1:(k1-1))]
    V.sm2.coef=V.sm.coef[c(k1:(k1+k2-2)),c(k1:(k1+k2-2))]
    
    GxE.stat=t(sm.coef) %*% solve(V.sm.coef) %*% (sm.coef)
    GxE1.stat=t(sm1.coef) %*% solve(V.sm1.coef) %*% (sm1.coef)
    GxE2.stat=t(sm2.coef) %*% solve(V.sm2.coef) %*% (sm2.coef)
    
    res = c(GxE.stat,GxE1.stat,GxE2.stat)
  }
  else {# if length(k)==1
    sm.coef = coef[-1] ## removing the intercept parameter
    ## Extracting the var-cov estimates obtained from
    ## fitting the original data
    V.sm.coef = var.cov[-1,-1]
    res=t(sm.coef) %*% solve(V.sm.coef) %*% (sm.coef)
  } 
  res = matrix(res, nrow=1)# row-matrix
  res
}

perm.test.stat <- function(data,Gp,cenv,pgenos,cgeno,penmod,k,sp,control){
  perm.dat <- data
  perm.obj <- NULL
  while(is.null(perm.obj)){
  if(any(Gp==1))
    perm.dat[Gp==1,cenv] = sample(data[Gp==1,cenv]) #mt1
  if(any(Gp==2))
    perm.dat[Gp==2,cenv] = sample(data[Gp==2,cenv]) #mt2
  if(any(Gp==3))
    perm.dat[Gp==3,cenv] = sample(data[Gp==3,cenv]) #mt3
  
  perm.obj <- trioGxE(data=perm.dat,pgenos=pgenos,cgeno=cgeno,cenv=cenv,
                      penmod=penmod,k=k,sp=sp,control=control,testGxE=TRUE)
  }## run until get the dataset with converged results
  
  perm.stat <- test.stat(coef=perm.obj$coef,var.cov=perm.obj$Vp,k=k)
  perm.stat
}
