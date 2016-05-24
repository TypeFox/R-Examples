qvSUDT <-
function(teststats,alternative="U",df=Inf,corr=0.5,corr.matrix=NA,mcs=1e+05){
  k <- length(teststats)
  names(teststats) <- paste("H",1:k,sep="")
  if(k > 16){
  	stop("The funtion is not applicable to situations where the number of tests exceeds 16 \n")
  }
  if(k >= 6 & k<=10){
    cat("The calculation may take up to several minutes, please be patient \n")
  }
  if(k > 10){
  	cat("The calculation will be time consuming \n")
  }
  
  #order the test statistcis
  if(alternative=="U"){od.index <- order(teststats,decreasing=FALSE)}
  if(alternative=="B"){od.index <- order(abs(teststats),decreasing=FALSE)}
  t.ordered <- teststats[od.index]
  if(is.matrix(corr.matrix)){corr.matrix <- corr.matrix[od.index,od.index]}
  
  pvSet <- NULL #to store the intermediate P-values
  qvSet <- NULL #to store the adjusted P-values (Q-values)
  z0 <- rnorm(mcs)
  if(df==Inf){u=1} else {u=sqrt(rchisq(mcs,df)/df)}
  
  #for different tail test, the function for calculate Psi(di)
  if(alternative=="U"){
  	psi.fun <- function(cj,corr){pnorm((cj*u+z0*sqrt(corr))/sqrt(1-corr))}
    }
  if(alternative=="B"){
    psi.fun <- function(cj,corr){pnorm((cj*u+z0*sqrt(corr))/sqrt(1-corr))-pnorm((-cj*u+z0*sqrt(corr))/sqrt(1-corr))}
    }

  #calculate the Q1
  if(alternative=="U"){
    P1 <- 1-pt(t.ordered[1],df)
  }
  if(alternative=="B"){
    P1 <- 2*(1-pt(abs(t.ordered[1]),df))
  }
  pvSet <- c(pvSet,P1)
  qvSet <- c(qvSet,P1)
  
list.J.Fun <-list(J2.fun,J3.fun,J4.fun,J5.fun,J6.fun,J7.fun,J8.fun,J9.fun,J10.fun,J11.fun,J12.fun,J13.fun,J14.fun,J15.fun,J16.fun)
#calculate the Q2-qvSet relying on the function Pj.fun given below
#Given a value of Pj and cj = t[j], evaluate c2-c(j-1), j=1,...,k, by numerical caculations
  Pj.fun <- function(Pj,j){
    #given Pj (act as alpha), the c1 is always evaluated by
    if(alternative=="U"){c1 <- qt(1-Pj,df)}
    if(alternative=="B"){c1 <- qt(1-Pj/2,df)}
    # find the corresponding averaged correlations between t1,...,tj 
    if(is.matrix(corr.matrix)){
    	subset <- 1:j
    	corr.subset <- corr.matrix[subset,subset]
    	corr <- sum(corr.subset[row(corr.subset)<col(corr.subset)])/(j*(j-1)/2)
    	}
    #function to solve for cj's given Pj
    cj.fun <- function(cj,J.fun,list.J,list.Psi){
    	Psi.dj <- psi.fun(cj,corr)
  	 	mean(J.fun(Psi.dj,list.J,list.Psi))-(1-Pj)
  	 	}    
  	# initialize 	
    J0 <- 1
    Psi.d1 <- psi.fun(c1,corr)
    J1 <-  Psi.d1
    list.J <- c(list(J0),list(J1))
    list.Psi <- list(Psi.d1)
    if(j==2){
        Psi.d2 <- psi.fun(t.ordered[2],corr)
        Error <- mean(J2.fun(Psi.d2,list.J,list.Psi))-(1-Pj)
      }
    if(j >= 3){
    	for(t in 2:(j-1)){
    		ct <- uniroot(cj.fun,J.fun=list.J.Fun[[t-1]],list.J=list.J,list.Psi=list.Psi,interval=c(-5,5),tol=1e-04)$root
    		Psi.dt <- psi.fun(ct,corr)
    		Jt <- list.J.Fun[[t-1]](Psi.dt,list.J,list.Psi)
    		list.J <- c(list.J,list(Jt))
    		list.Psi <- c(list.Psi,list(Psi.dt))
    		}
    	Psi.dj <- psi.fun(t.ordered[j],corr)
    	Error <- mean(list.J.Fun[[j-1]](Psi.dj,list.J,list.Psi))-(1-Pj)
    	}
    	 Error
   }
  
  #now solve for P2-pvSet by iterations, respectively
  for(j in 2:k){
    if(j>=6){
      cat(paste("In calculating Q",j,"\n",sep=""))}
    low.try<-1e-04
    up.try<-0.5
    low.err <- Pj.fun(low.try,j)
    up.err <- Pj.fun(up.try,j)
    err <- min(abs(c(low.err,up.err)))
    while(abs(err) > 1e-04){
      try.new <- (low.try+up.try)/2
      err <- Pj.fun(try.new,j)
      if(err > 1e-04){
        up.try <- try.new
      }
      if(err < -1e-04){
        low.try <- try.new
      }
    }
    pvSet <- c(pvSet,try.new)
    qvSet <- c(qvSet,min(pvSet))
  }
  Results <- list(t.ordered,as.vector(round(qvSet,3)))
  names(Results) <- c("ordered test statistics","Adjusted P-values of ordered test statistics")
  return(Results)
}
