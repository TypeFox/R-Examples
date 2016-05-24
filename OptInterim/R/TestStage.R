"TestStage" <-
function(tan,tstage,x,num.arm,num.stage,
         Y1,T1,Y0=NULL,T0=NULL,p0=NULL,
         C1L=NULL,C1U=NULL,C2L=NULL,C2U=NULL,C3U=NULL,
         printTest=TRUE,
         cen1=rep(1,length(T1)), cen0=rep(1,length(T0)))
{
  if( !identical(length(Y1),length(T1))|
    (identical(num.arm,2) && !identical(length(Y0),length(T0))) )
    stop("study times and failure times should be of equal length")

  if(isTRUE(x>=tan))
    stop("the survival time of interest should not exceed the interim study time")

  if(identical(num.arm,1) && !length(p0))
    stop("The null rate must be specified in single arm trials")

  if(identical(num.arm,2) && !length(Y0))
    stop("A control group must be specified in two arm trials")

  if(!identical(printTest,FALSE) && !identical(printTest,TRUE))
    stop("printTest should be a logical value")

  if(!length(C3U) && identical(tstage,3))
    stop("Stage 3 boundary missing")
    

    zout <- tst(tan,x,num.arm,Y1,T1,p0,Y0,T0,cen1,cen0)  
    z<-zout["z"]

    if(printTest){
      if(identical(tstage,1)){  
          cat(paste("interim test statistic Z1 =",z,"C1L =",C1L,"C1U =",C1U),"\n")
          if(z<=C1L)
            cat("Z1 <= C1L, stop the study and accept the null hypothesis \n") 
             
          if(z>=C1U)  
            cat("Z1 >= C1U, stop the study and reject the null hypothesis \n") 
          
          if(z>C1L && z<C1U)  
            cat("Z1 > C1L and Z1 < C1U, continue to the second stage \n")
      }else if(identical(tstage,2)){  
        if(num.stage==2){
          cat(paste("final test statistic Z2 =",z,"C2U =",C2U),"\n")    
          if(z>=C2U)
            cat("since Z2 >= C2U, reject the null hypothesis \n")
          if(z<C2U)
            cat("since Z2 < C2U, cannot reject the null hypothesis \n")
        }else if(num.stage==3){
          cat(paste("interim test statistic Z2 =",z,"C2L =",C2L,"C2U =",C2U),"\n") 
          if(z<=C2L)
            cat("Z2 < C2L, stop the study and accept the null hypothesis \n")          
          if(z>=C2U)  
            cat("Z2 >= C2U, stop the study and reject the null hypothesis \n")       
          if(z>C2L && z<C2U)  
            cat("Z2 > C2L and Z2 < C2U, continue to the third stage \n")      
        }
      }else{  ## third stage  
        cat(paste("final test statistic Z3 =",z,"C3U =",C3U),"\n")    
        if(z>=C3U)
          cat("since Z3 >= C3U, reject the null hypothesis \n")
        if(z<C3U)
          cat("since Z3 < C3U, cannot reject the null hypothesis \n")
      }  
    }
    return(zout)  
}  

  # Nelson-Aalen estimate of the cumulative hazard function
  cum <- function(Y,T,t,x,cen)
  {
    X <- pmin(T,x,pmax(0,t-Y))
    delta <- as.numeric(T<=pmin(x,pmax(0,t-Y)) & cen)
    ### set to censored status if explicitly specified
    n <- length(X)
    R <- as.numeric(n)
    for(i in 1:n)
      R[i] <- sum(X>=X[i])
    cum<-sum(delta[X<=x]/R[X<=x])
    varcum<-sum(delta[X<=x]/(R[X<=x]^2))
    return(c(cum=cum,varcum=varcum))
  }

  # test statistic Z(x;t) of equation 3
  tst <- function(t,x,num.arm,Y1,T1,p0,Y0,T0,cen1,cen0)
  {
    cumout1<-cum(Y1,T1,t,x,cen1)
    if(num.arm==1){cumout0<- c(cum=-log(1-p0),varcum=0)
                   names(cumout0)<-c("cum","varcum")
    }else cumout0<-cum(Y0,T0,t,x,cen0)
    
    se<-  sqrt( cumout0["varcum"]/(max(cumout0["cum"],1e-4)^2)+cumout1["varcum"]/(max(cumout1["cum"],1e-4)^2) )
    z<-(log(max(cumout0["cum"],1e-4))-log(max(cumout1["cum"],1e-4)))/se

	names(z)<-"z"
    names(se)<-"se"
    cumL<-log(c(max(cumout1["cum"],1e-4), max(cumout0["cum"],1e-4)))
    names(cumL)<-c("Log(cumL1)","Log(cumL0)")
    return(c(z,se,cumL))
  }
