"bounds" <-
function(t,t2=t,iuse=1,asf=NULL,alpha=0.05,phi=rep(1,length(alpha)),ztrun=rep(8,length(alpha))){
   if (length(t) != length(t2)){
      stop("Original and second time scales must be vectors of the same length.")
    }
   if ({min(t) <= 0}|{max(t) > 1}|{min(t2) <= 0}){
      stop("Analysis times must be in (0,1].  Second time scale values must be positive.")
    }
   t3 <- t2
   t2 <- t2/max(t2)
   if ({sum({t-c(0,t[-length(t)]) <= 0}) > 0}|{sum({t2-c(0,t2[-length(t)]) <= 0}) > 0}){
      stop("Analysis times must be ordered from smallest to largest.")
    }
   if ({sum(alpha <= 0) > 0}|{sum(alpha) > 1}){
      stop("Each component of alpha must be positive and their sum cannot exceed 1.")
    }
   if (length(iuse) != length(alpha)){
      stop("For two-sided bounds, the lengths of the iuse and alpha vectors must both be 2.")
    }
   if (length(asf)==2){
      if ({class(asf[[1]])!="function"}|{class(asf[[2]])!="function"}){
         stop("Alpha spending function must be of class 'function'.")
       }
    }
   if (length(asf)==1){
      if (class(asf)!="function"){
         stop("Alpha spending function must be of class 'function'.")
       }
    }
   if (sum(iuse==5)<length(asf)){
      stop("Can't specify 2 spending functions unless iuse=c(5,5).")
    }
   if ({sum(iuse==5)>0}&{length(asf)==0}){
      stop("If iuse=5, must specify spending function.")
    }
   if (sum({iuse==3}|{iuse==4}) > length(phi)){
      stop("Phi must be specified for each boundary that uses spending function 3 or 4.")
    }
   if (sum({iuse==3}&{phi <= 0}) > 0){
      stop("For power family (iuse=3), phi must be positive.")
    }
   if (sum({iuse==4}&{phi==0}) > 0){
      stop("For Hwang-Shih-DeCani family (iuse=4), phi can not be 0.")
    }
   if (length(phi)==1) phi <- rep(phi,2)
   if ({length(alpha)==1}|{{length(alpha)==2}&{alpha[1]==alpha[2]}&{iuse[1]==iuse[2]}&{length(asf)!=2}&{ztrun[1]==ztrun[2]}&{{length(phi)==1}|{phi[1]==phi[2]}}}){
      if (length(alpha)==1){
         type <- 1
         alph <- alpha
       }
      if (length(alpha)==2){
         type <- 2
         alph <- 2*alpha[1]
       }
      ld <- landem(t,t2,length(alpha),iuse[1],asf,alph,phi[1],ztrun[1])
      ubnd <- ld$upper.bounds
      lbnd <- ld$lower.bounds
      epr <- ld$exit.pr
      dpr <- ld$diff.pr
      spend <- ld$spend
    }
   else{
      type <- 3
      ld1 <- landem(t,t2,1,iuse[1],asf[[1]],alpha[1],phi[1],ztrun[1])
      ld2 <- landem(t,t2,1,iuse[2],asf[[2]],alpha[2],phi[2],ztrun[2])
      lbnd <- -ld1$upper.bounds
      ubnd <- ld2$upper.bounds
      epr <- ld1$exit.pr+ld2$exit.pr
      dpr <- ld1$diff.pr+ld2$diff.pr
      spend <- c(ld1$spend,ld2$spend)
    }
   ans <- list(bounds.type=type,spending.type=spend,time=t,time2=t3,alpha=alpha,overall.alpha=sum(alpha),lower.bounds=lbnd,upper.bounds=ubnd,exit.pr=epr,diff.pr=dpr)
   class(ans) <- "bounds"
   return(ans)
 }

