# Functions for calculating exact conditional independence tests for discrete data.
# Presently supported: deviance and wilcoxon


ciTest_ordinal <- function(x, set=NULL, statistic="dev", N=0, ...){

  statistic <- match.arg(statistic, c("deviance","wilcoxon","kruskal","jt"))
  if (inherits(x,"data.frame")){
    dataNames <- names(x)  
  } else {
    if (inherits(x,"table")){
      dataNames <- names(dimnames(x))
    } else {
      stop("'x' must be either a table or a dataframe")
    }
  }
  if (is.null(set)){
    set <- dataNames
    set.idx <- 1:length(set)
  } else {
    if (inherits(set,"numeric")){
      set.idx <- set
      set <- dataNames[set.idx]
    } else
    if (inherits(set,c("formula","character"))){
      set <- unlist(rhsFormula2list(set))
      set.idx <- match(set, dataNames)
    }
  }

  .CI.ordinal(set.idx, set, dataset=x, test=statistic, N=N)

}

.CI.ordinal <- function(set.idx, set, dataset, test="deviance", N=0) {  

  c1 <- set.idx[1]
  c2 <- set.idx[2]
  if (length(set.idx)>2)
    S <- set.idx[-(1:2)]
  else
    S <- NULL
  
#  cat(sprintf("CHK: c1: %s, c2: %s, S: %s test: %s\n", toString(c1), toString(c2), toString(S), test))
#  print(class(dataset))
#  str(dimnames(dataset))
  ## Calculates the deviance, degrees of freedom and asymptotic P given an ftable (m).   
  LRT <- function(m, d1, d2) {  
    oneslice <- function(t,d1,d2) {  
      dim(t) <- c(d1,d2)
      t1  <- addmargins(t)   
      cm  <- t1[d1+1,1:d2]
      rm  <- t1[1:d1,d2+1]
      N   <- t1[d1+1,d2+1]  
      df  <- (sum(cm>0)-1)*(sum(rm>0)-1)  
      dev <- 0  
      if (df>0) {
        fv <- (rm %o% cm)/N
        dev <- 2*sum(t*log(t/fv), na.rm=T)
      }  
      return(c(df=df, dev=dev))  
    }    
    ans <- apply(m, 1, oneslice, d1, d2)      
    obs.deviance <- sum(ans[2,])  
    df <- sum(ans[1,])  
    P <- 1 - pchisq(obs.deviance, df)  
    return(list(deviance=obs.deviance, df=df, P=P))  
  }  
  
  ## Calculates the wilcoxon test, its mean and asymptotic P given an ftable (m).            
  wilcoxon <- function(m, d1, d2) {  
    oneslice <- function(t,d1,d2) {  
      dim(t) <- c(d1,d2)
      t1 <- addmargins(t)   
      cm <- t1[d1+1,1:d2]
      rm <- t1[1:d1,d2+1]
      N  <- t1[d1+1,d2+1]  
      r  <- cumsum(c(0, cm[-d2]))+(1+cm)/2  
      W  <- sum(r*t[1,])  
      EW <- (rm[1]/N)*sum(r*cm)  
      VW <- (rm[1]*rm[2]/(N*(N-1)))*sum(((r-EW/rm[1])^2)*cm)  
      return(c(W, EW, VW))  
    }      
    ans <- apply(m, 1, oneslice, d1, d2)      
    W   <- sum(ans[1,])  
    EW  <- sum(ans[2,])  
    VW  <- sum(ans[3,])  
    P   <- 2*(1 - pnorm(abs(W-EW), sd=sqrt(VW)))  
    return(list(W=W, EW=EW, P=P))  
  }       
  
  ## Calculates the kruskal-wallis test, degrees of freedom and asymptotic P given an ftable (m).            
  kruskal <- function(m, d1, d2) {   
    oneslice <- function(t,d1,d2) {  
      dim(t) <- c(d1,d2)
      t1 <- addmargins(t)   
      cm <- t1[d1+1,1:d2]
      rm <- t1[1:d1,d2+1]
      N  <- t1[d1+1,d2+1]  
      r  <- cumsum(c(0, cm[-d2]))+(1+cm)/2  
      T  <- sum(cm[1:d2]^3-cm[1:d2])/(N^3-N)  
      f  <- 12*((N*(N+1)*(1-T))^(-1))  
      KW <- f*sum(((t%*%r-rm[1:d1]*(N+1)/2)^2)/rm[1:d1])  
      df <- (sum(rm>0)-1)  
      return(c(df=df, KW=KW))  
    }      
    ans    <- apply(m, 1, oneslice, d1, d2)      
    obs.KW <- sum(ans[2,])  
    df     <- sum(ans[1,])  
    P      <- 1 - pchisq(obs.KW,df)  
    return(list(KW=obs.KW, df=df, P=P))  
  }   
  
  ## Calculates the jonckheere-terpstra test, its mean and asymptotic P given an ftable (m).  
  jt <- function(m, d1, d2) {   
    oneslice<-function(t,d1,d2){  
      dim(t)<-c(d1,d2)
      t1 <-addmargins(t)  
      cm <-t1[d1+1,1:d2]
      rm <-t1[1:d1,d2+1]
      N  <-t1[d1+1,d2+1]  
      W  <-c()  
      T  <-c()  
      R  <-c()  
      for(i in c(2:d1)){  
        for(j in c(1:(i-1))){  
          if(i != j){  
            T<-c(t[i,],T)  
            R<-c((rm[i]*(rm[i]+1)/2),R)  
            W<-c(c(cumsum(c(0,t[i,-d2]+t[j,-d2])))+ c((((t[i,1:d2]+t[j,1:d2])+1)/2)),W)  
          }  
          W<-c(W)  
          T<-c(T)  
          R<-c(R)  
        }}  
      JT  <-sum(W*T)-sum(R)  
      EJT <-sum(N^2-sum(rm^2))/4  
      U1  <-N*(N-1)*(2*N+5)-sum(rm*(rm-1)*(2*rm+5))-sum(cm*(cm-1)*(2*cm+5))  
      U2  <-sum(rm*(rm-1)*(rm-2))*sum((cm)*(cm-1)*(cm-2))  
      U3  <-sum((rm)*(rm-1))*sum((cm)*(cm-1))  
      t1  <-72  
      t2  <-36*N*(N-1)*(N-2)  
      t3  <-8*N*(N-1)  
      VJT <-(U1/t1)+(U2/t2)+(U3/t3)   
      return(c(JT, EJT,VJT))  
    }     
    ans <- apply(m, 1, oneslice, d1, d2)      
    JT  <- sum(ans[1,])  
    EJT <- sum(ans[2,])  
    VJT <- sum(ans[3,])  
    P   <- 2*(1 - pnorm(abs(JT-EJT), sd=sqrt(VJT)))  
    return(list(JT=JT, EJT=EJT, P=P))  
  }   
  
  ## Returns row and column marginal totals for a given stratum.     
  rcsum <- function(t,d1,d2) {  
    dim(t) <- c(d1,d2)
    t1 <- addmargins(t)   
    cm <- t1[d1+1,1:d2]
    rm <- t1[1:d1,d2+1]   
    return(c(rm,cm))  
  }  
  
  ## Returns the deviances of Nsim random RxC tables with given margins  
  rdev <- function(tots, d1, d2, Nsim) {  
    rm <- tots[1:d1]
    cm <- tots[(d1+1):(d1+d2)]
    N  <- sum(rm)        
    fv <- (rm %o% cm)/N  
    tablist <- r2dtable(Nsim, rm, cm)  
    return(sapply(tablist, function(t) 2*sum(t*log(t/fv), na.rm=T)))  
  }  
  
  ## Returns the wilcoxon rank sum stats of Nsim random RxC tables with given margins  
  wdev <- function(tots, d1, d2, Nsim) {  
    rm <- tots[1:d1]
    cm <- tots[(d1+1):(d1+d2)]     
    r  <- cumsum(c(0, cm[-d2]))+(1+cm)/2    
    tablist <- r2dtable(Nsim, rm, cm)  
    return(sapply(tablist, function(t) sum(r*t[1,1:d2])))    
  }  
  
  ## Returns the kruskal stats of Nsim random RxC tables with given margins  
  kdev <- function(tots, d1, d2, Nsim) {  
    rm <- tots[1:d1]
    cm <- tots[(d1+1):(d1+d2)]
    N <- sum(rm)        
    r <- cumsum(c(0, cm[-d2]))+(1+cm)/2  
    T <- sum(cm[1:d2]^3-cm[1:d2])/(N^3-N)  
    f <- 12*((N*(N+1)*(1-T))^(-1))    
    tablist <- r2dtable(Nsim, rm, cm)  
    return(sapply(tablist, function(t) f*sum(((t[,1:d2]%*%r-rm[1:d1]*(N+1)/2)^2)/rm[1:d1])))  
  }  
  
  ## Returns the jonckheere terpstra rank sum stats of Nsim random RxC tables with given margins  
  jtdev <- function(tots, d1, d2, Nsim) {  
    rm <- tots[1:d1]
    cm <- tots[(d1+1):(d1+d2)]
    N <- sum(rm)  
    U<-function(t,d1,d2){  
      W<-c()  
      T<-c()  
      R<-c()  
      for(i in c(2:d1)){  
        for(j in c(1:(i-1))){  
          if(i != j){  
            T<-c(t[i,],T)  
            R<-c((rm[i]*(rm[i]+1)/2),R)  
            W<-c(c(cumsum(c(0,t[i,-d2]+t[j,-d2])))+ c((((t[i,1:d2]+t[j,1:d2])+1)/2)),W)  
          }  
          W<-c(W)  
          T<-c(T)  
          R<-c(R)  
        }}  
      return(sum(W*T)-sum(R))  
    }  
    tablist <- r2dtable(Nsim, rm, cm)  
    ans<-sapply(tablist, U,d1,d2)  
    return(c(ans)) 
  }  
  
  if (!(class(dataset) %in% c("data.frame", "table", "array", "matrix"))){
    stop("dataset incorrectly specified")
  }
  
  if (!(test %in% c("deviance", "wilcoxon", "kruskal", "jt"))){
    stop("test incorrectly specified")
  }
  if (class(dataset)=="data.frame") {
    d1 <- nlevels(dataset[,c1])
    d2 <- nlevels(dataset[,c2])
    ds <- dataset[,c(c1,c2,S)]
  } else {  
    d1 <- dim(dataset)[c1]
    d2 <- dim(dataset)[c2]
    ds <- apply(dataset, c(c1,c2,S), sum)
  }

  if ((d1<=1) | (d2<=1))
    stop("invalid factor(s)") 
  if (is.null(S))
    rv <- NULL
  else
    rv <- 3:(length(S)+2)
  
  ft      <- ftable(ds, col.vars=2:1, row.vars <- rv)    
  dim(ft) <- c(length(ft)/(d1*d2), d1*d2) 

  switch(test,
         "deviance"={obs <- LRT(ft, d1, d2)},
         "wilcoxon"={obs <- wilcoxon(ft, d1, d2)},
         "kruskal" ={obs <- kruskal(ft,d1,d2)},
         "jt"      ={obs <- jt(ft,d1,d2)})
  
  if (N>0){
    rcsums <- apply(ft, 1, rcsum, d1, d2)
    switch(test,
           "deviance"={
             strata.stats <- apply(rcsums, 2, rdev, d1, d2, N) 
             mc.P <- sum(rowSums(strata.stats) >= obs$deviance)/N
             obs <- c(obs, montecarlo.P=mc.P)   
           },
           "wilcoxon"={
             strata.stats <- apply(rcsums, 2, wdev, d1, d2, N)
             mc.P <- sum(abs(rowSums(strata.stats)-obs$EW) >= abs(obs$W-obs$EW))/N
             obs <- c(obs, montecarlo.P=mc.P)
           },
           "kruskal" ={
             strata.stats <- apply(rcsums, 2, kdev, d1, d2, N)
             mc.P <- sum((rowSums(strata.stats) >= obs$KW))/N
             obs <- c(obs, montecarlo.P=mc.P)             
           },
           "jt"      ={
             strata.stats <- apply(rcsums, 2, jtdev, d1, d2, N)
             mc.P <- sum((rowSums(strata.stats)-obs$EJT) >= abs(obs$JT-obs$EJT))/N 
             obs <- c(obs, montecarlo.P=mc.P)             
           })
  }
  
  obs$set <- set
  obs
}

##   if (test=="deviance")
##     obs <- LRT(ft, d1, d2)
##   else {
##     if (test=="wilcoxon")
##       obs <- wilcoxon(ft, d1, d2)
##     else {
##       if (test=="kruskal")
##         obs <- kruskal(ft,d1,d2)
##       else
##         obs <- jt(ft,d1,d2)
##     }
##   }


##   if (N>0){
##     rcsums <- apply(ft, 1, rcsum, d1, d2)
##     if (test=="deviance") {
##       strata.stats <- apply(rcsums, 2, rdev, d1, d2, N) 
##       mc.P <- sum(rowSums(strata.stats) >= obs$deviance)/N
##       obs <- c(obs, montecarlo.P=mc.P)   
##     } else {
##       if  (test=="wilcoxon") {
##         strata.stats <- apply(rcsums, 2, wdev, d1, d2, N)
##         mc.P <- sum(abs(rowSums(strata.stats)-obs$EW) >= abs(obs$W-obs$EW))/N
##         obs <- c(obs, montecarlo.P=mc.P)
##       } else {
##         if  (test=="kruskal") {
##           strata.stats <- apply(rcsums, 2, kdev, d1, d2, N)
##           mc.P <- sum((rowSums(strata.stats) >= obs$KW))/N
##           obs <- c(obs, montecarlo.P=mc.P)
##         } else {
##           strata.stats <- apply(rcsums, 2, jtdev, d1, d2, N)
##           mc.P <- sum((rowSums(strata.stats)-obs$EJT) >= abs(obs$JT-obs$EJT))/N 
##           obs <- c(obs, montecarlo.P=mc.P)
##         }}
##     }
##   }


# Functions for calculating exact conditional independence tests for discrete data.
# Presently supported: deviance and wilcoxon


.CI.exact <- function(c1,c2, S=NULL, dataset, test="deviance", N=0) {

       # Calculates the deviance, degrees of freedom and asymptotic P given an ftable (m).   
       LRT <- function(m, d1, d2) {  
         oneslice <- function(t,d1,d2) {  
             dim(t) <- c(d1,d2); t1 <- addmargins(t)   
             cm <- t1[d1+1,1:d2]; rm <- t1[1:d1,d2+1]; N<- t1[d1+1,d2+1]  
             df <- (sum(cm>0)-1)*(sum(rm>0)-1)  
             dev <- 0  
             if (df>0) {fv <- (rm %o% cm)/N; dev <- 2*sum(t*log(t/fv), na.rm=T)}  
             return(c(df=df, dev=dev))  
             }    
         ans <- apply(m, 1, oneslice, d1, d2)      
         obs.deviance <- sum(ans[2,])  
         df <- sum(ans[1,])  
         P <- 1 - pchisq(obs.deviance, df)  
         return(list(deviance=obs.deviance, df=df, P=P))  
       }  
                 
       # Calculates the wilcoxon test, its mean and asymptotic P given an ftable (m).            
       wilcoxon <- function(m, d1, d2) {  
         oneslice <- function(t,d1,d2) {  
             dim(t) <- c(d1,d2); t1 <- addmargins(t)   
             cm <- t1[d1+1,1:d2]; rm <- t1[1:d1,d2+1]; N<- t1[d1+1,d2+1]  
             r <- cumsum(c(0, cm[-d2]))+(1+cm)/2  
             W <- sum(r*t[1,])  
             EW <- (rm[1]/N)*sum(r*cm)  
             VW <- (rm[1]*rm[2]/(N*(N-1)))*sum(((r-EW/rm[1])^2)*cm)  
             return(c(W, EW, VW))  
             }      
         ans <- apply(m, 1, oneslice, d1, d2)      
         W <- sum(ans[1,])  
         EW <- sum(ans[2,])  
         VW <- sum(ans[3,])  
         P <- 2*(1 - pnorm(abs(W-EW), sd=sqrt(VW)))  
         return(list(W=W, EW=EW, P=P))  
       }       
         
       # Calculates the kruskal-wallis test, degrees of freedom and asymptotic P given an ftable (m).            
       kruskal <- function(m, d1, d2) {   
         oneslice <- function(t,d1,d2) {  
             dim(t) <- c(d1,d2); t1 <- addmargins(t)   
             cm <- t1[d1+1,1:d2]; rm <- t1[1:d1,d2+1]; N<- t1[d1+1,d2+1]  
             r <- cumsum(c(0, cm[-d2]))+(1+cm)/2  
           T <- sum(cm[1:d2]^3-cm[1:d2])/(N^3-N)  
           f <- 12*((N*(N+1)*(1-T))^(-1))  
             KW <- f*sum(((t%*%r-rm[1:d1]*(N+1)/2)^2)/rm[1:d1])  
             df <- (sum(rm>0)-1)  
             return(c(df=df, KW=KW))  
             }      
         ans <- apply(m, 1, oneslice, d1, d2)      
         obs.KW <- sum(ans[2,])  
         df<- sum(ans[1,])  
         P <- 1 - pchisq(obs.KW,df)  
         return(list(KW=obs.KW, df=df, P=P))  
       }   
         
       #Calculates the jonckheere-terpstra test, its mean and asymptotic P given an ftable (m).  
       jt <- function(m, d1, d2) {   
         oneslice<-function(t,d1,d2){  
           dim(t)<-c(d1,d2); t1<-addmargins(t)  
           cm<-t1[d1+1,1:d2]; rm<-t1[1:d1,d2+1]; N<-t1[d1+1,d2+1]  
           W<-c()  
           T<-c()  
           R<-c()  
               for(i in c(2:d1)){  
                   for(j in c(1:(i-1))){  
                       if(i != j){  
                       T<-c(t[i,],T)  
                       R<-c((rm[i]*(rm[i]+1)/2),R)  
                       W<-c(c(cumsum(c(0,t[i,-d2]+t[j,-d2])))+ c((((t[i,1:d2]+t[j,1:d2])+1)/2)),W)  
                       }  
               W<-c(W)  
               T<-c(T)  
               R<-c(R)  
               }}  
           JT<-sum(W*T)-sum(R)  
           EJT<-sum(N^2-sum(rm^2))/4  
           U1<-N*(N-1)*(2*N+5)-sum(rm*(rm-1)*(2*rm+5))-sum(cm*(cm-1)*(2*cm+5))  
           U2<-sum(rm*(rm-1)*(rm-2))*sum((cm)*(cm-1)*(cm-2))  
           U3<-sum((rm)*(rm-1))*sum((cm)*(cm-1))  
           t1<-72  
           t2<-36*N*(N-1)*(N-2)  
           t3<-8*N*(N-1)  
           VJT<-(U1/t1)+(U2/t2)+(U3/t3)   
           return(c(JT, EJT,VJT))  
           }     
         ans <- apply(m, 1, oneslice, d1, d2)      
         JT <- sum(ans[1,])  
         EJT <- sum(ans[2,])  
         VJT <- sum(ans[3,])  
         P<- 2*(1 - pnorm(abs(JT-EJT), sd=sqrt(VJT)))  
         return(list(JT=JT, EJT=EJT, P=P))  
       }   
         
         
         
       # Returns row and column marginal totals for a given stratum.     
       rcsum <- function(t,d1,d2) {  
             dim(t) <- c(d1,d2); t1 <- addmargins(t)   
             cm <- t1[d1+1,1:d2]; rm <- t1[1:d1,d2+1]   
             return(c(rm,cm))  
       }  
         
       # Returns the deviances of Nsim random RxC tables with given margins  
       rdev <- function(tots, d1, d2, Nsim) {  
         rm <- tots[1:d1]; cm <- tots[(d1+1):(d1+d2)]; N <- sum(rm)        
         fv <- (rm %o% cm)/N  
         tablist <- r2dtable(Nsim, rm, cm)  
         return(sapply(tablist, function(t) 2*sum(t*log(t/fv), na.rm=T)))  
       }  
         
       # Returns the wilcoxon rank sum stats of Nsim random RxC tables with given margins  
       wdev <- function(tots, d1, d2, Nsim) {  
         rm <- tots[1:d1]; cm <- tots[(d1+1):(d1+d2)]     
         r <- cumsum(c(0, cm[-d2]))+(1+cm)/2    
         tablist <- r2dtable(Nsim, rm, cm)  
         return(sapply(tablist, function(t) sum(r*t[1,1:d2])))    
       }  
         
       # Returns the kruskal stats of Nsim random RxC tables with given margins  
       kdev <- function(tots, d1, d2, Nsim) {  
         rm <- tots[1:d1]; cm <- tots[(d1+1):(d1+d2)]; N <- sum(rm)        
         r <- cumsum(c(0, cm[-d2]))+(1+cm)/2  
         T <- sum(cm[1:d2]^3-cm[1:d2])/(N^3-N)  
         f <- 12*((N*(N+1)*(1-T))^(-1))    
         tablist <- r2dtable(Nsim, rm, cm)  
         return(sapply(tablist, function(t) f*sum(((t[,1:d2]%*%r-rm[1:d1]*(N+1)/2)^2)/rm[1:d1])))  
       }  
         
       # Returns the jonckheere terpstra rank sum stats of Nsim random RxC tables with given margins  
       jtdev <- function(tots, d1, d2, Nsim) {  
         rm <- tots[1:d1]; cm <- tots[(d1+1):(d1+d2)]; N <- sum(rm)  
         U<-function(t,d1,d2){  
           W<-c()  
           T<-c()  
           R<-c()  
               for(i in c(2:d1)){  
                   for(j in c(1:(i-1))){  
                       if(i != j){  
                       T<-c(t[i,],T)  
                       R<-c((rm[i]*(rm[i]+1)/2),R)  
                       W<-c(c(cumsum(c(0,t[i,-d2]+t[j,-d2])))+ c((((t[i,1:d2]+t[j,1:d2])+1)/2)),W)  
                       }  
               W<-c(W)  
               T<-c(T)  
               R<-c(R)  
               }}  
         return(sum(W*T)-sum(R))  
         }  
         tablist <- r2dtable(Nsim, rm, cm)  
         ans<-sapply(tablist, U,d1,d2)  
         return(c(ans)) 
         }  

  if (!(class(dataset) %in% c("data.frame", "table"))) stop("dataset incorrectly specified")
  if (!(test %in% c("deviance", "wilcoxon", "kruskal", "jt"))) stop("test incorrectly specified")
  if (class(dataset)=="data.frame") {
       d1 <- nlevels(dataset[,c1]); d2 <- nlevels(dataset[,c2]); ds <- dataset[,c(c1,c2,S)]
     } else { d1 <- dim(dataset)[c1]; d2 <- dim(dataset)[c2]; ds <- apply(dataset, c(c1,c2,S), sum)}
  if ((d1<=1) | (d2<=1)) stop("invalid factor(s)") 
  if (is.null(S)) rv <- NULL else rv <- 3:(length(S)+2)
  ft <- ftable(ds, col.vars=2:1, row.vars <- rv)    
  dim(ft) <- c(length(ft)/(d1*d2), d1*d2) 
  if (test=="deviance") obs <- LRT(ft, d1, d2) else {
  if (test=="wilcoxon") obs <- wilcoxon(ft, d1, d2) else {
  if (test=="kruskal")  obs <- kruskal(ft,d1,d2) else obs <- jt(ft,d1,d2)}}
  if (N==0) return(obs) else {
    rcsums <- apply(ft, 1, rcsum, d1, d2)
    if (test=="deviance") {
       strata.stats <- apply(rcsums, 2, rdev, d1, d2, N) 
       mc.P <- sum(rowSums(strata.stats) >= obs$deviance)/N
       return(c(obs, montecarlo.P=mc.P))   
    } else {
    if  (test=="wilcoxon") {
        strata.stats <- apply(rcsums, 2, wdev, d1, d2, N)
        mc.P <- sum(abs(rowSums(strata.stats)-obs$EW) >= abs(obs$W-obs$EW))/N
        return(c(obs, montecarlo.P=mc.P))
    } else {
    if  (test=="kruskal") {
        strata.stats <- apply(rcsums, 2, kdev, d1, d2, N)
        mc.P <- sum((rowSums(strata.stats) >= obs$KW))/N
        return(c(obs, montecarlo.P=mc.P))
    } else {
        strata.stats <- apply(rcsums, 2, jtdev, d1, d2, N)
        mc.P <- sum((rowSums(strata.stats)-obs$EJT) >= abs(obs$JT-obs$EJT))/N 
        return(c(obs, montecarlo.P=mc.P))
    }}
    }}
}
