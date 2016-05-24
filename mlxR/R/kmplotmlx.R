#' Kaplan Meier plot
#' 
#' Plot empirical survival functions using the Kaplan Meier estimate.
#' 
#' See http://simulx.webpopix.org/mlxr/kmplotmlx/ for more details.
#' @param r a data frame with a column \samp{id}, a column \samp{time}, 
#' a column with values and possibly a column \samp{group}.
#' @param index an integer: \code{index=k} means that the survival function for the k-th event is displayed. 
#' Default is \code{index=1}.
#' @param level a number between 0 and 1:  confidence interval level. 
#' @examples
#' \dontrun{
#' tteModel1 <- inlineModel("
#'   [LONGITUDINAL]
#'   input = {beta,lambda}  
#'   EQUATION:
#'   h=(beta/lambda)*(t/lambda)^(beta-1)
#'   DEFINITION:
#'   e = {type=event, maxEventNumber=1, rightCensoringTime=70, hazard=h}
#'   ")
#'
#'   p1   <- c(beta=2.5,lambda=50)
#'   e    <- list(name='e', time=0)
#'   res1 <- simulx(model=tteModel1, parameter=p1, output=e, group=list(size=100))
#'   pl1  <- kmplotmlx(res1$e,level=0.95)
#'   print(pl1)
#' 
#'   p2   <- c(beta=2,lambda=45)
#'   g1   <- list(size=50, parameter=p1)
#'   g2   <- list(size=100, parameter=p2)
#'   res2 <- simulx(model=tteModel1, output=e, group=list(g1,g2))
#'   pl2  <- kmplotmlx(res2$e)
#'   print(pl2)
#' }
#' @importFrom ggplot2 ggplot geom_point theme aes geom_line xlab ylab
#' @importFrom stats qnorm
#' @export         
kmplotmlx  <-  function(r, index=1, level=NULL)
{ 
  r.name <- attr(r,"name")
  names(r)[names(r)==r.name] <- "y"
  
  N <- length(unique(r$id))
  r0 <- r1 <- NULL
  for (i in seq(1,N)){
    ri <- r[r$id==i,]
    cyi <- cumsum(ri$y)
    ri$y <- 1
    if (any(index<=cyi)){
      it <- min(which(index<=cyi))
      r1 <- rbind(r1,ri[it,]) 
    }else{
      di <- dim(ri)[1]
      r0 <- rbind(r0,ri[di,])  
    }
  }
  if (!is.null(r0)){
    names(r0)[names(r0)=="y"] <- "c"
    r0$d <- 0
  }
  if (!is.null(r1)){
    names(r1)[names(r1)=="y"] <- "d"
    r1$c <- 0
  }
  re <- rbind(r1,r0)
  re <- re[with(re, order(time)), ]
  
  if (any( "group" %in% names(re) )){
    g=as.numeric(levels(re$group))[re$group]
  }else{
    g=rep(1,length(re$id))
  }
  ng=max(g)
  re$id <- NULL
  
  S <- Se <- T <- G <- NULL
  S0 <- T0 <- G0 <- NULL
  t0=min(r$time)
  for (kg in seq(1,ng)){
    rk<-re[g==kg,]
    Nk <- dim(rk)[1]
    ut <- uniquemlx(rk$time)
    tu <- ut$uniqueValue
    iu <- ut$sortIndex
    nt <- length(tu)
    ru <- data.frame(time=c(t0,tu),d=0,c=0)
    nj <- Nk
    for (j in seq(1,nt)){
      ru$c[j+1] <- sum(rk$c[iu==j])
      ru$d[j+1] <- sum(rk$d[iu==j])
    }
    nj <- Nk
    sek <- vector(length=nt+1)
    Sk <- vector(length=nt+1)
    Sk[1] <- 1
    sek[1] <- 0
    V <- 0
    for (j in seq(2,nt+1)){
      nj <- nj - ru$c[j-1] - ru$d[j-1]
      pj <- (nj - ru$d[j])/nj
      Sk[j] <- Sk[j-1]*pj
      V <- V + ru$d[j]/nj/(nj - ru$d[j])
      #       sek[j] <- Sk[j]*sqrt(V)
      # sek[j] <- sqrt(V)/log(Sk[j])  # Kalbfleisch and Prentice (2002) 
      sek[j] <- sqrt(V)/(1-Sk[j])  # logit 
    }
    sek[which(is.nan(sek))] <- 0
    S <- c(S,rep(Sk,each=2))
    S <- S[-length(S)]
    Se <- c(Se,rep(sek,each=2))
    Se <- Se[-length(Se)]
    
    Tk<-c(t0,rep(tu,each=2))
    T<-c(T,Tk)
    G<-c(G,rep(kg,2*nt+1))
    
    i0 <- which(ru$c>0)
    if (length(i0)>0){
      T0 <- c(T0,ru$time[i0])
      S0 <- c(S0,Sk[i0])
      G0 <- c(G0, rep(kg, length(i0) ))
    }
    
  }
  
  group=factor(G)
  if (!is.null(level)){
    alpha <- (1-level)/2
    #  S1 <- pmax(S + Se*qnorm(alpha),0)
    #  S2 <- pmin(S + Se*qnorm(1-alpha),1)
    #     s1 <- log(-log(S)) + Se*qnorm(alpha)  # Kalbfleisch and Prentice (2002) 
    #     s2 <- log(-log(S)) + Se*qnorm(1-alpha)
    #     S1 <- exp(-exp(s1))
    #     S2 <- exp(-exp(s2))
    s1 <- log(S/(1-S)) + Se*qnorm(alpha)  # logit 
    s2 <- log(S/(1-S)) + Se*qnorm(1-alpha)
    S1 <- 1/(1+exp(-s1))
    S2 <- 1/(1+exp(-s2))
    D=data.frame(T,S,S1,S2,group)
  }else{
    D=data.frame(T,S,group)
  }
  
  
  plot1=ggplotmlx(data=D) +  geom_line(aes(x=T, y=S, colour=group), size=1)
  if (!is.null(level)){
    plot1=plot1+geom_line(aes(x=T, y=S1, colour=group), linetype="dotted", size=0.8) +
      geom_line(aes(x=T, y=S2, colour=group), linetype="dotted", size=0.8)
  }
  plot1 <- plot1 + xlab("time") + ylab("survival") 
  if (length(i0)>0){
    group <- factor(G0)
    D0 <- data.frame(T0,S0,group)
    plot1 <- plot1 + geom_point(data=D0, aes(x=T0,y=S0, colour=group), size=4)
  }
  if (ng>1){
    plot1 <- plot1 + theme(legend.position=c(0.1,0.15))
  }else{
    plot1 <- plot1 + theme(legend.position="none")
  }  
  return(plot1)
}



#--------------------------------------------------------
uniquemlx <- function(x) 
{ 
  d <- !duplicated(x) 
  u=list(uniqueValue=x[d], firstIndex=which(d), sortIndex=match(x,x[d])) 
  return(u)
}

