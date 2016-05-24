depth.projection <- function(x, data, method = "random", num.directions = 1000, seed = 0){
  if (seed!=0) set.seed(seed)
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.matrix(x) 
      && is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  if (method == "random"){
    dt <- as.vector(t(data))
    z <- as.vector(t(x))
    m <- nrow(x)
    d <- ncol(data)
    c <- nrow(data)
    q <- 1
    k <- num.directions
    newDirs <- 1
    rez <- .C("ProjectionDepth", as.double(dt), 
              as.double(z), 
              as.integer(m), 
              as.integer(d), 
              as.integer(c), 
              as.integer(q), 
              dirs=double(k*d), 
              prjs=double(k*c), 
              as.integer(k), 
              as.integer(1), 
              as.integer(seed),
              dps=double(m*q))
    return (rez$dps)
  }
  if (method == "linearize"){
    depths <- .zdepth(data, x)
    return (1/(1 + depths))
  }
}

depth.space.projection <- function(data, cardinalities, method = "random", num.directions = 1000, seed = 0){
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.vector(cardinalities, mode = "numeric") 
      || is.na(min(cardinalities)) 
      || sum(.is.wholenumber(cardinalities)) != length(cardinalities) 
      || min(cardinalities) <= 0 
      || sum(cardinalities) != nrow(data)){
    stop("Argument \"cardinalities\" should be a vector of cardinalities of the classes in \"data\" ")
  }
  if (sum(cardinalities < ncol(data) + 1) != 0){
    stop("Not in all classes sufficiently enough objetcs")
  }
  
  depth.space <- NULL
  for (i in 1:length(cardinalities)){
    pattern <- data[(1 + sum(cardinalities[0:(i - 1)])):sum(cardinalities[1:i]),]
    pattern.depths <- depth.projection (data, pattern, method, num.directions, seed)
    depth.space <- cbind(depth.space, pattern.depths, deparse.level = 0)
  }
  
  return (depth.space)
}


################################################################################
# R-codes of this function written by Subhajit Dutta,
# taken from http://www.isical.ac.in/~tijahbus/GAM/r_wilcox.txt.
################################################################################
.zdepth<-function(m,pts=m,zloc=median,zscale=mad){
  # 
  # Compute depth of points as in Zuo, Annals, 2003
  #
  if(!is.matrix(m))stop("argument m should be a matrix")
  if(!is.matrix(pts))stop("argument pts should be a matrix")
  if(ncol(m)!=ncol(pts))stop("Number of columns for m and pts are not equal")
  np<-ncol(m)
  val<-NA
  for(i in 1:nrow(pts)){
    pval<-pts[i,]
    START<-rep(1,np)/sqrt(np)
    temp<-.nelderv2(m,np,FN=.zdepth.sub,START=START,zloc=zloc,zscale=zscale,pts=pval)
    temp<-temp/sqrt(sum(temp^2))
    y<-t(t(m)*temp)
    y<-apply(y,1,sum)
    ppro<-sum(pval*temp)
    val[i]<-abs(ppro-zloc(y))/zscale(y)
  }
  val
}

################################################################################
# R-codes of this function written by Subhajit Dutta,
# taken from http://www.isical.ac.in/~tijahbus/GAM/r_wilcox.txt.
################################################################################
.zdepth.sub<-function(x,theta,zloc=median,zscale=mad,pts=NA){
  theta<-theta/sqrt(sum(theta^2))
  temp<-t(t(x)*theta)
  ppro<-sum(t(t(pts)*theta))
  yhat<-apply(temp,1,sum)
  val<-0-abs(ppro-zloc(yhat))/zscale(yhat)
  val
}

################################################################################
# R-codes of this function written by Subhajit Dutta,
# taken from http://www.isical.ac.in/~tijahbus/GAM/r_wilcox.txt.
################################################################################
.nelderv2<-function(x,N,FN,START=c(rep(1,N)),STEP=c(rep(1,N)),
                   XMIN=c(rep(0,N)),XSEC=c(rep(0,N)),...){
  #     NELDER-MEAD method for minimzing a function
  #
  #     TAKEN FROM OLSSON, J QUALITY TECHNOLOGY, 1974, 6, 56.
  #
  #     x= n by p matrix containing data; it is used by 
  #        function to be minimized.
  #     N= number of parameters 
  #
  #     FN=the function to be minimized
  #     FORM: FN(x,theta), theta is vector containing
  #     values for N parameters.
  #
  #     START = starting values.
  #     STEP=initial step.
  #     This function returns the N values for theta that minimize FN
  #
  ICOUNT<-500
  REQMIN<-.0000001
  NN<-N+1
  P<-matrix(NA,nrow=N,ncol=NN)
  P[,NN]<-START
  PBAR<-NA
  RCOEFF<-1
  ECOEFF<-2
  CCOEFF<-.5
  KCOUNT<-ICOUNT
  ICOUNT<-0
  DABIT<-2.04067e-35
  BIGNUM<-1.e38
  KONVGE<-5
  XN<-N
  DN<-N
  Y<-rep(0,NN)
  Y[NN]<-FN(x,START,...)
  ICOUNT<-ICOUNT+1
  for(J in 1:N){
    DCHK<-START[J]
    START[J]<-DCHK+STEP[J]
    for(I in 1:N){
      P[I,J]<-START[I]
    }
    Y[J]<-FN(x,START,...)
    ICOUNT<-ICOUNT+1
    START[J]<-DCHK
  }
  I1000<-T
  while(I1000){
    YLO<-Y[1]
    YNEWLO<-YLO
    ILO<-1
    IHI<-1
    for(I in 2:NN){
      if(Y[I] <  YLO){
        YLO<-Y[I]
        ILO<-I}
      if(Y[I] > YNEWLO){
        YNEWLO<-Y[I]
        IHI<-I}
    }
    DCHK<-(YNEWLO+DABIT)/(YLO+DABIT)-1
    if(abs(DCHK) < REQMIN){
      I1000<-F
      next
    }
    KONVGE<-KONVGE-1
    if(KONVGE == 0){
      KONVGE<-5
      for(I in 1:N){
        COORD1<-P[I,1]
        COORD2<-COORD1
        for(J in 2:NN){
          if(P[I,J] < COORD1)COORD1<-P[I,J]
          if(P[I,J] > COORD2)COORD2<-P[I,J]
        }     # 2010 CONTINUE
        DCHK<-(COORD2+DABIT)/(COORD1+DABIT)-1
        if(abs(DCHK) > REQMIN)break  
      } 
    }       
    if(ICOUNT >= KCOUNT){
      I1000<-F
      next
    }
    for(I in 1:N){
      Z<-0.0
      Z<-sum(P[I,1:NN]) # 6
      Z<-Z-P[I,IHI]
      PBAR[I]<-Z/DN
    }
    PSTAR<-(1.+RCOEFF)*PBAR-RCOEFF*P[,IHI]
    YSTAR<-FN(x,PSTAR,...)
    ICOUNT<-ICOUNT+1
    if(YSTAR < YLO && ICOUNT >= KCOUNT){
      P[,IHI]<-PSTAR
      Y[IHI]<-YSTAR
      next 
    }
    IFLAG<-T   
    if(YSTAR < YLO){
      P2STAR<-ECOEFF*PSTAR+(1-ECOEFF)*PBAR
      Y2STAR<-FN(x,P2STAR,...)
      ICOUNT<-ICOUNT+1
      if(Y2STAR >= YSTAR){ 
        P[,IHI]<-PSTAR
        Y[IHI]<-YSTAR
        next #In essence, go to 19 which goes to 1000
      }
      IFLAG<-T
      while(YSTAR < Y[IHI]){
        P[,IHI]<-P2STAR
        Y[IHI]<-Y2STAR
        IFLAG<-F
        break  
        L<-sum(Y[1:NN] > YSTAR)
        if(L > 1){
          P[,IHI]<-PSTAR
          Y[IHI]<-YSTAR
          IFLAG<-T
          break 
        }
        if(L > 1)break # go to 19
        if(L != 0){
          P[1:N,IHI]<-PSTAR[1:N]
          Y[IHI]<-YSTAR
        }
        I1000<-F
        break
        if(ICOUNT >= KCOUNT){
          I1000<-F
          next
        }
        P2STAR[1:N]<-CCOEFF*P[1:N,IHI]+(1-CCOEFF)*PBAR[1:N]
        Y2STAR<-FN(x,P2STAR,...)
        ICOUNT<-ICOUNT+1
      }   # END WHILE 
    }
    if(IFLAG){
      for(J in 1:NN){
        P[,J]=(P[,J]+P[,ILO])*.5
        XMIN<-P[,J]
        Y[J]<-FN(x,XMIN,...)
      }
      ICOUNT<-ICOUNT+NN
      if(ICOUNT < KCOUNT)next
      I1000<-F
      next
    }
    P[1:N,IHI]<-PSTAR[1:N]
    Y[IHI]<-YSTAR
  } 
  for(J in 1:NN){
    XMIN[1:N]<-P[1:N,J]
  }
  Y[J]<-FN(x,XMIN,...)
  YNEWLO<-BIGNUM
  for(J in 1:NN){
    if (Y[J] < YNEWLO){
      YNEWLO<-Y[J]
      IBEST<-J
    }}
  Y[IBEST]<-BIGNUM
  YSEC<-BIGNUM
  for(J in 1:NN){
    if(Y[J] < YSEC){
      YSEC<-Y[J]
      ISEC<-J
    }}
  XMIN[1:N]<-P[1:N,IBEST]
  XSEC[1:N]<-P[1:N,ISEC]
  XMIN
}