###Charles
### Playing with logcondens package and activeSetLogcon method
### with goal of converting it for the case where the location of mode is
### fixed.

##library(logcondens)


##returns (n-1)x1 matrix always;
##note if k==1, you don't subtract a convexity constraints whereas otherwise you do.
LocalConvexity.mode <- function (z, phi, k=NULL){
  ##k is index of a in z.
  n <- length(z)
  dphi <- diff(phi)
  dz <- diff(z) ##unnec
  deriv <- dphi/dz
  conv <- rep(0,n);
  conv[2:(n - 1)] <- diff(deriv[1:(n - 1)])
  conv <- matrix(conv,ncol=1);
  ##if (is.na(k)||is.null(k)) conv[2:(n-1)] <- diff(deriv[1:(n - 1)])##==diff(deriv); and why 'conv[2:n-1]'?
  ##else
  if (k>=2 && k<n){ ##diff(1)=numeric(0), so works for k==2, n==3
    ##conv <- c(diff(deriv[1:(k-1)]), dphi[(k-1):k], diff(deriv[k:(n-1)]));#diff(deriv[a:b]) checks a+1st,bth points for conv.
    mono <- c(-dphi[k-1], dphi[k]) ##<0 (if constraint holds)
  }
  else if (k==1){
    ##conv <- c(dphi[1], diff(deriv)); ##works for n==3 too.
    mono <- c(0,dphi[1]) ##<0 (if constraint holds)
  }
  else if (k==n){
    mono <- c(dphi[n-1],0)
  }
  shape <- list(conv=conv,mono=matrix(mono,ncol=1))
  return(shape)    
}


##AB is basis vectors; constant over each run.
LocalMLE.mode <- function (z, w, IsKnot, IsMIC, a, phi_o, prec, print=FALSE) 
{
  n <- length(z)
  k <- a$idx;
  r1 <- LocalCoarsen.mode(z, w, IsKnot, IsMIC,a);
  IsKnot <- syncIsKnot(IsKnot,IsMIC,k); ## REDUNDANT BUT J I C
  ##K <- ((1:n) * IsKnot)[IsKnot>0];
  K <- (1:n)[IsKnot>0];
  ## phi[as.logical(IsKnot)] == phi[K] right?
  res2 <- MLE.mode(r1$y,r1$constr, r1$w2, phi_o[K], print=print)
  ##phi <- LocalExtend(z, IsKnot, r1$y, res2$phi)
  phi <- LocalExtend(z, IsKnot, r1$y, res2$phi,
                     r1$constr)
  shape <- LocalConvexity.mode(z, phi,k=k)## * IsKnot
  shape$conv <- shape$conv * IsKnot
  shape$mono <- shape$mono * IsMIC
  H <- rep(0,n)
  HR <- rep(0,n)
  H.m <- matrix(nrow=2,ncol=1)
  ##AB <- getBasisVecs(z)
  ##DL <- getGrad.uncstr(z,w,phi)
  JJ <- (1:n) * IsKnot
  ##JJ[k] <- k ##we won't compute grad for k; OK b/c use H.m for that.
  JJ <- JJ[JJ > 0]
  p <- sum(JJ<k)  ##ks idx ; 
  p <- p+1 ##works for boundary.k==1,k==n. ##p \ge 2 UNLESS (k=1 AND z[1] is a knot)
  m <- length(JJ)
  if (is.null(a)){ ##why worry about this
    for (i in 1:(m - 1)) {
      if (JJ[i + 1] > JJ[i] + 1) {
        dtmp <- z[JJ[i + 1]] - z[JJ[i]]
        ind <- (JJ[i] + 1):(JJ[i + 1] - 1)
        mtmp <- length(ind)
        ztmp <- (z[ind] - z[JJ[i]])
        dstmp <- c(z[JJ[i]+1]-z[JJ[i]],diff(z[ind]),
                   z[JJ[i+1]]-z[JJ[i+1]-1])##len==mtmp+1
        wtmp <- w[ind]
        J01s <- J10(phi[ind],phi[ind-1]) * dstmp[1:mtmp]
        J10s <- J10(phi[ind],phi[ind+1]) * dstmp[2:(mtmp+1)]
        H[ind] <- cumsum(wtmp * ztmp) - ztmp * cumsum(wtmp) +
          ztmp * sum(wtmp * (1-ztmp/dtmp))
        jtmp <- - ztmp*cumsum(J01s+J10s) + cumsum(ztmp * (J01s + J10s)) +
          sum((J01s+J10s)*(1-ztmp/dtmp))*ztmp
        H[ind] <- H[ind] - jtmp
        H[IsKnot] <- 0; ##REDUNDANT but JNC.
      }
    }
  }
  else{
    ##H <- matrix(0,nrow=n,ncol=1); ##H[1,1] = 0  the constant vector.
    if (k!=1) { ##p != 1
      for (i in 1:(p-1)) {
        if (JJ[i + 1] > JJ[i] + 1) {
          ##This formula uses the fact that we're at the maximum over the knots.
          ## So outside of the nearest knots, the derivative is 0.
          ##deal with left-knot and right-knots at knots about the mode.
          cond1 <- r1$constr[1]!=r1$constr[2] && r1$constr[1]==1 && i==1
          cond2 <- r1$constr[1]!=r1$constr[2] && r1$constr[1]==i+1  && JJ[i+1]>JJ[i]+1 && JJ[i+2]>JJ[i+1]+1 #if cond2 false this won't be NULL
          if (cond1){
            LK <- 1 ## == JJ[1] also
            RK <- JJ[i+1]
          }
          else if (cond2){ ##m>2 here
            LK <- JJ[i+1]; ##Left Knot
            RK <- JJ[i+2]; ##Right Knot
          }
          else{ ##don't enter this if-block if p==1 (b/c then k==1)
            LK <- JJ[i+1]
            RK <- JJ[i+1]
          }
          ind <- (JJ[i] + 1):(RK - 1)
          mtmp <- length(ind)          
          ztmp <- (z[ind] - z[JJ[i]])
          dstmp <- c(z[JJ[i]+1]-z[JJ[i]],diff(z[ind]),
                     z[RK]-z[RK-1])##len==mtmp+1
          wtmp <- w[ind]
          J01s <- J10(phi[ind],phi[ind-1]) * dstmp[1:mtmp]
          J10s <- J10(phi[ind],phi[ind+1]) * dstmp[2:(mtmp+1)]
          H[ind] <- cumsum(wtmp * ztmp) - ztmp * cumsum(wtmp) 
          jtmp <- - ztmp*cumsum(J01s+J10s) + cumsum(ztmp * (J01s + J10s))
          if (cond1){## deals w/ "m==2" case.
            H0 <- -ztmp*w[1]
            J0 <- -ztmp*(J10(phi[1],phi[2])*(z[2]-z[1]))##z instead of ztmp etc            
          }
          else if (cond2){ ##may never happen, if constraint is on RHS of a
            ## m > 2 here.
            dtmp <- z[LK] - z[JJ[i]] ##not used if i==1 is constrained.
            tmptmp <- (JJ[i] + 1):(LK - 1)
            ind.lin <- (1:length(tmptmp))
            H0 <- ztmp*sum(wtmp[ind.lin] * (1-ztmp[ind.lin]/dtmp))
            J0 <- ztmp*sum((J01s[ind.lin]+J10s[ind.lin])*(1-ztmp[ind.lin]/dtmp))
            ##i <- i+1
            ## want to skip the next loop; i<-i+1 doesn't work tho
          }
          else{
            dtmp <- z[LK] - z[JJ[i]] ##not used if i==1 is constrained.
            H0 <- ztmp * sum(wtmp * (1-ztmp/dtmp))
            J0 <- sum((J01s+J10s)*(1-ztmp/dtmp))*ztmp
          }
          H[ind] <- H[ind] + H0  - jtmp - J0
          ## ## below used to be uncommented, now is commented to see if
          ## ## IsKnot is confusing L and R knots ...          
          H[IsKnot] <- 0; ##REDUNDANT but JNC.
          
          if (cond2)    break  ##would prefer i<-i+1; 
        } 
      }
      p2 <- p-1
      ##if (z[JJ[p]]==a$val) p2 <- p
    }
    else{
      p2 <- 1 ## 'a' might not be a knot.  So, p2 is firts index at or before a
    }
    for (i in p2:(m-1)) { ## care about x_j >= a; p-1 is 'softer' bound, if p==1, then p is hard bound
      if (JJ[i + 1] > JJ[i] + 1) {
        ##This formula uses the fact that we're at the maximum over the knots.
        ## So outside of the nearest knots, the derivative is 0.
        ##dtmp <- z[JJ[i + 1]] - z[JJ[i]]
        ##ind <- (JJ[i] + 1):(JJ[i + 1] - 1)

        cond1 <- r1$constr[1]!=r1$constr[2] && r1$constr[1]==m-1 && i==m-1 
        cond2 <- r1$constr[1]!=r1$constr[2] && r1$constr[1]==i  && JJ[i+1]>JJ[i]+1 && JJ[i+2]>JJ[i+1]
##########JUST rev the indices ! ! !         
        if (cond1){ ##case m==2
          LK <- JJ[i+1]
          RK <- JJ[i+1] ##both ==m
        }
        else if (cond2){ ## note m >2 here
          LK <- JJ[i+1]; ##Left Knot
          RK <- JJ[i+2];
        }
        else{ ##don't enter this if-block if p==1 (b/c then k==1)
          LK <- JJ[i+1]
          RK <- JJ[i+1]
        }
        ##ind <- (JJ[i] + 1):(JJ[i + 1] - 1)
        ind <- (JJ[i]+1):(RK-1)
        mtmp <- length(ind)

        
        ##ztmp <- (z[ind] - z[JJ[i]])
        ztmp <- rev(z[RK] - z[ind])
        dstmp <- c(z[JJ[i]+1]-z[JJ[i]],diff(z[ind]),
                   z[RK]-z[RK-1])##len==mtmp+1
        wtmp <- rev(w[ind])
        J01s <- rev(J10(phi[ind],phi[ind-1]) * dstmp[1:mtmp])
        J10s <- rev(J10(phi[ind],phi[ind+1]) * dstmp[2:(mtmp+1)])
        HR[ind] <- rev(cumsum(wtmp * ztmp) - ztmp * cumsum(wtmp) )
        jtmp <- rev( -ztmp*cumsum(J01s+J10s) + cumsum(ztmp * (J01s + J10s)) )

        if (cond1){##
          HR0 <- -rev(ztmp*w[n])
          JR0 <- -rev(ztmp*(J10(phi[n],phi[n-1])*(z[n]-z[n-1])))
          ##i <- i+1 ##no need, loop is already done.
        }
        else if (cond2){ ##may never happen, if constraint is on RHS of a
          ## m > 2 here.
          dtmp <- z[RK] - z[LK]
          ##tmptmp <- (LK + 1):(RK - 1)
          lentmp <- (RK-1) - (LK+1) +1
          ##translate by JJ[i]+1-1 b/c ws start at jj[i]+1
          ##ind.lin <- seq(from=LK+1-(JJ[i]+1-1), by=1, length=lentmp) ##"linear" indices
          ind.lin <- 1:lentmp  ## WE REVERSED EVERYTHING.
          HR0 <- rev(ztmp*sum(wtmp[ind.lin] * (1-ztmp[ind.lin]/dtmp)))
          JR0 <- rev(ztmp*sum((J01s[ind.lin]+J10s[ind.lin])*(1-ztmp[ind.lin]/dtmp)))
#### want to skip the next loop; this doesn't work though; 'for' will ignore it.
          ##i <- i+1
        }
        else{
          dtmp <- z[RK] - z[JJ[i]]
          HR0 <- rev( ztmp * sum(wtmp * (1-ztmp/dtmp)) )
          JR0 <- rev( sum((J01s+J10s)*(1-ztmp/dtmp))*ztmp )
        }
        HR[ind] <- HR[ind] + HR0  - jtmp - JR0
        HR[IsKnot] <- 0; ##REDUNDANT but JNC.
      }
    }
    if (k==n){
      H.m[1] <- H[k]
      H.m[2] <- 0
      ##H <- H  ## ie H = "HL"
    }
    else if (k==1){
      H.m[1] <- 0
      H.m[2] <- HR[k]
      H <- HR
    }
    else{
      H.m[1] <- H[k]
      H.m[2] <- HR[k]
      H[k:(n-1)] <- HR[k:(n-1)] 
    }
    ## ## HACK HACK HACK
    DL <- H.m.byhand <- NULL ## used only when mode is right next to another vector
### Del print statements
    ##print("p, p2")
    ##print(p) ## THOUGHT p could not be 1 but it CAN  ? 
    ##print(p2) 
    if (((p>1) && (JJ[p] == JJ[p-1]+1)) ||
        ((p2<n) && (JJ[p2+1] == JJ[p2]+1))){
      ## ## HACK HACK HACK
      ## BUG being fixed is: when the mode is a knot right next to a knot, the
      ## old code did not compute H.m correctly.  This is because it treats the
      ## mode as an inactive constraint (a knot), whereas really the modal
      ## constraint may be active. (i.e., if the mode is flat to the left, then
      ## the mode is a RK but is not a LK, but is treated as "a knot").  (This
      ## seems only to happen incorrectly when the mode is a data point?) So the
      ## following computes H.m by hand.  Overriding the above computations.
      ## with some error checking.

      ## ## names /vars slightly diff than above..
      ## just compute the formula by hand. note that because we don't make
      ## use of theoretical zeros, may result in more rounding errors
      dtmp <- diff(z)
      ind <- 1:n
      J10s <- dtmp[1:(n-1)] * J10(phi[1:(n-1)], phi[2:n])
      J01s <- dtmp[1:(n-1)] * J10(phi[2:n],phi[1:(n-1)]) ##J01(r,s) = J10(s,r)
      DL <- w ## derivative of L w.r.t. each coordinate
      DL[1:(n-1)] <- DL[1:(n-1)] - J10s
      DL[2:n] <- DL[2:n] - J01s

      delta.L <- pmin(z-z[k], 0)  ## left mode perturb
      delta.R <- pmin(z[k]-z,0)   ## right mode perturb
      ## automatically (approx) zero if k is 1 or n 
      H.m.byhand[1] <- sum(delta.L * DL)
      H.m.byhand[2] <- sum(delta.R * DL)
      ## but precision errors can cause crashes; will set phi[k] to -Inf,
      ## etc., if k==1||k==n.
      if (k==1 || k==n) H.m.byhand[H.m.byhand>0] <- 0
      if (print){
        print("Computing H.m.byhand. H.m and H.m.byhand are" )
        print(H.m)
        print(H.m.byhand)
      }
      ##now that we've printed H.m can override it with the byhand version
      H.m <- H.m.byhand * (1-IsMIC)
      ## do we want to enforce multiply by 1-IsMIC? We do implicitly for IsKnot
    }
    ## error check
    ## ## end HACK HACK HACK
    ##H[k:n] <- HR[k:n]
    H[k] <- H[n] <- H[1] <- 0
  }
  res <- list(phi = matrix(phi, ncol = 1), L = res2$L,
              shape=shape, H = matrix(H, ncol = 1),H.m=matrix(H.m,ncol=1))
  return(res)
}






#####z includes a.  a is the usual, idx,val,isx.
##### NOTE: if is.null(a) returns matrix.  otherwise returns list.!
##### Constraint vectors are columns
##### V refers to standard constraints; W to the 2 "modal" or "monotone" constraints.
##### Note W is always nx2, even when a is on the border.
##### keep in mind: < constraintvec_i , basisvec_i > <= 0 is the constraint.
getConstraintVecs <- function(z, a=NULL){
  dz <- diff(z);
  n <- length(z)
  invdz <- 1/dz;
  V <- matrix(0,nrow=n,ncol=n)
  for (i in 1:n){
    for (j in 2:(n-1)){
      if (i==j-1) V[i,j] <- invdz[i]
      else if (i==j) V[i,j] <- -(invdz[j-1]+invdz[j])
      else if (i==j+1) V[i,j] <- invdz[j]
    }
  }
  if (is.null(a)){
    ##return(V);
    return(list(V1=V, V2=V,W=NULL))
  }
  else{
    k <- a$idx;
    W <- matrix(0,nrow=n,ncol=2);
    W[k,1] <- W[k,2] <- -1; ## no reason to use +/-delta_k rather than +/- 1 i ithink
    W[k-1,1] <- W[k+1,2] <- 1;
    if (k==1) W[,1] <- rep(0,n)
    else if (k==n) W[,2] <- rep(0,n) ##keep in mind,if !a$isx, then n=n1+1
    return(list(V1=V[,1:(k-1)],V2=V[,(k+1):n],W=W))
  }
}

###again if a is not null, return a list; if it is null return just a matrix
### Columns as basis vectors.
### B is elbows nonzero (ie positive) to the left, A is elbows nonzero to the right
### right now it returns all of A and B so that the indexing is easy
### so esssentially ignores a.  could return A[k:n] and B[1:k] (overlap at k).
getBasisVecs <- function(z){
  n <- length(z)
  B <- -matrix(z,nrow=n,ncol=n,byrow=FALSE) +
    matrix(z,nrow=n,ncol=n,byrow=TRUE) ## [,i] column all equal z[i]
  A <- -B
  B[,1] <- rep(1,n)
  B <- apply(B, c(1,2), function(x){max(x,0)})
  A <- apply(A, c(1,2), function(x){max(x,0)})
  A[,n] <- rep(1,n)
  return(list(B=B,A=A))
  ##return(list(B=B[1:k],A=A[k:n]))
  ##  }
}




### y are knots
### w, phi0 same length as y
### constr is 2 (consecutive) idcs in y
MLE.mode <- function (y,constr, w = NA, phi_o = NA, prec = 10^(-7), print = FALSE) 
{
  n <- length(y)
  if (sum(y[2:n] <= y[1:n - 1]) > 0) {
    cat("We need strictly increasing numbers y(i)!\n")
    stop("Exiting because of bad ys")
  }
  if (max(is.na(w)) == 1) {
    w <- rep(1/n, n)
  }
  if (sum(w < 0) > 0) {
    cat("We need nonnegative weights w(i) !\n")
    stop("exiting because of bad ws")
  }
  ww <- w/sum(w)
  if (max(is.na(phi_o)) == 1) {
    m <- sum(ww * y)
    s2 <- sum(ww * (y - m)^2)
    phi <- LocalNormalize(y, -(y - m)^2/(2 * s2))
  }
  else {
    phi <- LocalNormalize(y, phi_o)
  }
  iter0 <- 0
  r1 <- LocalLLall.mode(y, ww, phi,constr)##comes up w new candidate (& corresponding dirderiv); ##returns 'knots' format
  L <- r1$ll
  phi_new <- r1$phi.new
  dirderiv <- r1$dirderiv
  ##while ((dirderiv >= prec) & (iter0 < 100)) {
  while ((dirderiv >= prec) && (iter0 < 100)) {## the "&&" is untested change from "&"
    iter0 <- iter0 + 1
    L_new <- Local_LL(y, ww, phi_new)##old version is fine here
    iter1 <- 0
    ##while ((L_new < L) & (iter1 < 20)) {
    if (print == TRUE) {
      print(paste("mle.mode:outer1: iter0=", iter0, " / iter1=", iter1, 
                  " / L_new=", round(L_new, 4),
                  " / L=", round(L,4),
                  " / L_new < L=", L_new<L,
                  sep = ""))
      ##           if (any(is.na(c(L_new,L))) || is.null(L_new) || is.null(L)){
      ##           }
    }
    while ((L_new < L) && (iter1 < 20)) { ## the "&&" is untested change from "&"
      ## NR moves in direction deriv is positive.
      ##Then we find the t such that the likelihood increases.
      ## 2^-20  approx 1e-
      iter1 <- iter1 + 1
      phi_new <- 0.5 * (phi + phi_new)
      L_new <- Local_LL(y, ww, phi_new)##old version is fine here
      dirderiv <- 0.5 * dirderiv
      if (print == TRUE) {
        print(paste("mle.mode:inner1: iter0=", iter0, " / iter1=", iter1, 
                    " / L_new=", round(L_new, 4), " / dirderiv=", round(dirderiv, 
                                                                        4), sep = ""))
      }
    }
##### don't understand this part.
#####  is this trying to provide a lower bound on the amount
##### of increase that we get in the LL (in terms of the dirderiv)?
#####  or ensuring that ActiveSet is increasing???
    if (L_new >= L) {
      tstar <- max((L_new - L)/dirderiv) ##max??
      if (tstar >= 0.5) {
        phi <- LocalNormalize(y, phi_new)
      }
      else {
        tstar <- max(0.5/(1 - tstar)) ##max??
        phi <- LocalNormalize(y, (1 - tstar) * phi + 
                              tstar * phi_new)
      }
      r1 <- LocalLLall.mode(y, ww, phi,constr)
      L <- r1$ll
      phi_new <- r1$phi.new
      dirderiv <- r1$dirderiv
    }
    else {
      dirderiv <- 0
    }
    if (print == TRUE) { ##print==true and print is a function...?
      print(paste("mle.mode:outer2: iter0=", iter0, " / iter1=",
                  iter1, 
                  ##" / L=", round(L, 4),
                  " / L=", round(L, 4),
                  " / dirderiv=", round(dirderiv, 4), sep = ""))
    }
  }
  r1 <- list(phi = matrix(phi, ncol = 1), L = L,
             Fhat = matrix(LocalF(y,  phi), ncol = 1))
  return(r1)
}







                                        #requires a not be NA
#######AGGH i think 'm' is called 'p' now
getm <- function(x,phi,a){
  print("getm: careful using this function: m is called p now i believe.")
  k <- a$idx;
  phim <- max(phi[k-1],phi[k])
  m <- k-2 + which(phim==phi[(k-1):k]);
}

###val can be a vector
### idx is single integer.
### values are inserted before the given index.
### so that in the result, at idx, you find val.
insert <- function(val,vec,idx){
  n <- length(vec)
  if (idx<=1)
    return(c(val,vec))
  else if (idx>=n+1)
    return(c(vec,val)) 
  else
    return(c(vec[1:(idx-1)],val,vec[idx:n]));
}
delete <- function(vec,idcs){ #works with NULL=idcs
  keep <- rep(TRUE,length(vec))
  keep[idcs] <- FALSE
  vec[keep];
}

###"fk2CK" = fine knots to coarse knots;
### 'fine' means includes extar point, 'coarse' means doesn't.
FK2CK <- function(y,w,phi,constr){
  ll <- length(y); ##"ell", "L" but too hard to distinguish "ell" from "one"
  o <- length(constr)
  if (o<=1 || constr[1]<1 || constr[o] > ll || constr[1]==constr[2]){ ##no mono constraint.
    res <- list(y=y,w=w,phi=phi,constr=constr)
  }
  else if (o == 2 && constr[1]==constr[2]-1){
    m <- min(constr); ## constr should be m and m+1
    wghtsum <- sum(w[m:(m+1)])
    w <- delete(w,m);
    w[m] <- wghtsum;
    v <- delete(y,m);
    phi <- delete(phi,m);
    res <- list(y=v, w=w,phi=phi,constr=constr)
  }
  else {
    print(paste("FK2CK ... : bad constraint vec passed. ",
                "len should be 2.  constr is"));
    print(constr);
  }
  return(res);
}

##returns 'p' which MAY NOT index z[k] in the knots!  ('a' may not be a knot!)
LocalCoarsen.mode <- function (z, w, IsKnot,IsMIC,a) 
{
  n <- length(z)
  ## Need to decide if z[k] will be a knot:
  idcs <- seq(from=1,by=1,length=a$idx-1) ##ie if k==1, get integer(0)
  KL <- idcs * IsKnot[idcs]
  KL <- KL[KL>0]
  p <- length(KL)+1 ## may or may not eventually correspond to a
  constr <- rep(p,length=2) ##default; including if k==1 or ==n; or IsMIC=c(1,1)
  aIsKnot <- NULL;
  if (identical(IsMIC,c(0,0))){
    aIsKnot <- FALSE
    if (a$idx==1) constr <- c(p,p+1)
    else constr <- c(p-1,p)
  }
  else{ ## all other scenarios, a is a knot
    ## get 2 consecutive idcs of y (ie _idcs_ of K) that are constrained
    aIsKnot <- TRUE
    if (a$idx != 1) constr[1] <- c(p-1,p)[IsMIC[1]+1] ##1, n are always knots
    if (a$idx != n) constr[2] <- c(p+1,p)[IsMIC[2]+1]
  }
##### this should be redundant here
######IsKnot <- syncIsKnot(IsKnot,IsMIC,a$idx)
  ##after we've decided on z[k]
  KR <- ((a$idx):n) * IsKnot[(a$idx):n]
  KR <- KR[KR>0]
  K <- c(KL,KR)
  x2 <- z[K]
  w2 <- w[K]
  for (k in 1:(length(K) - 1)){##at least 2 knots.
    if (K[k + 1] > (K[k] + 1)){
      ind <- (K[k] + 1):(K[k + 1] - 1)
      lambda <- (z[ind] - x2[k])/(x2[k + 1] - x2[k])
      w2[k] <- w2[k] + sum(w[ind] * (1 - lambda))
      w2[k + 1] <- w2[k + 1] + sum(w[ind] * lambda)
    }
  }
  w2 <- w2/sum(w2)
  return(list(y = matrix(x2, ncol = 1), w2 = matrix(w2, ncol = 1),
              constr=constr,  p=p, aIsKnot=aIsKnot)) ## p MAY NOT INDEX
}


####This should be run every time IsMIC is modified.
syncIsKnot <- function(IsKnot,IsMIC,k){
  if (max(is.na(IsMIC))==1)
    print("syncIsKnot:IsMIC should have negative or NULL if k==1 or n, not NA")
  if (k==1) IsKnot[k] <- 1 ##always
  else if (k==length(IsKnot)) IsKnot[k] <- 1
  else if (max(IsMIC)==1) IsKnot[k] <- 1 ##IsMIC may have NULL or negatives if k==1,n.
  else IsKnot[k] <- 0
  return(IsKnot)
}




###### Reinforces constancy of phi across constraint; REDUNDANT.
getGrad <- function(y,w,phi, constr=0){
### gradient NEEDs the full knot vector and the constraint.
### "y" is all knots.  this includes a and all modes regardless of constraints
### i.e. any data point that has a bend in the actual phi.
### Then "constr" is the binding constraint, integral w/ length 2
###  the indices that are constrained in y. One of these always pertains
### to the index of a in the knots.  at least one more pertains to which
### side is constrained.  3 if both are.
### at the end, basically: length(grad) == length(y) - (length(constr)-1)
### taking len(constr)=1 if there is no constraint.
  ll <- length(y); ##"ell", "L" but too hard to distinguish "ell" from "one"
  o <- length(constr)
  if (o<=1 || constr[1]<1 || constr[o] > ll || constr[1]==constr[2]){ ##no mono constraint.
    grad <- matrix(0,ncol=1,nrow=ll);
    grad <- t(getGrad.uncstr(x=y,w=w,phi=phi)); ##length == ll
  }
  else if (o==2 && constr[1]==constr[2]-1){
### the 'om1' is holdover code from when o could be 3 or 2
    om1 <- o-1
    grad <- matrix(0,ncol=1,nrow=ll-om1);
    m <- min(constr); ##ie constr[1]
    if (m==1){ ##could do this w/ similar matrix as one below...
      endseq <- seq(from=m+o, by=1, length=ll-m-o+1)
      phi <- c(rep(phi[m],o), phi[endseq]); ##redundant perhaps
      ##phi <- c(rep(phi[m],o), phi[(m+o):ll]); ##redundant perhaps
    }
    else if (m==ll-1)
      {phi <- c(phi[1:(m-1)], rep(phi[m],o));} ##redundant perhaps
    else     ## m != ll ever when o==2.
      {phi <- c(phi[1:(m-1)], rep(phi[m],o), phi[(m+o):ll]);} ##redundant perhaps
    xtraRow <- rep(0,ll-1); xtraRow[m] <- 1;
    projmat <- diag(rep(1,ll-1)) ##this works even if m==1
    if (m<ll-1) {projmat <- rbind(projmat[1:m,],xtraRow,projmat[(m+1):(ll-1),])}
    else {projmat <- rbind(projmat[1:m,],xtraRow)}
    grad.tmp <- getGrad.uncstr(y,w,phi);
    ##     grad[m] <- sum(grad.tmp[m:(m+om1)]);
    ##     grad[1:(m-1)] <- grad.tmp[1:(m-1)];
    ##     grad[(m+1):(ll-om1)] <- grad.tmp[(m+o):ll]
    grad <- t(grad.tmp) %*% projmat
  }
  else {
    print(paste("getGrad ... : bad constraint vec passed. ",
                "len should be 2.  constr is"));
    print(constr);
  }
  grad;
}


reducePhi <- function(phi,constr){
  ll <- length(phi)
  o <- length(constr)
  if (o<=1 || constr[1]<1 || constr[o] > ll || constr[1]==constr[2]){ ##no mono constraint.
    return(phi)
  }
  else if (o==2 && constr[1]==constr[2]-1){
    m <- min(constr); ##ie constr[1]
    ##     ##om1 is holdover from when o could be 3 or 2.
    ##     om1 <- o-1
    ##     m <- min(constr); ##ie constr[1]
    ##     if (m==1){ ##could do this w/ similar matrix as one below...
    ##       endseq <- seq(from=m+o, by=1, length=ll-m-o+1)
    ##       phi <- c(phi[m], phi[endseq])
    ##     }
    ##     else if (m==ll-1)
    ##       phi <- c(phi[1:(m-1)], phi[m])
    ##     else     ## m != ll ever when o==2.
    ##       phi <- c(phi[1:(m-1)], phi[m], phi[(m+o):ll])
    phi <- delete(phi,m+1) ##m+1=constr[2], should be allowable index for phi.
  }
  else {
    print(paste("reducePhi ... : bad constraint vec passed. ",
                "len should be 2.  constr is"));
    print(constr);
  }
  return(phi)
}

###note: constr refers to phi NOT reduced, not to phi.red(uced)!!!
unreducePhi <- function(phi.red,constr){
  ll <- length(phi.red)
  o <- length(constr)
  if (o<=1 || constr[1]<1 || constr[1]==constr[2] || constr[o] > ll+1){ ##no mono constraint.
    return(phi.red)
  }
  else if (o==2 && constr[1]==constr[2]-1){
    m <- min(constr); ##ie constr[1]
    if (m==1){ 
      endseq <- seq(from=m+o-1, by=1, length=ll+1-m-o+1)
      phi.red <- c(rep(phi.red[m],o), phi.red[endseq]); 
    }
    else if (m==ll)
      {phi.red <- c(phi.red[1:(m-1)], rep(phi.red[m],o));}
    else     ## m != ll-1 ever when o==2.
      {phi.red <- c(phi.red[1:(m-1)], rep(phi.red[m],o), phi.red[(m+1):ll]);}
  }
  else {
    print(paste("unreducePhi ... : bad constraint vec passed. ",
                "len should be 2.  constr is"));
    print(constr);
  }
  return(phi.red)
}

##constr as usual is either 0 ie no constraint or length 2.
## could specify it to be length 1 but if its 2 it's clear exactly what
## most of this function is redundant
getHess <- function(y,w,phi,constr=0){
  ll <- length(y); ##"ell", "L" , not "one"!
  o <- length(constr)
  if (o<1 || constr[1]<1 || constr[o] > ll || constr[1]==constr[2]){ ##no mono constraint.
    hess <- matrix(0,ncol=ll,nrow=ll);
    hess <- getHess.uncstr(x=y,w=w,phi=phi); ##length == ll
  }
  else if (o==2 && constr[1]==constr[2]-1){
### NOTE: because 'a' is passed always as part of y,
### o might be 3.  could change this so that a isn't
### always passed in.  then, o would always be 2.
### the 'om1' is holdover code from when o could be 3 or 2
    om1 <- o-1
    hess <- matrix(0,ncol=ll-om1,nrow=ll-om1);
    m <- min(constr);
#######
####### REDUNDANT: reinforcing constancy of phi on constraint.
    if (m==1) {##could do this w/ similar matrix as one below...
      endseq <- seq(from=m+o, by=1, length=ll-m-o+1)
      phi <- c(rep(phi[m],o), phi[endseq]);
    }
    else if (m==ll-1)
      {phi <- c(phi[1:(m-1)], rep(phi[m],o));}
    else
      {phi <- c(phi[1:(m-1)], rep(phi[m],o), phi[(m+o):ll]);} 
    ##else if (m==ll){ ##cant happne
    ##}
######
    xtraRow <- rep(0,ll-1); xtraRow[m] <- 1;
    projmat <- diag(rep(1,ll-1)) ##this works even if m==1
    if (m<ll-1) {projmat <- rbind(projmat[1:m,],xtraRow,projmat[(m+1):(ll-1),])}
    else {projmat <- rbind(projmat[1:m,],xtraRow)}
    hess.tmp <- getHess.uncstr(y,w,phi);
    hess <- t(projmat) %*% hess.tmp %*% projmat;         
  }
  else {
    print(paste("getHess ... : bad constraint vec passed. ",
                "len should be 2.  constr is"));
    print(constr);
  }
  hess;
}

##This function is generally useless but useful for
## testing whether gradient /  hessian are correct
## constr is 2 (consecutive) constrained indices in y
Local_LL.mode <- function(y,w,phi,constr=0){
  ll <- length(y); ##"ell", "L" , not "one"!
  o <- length(constr)
  if (o<1 || constr[1]<1 || constr[o] > ll || constr[1]==constr[2]){ ##no mono constraint.
    return(Local_LL(y,w,phi));
  }
  else if (o == 2 && constr[1]==constr[2]-1){
    ##r <- FK2CK(y,w,phi,constr);
    ##return(Local_LL(r$y,r$w,r$phi));
    m <- min(constr);
    phi[m+1] <- phi[m]
    return(Local_LL(y,w,phi));
  }
  else {
    print(paste("Local_LL.mode ... : bad constraint vec passed. ",
                "len should be 2.  constr is"));
    print(constr);
  }
  -1;
}




##uncstr = unconstrained
getGrad.uncstr <- function(x,w,phi){
  grad <- matrix(w, ncol = 1)
  n <- length(phi);
  dx <- diff(x)
  J10s <- J10(phi[1:(n - 1)], phi[2:n])
  J01s <- J10(phi[2:n], phi[1:(n - 1)])
  grad[1:(n - 1)] <- grad[1:(n - 1)] - (dx * J10s)
  grad[2:n] <- grad[2:n] - (dx * J01s);
  grad
}

### y is vector of knots.
### constr is first index of mono constrained.
### (there are always 2 constrained togethre or 0)
LocalLLall.mode <- function(y,w,phi,constr=NULL){
  ll <- Local_LL.mode(y=y,w=w,phi=phi,constr=constr);
  LIs <- (1:length(phi))[exp(phi)==Inf] ##Large Indices
  if (length(LIs) > 0 ) stop("LocalLLall.mode Error: We have not yet accounted for extraordinarily steep/large phi.  This is probably an error.")
  grad <- getGrad(y=y,w=w,phi=phi,constr=constr)
  hess <- getHess(y=y,w=w,phi=phi,constr=constr)
  ##note, i changed hessian to be the actual hessian, not its negative.
  phi.red <- reducePhi(phi,constr)
  phi.red.new <- phi.red - solve(hess) %*% t(grad) ## one newton-raphson step
  dirderiv <- grad %*% (phi.red.new - phi.red)
  phi.new <- unreducePhi(phi.red.new,constr)
  return(list(ll=ll, phi.new=phi.new, dirderiv=dirderiv));
}


getHess.uncstr <- function(x,w,phi){
  dx <- diff(x);
  n <- length(x);
  tmp <- c(dx * J20(phi[1:(n - 1)], phi[2:n]), 0) +
    c(0, dx * J20(phi[2:n], phi[1:(n - 1)]))
  tmp <- tmp + mean(tmp) * 10^(-12)
  mhess2 <- matrix(0, nrow = n, ncol = n)
  mhess3 <- mhess2
  mhess1 <- tmp
  tmp <- c(0, dx * J11(phi[1:(n - 1)], phi[2:n]))
  tmp.up <- diag(tmp[2:n], nrow = n - 1, ncol = n - 1)
  mhess2[1:(n - 1), 2:n] <- tmp.up
  mhess3[2:n, 1:(n - 1)] <- diag(tmp[2:n], nrow=n-1, ncol=n-1)
  mhess <- diag(mhess1) + mhess2 + mhess3;
  ##mhess;
  -mhess; ## prefer actual hessian not its negative.
}




##k=NULL just get standard lambdaaaaa
## NOTE: k indexes z, and conv is the same indexing as z.
getLambda <- function(shape.new,shape,IsKnot,IsMIC, k){
  conv.new <- shape.new$conv; conv <- shape$conv
  mono.new <- shape.new$mono; mono <- shape$mono
  n <- length(conv);
  JJ1 <- (1:n) * (conv.new > 0)
  JJ1 <- JJ1[JJ1 > 0]
  JJ2 <- (1:2) * (mono.new>0)
  JJ2 <- JJ2[JJ2>0]
  tmp1 <- conv[JJ1]/(conv[JJ1] - conv.new[JJ1])
  tmp2 <- mono[JJ2]/(mono[JJ2] - mono.new[JJ2])
  lambda <- min(c(tmp1,tmp2))
  if (!is.null(IsMIC) && !is.null(k)){
    KK1 <- (1:length(JJ1)) * (tmp1 == lambda)
    KK1 <- KK1[KK1 > 0]
    IsKnot[JJ1[KK1]] <- 0
    KK2 <- (1:length(JJ2)) * (tmp2 ==lambda)
    KK2 <- KK2[KK2>0]
    IsMIC[JJ2[KK2]] <- 0
    if (k==1) IsMIC[1] <- -1 ##redundant i think.
    else if (k==n) IsMIC[2] <- -1
    IsKnot <- syncIsKnot(IsKnot,IsMIC,k)
  }
  else{ ##what's the point of dealing with this case?
    ## also, not sure if i deal
    ## with JJ1 correctly below, changed without super careful inspection.
    KK <- (1:length(JJ1)) * (tmp1 == lambda)
    KK <- KK[KK > 0]
    ## IsKnot[JJ[KK]] <- 0 ## HERE HERE ERROR: WHAT IS JJ
    ## changed "JJ" to "JJ1" here. is this correct?
    IsKnot[JJ1[KK]] <- 0 ## HERE HERE: IS JJ1 correct
  }
  return(list(lambda=lambda, IsKnot=IsKnot, IsMIC=IsMIC));
}


##k is sort of unnecessary. good programming practice to call syncIsKnot tho.
inactivate <- function(H,H.m, IsKnot, IsMIC, k){
  tmp <- max(c(H,H.m))
  j1 <- (1:length(H)) * (H == tmp)
  j2 <- (1:2) * (H.m==tmp)
  j1 <- j1[j1>0]
  j2 <- j2[j2>0]
  ##j1 <- min(j1[j1>0])
  ##j2 <- min(j2[j2>0])
  if ( length(j2) > 0 )
    IsMIC[j2[1]] <- 1
  else
    IsKnot[j1[1]] <- 1
  IsKnot <- syncIsKnot(IsKnot,IsMIC,k=k) ##the way it is now, IsKnot won't change.
  return(list(IsKnot=IsKnot, IsMIC=IsMIC))
}

## activeSetLogCon.mode <- 
## function (x, aval = x[1], w = NA, print = FALSE, logfile = NULL, 
##     prec = 10^-10) 
## {
##     if (!is.null(logfile) && !is.na(logfile)) 
##         sink(logfile)
##     n1 <- length(x)
##     if (sum(x[2:n1] <= x[1:(n1 - 1)]) > 0) {
##         cat("We need strictly increasing numbers x(i)!\n")
##         print("x is ")
##         print(x)
##         save(x, file = "ASLCMbadxs.rsav")
##         stop("Exiting because of bad xs")
##     }
##     if (max(is.na(w)) == 1) {
##         w <- rep(1/n1, n1)
##     }
##     if (sum(w < 0) > 0) {
##         cat("We need nonnegative weights w(i)! \n")
##         save(w, file = "ASLCMbadws.rsav")
##         stop("Exiting because of bad ws")
##     }
##     w <- w/sum(w)
##     phi <- LocalNormalize(x, 1:n1 * 0)
##     r1 <- x2z(x = x, w = w, phi = phi, aval = aval)
##     z <- r1$z
##     w <- r1$w.a
##     phi <- r1$phi
##     a <- r1$a
##     n <- length(z)
##     IsKnot <- 1:n * 0
##     IsMIC <- c(0, 0)
##     IsKnot[c(1, n)] <- 1
##     IsKnot <- syncIsKnot(IsKnot, IsMIC, a$idx)
##     res1 <- LocalMLE.mode(z, w, IsKnot, IsMIC, a = a, phi_o = phi, 
##         prec, print = print)
##     phi <- res1$phi
##     L <- res1$L
##     shape <- res1$shape
##     H <- res1$H
##     H.m <- res1$H.m
##     iter1 <- 1
##     while ((iter1 < 500) && ((max(H) > prec * mean(abs(H))) || 
##         (max(H.m) > prec * mean(abs(H.m))))) {
##         IsKnot_old <- IsKnot
##         IsMIC_old <- IsMIC
##         iter1 <- iter1 + 1
##         IC <- inactivate(H, H.m, IsKnot, IsMIC, k = a$idx)
##         IsKnot <- IC$IsKnot
##         IsMIC <- IC$IsMIC
##         res2 <- LocalMLE.mode(z, w, IsKnot, IsMIC, a = a, phi_o = phi, 
##             prec, print = print)
##         phi_new <- res2$phi
##         L <- res2$L
##         shape_new <- res2$shape
##         H <- res2$H
##         H.m <- res2$H.m
##         if (print == TRUE) {
##             Kidx <- (1:n)[as.logical(IsKnot - IsKnot_old)]
##             print(paste("ASLCM:Proc2: ", "iter1=", iter1, " / L=", 
##                 round(L, 4), " / max(H)=", round(max(H), 4), 
##                 " / max(Hm)=", round(max(H.m), 4), " / IsKnot idx=", 
##                 Kidx, sep = ""))
##             save(H, H.m, file = "activesetlogconmode.print.rsav")
##         }
##         while (max(shape_new$conv) > prec * max(abs(shape_new$conv)) || 
##             max(shape_new$mono) > prec * max(abs(shape_new$mono))) {
##             IsKnot_old2 <- IsKnot
##             IsMIC_old2 <- IsMIC
##             s1 <- getLambda(shape_new, shape, IsKnot, IsMIC, 
##                 k = a$idx)
##             lambda <- s1$lambda
##             IsKnot <- s1$IsKnot
##             IsMIC <- s1$IsMIC
##             phi <- (1 - lambda) * phi + lambda * phi_new
##             shape <- LocalConvexity.mode(z = z, phi, k = a$idx)
##             shape$conv <- pmax(shape$conv, 0)
##             shape$mono <- pmax(shape$mono, 0)
##             res3 <- LocalMLE.mode(z, w, IsKnot, IsMIC, a = a, 
##                 phi, prec, print = print)
##             phi_new <- res3$phi
##             L <- res3$L
##             shape_new <- res3$shape
##             H <- res3$H
##             H.m <- res3$H.m
##             H <- H * (1 - IsKnot)
##             H.m <- H.m * (1 - IsMIC)
##             if (print == TRUE) {
##                 Kidx <- (1:n)[as.logical(IsKnot - IsKnot_old2)]
##                 print(paste("ASLCM:Proc1: ", "iter1=", iter1, 
##                   " / L=", round(L, 4), " / max(H)=", round(max(H), 
##                     4), " / max(Hm)=", round(max(H.m), 4), " / max(convnew)=", 
##                   round(max(shape_new$conv), 4), " / max(mononew)=", 
##                   round(max(shape_new$mono), 4), " / IsKnot idx=", 
##                   Kidx, sep = ""))
##                 print(paste("IsMic"))
##                 print(IsMIC)
##             }
##         }
##         phi <- phi_new
##         shape <- shape_new
##         if (sum(IsKnot != IsKnot_old) == 0 && sum(IsMIC != IsMIC_old) == 
##             0) {
##             if (print == TRUE) {
##                 print("No change in constraints. Ending.")
##             }
##             break
##         }
##     }
##     Fhat <- LocalF(z, phi)
##     phi.f <- function(x0) {
##         apply(matrix(x0), 1, function(x00) {
##             evaluateLogConDens(x00, x = z, phi, Fhat, IsKnot)[1]
##         })
##     }
##     fhat.f <- function(x0) {
##         apply(matrix(x0), 1, function(x00) {
##             evaluateLogConDens(x00, x = z, phi, Fhat, IsKnot)[2]
##         })
##     }
##     Fhat.f <- function(x0) {
##         apply(matrix(x0), 1, function(x00) {
##             evaluateLogConDens(x00, x = z, phi, Fhat, IsKnot)[3]
##         })
##     }
##     E.f <- intFfn(z, phi, Fhat)
##     tmp <- LocalCoarsen.mode(z, w, IsKnot, IsMIC, a)
##     MI <- z[IsKnot > 0][tmp$constr]
##     {
##         KK <- (1:n)[as.logical(IsKnot)]
##         phiK <- phi[as.logical(IsKnot)]
##         slopes <- diff(phiK)/diff(z[KK])
##         phiPL.f <- stepfun(x = z[KK], c(slopes[1], slopes, tail(slopes, 
##             1)), right = TRUE, f = 1)
##         phiPR.f <- stepfun(x = z[KK], c(slopes[1], slopes, tail(slopes, 
##             1)), right = FALSE, f = 0)
##         phiPL <- phiPL.f(z)
##         phiPR <- phiPR.f(z)
##     }
##     if (!is.null(logfile) && !is.na(logfile)) 
##         sink(NULL)
##     return(list(x = x, z = z, w = w, a = a, MI = MI, phi = phi, 
##         IsKnot = IsKnot, n1 = n1, n = n, knots = KK, IsMIC = IsMIC, 
##         constr = tmp$constr, L = L, fhat = exp(phi), Fhat = Fhat, 
##         H = H, H.m = H.m, phi.f = phi.f, fhat.f = fhat.f, Fhat.f = Fhat.f, 
##         E.f = E.f, phiPL = phiPL, phiPR = phiPR, phiPL.f = phiPL.f, 
##         phiPR.f = phiPR.f))
## }


### CAREFUL: for large length(x), prec really matters.
##  prec=1e-10 seems good; 1e-12  is too small!
activeSetLogCon.mode <- function (x,xgrid=NULL,
                                  ##aval=x[1],
                                  mode=x[1],
                                  print = FALSE,
                                  w = NA,
                                  ## logfile=NULL, ## just use sink() outside of func call
                                  prec=10^-10) {
  ##hopefully code works for boundary case, k=1; k=1 iff a=x[1].
  ##if (!is.null(logfile) && !is.na(logfile)) sink(logfile)
  aval <- mode; ## last minute renaming 
  if (print){ print("ASLCM: Beginning")}
  xn <- sort(x)
  if ((!identical(xgrid, NULL) & (!identical(w, NA)))) {
    stop("If w != NA then xgrid must be NULL!\n")
  }
  if (identical(w,NA)){
    tmp <- preProcess(x,xgrid=xgrid)
    x <- tmp$x
    w <- tmp$w
    sig <- tmp$sig ##nonpara est of sd.
  }
  else { ## if (!identical(w, NA)) {
    if (abs(sum(w) - 1) > prec) stop("activeSetLogCon.mode Error: weights w do not sum to 1.")
    tmp <- cbind(x, w)
    tmp <- tmp[order(x), ]
    x <- tmp[, 1]
    w <- tmp[, 2]
    est.m <- sum(w * x)
    est.sd <- sum(w * (x - est.m)^2)
    est.sd <- sqrt(est.sd * length(x)/(length(x) - 1))
    sig <- est.sd
  }
  n1 <- length(x)
  phi <- LocalNormalize(x, 1:n1 * 0)
  r1 <- x2z(x=x,w=w,phi=phi,aval=aval);
  z <- r1$z;
  w <- r1$w.a;
  phi <- r1$phi;
  a <- r1$a
  class(a) <- "dlc.mode"
  n <- length(z); 
  IsKnot <- 1:n * 0
  IsMIC <- c(0,0) ##not sure how to start ... flat for now.
  ##IsMIC <- c(1,1) ##not sure how to start ... flat for now.  
  IsKnot[c(1, n)] <- 1
  IsKnot <- syncIsKnot(IsKnot,IsMIC,a$idx); ##good practice.
  ##res1 <- LocalMLE.mode(z, w, IsKnot, IsMIC, a=a, phi_o=phi, prec,AB=AB)
  res1 <- LocalMLE.mode(z, w, IsKnot, IsMIC, a=a, phi_o=phi, prec, print=print)
  phi <- res1$phi
  L <- res1$L
  shape <- res1$shape ##includes 'modal' type constraints too
  H <- res1$H
  H.m <- res1$H.m
  iter1 <- 1
  ## don't understand why start w/ proc2.   res1 is (probably) not in K.
  while ((iter1 < 500) &&
         ((max(H) > prec * mean(abs(H))) || (max(H.m) > prec * mean(abs(H.m))))) { ##procedure 2
    IsKnot_old <- IsKnot
    IsMIC_old <- IsMIC
    iter1 <- iter1 + 1
    IC <- inactivate(H,H.m,IsKnot,IsMIC, k=a$idx)
    IsKnot <- IC$IsKnot;
    IsMIC <- IC$IsMIC
    res2 <- LocalMLE.mode(z, w, IsKnot,IsMIC, a=a, phi_o=phi, prec,print=print) ##not convex/modal
    phi_new <- res2$phi
    L <- res2$L
    shape_new <- res2$shape
    H <- res2$H
    H.m <- res2$H.m
    if (print == TRUE) {
      Kidx <- (1:n)[as.logical(IsKnot-IsKnot_old)] ##want to check if going back 'n forth
      ##Midx <- (1:n)[as.logical(IsMIC-IsMIC_old)] ##want to check if going back 'n forth
      ##Midx is broken!      
      print(paste("ASLCM:Proc2: ", "iter1=", iter1, " / L=", round(L, 4),
                  ## " / max(H)=", round(max(H), 4),
                  ## " / max(Hm)=", round(max(H.m), 4),
                  " / max(H)=", max(H),
                  " / max(Hm)=", max(H.m),
                  " / IsKnot idx=", Kidx,
                  ##" / IsMIC idx=", Midx,
                  sep = ""))
      save(H,H.m, file="activesetlogconmode.print.rsav"); ## too long  to print
    }
    while (max(shape_new$conv) > prec * max(abs(shape_new$conv)) ||
           max(shape_new$mono) > prec * max(abs(shape_new$mono))   ){
      IsKnot_old2 <- IsKnot
      IsMIC_old2 <- IsMIC
      s1 <- getLambda(shape_new, shape, IsKnot, IsMIC, k=a$idx)
      lambda <- s1$lambda;
      IsKnot <- s1$IsKnot;
      IsMIC <- s1$IsMIC;
      phi <- (1 - lambda) * phi + lambda * phi_new
      shape <- LocalConvexity.mode(z=z, phi, k=a$idx)
      shape$conv <- pmax(shape$conv,0) ##positive gets 0s
      shape$mono <- pmax(shape$mono,0)
      ##res3 <- LocalMLE.mode(z, w, IsKnot,IsMIC, a=a,phi, prec,AB=AB)  ##not convex/modal
      res3 <- LocalMLE.mode(z, w, IsKnot,IsMIC, a=a,phi, prec,print=print)  ##not convex/modal
      phi_new <- res3$phi
      L <- res3$L
      shape_new <- res3$shape
      H <- res3$H
      H.m <- res3$H.m
      H <- H * (1-IsKnot)
      H.m <- H.m * (1-IsMIC)
      if (print == TRUE) {
        Kidx <- (1:n)[as.logical(IsKnot-IsKnot_old2)] ##want to check if going back 'n forth
        ##Midx <- (1:n)[as.logical(IsMIC-IsMIC_old2)] ##want to check if going back 'n forth
        print(paste("ASLCM:Proc1: ","iter1=", iter1, " / L=", round(L, 4),
                    " / max(H)=", round(max(H), 4),
                    " / max(Hm)=", round(max(H.m), 4),
                    " / max(convnew)=", round(max(shape_new$conv), 4),
                    " / max(mononew)=", round(max(shape_new$mono), 4),
                    " / IsKnot idx=", Kidx,
                    ##" / IsMIC idx=", Midx,
                    sep = ""))
        print(paste("IsMIC")); print(IsMIC)
      }
    }
    ##note at end of procedure 1 (inner loop) we have a max over the give knots.
    ## so for all knot indices j, H[j] should be 0.
    phi <- phi_new
    shape <- shape_new;
    ##prevent infinite loops (bugs could cause lambda and H to work opposite each other)
    if (sum(IsKnot != IsKnot_old) == 0 && sum(IsMIC!=IsMIC_old)==0) { ## identical(isKnot,isKnot_old)
      if (print==TRUE){
        print("No change in constraints. Ending.")
      }
      break
    }
    ##     if (print == TRUE) {
    ##       print(paste("iter1=", iter1, " / L=", round(L, 4), 
    ##                   " / max(H)=", round(max(H), 4),
    ##                   " / max(Hm)=", round(max(H.m), 4),
    ##                   " / max(convnew)=", round(max(shape_new$conv), 4),
    ##                   " / max(mononew)=", round(max(shape_new$mono), 4),
    ##                   " / IsKnot idx=", Kidx,
    ##                   " / IsMIC idx=", Midx,
    ##                   sep = ""))
    ##       ##print("H")
    ##       ##print(H)
    ##       ##print(H.m)
    ##     }
  }
  Fhat <- LocalF(z, phi)
  tmp <- LocalCoarsen.mode(z,w,IsKnot,IsMIC,a)
  MI <- z[IsKnot>0][tmp$constr]#modal interval
  KK <- (1:n)[as.logical(IsKnot)] ## Z-index-knots
  ## create the various representations of the results.
  res1 <- list(xn=xn, ##duplicates may exist
               ##x=x,## no duplicates
               ## there is currently no "no duplicates, no 'a' "
               x=z,  ## no duplicates, includes 'a'
               ##z=z, ##no duplicates, includes 'a'
               w=w,
               L=L,
               MI=MI,
               IsKnot=IsKnot,
               IsMIC=IsMIC,
               constr=tmp$constr,
               knots= z[KK],
               phi=as.vector(phi),
               fhat=as.vector(exp(phi)),
               Fhat=as.vector(Fhat),
               H=as.vector(H),
               H.m=as.vector(H.m),
               n=length(xn), # n \ge m. equal iff x has no repeats.
               m=length(z), ##either =m1 or =m1+1. notation DIFFERS FROM MANUSCRIPT!
               m1=n1,  # = length(x) or number of unique x's
               ##a=a, ## REMOVED "$a", replaced with $dlcMode!
               dlcMode=a,
               sig=sig);
  ## res1.tmp <- list(xn=xn, ##duplicates may exist
  ##                  x=z,##  X=Z here, for evaluatelogcondens!
  ##            z=z, ##no duplicates, includes 'a'
  ##            w=w,
  ##            L=L,
  ##            MI=MI,
  ##            IsKnot=IsKnot,
  ##            IsMIC=IsMIC,
  ##            constr=tmp$constr,
  ##            knots= z[KK],
  ##            phi=as.vector(phi),
  ##            fhat=as.vector(exp(phi)),
  ##            Fhat=as.vector(Fhat),
  ##            H=as.vector(H),
  ##            H.m=as.vector(H.m),
  ##            n=length(xn), # n \ge m. equal iff x has no repeats.
  ##            m=length(z), ##either =m1 or =m1+1. notation DIFFERS FROM MANUSCRIPT!
  ##            m1=n1,  # = length(x)
  ##            ##a=a, ## REMOVED "$a", replaced with $dlcMode!
  ##            dlcMode=a,
  ##            sig=sig);
  
  class(res1) <- "dlc"; ## could make a subclass ... ? "dlc.mc"
  ##class(res1.tmp) <- "dlc";
  phi.f <- function(x0){
    ##apply(matrix(x0),1, function(x00){evaluateLogConDens(x00,x=z,phi,Fhat,IsKnot)[1]})
    evaluateLogConDens(xs=x0, res=res1, which=1)[,2]
  }
  fhat.f <- function(x0){
    ##apply(matrix(x0),1,function(x00){evaluateLogConDens(x00,x=z,phi,Fhat, IsKnot)[2]})
    evaluateLogConDens(xs=x0, res=res1, which=2)[,3]
  }
  Fhat.f <- function(x0){
    ##apply(matrix(x0),1,function(x00){evaluateLogConDens(x00,x=z,phi,Fhat, IsKnot)[3]})
    evaluateLogConDens(xs=x0, res=res1, which=3)[,4]
  }
  E.f <- intFfn(z,phi,Fhat)

  ##   E.MC.L <- myE.MC;  
  ##    E.MC.R <- function(t){ ##"HR"
  ##      Xn <- tail(myxx,1)  
  ##      (Xn -t) - (E.MC.L(Xn)-E.MC.L(t));  
  ##    }

  
  EL.f <- E.f
  ER.f <- intFfn(z,phi,Fhat,side="right") ##fixing intffn here here here 
  
  {##get phiP=phi prime= deriv of phi; left and right derivs are L and R.
    phiK <- phi[as.logical(IsKnot)]
    slopes <- diff(phiK)/diff(z[KK])
    phiPL.f <- stepfun(x=z[KK], c(slopes[1],slopes,tail(slopes,1)), right=TRUE,f=1) ##DONT TRUST stepfun() DOCUMENTATION!
    phiPR.f <- stepfun(x=z[KK], c(slopes[1],slopes,tail(slopes,1)), right=FALSE,f=0)
    phiPL <- phiPL.f(z)
    phiPR <- phiPR.f(z)
    ##     phiPL2 <- c(slopes[1],rep(slopes, diff(KK)))
    ##     phiPR2 <- c(rep(slopes, diff(KK)),tail(slopes,1))
    ##     if ( !all(phiPL==phiPL2 & phiPR==phiPR2)) print("phiP error LCmode")
  }
  ##if (!is.null(logfile) && !is.na(logfile)) sink(NULL)
  if (print) {
    print("ASLCM: returning.")
  }
  return(c(res1,
           list(phi.f=phi.f, fhat.f=fhat.f, Fhat.f=Fhat.f,
                ## E.f=E.f,
                EL.f=EL.f,ER.f=ER.f, ##note that EL.f was previously known as E.f, HL.f...
                phiPL=phiPL,phiPR=phiPR,
                phiPL.f=phiPL.f,phiPR.f=phiPR.f)))
}




