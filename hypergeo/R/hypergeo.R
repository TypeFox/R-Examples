`.f3` <- function(a1,a2,a3){
  exp(complex_gamma(a1,log=TRUE) - complex_gamma(a2,log=TRUE) - complex_gamma(a3,log=TRUE))
}

`.f4` <- function(a1,a2,a3,a4){
  exp(complex_gamma(a1,log=TRUE) + complex_gamma(a2,log=TRUE) - complex_gamma(a3,log=TRUE) - complex_gamma(a4,log=TRUE))
}

"f15.1.1" <- function(A, B, C, z, tol=0, maxiter=2000){
    if(!is.null(getOption("showHGcalls"))){print(match.call())}
    genhypergeo(U=c(A,B), L=C, z=z, tol=tol, maxiter=maxiter)
}

"f15.3.1" <- function(A,B,C,z,h=0){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
# mult <- exp(lgamma(C)-lgamma(B)-lgamma(C-B))
  mult <- .f3(C,B,C-B)

  f <- function(t){t^(B-1)*(1-t)^(C-B-1)*(1-t*z)^(-A)}
  if(length(h)==1){
    if(h==0){
      return(mult * myintegrate(f,lower=0,upper=1))
    } else {
      if(is.double(h)){
        h <- 0.5 + h*1i
      }
    }
  } 
  return(mult * integrate.segments(f,c(0,h,1),close=FALSE))
}

"f15.3.3" <- function(A,B,C,z,tol=0,maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  (1-z)^(C-A-B)*genhypergeo(U=c(C-A,C-B),L=C,z=z,tol=tol,maxiter=maxiter)
}

"f15.3.4" <- function(A,B,C,z,tol=0,maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  (1-z)^(-A)*genhypergeo(U=c(A,C-B),L=C,z=z/(z-1),tol=tol,maxiter=maxiter)
}

"f15.3.5" <- function(A,B,C,z,tol=0,maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  (1-z)^(-B)*genhypergeo(U=c(B,C-A),L=C,z=z/(z-1),tol=tol,maxiter=maxiter)
}

"i15.3.6" <- function(A,B,C){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  c(
     ifelse(is.nonpos(C-A) | is.nonpos(C-B), 0, .f4(C, C-A-B, C-A,C-B)),
     ifelse(is.nonpos(A  ) | is.nonpos(B  ), 0, .f4(C, A+B-C, A  , B ))
    )
}

"j15.3.6" <- function(A,B,C){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  is.nonpos(c(
            C , C-A-B ,
            C , A+B-C 
            ))
}

"f15.3.6" <- function(A,B,C,z,tol=0,maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  if(length(z)==0){
    return(z)
  }
  jj <- i15.3.6(A,B,C)
      jj[1] * genhypergeo(U=c(  A,  B),L=A+B-C+1,z=1-z,tol=tol,maxiter=maxiter) +
      jj[2] * genhypergeo(U=c(C-A,C-B),L=C-A-B+1,z=1-z,tol=tol,maxiter=maxiter) * (1-z)^(C-A-B)
}

"i15.3.7" <- function(A,B,C){
    if(!is.null(getOption("showHGcalls"))){print(match.call())}

    c(
        ifelse(is.nonpos(B) | is.nonpos(C-A), 0, .f4(C,B-A,B,C-A)),
        ifelse(is.nonpos(A) | is.nonpos(C-B), 0, .f4(C,A-B,A,C-B))
    )
}

"j15.3.7" <- function(A,B,C){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  is.nonpos(c(
            C , B-A,
            C , A-B 
            ))
}

"f15.3.7" <- function(A,B,C,z,tol=0,maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
    if(length(z)==0){
    return(z)
  }
  jj <- i15.3.7(A,B,C)
    jj[1] * (-z)^(-A) * genhypergeo(U=c(A,1-C+A),L=1-B+A,z=1/z,tol=tol,maxiter=maxiter) +
    jj[2] * (-z)^(-B) * genhypergeo(U=c(B,1-C+B),L=1-A+B,z=1/z,tol=tol,maxiter=maxiter)
}

"i15.3.8" <- function(A,B,C){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}

   c(
      ifelse(is.nonpos(B) | is.nonpos(C-A), 0, .f4(C,B-A,B,C-A)),
      ifelse(is.nonpos(A) | is.nonpos(C-B), 0, .f4(C,A-B,A,C-B))
    )

}

"j15.3.8" <- function(A,B,C){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  is.nonpos(c(
              C , B-A , 
              C , A-B 
              ))
}

"f15.3.8" <- function(A,B,C,z,tol=0,maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  if(length(z)==0){
    return(z)
  }
  jj <- i15.3.8(A,B,C)
  return(
      jj[1] * (1-z)^(-A) * genhypergeo(U=c(A,C-B),L=A-B+1,z=1/(1-z),tol=tol,maxiter=maxiter) + 
      jj[2] * (1-z)^(-B) * genhypergeo(U=c(B,C-A),L=B-A+1,z=1/(1-z),tol=tol,maxiter=maxiter)
      )
}

"i15.3.9" <- function(A,B,C){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  return(c(
      ifelse(is.nonpos(C-A)|is.nonpos(C-B), 0, .f4(C,C-A-B,C-A,C-B)),
      ifelse(is.nonpos(  A)|is.nonpos(B)  , 0, .f4(C,A+B-C,A,    B))
    ))
} 
"j15.3.9" <- function(A,B,C){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  is.nonpos(c(
              C , C-A-B , 
              C , A+B-C 
              ))
}              

"f15.3.9" <- function(A,B,C,z,tol=0,maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  if(length(z)==0){
    return(z)
  }
  jj <- i15.3.9(A,B,C)
    jj[1] *               z^( -A)*genhypergeo(U=c(A,A-C+1),L=A+B-C+1,z=1-1/z,tol=tol,maxiter=maxiter) +
    jj[2] * (1-z)^(C-A-B)*z^(A-C)*genhypergeo(U=c(C-A,1-A),L=C-A-B+1,z=1-1/z,tol=tol,maxiter=maxiter)
  }

"isgood" <- function(x,tol){ all(abs(x[!is.na(x)]) <= tol)}

"genhypergeo" <- function (U, L, z, tol = 0, maxiter=2000, check_mod=TRUE, polynomial=FALSE, debug=FALSE, series=TRUE)
{
  if(series){
    return(genhypergeo_series(U, L, z, tol = tol, maxiter=maxiter, check_mod=check_mod, polynomial=polynomial, debug=debug))
  } else {
    return(genhypergeo_contfrac(U, L, z, maxiter=maxiter))
  }
}

"genhypergeo_series" <-
function (U, L, z, tol = 0, maxiter=2000, check_mod=TRUE, polynomial=FALSE, debug=FALSE) 
{
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  if(debug){
    stopifnot(length(z)==1)
    out <- NULL
  }
  
  if(check_mod){
    lU <- length(U)
    lL <- length(L)
    
    if(lU > lL+1){
      greater <- Mod(z)>0
    } else if(lU > lL) {
      greater <- Mod(z)>1
    } else {
      greater  <- Mod(z)<0
    }
    if(all(greater)){
      return(z*NA)
    } else {
      z[greater] <- NA
    }
  }
  
  fac <- 1
  temp <- fac
  if(debug){out <- temp}
  
  if(maxiter==0){
    return(z*0+fac)
  }
  for (n in seq_len(maxiter)) {
    fac <- fac * (prod(U)/prod(L)) * (z/n)
    series <- temp + fac
    if(debug){out <- c(out,fac)}
    if (isgood(series-temp,tol)){
      if(debug){
        return(list(series,out))
      } else {

        return(series)
      }
    }
    temp <- series
    
    U <- U + 1

      
    L <- L + 1
  }

  if(debug){
    return(list(series,out))
  }
  
  if(polynomial){
    return(series)
  } else {
    warning("series not converged")
    return(z*NA)
  }
}

"hypergeo_taylor" <- function(A, B, C, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  genhypergeo(U=c(A,B), L=C, z=z, tol=tol, maxiter=maxiter, check_mod=FALSE, polynomial=TRUE)
}

"is.near_integer" <- function(i , tol=getOption("tolerance")){
  if(is.null(tol)){
    tol <- 1e-11
  }
  abs(i-round(Re(i))) <= tol
}

"is.nonpos" <- function(i){
  is.near_integer(i) & (Re(i) < 0.5)
}

"is.zero" <- function(i){
  is.near_integer(i) & (abs(i) < 0.5)
}

"hypergeo_A_nonpos_int" <- function(A, B, C, z, tol=0){
                                        # Assumed: A integer <=0 , B
                                        # either non-integer or (if an
                                        # integer) <= A. (for example:
                                        # A = -2, B = -5).  The
                                        # hypergeometric series is a
                                        # polynomial.
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(is.nonpos(A))
  if(( is.near_integer(C) ) & is.near_integer(C) & (abs(C-A) < 0.5) ){  # A==C==integer
    warning("this case is not uniquely defined: proceed, assuming both A and C approach the same nonpositive integer at the same speed [that is, (a)_n cancels (c)_n for all 'n']")
    return(genhypergeo(U=B,L=NULL,z,tol=tol,check_mod = FALSE))
  } else {
    return(hypergeo_taylor(A,B,C,z,tol=tol,maxiter = -A))
  }
}

"hypergeo_AorB_nonpos_int" <- function(A, B, C, z, tol=0){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(is.nonpos(A) | is.nonpos(B))
  if(is.nonpos(A) & is.nonpos(B)){
    if(A>B){  # eg A = -2,  B = -5
      return(hypergeo_A_nonpos_int(A,B,C,z,tol=tol)) # Note A,B not swapped over
    } else {
        return(hypergeo_A_nonpos_int(B,A,C,z,tol=tol)) # Note A,B swapped over
    }
  }
  
  ## Thus from here on, A is a nonpositive integer and B is not an
  ## integer.
  if(is.nonpos(A)){
    return(hypergeo_A_nonpos_int(A,B,C,z,tol=tol))
  } else {  # Former bug!
    return(hypergeo_A_nonpos_int(B,A,C,z,tol=tol))
  }
}


".process_args" <- function(...){  # slight modification of process.args() of package gsl...
  a <- list(...)
  attr <- attributes(a[[which.max(unlist(lapply(a,length)))]])
  a <- lapply(a,as.vector)
  out <- do.call("cbind",a)
  
  return(list(out=out, attr = attr))
}

"crit" <- function(...){
    c(
        1/2 + 1i*sqrt(3)/2,
        1/2 - 1i*sqrt(3)/2
        )
}
   
"hypergeo" <- function(A, B, C, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  
  if(length(A)>1 | length(B)>1 | length(C)>1){
    jj <- .process_args(A,B,C,z)
    f <- function(x){hypergeo(A=Re(x[1]), B=Re(x[2]),C=Re(x[3]),z=x[4],tol=tol,maxiter=maxiter)}
    out <- apply(jj$out , 1, f)
    attributes(out) <- jj$attr
    return(out)
  }


  # if you are here, length(A)=length(B)=length(C)=1.

  jj <- crit()
  c1 <- jj[1]
  c2 <- jj[2]
  
  close_to_crit <- (abs(z-c1) < 0.1) | (abs(z-c2) < 0.1)

  ## following lines commented out because ifelse() evaluates both
  ## functions for *every* value of z, irregardless of the value of
  ## close_to_crit.  So both hypergeo_residue_close_to_crit() *and*
  ## hypergeo_powerseries() return errors [and there is also the risk
  ## of an infinite regress].
  ## out <- ifelse(close_to_crit,
  ##              hypergeo_residue_close_to_crit_multiple(A,B,C,z, tol=tol, maxiter=maxiter),
  ##              hypergeo_powerseries                   (A,B,C,z, tol=tol, maxiter=maxiter)
  ##              )

  out <- z*NA
#  if(any( close_to_crit)){out[ close_to_crit] <- hypergeo_residue_close_to_crit_multiple(A,B,C,z[ close_to_crit], tol=tol, maxiter=maxiter)}
#  if(any(!close_to_crit)){out[!close_to_crit] <- hypergeo_powerseries                   (A,B,C,z[!close_to_crit], tol=tol, maxiter=maxiter)}
  
  if(any( close_to_crit)){out[ close_to_crit] <- hypergeo_gosper      (A,B,C,z[ close_to_crit], tol=tol, maxiter=maxiter)}
  if(any(!close_to_crit)){out[!close_to_crit] <- hypergeo_powerseries (A,B,C,z[!close_to_crit], tol=tol, maxiter=maxiter)}
  
  do_with_cf <- !is.na(z) & is.na(out)   # ie failures to converge; do_with_cf == "do with Continued Fraction"
  if(any(do_with_cf)){
    out[do_with_cf] <- hypergeo_contfrac(A=A, B=B, C=C, z=z[do_with_cf], maxiter=maxiter)
  }
  do_with_integration <- !is.na(z) & is.na(out)
  if(any(do_with_integration)){
    g <- function(z){f15.3.1(A=A, B=B, C=C, z=z)}
    out[do_with_integration] <- sapply(z[do_with_integration] , g)
  }
  
  return(out)
}

"hypergeo_powerseries" <- function(A, B, C, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}

  z <- z+0i

  if(is.zero(A) | is.zero(B)){
    if(is.zero(C)){
      return(z*NA)
    } else {
      return(z*0+1)
    }
  }
  
  if(is.zero(C)){
    return(z*Inf)
  }
  
  if(is.zero(A-C)){
    return( (1-z)^(-B) )
  } else if (is.zero(B-C)){
    return( (1-z)^(-A) )
  }
    
  if(is.nonpos(A) | is.nonpos(B)){
    return(hypergeo_AorB_nonpos_int(A,B,C,z,tol=tol))
  }

  if(is.nonpos(C)){    # C is a nonpositive integer; series not defined [unless it terminates in which case a limit is used]
      return(z*NA)

  }

  ## So from here on, A, B, C are either non-integer, or integers >0.

  if(Re(A) > Re(B)){
    swap <- A
    A    <- B
    B    <- swap
  }  # So from here on,  A <= B
    
    
  m <- C-A
  n <- B-A   #  remember: 'n' must be >= 0 because of the 'swap' above.
  
  if(is.near_integer(m)){
    if(m <= 0){
      return( (1-z)^(C-A-B)*Recall(C-A,C-B,C,z=z,tol=tol,maxiter=maxiter) )  # This is 15.3.3, but do not call f15.3.3(), because this leads to an infinite recursion
    } else {
      if(is.near_integer(n)){ # This means B-A and C-A are both integers; the "limiting process" on p560 [just after 15.3.4] needs hypergeo_cover3()
        return(hypergeo_cover3(A,n,m,z,tol=tol,maxiter=maxiter))
      } 
    }
  } 
  
  m <- -(A+B-C)   # Former bug!

  if(is.near_integer(m)){  # This is the "Each term of 15.3.6 has a pole..." on p559
    return(hypergeo_cover1(A,B,m,z,tol=tol,maxiter=maxiter))
  }
  
  m <- B-A
  if(is.near_integer(m)){  # This is the "Similarly each term of 15.3.7..." on p560
    return(hypergeo_cover2(A,C,m,z,tol=tol,maxiter=maxiter))
  }
  
  return(hypergeo_general(A,B,C,z,tol=tol,maxiter=maxiter))
}

"hypergeo_general" <- function(A, B, C, z, tol=0, maxiter=2000,  give=FALSE){

  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  attr <- attributes(z)
  z <- as.vector(as.complex(z))
  
  things <- thingfun(z)
  choice <- apply(things,1,which.min)
  if(!is.null(getOption("showHGcalls"))){
    print("choice: ")
    print(choice)
  }
  
  u15.1.1 <- choice==1
  u15.3.4 <- choice==2
  u15.3.6 <- choice==3
  u15.3.7 <- choice==4
  u15.3.8 <- choice==5
  u15.3.9 <- choice==6

  out <- z*NA
  if(any(u15.1.1)){ out[u15.1.1] <- f15.1.1(A=A,B=B,C=C, z[u15.1.1], tol=tol,maxiter=maxiter) }  # 1
  if(any(u15.3.4)){ out[u15.3.4] <- f15.3.4(A=A,B=B,C=C, z[u15.3.4], tol=tol,maxiter=maxiter) }  # 2
  if(any(u15.3.6)){ out[u15.3.6] <- f15.3.6(A=A,B=B,C=C, z[u15.3.6], tol=tol,maxiter=maxiter) }  # 3 
  if(any(u15.3.7)){ out[u15.3.7] <- f15.3.7(A=A,B=B,C=C, z[u15.3.7], tol=tol,maxiter=maxiter) }  # 4
  if(any(u15.3.8)){ out[u15.3.8] <- f15.3.8(A=A,B=B,C=C, z[u15.3.8], tol=tol,maxiter=maxiter) }  # 5 
  if(any(u15.3.9)){ out[u15.3.9] <- f15.3.9(A=A,B=B,C=C, z[u15.3.9], tol=tol,maxiter=maxiter) }  # 6

  attributes(out) <- attr
  if(give){
    return(list(choice,out))
  } else {
    return(out)
  }
}

"thingfun" <- function(z,complex=FALSE){
    things <- cbind("z"       = z,        # 1
                    "z/(z-1)" = z/(z-1),  # 2
                    "1-z"     = 1-z,      # 3
                    "1/z"     = 1/z,      # 4
                    "1/(1-z)" = 1/(1-z),  # 5
                    "1-1/z"   = 1-1/z     # 6
                    )

    if(complex){return(things)}

    things <- Mod(things)

  if(any(apply(things,1,min, na.rm=TRUE)>1)){ # Thanks to Igor Kojanov for fixing this
    stop("odd: none of the transformations take the argument inside the unit disk.  Contact the package maintainer")
  }
  return(things)
}

"hypergeo_cover1" <- function(A, B, m, z, tol=0, maxiter=2000,  method="a", give=FALSE){
  
  ## use equation 15.3.3 - 15.3.9 EXCEPT 15.3.6, which has a pole when
  ## a+b-c is an integer.  See the bit between 15.3.9 and 15.3.10,
  ## p559.
  
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  
  stopifnot(is.near_integer(m))
  C <- A+B+m

  attr <- attributes(z)
  z <- as.vector(as.complex(z))

  things <- thingfun(z)

  ## Now to discourage bad ones:
  if(any(j15.3.7(A,B,C))){ things[,4] <- Inf }
  if(any(j15.3.8(A,B,C))){ things[,5] <- Inf }
  if(any(j15.3.9(A,B,C))){ things[,6] <- Inf }
  ## thus we take the minimum modulus of non-forbidden options.
  ## Compare similar lines in hypergeo_cover2(): here the functions
  ## are 7,8,9; there they are 6,8,9
  
  choice <- apply(things,1,which.min)
  u15.1.1 <- choice==1
  u15.3.4 <- choice==2
  u15.3.x <- choice==3  # This one!  [corresponds to u15.3.6()]
  u15.3.7 <- choice==4
  u15.3.8 <- choice==5
  u15.3.9 <- choice==6


  out <- z*NA
  if(any(u15.1.1)){ out[u15.1.1] <- f15.1.1       (A=A,B=B,C=C, z[u15.1.1], tol=tol,maxiter=maxiter) }
  if(any(u15.3.4)){ out[u15.3.4] <- f15.3.4       (A=A,B=B,C=C, z[u15.3.4], tol=tol,maxiter=maxiter) }
  if(any(u15.3.x)){ out[u15.3.x] <- f15.3.10_11_12(A=A,B=B,m=m, z[u15.3.x], tol=tol,maxiter=maxiter, method=method) }
  if(any(u15.3.7)){ out[u15.3.7] <- f15.3.7       (A=A,B=B,C=C, z[u15.3.7], tol=tol,maxiter=maxiter) }
  if(any(u15.3.8)){ out[u15.3.8] <- f15.3.8       (A=A,B=B,C=C, z[u15.3.8], tol=tol,maxiter=maxiter) }
  if(any(u15.3.9)){ out[u15.3.9] <- f15.3.9       (A=A,B=B,C=C, z[u15.3.9], tol=tol,maxiter=maxiter) }
  
  attributes(out) <- attr
  if(give){
    return(list(choice,out))
  } else {
    return(out)
  }
}

"hypergeo_cover2" <- function(A, C, m, z, tol=0, maxiter=2000, method="a", give=FALSE){

  if(!is.null(getOption("showHGcalls"))){print(match.call())}

  ## use equation 15.3.3 - 15.3.9 EXCEPT 15.3.7, which has a pole when
  ## a+b-c is an integer.  See the bit between 15.3.13 and 15.3.15,
  ## p559.

  stopifnot(is.near_integer(m))

  B <- A+m
  
  attr <- attributes(z)
  z <- as.vector(as.complex(z))

  things <- thingfun(z)


  ## Now to discourage bad ones:
  if(any(j15.3.6(A,B,C))){ things[,3] <- Inf }
  if(any(j15.3.8(A,B,C))){ things[,5] <- Inf }
  if(any(j15.3.9(A,B,C))){ things[,6] <- Inf }  
  
  choice <- apply(things,1,which.min)
  u15.1.1 <- choice==1
  u15.3.4 <- choice==2
  u15.3.6 <- choice==3
  u15.3.x <- choice==4 #  This one!  [corresponds to u15.3.7()]
  u15.3.8 <- choice==5
  u15.3.9 <- choice==6

  out <- z*NA
  if(any(u15.1.1)){ out[u15.1.1] <- f15.1.1    (A=A,B=B,C=C, z[u15.1.1], tol=tol,maxiter=maxiter) }
  if(any(u15.3.4)){ out[u15.3.4] <- f15.3.4    (A=A,B=B,C=C, z[u15.3.4], tol=tol,maxiter=maxiter) }
  if(any(u15.3.6)){ out[u15.3.6] <- f15.3.6    (A=A,B=B,C=C, z[u15.3.6], tol=tol,maxiter=maxiter) }
  if(any(u15.3.x)){ out[u15.3.x] <- f15.3.13_14(A=A,C=C,m=m, z[u15.3.x], tol=tol,maxiter=maxiter, method=method) }
  if(any(u15.3.8)){ out[u15.3.8] <- f15.3.8    (A=A,B=B,C=C, z[u15.3.8], tol=tol,maxiter=maxiter) }
  if(any(u15.3.9)){ out[u15.3.9] <- f15.3.9    (A=A,B=B,C=C, z[u15.3.9], tol=tol,maxiter=maxiter) }
  
  attributes(out) <- attr
  if(give){
    return(list(choice,out))
  } else {
    return(out)
  }
}

"hypergeo_cover3" <- function(A, n, m, z, tol=0, maxiter=2000, method="a", give=FALSE){

  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(is.near_integer(n))
  stopifnot(is.near_integer(m))
  
  attr <- attributes(z)
  z <- as.vector(as.complex(z))
  
  ## following is a cut-down version of thingfun(), tailored for the Wolfram functions:
  things <- Mod(cbind(
                      "z"       = z,        # 1
                      "1/z"     = 1/z       # 4
                      )
                )

  if(any(apply(things,1,min,na.rm=TRUE)>1)){
    stop("odd: none of the transformations take the argument inside the unit disk.  Contact the package maintainer")
  }
  
  choice <- apply(things,1,which.min)
  
  
  u15.1.1           <- choice==1
  u07.23.06.0026.01 <- (choice==2) & (m >  n)
  u07.23.06.0031.01 <- (choice==2) & (m <= n)
  
  out <- z*NA
  if(any(u15.1.1)){ out[u15.1.1] <- f15.1.1(A=A,B=A+n,C=A+m, z[u15.1.1], tol=tol,maxiter=maxiter) }
  if(any(u07.23.06.0026.01)){
    out[u07.23.06.0026.01] <- w07.23.06.0026.01(A=A,n,m, z[u07.23.06.0026.01], tol=tol, maxiter=maxiter, method=method)
  }
  if(any(u07.23.06.0031.01)){
    out[u07.23.06.0031.01] <- w07.23.06.0031.01(A=A,n,m, z[u07.23.06.0031.01], tol=tol, maxiter=maxiter)
  }
  
  attributes(out) <- attr
  if(give){
    return(list(choice,out))
  } else {
    return(out)
  }
}

"f15.3.10_a" <- function(A, B, z, tol=0, maxiter=2000){ #"_a" means use psigamma, "_b" means use 6.3.5, p258
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  U <- c(A,B)
  z[Mod(1-z) >= 1]  <- NA
  fac <- 1
  l1mz <- log(1+0i-z)
  
  temp <- 2*psigamma(0+1)-psigamma(A+0)-psigamma(B+0)-l1mz  # n=0
  for(n in seq_len(maxiter)){
    fac <- fac * prod(U) * ((1-z)/n^2)
    series <-
      temp + fac * (2*psigamma(n+1)- psigamma(A+n) - psigamma(B+n) - l1mz)
    if(isgood(series-temp,tol)){
      return(series/beta(A,B))
    }
    temp <- series
    U <- U+1
  }
  warning("series not converged")
  return(z*NA)
}

"f15.3.10_b" <- function(A, B, z, tol=0, maxiter=2000){ #"_a" means use psigamma, "_b" means use 6.3.5, p258
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  U <- c(A,B)
  z[Mod(1-z) >= 1]  <- NA
  fac <- 1
  pn <- psigamma(1)
  pa <- psigamma(A)
  pb <- psigamma(B)

  l1mz <- log(1+0i-z)
  
  temp <- 2*pn-pa-pb-l1mz  # n=0
  for(n in seq_len(maxiter)){
    fac <- fac * prod(U) * ((1-z)/n^2)
    pn <- pn + 1/n
    pa <- pa + 1/(A+n-1)  # no repeated psigamma() calls; cf 6.3.2, 6.3.5, p258
    pb <- pb + 1/(B+n-1)
    series <-
      temp + fac * (2*pn - pa - pb - l1mz)
    if(isgood(series-temp,tol)){
      return(series/beta(A,B))
    }
    temp <- series
    U <- U+1
  }
  warning("series not converged")
  return(z*NA)
}

"f15.3.10" <- function(A, B, z, tol=0, maxiter=2000, method="a"){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  switch(method,
         a = f15.3.10_a(A,B,z,tol=tol,maxiter=2000),
         b = f15.3.10_b(A,B,z,tol=tol,maxiter=2000),
         stop("method must be either 'a' or 'b'")
         )
}

"f15.3.11_bit1" <- function(A, B, m, z, tol=0){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(is.near_integer(m))
  stopifnot(m>0)
  m <- round(m)
  
  U <- c(A,B)
  L <- 1-m
# mult <- exp(lgamma(m)+lgamma(A+B+m)-lgamma(A+m)-lgamma(B+m))
  mult <- .f4(m,A+B+m,A+m,B+m)
  series <- z*0+1
  z[Mod(1-z)>1] <- NA
  fac <- 1
  temp <- fac
  for (n in seq_len(m-1)) {
    fac <- fac * (prod(U)/prod(L)) * (1-z)/n
    series <- temp + fac
    if (isgood(series-temp,tol)){
      return(series * mult)
    }
    temp <- series
    
    U <- U + 1
    L <- L + 1
  }
  return(series*mult)
}

"f15.3.11_bit2_a" <- function(A, B, m, z, tol=0, maxiter=2000){  #"_a" means use psigamma, "_b" means use 6.3.5.
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(is.near_integer(m))
  stopifnot(m>0)

  U <- c(A+m , B+m)  # sic

  z[Mod(1-z) >= 1]  <- NA
  fac <- 1/factorial(m)
  l1mz <- log(1+0i-z)
  temp <- (l1mz-psigamma(0+1)-psigamma(0+m+1) + psigamma(A+0+m) + psigamma(B+0+m) ) * fac
  for(n in seq_len(maxiter)){
    fac <- fac * prod(U) * (1-z)/(n*(n+m))
    series <-
      temp + fac * (l1mz - psigamma(n+1) - psigamma(n+m+1) + psigamma(A+n+m) + psigamma(B+n+m))
    if(isgood(series-temp,tol)){
#      return((z-1)^m * exp(lgamma(A+B+m)-lgamma(A)-lgamma(B)) * series)
       return((z-1)^m * .f3(A+B+m,A,B) * series)
    }
    temp <- series
    U <- U+1
  }
  warning("series not converged")
  return(z*NA)
}

"f15.3.11_bit2_b" <- function(A, B, m, z, tol=0, maxiter=2000){  # "_a" means use psigamma, "_b" means use 6.3.5.
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(is.near_integer(m))
  stopifnot(m>0)

  U <- c(A+m , B+m)  # sic
  z[Mod(1-z) >= 1]  <- NA
  fac <- 1/factorial(m)
  pn <- psigamma(  1)
  pm <- psigamma(m+1)
  pa <- psigamma(m+A)
  pb <- psigamma(m+B)

  l1mz <- log(1+0i-z)
  
  temp <- (l1mz - pn - pm + pa + pb ) * fac
  for(n in seq_len(maxiter)){
    fac <- fac * prod(U) * (1-z)/(n*(n+m))
    pn <- pn + 1/n
    
    pm <- pm + 1/(n+m)
    pa <- pa + 1/(A+n+m-1)
    pb <- pb + 1/(B+n+m-1)
    series <-
      temp + fac * (l1mz - pn - pm + pa + pb)
    if(isgood(series-temp,tol)){
#      return((z-1)^m * exp(lgamma(A+B+m)-lgamma(A)-lgamma(B)) * series)
       return((z-1)^m * .f3(A+B+m,A,B) * series)
    }
    temp <- series
    U <- U+1
  }
  warning("series not converged")
  return(z*NA)
}

"f15.3.11" <- function(A,B,m,z,tol=0, maxiter=2000,method="a"){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  switch(method,
         a = f15.3.11_bit1(A,B,m,z,tol=tol) -  f15.3.11_bit2_a(A,B,m,z,tol=tol, maxiter=maxiter),
         b = f15.3.11_bit1(A,B,m,z,tol=tol) -  f15.3.11_bit2_a(A,B,m,z,tol=tol, maxiter=maxiter),
         stop("method must be either 'a' or 'b'")
       )
}

"f15.3.12_bit1" <- function(A, B, m, z, tol=0){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  m <- round(m)
  U <- c(A-m,B-m)
  L <- 1-m
# mult <- exp(lgamma(m)+lgamma(A+B-m)-lgamma(A)-lgamma(B)) / (1-z)^m
  mult <- .f4(m,A+B-m,A,B) / (1-z)^m
  z[Mod(1-z)>1] <- NA
  fac <- 1
  temp <- fac
  series <- z*0+1
  for (n in seq_len(m-1)) {
    fac <- fac * (prod(U)/prod(L)) * (1-z)/n
    series <- temp + fac
    if (isgood(series-temp,tol)){
      return(series * mult)
    }
    temp <- series
    
    U <- U + 1
    L <- L + 1
  }
  return(series*mult)
}

"f15.3.12_bit2_a" <- function(A, B, m, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  m <- round(m)
  if(is.nonpos(A-m)|is.nonpos(B-m)){return(z*0)}
# mult <- (-1)^m * exp(lgamma(A+B-m)-lgamma(A-m)-lgamma(B-m))
  mult <- (-1)^m * .f3(A+B-m,A-m,B-m)

  U <- c(A , B)  # sic
  z[Mod(1-z) >= 1]  <- NA
  fac <- 1/factorial(m)
  l1mz <- log(1+0i-z)
  temp <- (l1mz-psigamma(1)-psigamma(m+1) + psigamma(A) + psigamma(B) ) * fac
  for(n in seq_len(maxiter)){

    fac <- fac * prod(U) * (1-z)/(n*(n+m))
    series <-
      temp + fac * (l1mz - psigamma(n+1) - psigamma(n+m+1) + psigamma(A+n) + psigamma(B+n))
    if(isgood(series-temp,tol)){
      return(mult * series)
    }
    temp <- series
    U <- U+1
  }
  warning("series not converged")
  return(z*NA)
}

"f15.3.12_bit2_b" <- function(A, B, m, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  m <- round(m)
  if(is.nonpos(A-m)|is.nonpos(B-m)){return(z*0)}
#  mult <- (-1)^m * exp(lgamma(A+B-m)-lgamma(A-m)-lgamma(B-m))
   mult <- (-1)^m * .f3(A+B-m,A-m,B-m)

  U <- c(A , B)  # sic
  z[Mod(1-z) >= 1]  <- NA
  fac <- 1/factorial(m)
  pn <- psigamma(1)
  pm <- psigamma(m+1)
  pa <- psigamma(A)
  pb <- psigamma(B)
  l1mz <- log(1+0i-z)
  temp <- (l1mz-pn - pm + pa + pb ) * fac
  for(n in seq_len(maxiter)){
    fac <- fac * prod(U) * (1-z)/(n*(n+m))
    pn <- pn + 1/n
    pm <- pm + 1/(n+m)
    pa <- pa + 1/(A+n-1)
    pb <- pb + 1/(B+n-1)
    series <-
      temp + fac * (l1mz - pn - pm + pa + pb)
    if(isgood(series-temp,tol)){
      return(mult * series)
    }
    temp <- series
    U <- U+1
  }
  warning("series not converged")
  return(z*NA)
}

"f15.3.12" <- function(A, B, m, z, tol=0, maxiter=2000, method = "a"){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  switch(method,
                a =  f15.3.12_bit1(A,B,m,z,tol=tol) -  f15.3.12_bit2_a(A,B,m,z,tol=tol, maxiter=maxiter),
                b =  f15.3.12_bit1(A,B,m,z,tol=tol) -  f15.3.12_bit2_b(A,B,m,z,tol=tol, maxiter=maxiter),
                stop("method must be one of 'a' or 'b'")
                )
}

"f15.3.13" <- function(A, C, z, tol=0, maxiter=2000, method = "a"){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  switch(method,
         a = f15.3.13_a(A,C,z,tol=tol,maxiter=maxiter),
         b = f15.3.13_b(A,C,z,tol=tol,maxiter=maxiter)
         )
}

"f15.3.13_a" <- function(A, C, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  U <- c(A,1-C+A)
  z[Mod(z) < 1]  <- NA
  fac <- 1
  pn <- psigamma(1)
  pa <- psigamma(A)
  pc <- psigamma(C-A)
  lmz <- log(0i-z)
  temp <- lmz + 2*psigamma(1) - psigamma(A) - psigamma(C-A)  # n=0
  for(n in seq_len(maxiter)){
    fac <- fac * prod(U) / (z*n^2)
    series <- temp + fac * (lmz + 2*psigamma(n+1) - psigamma(A+n) - psigamma(C-A-n))
    if(isgood(series-temp,tol)){
#      return(series * exp(lgamma(C)-lgamma(A)-lgamma(C-A)) * (0i-z)^(-A))
       return(series * .f3(C,A,C-A) * (0i-z)^(-A))
    }
    temp <- series
    U <- U+1
  }
  warning("series not converged")
  return(z*NA)
}

"f15.3.13_b" <- function(A, C, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  U <- c(A,1-C+A)
  z[Mod(z) < 1]  <- NA
  fac <- 1
  pn <- psigamma(1)
  pa <- psigamma(A)
  pc <- psigamma(C-A)
  lmz <- log(0i-z)
  temp <- lmz + 2*pn - pa - pc
  for(n in seq_len(maxiter)){
    fac <- fac * prod(U) / (z*n^2)
    pn <- pn + 1/n
    pa <- pa + 1/(A+n-1)  
    pc <- pc - 1/(C-A-n)  # The term is psi(c-a-n), not psi(c-a+n)
    series <- temp + fac * (lmz + 2*pn - pa - pc)
    if(isgood(series-temp,tol)){
#      return(series * exp(lgamma(C)-lgamma(A)-lgamma(C-A)) * (0i-z)^(-A))
       return(series * .f3(C,A,C-A) * (0i-z)^(-A))
    }
    temp <- series
    U <- U+1
  }
  warning("series not converged")
  return(z*NA)
}

"f15.3.14_bit1_a" <- function(A, C, m, z, tol=0, maxiter=2000){ # "_a" means use psigamma, "_b" means use 6.3.5.
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  m <- round(m)
  U <- c(A+m, 1-C+A+m)
  z[Mod(z) < 1]  <- NA
#  fac <-  exp(lgamma(A+m)-lgamma(A)+ lgamma(1-C+A+m)-lgamma(1-C+A)  - lfactorial(m))
  fac <-  exp(
      +complex_gamma(A+m,log=TRUE)
      -complex_gamma(A,log=TRUE)
      +complex_gamma(1-C+A+m,log=TRUE)
      -complex_gamma(1-C+A,log=TRUE)
      -complex_factorial(m,log=TRUE)
      )
  lmz <- log(0i-z) 
  temp <- (lmz + psigamma(1+m) + psigamma(1) - psigamma(A+m) - psigamma(C-A-m)) * fac
  for(n in seq_len(maxiter)){
    fac <- fac * prod(U) / (z*n*(n+m))
    series <- 
      temp + fac * (lmz + psigamma(1+m+n) + psigamma(1+n) - psigamma(A+m+n) - psigamma(C-A-m-n))
    if(isgood(series-temp,tol)){
#     return(  (0i-z)^(-A-m) * exp(lgamma(C) -lgamma(A+m)-lgamma(C-A)) * series)
      return(  (0i-z)^(-A-m) * .f3(C,A+m,C-A) * series)
    }
    temp <- series
    U <- U+1
  }
  warning("series not converged")
  return(z*NA)
}

"f15.3.14_bit1_b" <- function(A, C, m, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  m <- round(m)
  U <- c(A+m, 1-C+A+m)
  z[Mod(z) < 1]  <- NA
#  fac <-  exp(lgamma(A+m)-gamma(A) + lgamma(1-C+A+m)-lgamma(1-C+A)- factorial(m))  #NB gamma(A) should be lgamma(A)
  fac <-  exp(
      +complex_gamma(A+m,log=TRUE)
      -complex_gamma(A,log=TRUE)
      +complex_gamma(1-C+A+m,log=TRUE)
      -complex_gamma(1-C+A)
      -complex_factorial(m,log=TRUE)
      )
  pm <- psigamma(m+1)  
  pn <- psigamma(1)
  pa <- psigamma(m+A)
  pc <- psigamma(C-A-m)
  lmz <- log(0i-z) 
  temp <- (lmz + pm + pn - pa - pc) * fac
  for(n in seq_len(maxiter)){
    fac <- fac * prod(U) / (z*n*(n+m))
    pm <- pm + 1/(n+m)
    pn <- pn + 1/n
    pa <- pa + 1/(m+A+n-1)
    pc <- pc - 1/(C-A-m-n)
    series <-
      temp + fac * (lmz + pm + pn - pa - pc)
    if(isgood(series-temp,tol)){
#      return(  (0i-z)^(-A-m) * exp(lgamma(C)-lgamma(A+m)-gamma(C-A)) * series)
      return(  (0i-z)^(-A-m) * .f3(C,A+m,C-A) * series)
    }
    temp <- series
    U <- U+1
  }

  warning("series not converged")
  return(z*NA)
}
  
"f15.3.14_bit2" <- function(A, C, m, z, tol=0){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  m <- round(m)
  stopifnot(m>0)
  stopifnot(is.near_integer(m))
  U <- c(A)
#  mult <- (0i-z)^(-A) * exp(lgamma(C) - lgamma(A+m))
  mult <- (0i-z)^(-A) * .f3(C,A+m,1)  # NB log(gamma(1))=0
  z[Mod(z)<1] <- NA
  fac <- 1
  temp <- gamma(m)/gamma(C-A)
  series <- z*0+temp
  for (n in seq_len(m-1)) {
    fac <- fac * prod(U) / (z*n)
#    series <- temp + fac * exp(lgamma(m-n)-lgamma(C-A-n))
    series <- temp + fac * .f3(m-n,C-A-n,1)
    if (isgood(series-temp,tol)){
      return(series * mult)
    }
    temp <- series
    U <- U + 1
  }
  return(series*mult)
}

"f15.3.14" <- function(A, C, m, z, tol=0, maxiter=2000, method="a"){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  a1 <- f15.3.14_bit1_a(A,C,m,z,tol=tol,maxiter=maxiter)
  a2 <- f15.3.14_bit2(A,C,m,z,tol=tol)
  switch(method,
         a=f15.3.14_bit1_a(A,C,m,z,tol=tol,maxiter=maxiter) + f15.3.14_bit2(A,C,m,z,tol=tol),
         b=f15.3.14_bit1_b(A,C,m,z,tol=tol,maxiter=maxiter) + f15.3.14_bit2(A,C,m,z,tol=tol),
         stop("method must be one of 'a' or 'b'")
         )
}

"f15.3.10_11_12" <- function(A,B,m,z,tol=0,maxiter=2000,method="a"){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(is.near_integer(m))
  m <- round(m)
  if(is.zero(m)){
    return(f15.3.10(A,B,   z,tol=tol,maxiter=maxiter,method=method))
  } else if (m>0){
    return(f15.3.11(A,B, m,z,tol=tol,maxiter=maxiter,method=method))
  } else if (m<0){
    return(f15.3.12(A,B,-m,z,tol=tol,maxiter=maxiter,method=method))
  } else {
    stop("this cannot happen")
  }
}

"f15.3.13_14" <- function(A, C, m, z, tol=0, maxiter=2000, method="a"){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(is.near_integer(m))
  m <- round(m)
  if(is.zero(m)){
    return(f15.3.13(A  ,C   ,z,tol=tol,maxiter=maxiter,method=method))
  } else if (m>0){
    return(f15.3.14(A  ,C, m,z,tol=tol,maxiter=maxiter,method=method))
  } else if (m<0){
    return(f15.3.14(A+m,C,-m,z,tol=tol,maxiter=maxiter,method=method)) #F(a,b,c;z)==F(b,a,c;z)
  } else {
    stop("this cannot happen")
  }
}

"w07.23.06.0029.01" <- function(A, n, m, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  ((-1)^m*gamma(A-m)*factorial(m+n)*(0i-z)^(-A-n)/(gamma(A)*factorial(n)))*
    hypergeo(A+n , m+n+1, n+1, 1/z,tol=tol,maxiter=maxiter)
}

"w07.23.06.0031.01_bit1" <- function(A, n, m, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(is.near_integer(m))
  stopifnot(m>0)
  
  U <- c(A,1-m)
  L <- 1-n

# mult <- exp(lgamma(A+m)+lgamma(n) - lgamma(m)-lgamma(A+n)) * (0i-z)^(-A)
  mult <- .f4(A+m,n,m,A+n) * (0i-z)^(-A)
  series <- z*0+1
  z[Mod(z) < 1] <- NA
  fac <- 1
  temp <- fac
  for (k in seq_len(m-1)) {  # Note iteration is over "k", not "n", as per 07.23.06.0031.01
    fac <- fac * (prod(U)/prod(L)) / (k*z)
    series <- temp + fac
    if (isgood(series-temp,tol)){
      return(series * mult)
    }
    temp <- series
    
    U <- U + 1
    L <- L + 1
  }
  return(series*mult)
}

"w07.23.06.0031.01_bit2" <- function(A, n, m, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
# (-1)^m* (gamma(A+m)/gamma(A)) * factorial(n-m) *(0i-z)^(-A-n) / factorial(n) *
  (-1)^m*(0i-z)^(-A-n) * .f4(A+m,n-m+1,A,n+1)*
    hypergeo(A+n , 1-m+n , n+1 , 1/z , tol=tol , maxiter=maxiter)
}

"w07.23.06.0031.01" <- function(A, n, m, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(m <= n)
    w07.23.06.0031.01_bit1(A, n, m, z, tol=tol, maxiter=maxiter) + 
    w07.23.06.0031.01_bit2(A, n, m, z, tol=tol, maxiter=maxiter)
  }

"w07.23.06.0026.01" <- function(A, n, m, z, tol=0, maxiter=2000, method="a"){  # checks out with maple.
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  stopifnot(m >= n)
  stopifnot(m >= 0)
  stopifnot(n >= 0)
  stopifnot(is.near_integer(n))
  stopifnot(is.near_integer(m))

  m <- round(m)
  n <- round(n)
  
  z <- z+0i

  bit1 <- w07.23.06.0026.01_bit1(A, n, m, z, tol=tol)  
  bit2 <- w07.23.06.0026.01_bit2(A, n, m, z, tol=tol, maxiter=maxiter) 
  bit3 <- switch(method,
                 a = w07.23.06.0026.01_bit3_a(A, n, m, z, tol=tol),
                 b = w07.23.06.0026.01_bit3_b(A, n, m, z, tol=tol),
                 c = w07.23.06.0026.01_bit3_c(A, n, m, z, tol=tol),
                 stop("method must be 'a' or 'b' or 'c'")
                 )
  return(bit1 + bit2 + bit3)
}

"w07.23.06.0026.01_bit1" <- function(A, n, m, z, tol=0){  # Checks with Maple
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  if(length(z)==0){return(z)}
  if(is.zero(n)){ return(0) }
  
# mult <- gamma(n)*gamma(A+m)*(-z)^(-A) / (gamma(m)*gamma(A+n))
  mult <- (0i-z)^(-A) * .f4(n,A+m,m,A+n)

  U <- c(A,1-m)
  L <- 1-n
  
  series <- z*0+1
  z[Mod(z) < 1] <- NA
  fac <- 1   # k=0
  temp <- fac
  for (k in seq_len(n-1)) {
    fac <- fac * (prod(U)/prod(L)) /(z*k)
    series <- temp + fac
    if (isgood(series-temp,tol)){
      return(series * mult)
    }
    temp <- series
    
    U <- U + 1
    L <- L + 1
  }
  return(series*mult)
}

"w07.23.06.0026.01_bit2" <- function(A, n, m, z, tol=0, maxiter = 2000){  # checks with Maple
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  mult <-
    (-1)^n * (-z)^(-A-m) * .f3(A+m,A,A+n) * .f3(A+m,m+1,m-n+1)
#      (gamma(A)*gamma(A+n)*factorial(m)*factorial(m-n))
  return(mult * genhypergeo(U=c(1,1,A+m),L=c(m+1,m-n+1), z=1/z, tol=tol, maxiter=maxiter))
}

"w07.23.06.0026.01_bit3_a" <- function(A, n, m, z, tol=0){ #"_a" means use psigamma, "_b" means use 6.3.5.
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  U <- c(A+n , 1-m+n)
#  mult <- (-1)^n * exp(lgamma(A+m)-lgamma(A)-lfactorial(m-n-1)) * (-z)^(-A-n)
   mult <- (-1)^n * .f3(A+m,A,m-n) * (-z)^(-A-n)
  
  fac <- 1/factorial(n)
  lmz <- log(0i-z)
  
  temp <- (lmz - psigamma(m-n-0) +psigamma(0+1) + psigamma(0+n+1) - psigamma(A+0+n)) * fac    #k=0
  series <- temp
  for(k in seq_len(m-n-1)){
    fac <- fac * prod(U) / (z * k * (k+n))
    series <-
      temp + fac * (lmz - psigamma(m-n-k) + psigamma(k+1) + psigamma(k+n+1) - psigamma(A+k+n))
    if(isgood(series-temp,tol)){
      return(series*mult)
    }
    temp <- series
    U <- U+1
  }
  return(series*mult)
}

"w07.23.06.0026.01_bit3_b" <- function(A, n, m, z, tol=0){ #"_a" means use psigamma, "_b" means use 6.3.5.
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  U <- c(A+n , 1-m+n)
# mult <- (-1)^n*exp(lgamma(A+m)-lgamma(A)-lfactorial(m-n-1)) * (-z)^(-A-n)
  mult <- (-1)^n*.f3(A+m,A,m-n) * (-z)^(-A-n)
  
  fac <- 1/factorial(n)
  lmz <- log(0i-z)

  p1 <- psigamma(m-n)
  p2 <- psigamma(1)
  p3 <- psigamma(n+1)
  p4 <- psigamma(A+n)
  
  temp <- (lmz - p1 + p2 + p3 - p4) * fac
  series <- temp
  for(k in seq_len(m-n-1)){
    fac <- fac * prod(U) / (z * k * (k+n) )

    p1 <- p1 - 1/(m-n-k)
    p2 <- p2 + 1/k
    p3 <- p3 + 1/(k+n)
    p4 <- p4 + 1/(A+k+n-1)
      
    series <-
      temp + fac * (lmz - p1 + p2 + p3 - p4)
    if(isgood(series-temp,tol)){
      return(series*mult)
    }
    temp <- series
    U <- U+1
  }
  return(series*mult)
}

"w07.23.06.0026.01_bit3_c" <- function(A, n, m, z, tol=0){ #"_a" means
  # use psigamma, "_b" means use 6.3.5; here "_c" means use a totally
  # dull, slow, direct (but clearly correct) summation, for the
  # purposes of debugging.
  
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  poch <- function(x,j){ prod(x + (seq_len(j)-1)) }
  
#  mult <- ((-1)^n*gamma(A+m)/(gamma(A)*factorial(m-n-1)))*(-z)^(-A-n)
   mult <- ((-1)^n*.f3(A+m,A,m-n))*(-z)^(-A-n)
  out <- 0

  for(k in 0:(m-n-1)){
    out <- out +
      (
       (poch(A+n,k) * poch(1-m+n,k))/(factorial(k)*factorial(k+n))
       ) *
         (log(-z) - psigamma(m-n-k)+psigamma(k+1)+psigamma(k+n+1)-psigamma(A+k+n))*z^(-k)
  }
  return(out * mult)
}

"genhypergeo_contfrac_single" <- function(U, L, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  f <- function(k){prod(U+k)/prod(k+c(1,L))}
  alpha <- z*sapply(seq_len(maxiter), f)
  1+z*prod(U)/(prod(L)*(1+GCF(a = -alpha , b = 1+alpha, tol=tol)))
}

"genhypergeo_contfrac" <- function(U, L, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  attr <- attributes(z)
  f <- function(z){genhypergeo_contfrac_single(U, L, z=z, tol=tol, maxiter=maxiter)}
  out <- sapply(z,f)
  attributes(out) <- attr
  return(out)
}

"hypergeo_contfrac" <- function(A, B, C, z, tol=0, maxiter=2000){
  if(!is.null(getOption("showHGcalls"))){print(match.call())}
  attr <- attributes(z)
  f <- function(z){genhypergeo_contfrac_single(U=c(A, B), L=C, z=z, tol=tol, maxiter=maxiter)}
  out <- sapply(z,f)
  attributes(out) <- attr
  return(out)
}


"hypergeo_residue_general" <- function(A, B, C, z, r, O=z, tol=0, maxiter=2000){
    if(!is.null(getOption("showHGcalls"))){print(match.call())}
    stopifnot(length(z)==1)
    residue(f=function(z){hypergeo(A,B,C,z,tol=tol,maxiter=maxiter)}, z0=z, r=0.15, O=O) # NB: residue() is defined in the elliptic package
}

"hypergeo_residue_close_to_crit_single" <- function(A, B, C, z, strategy='A', tol=0, maxiter=2000){

    if(!is.null(getOption("showHGcalls"))){print(match.call())}

    jj <- crit()
    c1 <- jj[1]
    c2 <- jj[2]
    
    if(
        (abs(z-c1) <= 0.1)   &
        (abs(z-c2) <= 0.1)
        ) {stop("this cannot happen")}

    
    stopifnot(
        (abs(z-c1) <= 0.1)   |
        (abs(z-c2) <= 0.1)
        )

    if(abs(z-c1) <= 0.1){
        crit <- c1
    } else {
        crit <- c2
    }

    O <- switch(
        strategy,
        A = crit,
        B = z,
        stop('strategy must be A or B')
        )
        
    hypergeo_residue_general(A=A,B=B,C=C, z=z, r=0.15, O=O, tol=tol, maxiter=maxiter)

}

"hypergeo_residue_close_to_crit_multiple" <- function(A, B, C, z, strategy='A', tol=0, maxiter=2000){
    if(!is.null(getOption("showHGcalls"))){print(match.call())}
    sapply(z, function(z){
        hypergeo_residue_close_to_crit_single(A,B,C,z,strategy=strategy,tol=tol,maxiter=maxiter)
    } )
}    

"lpham" <- function(x,n){lgamma(x+n)-lgamma(x)}

"buhring_eqn11" <- function(n,S,A,B,C,z0=1/2){  #NB no z
    stopifnot(length(z0)==1)
    if(length(n)>1) {return(sapply(n,function(nn){buhring_eqn11(n=nn,S,A,B,C,z0=z0)}))}
    return(
        exp(
            +lpham(S,n)
            +lpham(1+S-C,n)
            -lpham(1+2*S-A-B,n)
            -lfactorial(n)
            ) * hypergeo(-n, A+B-2*S-n, C-S-n, z=z0)
        )
}

"buhring_eqn12" <- function(n,S,A,B,C,z0=1/2){
    stopifnot(length(z0)==1)
    if(length(n)>1) {return(sapply(n,function(nn){buhring_eqn12(n=nn,S,A,B,C,z0=z0)}))}
    
    return(
        (-1)^n*
        exp(
            +lpham(S,n)
            +lpham(S+C-A-B,n)
            -lpham(1+2*S-A-B,n)
            -lfactorial(n)
            ) * hypergeo(-n,A+B-2*S-n, 1+A+B-S-C-n, z=1-z0)
        )
}
    
"buhring_eqn5_factors" <- function(A,B,C,z,z0=1/2){
    c(
        exp(
            +complex_gamma(C,log=TRUE)
            +complex_gamma(B-A,log=TRUE)
            -complex_gamma(B,log=TRUE)
            -complex_gamma(C-A,log=TRUE)
            -A*log(z0-z)
            ),
        exp(
            +complex_gamma(C,log=TRUE)
            +complex_gamma(A-B,log=TRUE)
            -complex_gamma(A,log=TRUE)
            -complex_gamma(C-B,log=TRUE)
            -B*log(z0-z)
            )
        )    
}

"buhring_eqn5_series" <- function(S,A,B,C,z,z0=1/2,use11=FALSE,tol=0,maxiter=2000){  # sum
    if(!is.null(getOption("showHGcalls"))){print(match.call())}
    if(length(z)==0){return(z)}
    
    if(use11){
        f <- buhring_eqn11
    } else {
        f <- buhring_eqn12
    }
    temp <- 1
    n <- 1
    while(n < maxiter){
        out <- temp + f(n,S=S,A=A,B=B,C=C,z0=z0)/(z-z0)^n
        if(isgood(out-temp,tol)){return(out)}
        temp <- out
        n <- n+1
    }
    warning("series not converged")
    return(out)
}

"hypergeo_buhring" <- function(A,B,C,z,z0=1/2,tol=0,maxiter=2000,use11=TRUE){
    jj <- buhring_eqn5_factors(A,B,C,z,z0)
    return(
        jj[1]*buhring_eqn5_series(S=A,A,B,C,z,z0=1/2,use11=use11,tol=tol,maxiter=maxiter)+
        jj[2]*buhring_eqn5_series(S=B,A,B,C,z,z0=1/2,use11=use11,tol=tol,maxiter=maxiter)
        )
}

"shanks" <- function(Last,This,Next){
    if(identical(Next,This)){return(Next)}
    num <- Next*Last - This^2
    den <- Next-2*This+Last

    if(den==0){
        return(Next)
    } else {
        return(num/den)
    }
}
 
"genhypergeo_shanks" <-
function (U, L, z,  maxiter=20){
    if(!is.null(getOption("showHGcalls"))){print(match.call())}

    fac <- 1
    temp <- fac
    
    if(maxiter==0){ return(z*0+fac) }
    
    Last <- 0
    This <- 1
    Next <- 2
    
    Shanks <- shanks(Last,This,Next)
    
    for (n in seq_len(maxiter)) {
        
        fac.old <- fac
        fac <- fac * (prod(U)/prod(L)) * (z/n)
        fac.new <- fac
        
        series <- temp + fac
        ## following three lines a "conveyor belt" Next -> This -> Last
        Last <- This
        This <- Next
        Next <- series
        
        Shanks.old <- Shanks
        Shanks <- shanks(Last,This,Next)
        
        temp <- series
        U <- U + 1      
        L <- L + 1
    }
    return(series)
} 

"hypergeo_shanks" <- function (A, B, C, z, maxiter = 20){
    genhypergeo_shanks(U=c(A,B), L=C, z=z,maxiter=maxiter)
}

"hypergeo_gosper" <- function(A, B, C, z, tol=0, maxiter=2000){
    d <- 0
    e <- 1
    f <- 0

    for(k in 0:maxiter){
        dnew <- (k+A)*(k+B)*z*(e-(k+C-B-A)*d*z/(1-z))  /(4*(k+1)*(k+C/2)*(k+(C+1)/2))
        enew <- (k+A)*(k+B)*z*(A*B*d*z/(1-z) + (k+C)*e)/(4*(k+1)*(k+C/2)*(k+(C+1)/2))
        fnew <- f-d*(k*((C-B-A)*z+k*(z-2)-C)-A*B*z)    /(2*      (k+C/2)*(1-z)      )+e
        
        if(isgood(f-fnew,tol)){return(f)}
        d <- dnew
        e <- enew
        f <- fnew
    }
    warning("not converged")
    return(f)
}
    
    
