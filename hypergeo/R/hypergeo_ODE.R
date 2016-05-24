"to_real" <- function(o){
    out <- c(rbind(Re(o),Im(o)))
    if(!is.null(names(o))){
        names(out) <-   # pathetic kludge
            apply(expand.grid(c("_real","_imag"),names(o))[,2:1],1,paste,collapse="")
    } else {
        names(out) <- NULL
    }
    return(out)
}

"to_complex" <- function(p){
    if(is.vector(p)){
        jj <- Recall(t(p))
        out <- c(jj)
        names(out) <- colnames(jj)
        return(out)
    }  # if not a vector, assumed to be a matrix
    out <- (
        p[,seq(from=1,by=2,to=ncol(p)),drop=FALSE] +
            1i*p[,seq(from=2,by=2,to=ncol(p)),drop=FALSE]
        )
    f <- function(string){sub("_real","",string)}
    colnames(out) <- sapply(colnames(out),f)
    return(out)
}

"complex_ode" <- function(y, times, func, parms=NA, method=NULL, u, udash, ...){
    out <-
        ode(y=to_real(y), times=times, func=func, parms=to_real(parms), method, u=u, udash=udash, ...)
    
    out <- cbind(z=u(out[,1]),to_complex(out[,-1]))
    class(out) <- c("deSolve", "matrix")
    return(out)
}

hypergeo_press <- function(A,B,C,z, ...){    # Press et al, 5.14
    
    if(Re(z)<=0){
        startz <- -0.5
    } else if( (Re(z)<=0.5)){
        startz <- 0.5
    } else if(Im(z)>=0){
        startz <- 0.5i
    } else if(Im(z)<0){
        startz <- -0.5i
    }
    
    initial_value <- hypergeo(A,B,C,z=startz)
    initial_deriv <- (A*B)/C*hypergeo(A+1,B+1,C+1,z=startz)   # 15.2.1
    
    complex_ode(y     = c(F=initial_value, Fdash=initial_deriv),
                times = seq(0,1,by=0.05),
                func  = hypergeo_func,
                parms = c(A=A, B=B, C=C)+0i,
                u     = function(u){startz + (z-startz)*u}, # path
                udash = function(u){z-startz},              # derivative of path
                ...)
}

"hypergeo_func" <- function(Time, State, Pars, u, udash) {
    with(as.list(c(to_complex(State), to_complex(Pars))), {
             z <- u(Time)
             dz <- udash(Time)
             
             ## 'meat' of function: AMS-55 15.5.1; w -> F
             dF     <- dz * Fdash
             dFdash <- dz * (A*B*F -(C-(A+B+1)*z)*Fdash)/(z*(1-z)) 
             
             ## Now coerce back to real:
             out <- to_real(c(dF,dFdash))
             names(out) <- names(State)
             return(list(out))
         })
}

f15.5.1 <- function(A,B,C,z,startz,u,udash,give=FALSE, ...){   # solves the ODE, 15.5.1, directly.

    out <-
    complex_ode(y     = c(F=hypergeo(A,B,C,startz), Fdash=hypergeo(A+1,B+1,C+1,startz)*A*B/C),
                times = seq(0,1,by=0.1),
                func  = hypergeo_func,
                parms = c(A=A, B=B, C=C)+0i,
                u     = u,
                udash = udash,
                ...)
    if(give){
       return(out)
    } else {
       return(unname(out[11,2]))
    }
}

"semicircle" <- function(t,z0,z1,clockwise=TRUE){
    if(clockwise){m <- -1} else {m <- 1}
    center <- (z0+z1)/2
    center + (z0-center)*exp(1i*t*pi*m)
}

"semidash" <- function(t,z0,z1,clockwise=TRUE){
    if(clockwise){m <- -1} else {m <- 1}
    center <- (z0+z1)/2
    (z0-center)*(1i*pi*m)*exp(1i*t*pi*m)
}

"straight" <- function(t,z0,z1){ z0 + t*(z1-z0) }

"straightdash" <- function(t,z0,z1){ (z1-z0) }
