
testqz <- function(A,B,z) {
    VL <- z$Q
    VR <- z$Z

    tol <- 100 * sqrt(.Machine$double.eps)
    ret <- logical(4)

    if( is.complex(A) || is.complex(B) ) {
        tmp <- A - VL %*% z$S %*% Conj(t(VR))
        ret[1] <- all(abs(tmp)<=tol)
        tmp <- B - VL %*% z$T %*% Conj(t(VR))    
        ret[2] <- all(abs(tmp)<=tol)
        
        ret[3] <- all(abs(z$Q %*% Conj(t(z$Q))-diag(1.0,nrow(z$Q),ncol(z$Q)))<=tol)
        ret[4] <- all(abs(z$Z %*% Conj(t(z$Z))-diag(1.0,nrow(z$Z),ncol(z$Z)))<=tol)

    }
    else {
        tmp <- A - VL %*% z$S %*% t(VR)
        ret[1] <- all(abs(tmp)<=tol)
        tmp <- B - VL %*% z$T %*% t(VR)
        ret[2] <- all(abs(tmp)<=tol)
        ret[3] <- all(abs(z$Q %*% t(z$Q)-diag(1.0,nrow(z$Q),ncol(z$Q)))<=tol)
        ret[4] <- all(abs(z$Z %*% t(z$Z)-diag(1.0,nrow(z$Z),ncol(z$Z)))<=tol)
    }

    ret
}
