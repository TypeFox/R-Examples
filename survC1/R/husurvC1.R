#-- R wrapper (husurvC1v1.R - specify package name)--

#=============================================
# cf.unoCWv1.R
#=============================================
unoCW <- function(X, D, W, Wstar, RS, RSstar, REX, Chat) {
	#----------------------
	#-- X: observed time
	#-- D: event indicator (1:event, 0:no)
	#-- W: weight (original)
	#-- Wstar: weight (perturbed)
	#-- RS: risk score (original)
	#-- RSstar: risk score (perturbed)
	#-- REX: random number by perturbation
	#----------------------
	    out <- .Fortran("unoCW",
                n=as.integer(length(X)),
                time=as.integer(X*10000),
                status=as.integer(D),
                weight=as.double(W),
                wtstar=as.double(Wstar),
                rs=as.integer(RS*100000),
                rsstar=as.integer(RSstar*100000),
                rex=as.double(REX),
                Dhat=as.double(Chat),
                Wa=as.double(0),
                Wb=as.double(0),
                Wg=as.double(0),
                WK1=as.double(0),
                WK3A=as.double(0),
                WK3B=as.double(0),
                USEP=as.double(0),
                USEPX=as.double(0),
                PACKAGE="survC1")
          return(out)
}

#=============================================
# cf.unoU2pv1.R
#=============================================
unoU2P <- function(A, B) {
	#----------------------
	#-- A: matrix (nxp)
	#-- B: vector (nx1)
	#-- retrun -- sum_{i<j}(A[i,]+A[j,])*B[i]*B[j]/2/{n choose 2}
	#----------------------
	    out <- .Fortran("unoU2P",
                n=as.integer(nrow(A)),
                p=as.integer(ncol(A)),
                A=as.double(A),
                B=as.double(B),
                OUT=rep(0,ncol(A)),
                UP=as.integer(0),
                PACKAGE="survC1")
        return(out$OUT)
   }

#=============================================
# cf.conc4.R
#=============================================
conc <- function(X,D,W,R) {
	#----------------------
	#-- X: observed time
	#-- D: event indicator (1:event, 0:no)
	#-- W: weight (IPCW)
	#-- R: risk score
	#----------------------
	    out <- .Fortran("conc",
                n=as.integer(length(X)),
                time=as.integer(X*1000),
                status=as.integer(D),
                rs=as.integer(R*100000),
                weight=as.double(W),
                hfwt=as.double(W/2),
                WK1=as.double(0),
                CSTAT=as.double(0),
                PACKAGE="survC1")
        return(out$CSTAT)
   }


