############################################################
## Functions to generate a MA process
############################################################


############################## Generating the filter 'a'
#' Generate the Filter of a multivariate MA process 
#' 
#' @param d.ts dimension of the (output) time series
#' @param d.n    dimension of the noise that is filtered
#' @param MA.len Length of the filter. Set to 3 by default.
#' @param ma.scale   scaling factor of each lag matrix. See details.
#' @param a.smooth.coef A coefficient to shrink coefficients of filter. Set to 0 by default.
#' @param seed The random seed used to generate the filter. Set to 1 by default.
#' 
#' @export
#'
#' @return A \code{d.ts x d.n x MA.len} array
#'
#' @section Details:
#'
#' Generates a filter (i.e. a \code{d.ts x d.n x MA.len} array) for a moving
#' average process. The entries of the filter are generate randomly, but can be
#' reproduced by specifying the random seed \code{seed}. 
#' 
#' The \code{ma.scale} parameter should be a vector of length \code{MA.len},
#' and corresponds to a scaling factor applied to each lag of the filter of the
#' MA process that is generated.
#'
#' @examples
#' ma.scale1=c(-1.4,2.3,-2)
#' a1=Generate_filterMA(10, 10, MA.len=3, ma.scale=ma.scale1, seed=10)
#' str(a1)
#' rm(a1)

Generate_filterMA = function(d.ts, d.n, MA.len=3, ma.scale=rep(1, MA.len), a.smooth.coef=0,  seed=1){
    if( MA.len != length(ma.scale) ) stop("ma.scale must be of length equal to MA.len")
    set.seed(seed)
    a<-array(stats::rnorm(d.ts*d.n*MA.len, mean=1, sd=0.5), c(d.ts, d.n, MA.len) ) 
    a<- sapply(1:MA.len, function(j){ diag( (1:d.ts)^(-a.smooth.coef) ) %*% a[,,j] }, simplify="array")
    rm(list=".Random.seed", envir=globalenv())
#
    for(j in 1:dim(a)[3]){
            a[,,j]=ma.scale[j]*a[,,j]
    }
    a
} # Generate_filterMA


#' Get the square root of the covariance matrix associated to a noise type
#' 
#' @param noise.type the type of noise that is driving the MA process. See Details section.
#' @param d.n    dimension of the noise that is filtered

Get_noise_sd<-function(noise.type, d.n){
    if(noise.type=="white-noise"){ 
        noise.sd <- diag(rep(1,d.n), d.n, d.n)
    } else if(noise.type=="wiener") {
        noise.sd <- diag( ( ((1:d.n) - 0.5 )*pi)^(-1), d.n, d.n )  # Wiener process type noise
    } else if (substring(noise.type, 1, 7) == "student") {
        df.stud=eval(parse(text=substring(noise.type, 8)))
        if(df.stud <= 0){
            stop("wrong df for student noise")
        }
        noise.sd<- diag(rep(sqrt((df.stud-2)/df.stud),d.n), d.n, d.n)
    } else { 
        stop("This noise type is unknown"); 
    } 
    return(noise.sd)
} # Get_noise_sd




#' Simulate a new Moving Average (MA) vector time series and return the time series 
#' 
#' 
#' @param a Array, returned by \code{Generate_filterMA}, containing the filter of the MA process
#' @param T.len Numeric, the length of the time series to generate
#' @param noise.type the type of noise that is driving the MA process. See Details section.
#' @param DEBUG Logical, for outputting information on the progress of the function
#' 
#' @section Details:
#'
#' The function simulates a moving average process of dimension
#' \code{dim(a)[1]}, defined by \deqn{X[t,] =  a[,,1]  * epsilon[,t-1] + a[,,2]
#' * epsilon[,t-2] + ... + a[,,dim(a)[3]] * epsilon[t-dim(a)[3]]  } 
#' 
#' \code{noise.type} specifies the nature and internal correlation of the noise
#' that is driving the MA process. It can take the values 
#' \describe{
#' \item{\code{white-noise}}{ the noise is Gaussian with covariance matrix identity}
#' \item{\code{white-noise}}{the noise is Gaussian with diagonal covariance
#' matrix, whose j-th diagonal entry is \eqn{((j - 0.5 )*pi)^(-1)}}
#' \item{\code{studentk}}{the coordinates  of the noise are independent and
#' have a student t distribution with 'k' degrees of freedom, standardized to
#' have variance 1}
#' }
#' 
#' @return A \code{T.len x dim(a)[1]} matrix, where each column corresponds to a
#' coordinate of the vector time series
#' 
#' @export
#' 
#' @examples 
#' ma.scale1=c(-1.4,2.3,-2)
#' a1=Generate_filterMA(6, 6, MA.len=3, ma.scale=ma.scale1)
#' X=Simulate_new_MA(a1, T.len=512, noise.type='wiener')
#' plot.ts(X)

Simulate_new_MA <- function(a, T.len, noise.type, DEBUG=FALSE){
#
  tmp <- dim(a)
  d.ts = tmp[1] 
  d.n = tmp[2]
  MA.length = tmp[3]
  Epsilon.length= T.len + MA.length - 1
  noise.sd=Get_noise_sd(noise.type, d.n)
  rm(tmp) # clean temporary variable
# 
  #simulate epsilon
  if(noise.type %in% c("wiener","white-noise", "exponential")){
    epsilon<- noise.sd %*% array( stats::rnorm(d.n*Epsilon.length, sd=1) , c(d.n, Epsilon.length)) #(gaussian noise)
  } else if (substring(noise.type, 1, 7) == "student"){
      df.stud=eval(parse(text=substring(noise.type, 8)))
      if(df.stud <= 0){
          stop("wrong df for student noise")
      }
    epsilon<- noise.sd %*% array( stats::rt(d.n*Epsilon.length, df=df.stud) , c(d.n, Epsilon.length)) #(Student noise)
  } else {
      stop("unknown noise type")
  }
  
  X<- array(numeric(1), c(T.len, d.ts))
  
  
  #this is quite slow
  for(j in 1:T.len) # careful: t=j-1
  {
    if(DEBUG){     if((j %% floor(T.len/10))==0){print(paste(j, "of", T.len, "in computation of X"))} }
    for(s in 1:MA.length)
    {
      X[j,] <- X[j,] + a[,, s ] %*% epsilon[,j - (s-1) + (MA.length-1)]
    }
  }
  return(X)
} # Simulate_new_MA




#' 'Spectral density operator of a MA vector process' Object
#' 
#' @param a the filter of the moving average
#' @param nfreq the number of frequencies between 0 and pi at which the
#' spectral density has to be computed 
#' @param noise.type the type of noise that is driving the MA process. See
#' \code{\link{Simulate_new_MA}} 
#' 
#' @examples
#' ma.scale1=c(-1.4,2.3,-2)
#' a1=Generate_filterMA(6, 6, MA.len=3, ma.scale=ma.scale1)
#' a1.spec=SpecMA(a1, nfreq=512, noise.type='wiener')
#' plot(a1.spec)
#' rm(a1, a1.spec)
#' @export
#' @importFrom stats fft kernapply kernel mvfft p.adjust pchisq rnorm rt

SpecMA<-function(a, nfreq = 2^9, noise.type){
  T.len  <- 2*nfreq
  d.ts<-dim(a)[1]
  d.n<-dim(a)[2]
  noise.sd=Get_noise_sd(noise.type, d.n)
  noise.cov=noise.sd %*% noise.sd
  
  a.padded<-array(0, c(d.ts, d.n, T.len)) 
  a.padded[,,1:dim(a)[3] ] <- a
  if(d.ts == 1) {
    a.tilde <- array(stats::fft(a.padded), c(1,d.n,T.len))
    specdK<- sapply(1:T.len, function(j){ a.tilde[,,j]  %*% noise.cov %*% Conj((a.tilde[,,j])) }, simplify="array")/(2*pi)
    
  } else {
    a.tilde<-array(NA, c(d.ts, d.n, T.len))
    
    for(i in 1L:d.ts){
      a.tilde[i,,]<- t(stats::mvfft(t(a.padded[i,,])))
    }
    specdK<- sapply(1:T.len, function(j){ a.tilde[,,j]  %*% noise.cov %*% Conj(t(a.tilde[,,j])) }, simplify="array")/(2*pi)
  }
  
  spec.density <- list(a=a, omega=seq(0, pi, len=nfreq + 1), spec=array(specdK, c(d.ts, d.ts, T.len))[,,1:(floor(nfreq)+1)])
  class(spec.density) <- "SpecMA"
  return(spec.density)
} 

#' Plotting method for object inheriting from class \code{SpecMA}
#' @param x A object of the class \code{SpecMA}
#' @param ... additional parameters to be passed to plot()
#' @export
plot.SpecMA <- function(x, ...){
    graphics::plot(x$omega, apply(x$spec, 3, function(x) sum(diag(Re(x)))), type = 'l', main = 'Trace
      of spectral density of MA process', xlab = 'omega', ylab = 'intensity', ...)
}



############################################################
## (END) Functions to generate a MA process
############################################################


#' The Epanechnikov weight function, with support in \eqn{[-1,1]}
#'
#' @export
#' @param x argument at which the function is evaluated
#'
#'

Epanechnikov_kernel<- function(x){ 
    ifelse( abs(x)<= 1,    (3/4)*(1-x^2), 0)	
}


#' Compute Spectral Density of Functional Time Series
#' 
#' This function estimates the spectral density operator of a Functional Time Series (FTS) 
#' 
#' @param X A \eqn{T \times nbasis} matrix of containing the coordinates of the FTS
#' expressed in a basis. Each row corresponds to a time point, and each column
#' corresponds to the coefficient of the corresponding basis function of the  FTS.
#' @param W The weight function used to smooth the periodogram operator. Set by
#' default to be the Epanechnikov kernel
#' @param B.T The bandwidth of frequencies over which the periodogram operator
#' is smoothed. If \code{B.T=0},  the periodogram operator is returned.
#' @param only.diag A logical variable to choose if the function only computes
#' the marginal spectral density of each basis coordinate
#' (\code{only.diag=TRUE}). \code{only.diag=FALSE} by default, the full spectral
#' density operator is computed .
#' @param trace A logical variable to choose if only the trace of the spectral
#' density operator is computed. \code{trace=FALSE} by default.
#' @param demean A logical variable to choose if the FTS is centered before
#' computing its spectral density operator.
#' @param subgrid A logical variable to choose if the spectral density operator
#' is only returned for a subgrid of the Fourier frequencies, which can be
#' useful in large datasets to reduce memory usage. \code{subgrid=FALSE} by
#' default.
#' @param subgrid.density Only used if \code{subgrid=TRUE}. Specifies the
#' approximate number of frequencies within the bandwidth over which the
#' periodogram operator is smoothed. 
#' @param subgrid.density.relative.to.bandwidth logical parameter to specify if
#' \code{subgrid.density} is specified relative to the bandwidth parameter \code{B.T}
#' @param verbose A  variable to show the progress of the computations. By
#' default, \code{verbose=0}.
#' 
#' @export
#'
#' @return A list containing the following elements:
#' \describe{
#'  \item{spec}{The estimated spectral density operator. The first dimension corresponds to the different frequencies over which the spectral density operators are estimated.}
#'  \item{omega}{The frequencies over which the spectral density is estimated.}
#'  \item{m}{The number of Fourier frequencies over which the periodogram operator was smoothed.}
#'  \item{bw}{The equivalent Bandwidth used in the weight function W(), as defined in Bloomfield (1976, p.201).}
#'  \item{weight}{The weight function used to smooth the periodogram operator.}
#'  \item{kappa.square}{The L2 norm of the weight function W.}
#'  }
#'
#' @examples 
#' ma.scale1=c(-1.4,2.3,-2)
#' a1=Generate_filterMA(10, 10, MA.len=3, ma.scale=ma.scale1)
#' X=Simulate_new_MA(a1, T.len=512, noise.type='wiener')
#' ans=Spec(X, trace=FALSE, only.diag=FALSE)
#' plot(ans)
#' plot(Spec(X, trace=FALSE, only.diag=FALSE, subgrid=TRUE, subgrid.density=10,
#' subgrid.density.relative.to.bandwidth=FALSE))
#' rm(ans)
#'
#' @references \link{spec.pgram} function of R.
#' @references \cite{ Bloomfield, P. (1976) "Fourier Analysis of Time Series:
#' An Introduction", Wiley.}
#' @references \cite{Panaretos, V. M. and Tavakoli, S., "Fourier Analysis of Functional Time Series", Ann. Statist. Volume 41, Number 2 (2013), 568-603.}
#' @importFrom utils setTxtProgressBar txtProgressBar

Spec<-function(X, W=Epanechnikov_kernel, B.T=(dim(X)[1])^(-1/5), only.diag=FALSE, trace=FALSE, demean=TRUE, subgrid=FALSE, subgrid.density=10, verbose=0, subgrid.density.relative.to.bandwidth=TRUE){
    X<-as.matrix(X)
    T<-dim(X)[1]
    res<-dim(X)[2]
    M=1e4; xtmp=(0:M)/M 
    sigma.W=sqrt(2*sum((xtmp^2)*W(xtmp))/M)  ## standard deviation of the weight function 
    kappa.square=2*sum(W(xtmp)^2)/M ## L2 norm of Weight function 
    rm(M, xtmp) 
    m=floor( T*B.T/(2*pi) ) 
    bw=B.T*sigma.W

    sample.spec=list(spec=NA, omega=NA, T.len=T,  m=m, B.T=B.T, bw=bw,
                     W=W, kappa.square=kappa.square, sigma.W=sigma.W ,
                     subgrid.density.relative.to.bandwidth=subgrid.density.relative.to.bandwidth)
    class(sample.spec) <- 'SampleSpec' 

    omega<-2*pi*(0:floor(T/2))/T
    ###  allow subgrid.density only between 1 and 2*m if subgrid.density.relative.to.bandwidth=TRUE
    ###  allow subgrid.density only between 1 and 100 if subgrid.density.relative.to.bandwidth=FALSE
    omega.i.subseq<-1:length(omega)
    if(subgrid){ # take a subgrid of omega to lower storage size
        ## if omega is too fine, take a subsequence of it for faster computations
        if(subgrid.density.relative.to.bandwidth){ ### the 'subgrid.density' is specified relative to the bandwidth parameter B.T
            if((subgrid.density < 1) || (subgrid.density > 2*m)){
                stop(paste0("subgrid.density must be between 1 and 2*m (", 2*m, ") since subgrid.density.relative.to.bandwidth=TRUE"))
            }
            min.length.omega<- floor(subgrid.density*pi/(2*B.T)) ## you want to have a frequency grid dense enough compared to the bandwidth ## value of subgrid.density between 1 and 2*m 
        } else {
            if((subgrid.density < 1) || (subgrid.density > 100)){
                stop(paste0("subgrid.density must be between 1 and 100 since subgrid.density.relative.to.bandwidth=FALSE"))
            }
            min.length.omega<- floor( (subgrid.density - 1)/99*(length(omega) -
                                                                pi/(2*B.T)) +
                                     pi/(2*B.T) ) 
            ## linear interpolation between min.length=pi/(2 B.T) for
            ## subgrid.density=1 and min.length=length(omega) for subgrid.density=100 
        }
        if(length(omega)> min.length.omega){
            omega.i.subseq<-floor(seq(1, to=length(omega),
                                      by=round(length(omega)/min.length.omega)))
            # this yields uniform frequency grids of size asymptotically min.length.omega
            omega.i.subseq[length(omega.i.subseq)]=length(omega) ## We want to have the frequency pi in the subgrid
        }
    }
    if(demean){ ## Recenter each column of X
        X<-scale(X, center=TRUE, scale=FALSE)
    }
    if(m == 0) {
        kern<-stats::kernel(1)
    }  else {
        coef<-sapply((-m:m)/m, W) / m
        coef<-(coef/sum(coef))[-c(1:(m))]
        kern<-stats::kernel(coef)
    }
    if(trace){
        Pgram<-Mod(stats::mvfft(X))^2/(2*pi*T)
        spec.est<-stats::kernapply(rowSums(Pgram), kern, circular = TRUE)
        sample.spec$spec=spec.est[omega.i.subseq]
        sample.spec$omega=omega[omega.i.subseq]
        return(sample.spec)
    } else
        if(only.diag){
            Pgram<-Mod(stats::mvfft(X))^2/(2*pi*T)
            spec.est<-as.matrix(stats::kernapply(Pgram, kern, circular = TRUE))
            #
            sample.spec$spec=as.matrix(spec.est[omega.i.subseq,])
            sample.spec$omega=omega[omega.i.subseq]
            return(sample.spec)
        }
        else{
            Xfft<-stats::mvfft(X)/sqrt(2*pi*T)
            Pgram<-array(NA, c(length(omega.i.subseq), res,res))
            if(verbose==1){ 
                pb<-txtProgressBar(0, res , style=2)
            }
            for(i in 1L:res) {
                if(verbose==1){ 
                    cat("- Computation of Spectral Density\n")
                    setTxtProgressBar(pb, i)
                } else if(verbose==2){
                    cat("- Computation of Spectral Density: ")
                    cat("Smooth Pgram step ", i, " of ", res, "\n", sep='') 
                }
                for(j in i:res){
                    #         cat(paste(j, "/", res, "- "));
                    tmp<- Xfft[,i] * Conj(Xfft[,j])
                    Pgram[,i,j]<- stats::kernapply(tmp, kern, circular = TRUE)[omega.i.subseq] 
                    if(j != i) { Pgram[,j,i] = Conj(Pgram[,i,j]) }
                    rm(tmp)
                }
            }
            if(verbose==1){ 
                close(pb)
            }
            sample.spec$spec=Pgram
            sample.spec$omega=omega[omega.i.subseq]
            return(sample.spec)
        }
}


#' Plotting method for object inheriting from class \code{SampleSpec}
#' 
#' @param x An object of the class \code{SampleSpec}
#' @param ... additional parameters to be passed to plot()
#' @export
#' @importFrom graphics image lines matlines matplot mtext par plot rug
plot.SampleSpec  <- function(x,  ...){
    if(!exists('bw.col')) bw.col='blue'
    if( is.array(x$spec) ){
        if( length(dim(x$spec)) == 2L ){
    trace.spec=rowSums(x$spec)
        } else if (length(dim(x$spec)) == 3L){
       trace.spec  <-   apply(x$spec, 1, function(x) Re(sum(diag(x))) ) 
        }
    } else {
       trace.spec <- x$spec
    }
    graphics::plot(x$omega, trace.spec, type='l', main='Trace spectrum',  ylab =
         'intensity', ...)
    graphics::lines(c(pi-x$bw, pi), rep(max(trace.spec), 2), col=bw.col)
    graphics::rug(x$omega, ticksize=.01 )
}






#' Test if two spectral density operators at  some  fixed frequency are equal.
#'
#' A test for the null hypothesis that two spectral density operators (at the
#' same frequency \eqn{\omega}) are equal, using a pseudo-AIC criterion for the
#' choice of the truncation parameter. (used in
#' \code{\link{Spec_compare_localize_freq}})
#'
#' @param spec1,spec2 The two sample spectral densities (at the same frequency \eqn{\omega}) to be compared.
#' @param is.pi.multiple A logical variable, to specify if \eqn{\omega = {0, \pi}} or not. 
#' @param m The number of Fourier frequencies over which the periodogram
#' operator was smoothed.
#' @param kappa.square the L2-norm of the weight function used to estimate the
#' spectral density operator
#' @param autok A variable used to specify if (and which) pseudo-AIC criterion
#' is used to select the truncation parameter \eqn{K}. 
#' @param K.fixed The value of K used if \code{autok=0}.
#' 
#'
#' @examples 
#' 
#' ma.scale2=ma.scale1=c(-1.4,2.3,-2)
#' ma.scale2[3] = ma.scale1[3]+.3
#' a1=Generate_filterMA(10, 10, MA.len=3, ma.scale=ma.scale1)
#' a2=Generate_filterMA(10, 10, MA.len=3, ma.scale=ma.scale2)
#' X=Simulate_new_MA(a1, T.len=512, noise.type='wiener')
#' Y=Simulate_new_MA(a2, T.len=512, noise.type='wiener') 
#' spec.X = Spec(X)
#' spec.Y = Spec(Y)
#' Spec_compare_fixed_freq(spec.X$spec[1,,], spec.Y$spec[1,,],
#' is.pi.multiple=TRUE, spec.X$m, spec.X$kappa.square)
#' 
#' @seealso \code{\link{Spec_compare_localize_freq}}
#' @export
#'
#' @references \cite{Tavakoli, Shahin and Panaretos, Victor M. "Detecting and Localizing Differences in Functional Time Series Dynamics: A Case Study in Molecular Biophysics", 2014, under revision}
#' @references \cite{Panaretos, Victor M., David Kraus, and John H. Maddocks.
#' "Second-order comparison of Gaussian random functions and the geometry of
#' DNA minicircles." Journal of the American Statistical Association 105.490
#' (2010): 670-682.}

Spec_compare_fixed_freq=function(spec1, spec2, is.pi.multiple, m, kappa.square, autok=2, K.fixed=NA){
    n=dim(spec1)[2]
    #
    ##### To bypass some problems with svd: "Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'"
    # svd1 <- svd(spec1)
    svd1=tryCatch({svd(spec1)},
        error = function(err){
            svd.tmp=(eigen(spec1, symmetric=TRUE))
            names(svd.tmp)=c("d", "u")
            svd.tmp$d=pmax(svd.tmp$d, 0) ## correct the eigenvalues if they are negative
            return(svd.tmp)
        },
        finally = {
        })
    # svd2 <- svd(spec2)
    svd2=tryCatch({svd(spec2)},
        error = function(err){
            svd.tmp=(eigen(spec2, symmetric=TRUE))
            names(svd.tmp)=c("d", "u")
            svd.tmp$d=pmax(svd.tmp$d, 0) ## correct the eigenvalues if they are negative
            return(svd.tmp)
        },
        finally = {
        })
    #svd.pooled=svd((spec1+spec2)/2)
    svd.pooled=tryCatch({svd((spec1+spec2)/2)},
        error = function(err){
            svd.tmp=(eigen((spec1+spec2)/2, symmetric=TRUE))
            names(svd.tmp)=c("d", "u")
            svd.tmp$d=pmax(svd.tmp$d, 0) ## correct the eigenvalues if they are negative
            return(svd.tmp)
        },
        finally = {
        })
    #
    u1=svd1$u
    d1=svd1$d
    rm(svd1)
    #
    u2=svd2$u
    d2=svd2$d
    rm(svd2)
    #
    u=svd.pooled$u
    d=svd.pooled$d
    #
    if(autok == 0){ ## no automatic choice of K
        if( is.na(K.fixed) ){
            stop("K.fixed should be specified if autok=0")
        } else {
            K = K.fixed 
        }
    } else if(autok == 1){ ## AIC
        m.loc= m/kappa.square
        gof1=m.loc* Re(c(rev(cumsum(rev(diag(Conj(t(u))%*% spec1 %*% u))[-n])),0))
        gof2=m.loc* Re(c(rev(cumsum(rev(diag(Conj(t(u))%*% spec2 %*% u))[-n])),0))
        #
        ## Victor's version
        pen1=numeric(n)
        pen2=numeric(n)
        for(k in 1:n){
            pen1[k] = sum((diag(Re(Conj(t(u1)) %*% u[, 1:k] %*% Conj(t(u[,1:k])) %*% spec1 %*% u[, 1:k] %*% Conj(t(u[,1:k])) %*% u1))/d1))
            pen2[k] = sum((diag(Re(Conj(t(u2)) %*% u[, 1:k] %*% Conj(t(u[,1:k])) %*% spec2 %*% u[, 1:k] %*% Conj(t(u[,1:k])) %*% u2))/d2))
        }
        pen1=  sum(d) * pen1
        pen2=  sum(d) * pen2
        #
        K=which.min(pen1 + pen2 + gof1 +gof2) 
        #print(K)
    } else if (autok == 2){ ## AIC*
        m.loc= m/kappa.square
        #
        tail.d1=rev(cumsum(rev(d1)[-n]))
        tail.d2=rev(cumsum(rev(d2)[-n]))
        gof1=m.loc* Re(c(rev(cumsum(rev(diag(Conj(t(u))%*% spec1 %*% u))[-n])))) 
        gof2=m.loc* Re(c(rev(cumsum(rev(diag(Conj(t(u))%*% spec2 %*% u))[-n]))))
        #
        tmp.invdiff=1/(d[-n]-d[-1])
        invdiff=apply(cbind(c(tmp.invdiff,0), c(0, tmp.invdiff)), 1, max)
        rm(tmp.invdiff)
        tmp.invdiff=1/(d1[-n]-d1[-1])
        invdiff1=apply(cbind(c(tmp.invdiff,0), c(0, tmp.invdiff)), 1, max)
        rm(tmp.invdiff)
        tmp.invdiff=1/(d2[-n]-d2[-1])
        invdiff2=apply(cbind(c(tmp.invdiff,0), c(0, tmp.invdiff)), 1, max)
        #
        pen1=numeric(n)
        pen2=numeric(n)
        #
        for(k in 1:(n)){
            pen1[k] = sum((diag(Re(Conj(t(u1)) %*% u[, 1:k] %*% Conj(t(u[,1:k])) %*% spec1 %*% u[, 1:k] %*% Conj(t(u[,1:k])) %*% u1)))*sqrt(invdiff1/d1))
            pen2[k] = sum((diag(Re(Conj(t(u2)) %*% u[, 1:k] %*% Conj(t(u[,1:k])) %*% spec2 %*% u[, 1:k] %*% Conj(t(u[,1:k])) %*% u2)))*sqrt(invdiff2/d2))
        }
        pen1=  cumsum(d) * pen1
        pen2=  cumsum(d) * pen2
        pen1=pen1[-n]
        pen2=pen2[-n]
        #
        K=which.min(pen1 + pen2 + gof1 + gof2)
        #print(K)
    } else {
        stop("please set value of 'autok' to either 1 or 2")
    }
    #
    Delta=spec1-spec2
    #
    svd.glob=svd.pooled
    eigvect.glob=svd.glob$u[,1:K]
    eigval.glob=svd.glob$d[1:K]
    test.res=2*pi*m*({Mod(Conj(t(eigvect.glob)) %*% Delta[,] %*% eigvect.glob)^2 }/ ((4*pi*kappa.square*(eigval.glob %o% eigval.glob ))))
    # 
    if(is.pi.multiple){## Adjust formula of test statistic and degrees of freedom of the chi2
        test.res=test.res/2 
        if(autok==0){ ## if no automatic choice of K, then return values of the statistic for k=1,2, ..., K
            pval=numeric(K)
            for(k in 1:K){
                pval[k]=stats::pchisq(sum(test.res[1:k,1:k]), df=k*(k+1)/2, lower.tail=F)
            }
        } else {
            pval=stats::pchisq(sum(test.res[1:K,1:K]), df=K*(K+1)/2, lower.tail=F)
        }
    } else {
        if(autok==0){ ## if no automatic choice of K, then return values of the statistic for k=1,2, ..., K
            pval=numeric(K)
            for(k in 1:K){
                pval[k]=stats::pchisq(sum(test.res[1:k,1:k]), df=k^2, lower.tail=F)
            }
        } else {
            pval=stats::pchisq(sum(test.res[1:K,1:K]), df=K^2, lower.tail=F)
        }
    }    
    return(list(pval=pval, K=K))
} 

#' Generic function to adjust pvalues
#'
#' @param sample.spec.diff Object of the class
#' \code{SampleSpecDiffFreq}
#' @param method method used to adjust p-values
#'
#' @export
PvalAdjust <- function(sample.spec.diff, method) UseMethod("PvalAdjust", sample.spec.diff)

#' function to adjust pvalues for class \code{SampleSpecDiffFreq}
#'
#' @rdname PvalAdjust
#'
#' @seealso \code{\link{Spec_compare_localize_freq}}
#' @export
PvalAdjust.SampleSpecDiffFreq <- function(sample.spec.diff, method){
    ## distinguish autok
    if(sample.spec.diff$autok == 0){
        return( apply(sample.spec.diff$pval.raw, 2, function(x) stats::p.adjust(x, method=method)))
    } else {
        return(stats::p.adjust(sample.spec.diff$pval.raw, method=method))
    }
}



#' Plotting function for \code{SampleSpecDiffFreq} class
#'
#' @param x object of the class \code{SampleSpecDiffFreq}
#' @param method method used to adjust p-values
#' @param Kmax maximum number of levels K for which the pvalues are plotted
#' (used only if autok==0)
#' @param pch the plot character to be used
#' @param ... additional parameters to be passed to plot()
#'
#' @seealso \code{\link{Spec_compare_localize_freq}}
#' @export

plot.SampleSpecDiffFreq <- 
    function(x, method=NA, Kmax=4, pch=20, ...)
    {
        if( is.na(method) ) {
            pvals  <- x$pval.raw
            main <- 'raw pvalues'
        } else {
            pvals <- PvalAdjust(x, method=method)
            main <- paste0('adjusted pvalues (', method, ')')
        }
        #
        if(x$autok == 0){
            Kmax=min( ncol(pvals), Kmax )
            graphics::matplot(x$omega, pvals[,1:Kmax], ylim=c(0, 1.01), type='o', col=1, pch=pch, main=main, ...)
        } else {
            graphics::plot(x$omega, pvals, ylim=c(0, 1.01), type='o', col=1, pch=pch, main=main, ...)
        }
        rm(pvals)
    }


#' Plotting function for \code{SampleSpecDiffFreq} class
#'
#' @inheritParams plot.SampleSpecDiffFreq
#'
#' @seealso \code{\link{Spec_compare_localize_freq}}
#' @export

lines.SampleSpecDiffFreq <- 
    function(x, method=NA, Kmax=4, pch=20, ...)
    {
        if( is.na(method) ) {
            pvals  <- x$pval.raw
            main <- 'raw pvalues'
        } else {
            pvals <- PvalAdjust(x, method=method)
            main <- paste0('adjusted pvalues (', method, ')')
        }
        #
        if(x$autok == 0){
            Kmax=min( ncol(pvals), Kmax )
            graphics::matlines(x$omega, pvals[,1:Kmax], ylim=c(0, 1.01), type='o', pch=pch, main=main, ...)
        } else {
            graphics::lines(x$omega, pvals, ylim=c(0, 1.01), type='o',  pch=pch, main=main, ...)
        }
        rm(pvals)
    }




#' Compare the spectral density operator of two Functional Time Series and localize frequencies at which they differ.
#'
#' @param X,Y The \eqn{T \times nbasis} matrices of containing the coordinates, expressed in some functional basis,  of the two FTS that to be compared.  
#' expressed in a basis. 
#' @inheritParams Spec 
#' @inheritParams Spec_compare_fixed_freq 
#' 
#' @section Details:
#' 
#' \code{X,Y} must be of equal size \eqn{T.len \times d}, where T.len is the length of the time series, and \eqn{d} is the number of basis functions. Each row corresponds to a time point, and each column
#' corresponds to the coefficient of the corresponding basis function of the  FTS.
#' 
#' \code{autok=0} returns the p-values for \eqn{K=1, \ldots, \code{K.fixed}}.
#' \code{autok=1} uses the \code{AIC} criterion of Tavakoli \& Panaretos
#' (2015), which is a generalization of the pseudo-AIC introduced in Panaretos
#' et al (2010).
#' \code{autok=2} uses the \code{AIC*} criterion of Tavakoli \& Panaretos
#' (2015), which is an extension of the \code{AIC} criterion that takes into
#' account the difficulty associated with the estimation of eigenvalues of a
#' compact operator. 
#' 
#' @examples 
#' 
#' ma.scale2=ma.scale1=c(-1.4,2.3,-2)
#' ma.scale2[3] = ma.scale1[3]+.0
#' a1=Generate_filterMA(10, 10, MA.len=3, ma.scale=ma.scale1)
#' a2=Generate_filterMA(10, 10, MA.len=3, ma.scale=ma.scale2)
#' X=Simulate_new_MA(a1, T.len=512, noise.type='wiener')
#' Y=Simulate_new_MA(a2, T.len=512, noise.type='wiener') 
#' ans0=Spec_compare_localize_freq(X, Y, W=Epanechnikov_kernel, autok=2,
#' subgrid.density=10, verbose=0, demean=FALSE,
#' subgrid.density.relative.to.bandwidth=TRUE)
#' plot(ans0)
#' plot(ans0, method='fdr')
#' PvalAdjust(ans0, method='fdr') ## print FDR adjusted p-values
#' abline(h=.05, lty=3)
#' ans0=Spec_compare_localize_freq(X, Y, W=Epanechnikov_kernel, autok=0,
#' subgrid.density=10, verbose=0, demean=FALSE,
#' subgrid.density.relative.to.bandwidth=TRUE, K.fixed=4) ## fixed values of K
#' plot(ans0)
#' plot(ans0, 'fdr')
#' plot(ans0, 'holm')
#' PvalAdjust(ans0, method='fdr')
#' rm(ans0)
#' 
#' @export
#'
#' @references \cite{Tavakoli, Shahin and Panaretos, Victor M. "Detecting and Localizing Differences in Functional Time Series Dynamics: A Case Study in Molecular Biophysics", 2014, under revision}
#' @references \cite{Panaretos, Victor M., David Kraus, and John H. Maddocks.
#' "Second-order comparison of Gaussian random functions and the geometry of
#' DNA minicircles." Journal of the American Statistical Association 105.490
#' (2010): 670-682.}

Spec_compare_localize_freq=function(X, Y, B.T=(dim(X)[1])^(-1/5), W, autok=2, subgrid.density, verbose=0, demean=FALSE, K.fixed=NA, subgrid.density.relative.to.bandwidth){
#
    X<-as.matrix(X)
    Y<-as.matrix(Y)
    T<-dim(X)[1]
#
    #### perform some checks
    if(!isTRUE(all.equal(dim(X), dim(Y)))) stop("Dimensions of 'X' and 'Y' must be equal")
    #### compute sample spectral density operators for X and Y
    X.spec=Spec(X, W=W, B.T=B.T, only.diag=FALSE, trace=FALSE,  subgrid=TRUE, subgrid.density=subgrid.density, verbose=verbose, demean=demean, subgrid.density.relative.to.bandwidth=subgrid.density.relative.to.bandwidth)
    Y.spec=Spec(Y, W=W, B.T=B.T, only.diag=FALSE, trace=FALSE,  subgrid=TRUE, subgrid.density=subgrid.density, verbose=verbose, demean=demean, subgrid.density.relative.to.bandwidth=subgrid.density.relative.to.bandwidth)
    m=X.spec$m
    kappa.square=X.spec$kappa.square
    omega=X.spec$omega
    sigma.W=X.spec$sigma.W
    bw=X.spec$bw
#
    #### compute raw pvalue at each frequency for the test of equality of the spectral density operators
    #
    ## remove frequencies too close to 0,pi than B.T=bw/sigma.W
    low.cut=min(which( omega>= B.T))
    up.cut=max(which( omega<= pi-B.T)) 
    testgrid.i=c(1,low.cut:up.cut, length(omega)) ## the subgrid of omega; contains the indices of the frequencies omega that will be used for inference
    #
    if(autok==0){ ## fixed choice of K
        if(is.na(K.fixed)) stop("K.fixed must be specified if autok==0")
    pval.raw=array(NA, c(length(testgrid.i), K.fixed ), dimnames=list(c(), paste0("K=", 1:K.fixed)))
#
    for(freqi in 1:length(testgrid.i)){
    if(verbose)    cat("Frequency : ", freqi, "/", length(testgrid.i), "\n", sep="")
        if(freqi %in% c(1, length(testgrid.i))){
            is.pi.multiple=TRUE
        }else{
            is.pi.multiple=FALSE
        }
        ans.tmp=Spec_compare_fixed_freq(X.spec$spec[testgrid.i[freqi],,], Y.spec$spec[testgrid.i[freqi],,], m=m, is.pi.multiple=is.pi.multiple, kappa.square=kappa.square, autok=autok, K.fixed=K.fixed)
        pval.raw[freqi,]=ans.tmp$pval
    }
#
    #### return results
    spec.diff.freq=list(omega=omega[testgrid.i], pval.raw=pval.raw, K=K.fixed, B.T=B.T, bw=bw, m=m, autok=autok)
    class(spec.diff.freq) <- 'SampleSpecDiffFreq'
    return(spec.diff.freq)
    } else {  ## automatic choice of K
    pval.raw=array(NA, length(testgrid.i))
    K=array(NA, length(testgrid.i))
    for(freqi in 1:length(testgrid.i)){
    if(verbose)    cat("Frequency : ", freqi, "/", length(testgrid.i), "\n", sep="")
        if(freqi %in% c(1, length(testgrid.i))){
            is.pi.multiple=TRUE
        }else{
            is.pi.multiple=FALSE
        }
        ans.tmp=Spec_compare_fixed_freq(X.spec$spec[testgrid.i[freqi],,], Y.spec$spec[testgrid.i[freqi],,], m=m, is.pi.multiple=is.pi.multiple, kappa.square=kappa.square, autok=autok, K.fixed=K.fixed)
        pval.raw[freqi]=ans.tmp$pval
        K[freqi]=ans.tmp$K
    }
    #### return results
    spec.diff.freq=list(omega=omega[testgrid.i], pval.raw=pval.raw, K=K, B.T=B.T, bw=bw, m=m, autok=autok)
    class(spec.diff.freq) <- 'SampleSpecDiffFreq'
    return(spec.diff.freq)
    }
}



## #############################################################################
## Functions for localization of differences in frequency and along curve length
## #############################################################################


#' Compute the marginal  p-values at each basis coefficients of for testing the equality of two spectral density kernels
#' 
#' @inheritParams Spec_compare_fixed_freq
#' 
#' @export
#' 
Marginal_basis_pval=function(spec1, spec2, m, kappa.square, is.pi.multiple){
    data=spec1-spec2
    var.est=(spec1+spec2)/2
    kappa=sqrt(kappa.square)
#
    if(!is.pi.multiple){
        tmp.var=Re(diag(var.est)) 
        normalizing.factor.diag2=2*kappa.square*tmp.var^2/m
        specdiff.basis=array(NA, dim(spec1))
       # 
        P=Re(2*kappa^2*(outer(tmp.var,tmp.var) - Mod(var.est)^4/outer(tmp.var,tmp.var))/m)
        diag(P)=1 ## avoid division by zero on the diagonal; we shall correct the diagonal later on
        R=Conj(var.est^2)/outer(tmp.var,tmp.var)
        specdiff.basis= 2*( Mod(data)^2/Conj(P) - Re( (data)^2*R/Conj(P)) )
        ### correcting the diagonal
        diag(specdiff.basis)=NA
        diag(specdiff.basis)=Mod(diag(data))^2/normalizing.factor.diag2
        #
        df.tmp=array(2, dim(specdiff.basis))
        diag(df.tmp)=1
        pval.basis<-stats::pchisq(specdiff.basis, df=df.tmp, lower.tail=FALSE) 
        rm(tmp.var, df.tmp, normalizing.factor.diag2)
    }else{
        tmp.var=diag(Re(var.est)) 
        normalizing.factor2=2*kappa^2*(outer((tmp.var), (tmp.var)) + Mod(var.est)^2)/m   
        specdiff.basis=(Mod(data)^2)/normalizing.factor2 
        pval.basis<-stats::pchisq(specdiff.basis, df=1, lower.tail=FALSE) 
        rm(tmp.var, normalizing.factor2)
    }
    return(pval.basis)
}



#' Compare the spectral density operator of two Functional Time Series and localize frequencies at which they differ, and (spatial) regions where they differ
#'
#' @inheritParams Spec_compare_localize_freq
#' @param alpha level for the test
#' @param accept,reject values for accepted, rejected regions
#' 
#' @examples
#' ma.scale2=ma.scale1=c(-1.4,2.3,-2)
#' ma.scale2[3] = ma.scale1[3]+.4
#' a1=Generate_filterMA(10, 10, MA.len=3, ma.scale=ma.scale1)
#' a2=Generate_filterMA(10, 10, MA.len=3, ma.scale=ma.scale2)
#' X=Simulate_new_MA(a1, T.len=2^9, noise.type='wiener')
#' Y=Simulate_new_MA(a2, T.len=2^9, noise.type='wiener') 
#' ans0=Spec_compare_localize_freq_curvelength(X, Y, W=Epanechnikov_kernel, alpha=.01, demean=TRUE)
#' print(ans0)
#' plot(ans0)
#' rm(ma.scale1, ma.scale2, a1, a2, X, Y, ans0)
#' @export

Spec_compare_localize_freq_curvelength = function(X, Y, B.T=(dim(X)[1])^(-1/5), W,  alpha=0.05, accept=0, reject=1, verbose=0, demean=FALSE){
    #
    ####### local functions
    # Adjust the p-values of a symmetric matrix of p-values for multiplicity 
    Upper_tri_p_adjust=function(pval.matrix, method){
        #require(sna)
        M=pval.matrix
        up=upper.tri(M, diag=TRUE)
        M[!up]=0
        M[up]=stats::p.adjust(M[up], method=method)
        return(sna::symmetrize(M, "upper"))
    }
    # Adjust for multiplicity the p-values of an array of symmetric p-value matrices,  marginally for each matrix
    P_array_upper_tri_adjust=function(pval.array, method){
        aperm(sapply(1:dim(pval.array)[1], function(freqi){ Upper_tri_p_adjust(pval.array[freqi,,], method=method)}, simplify="array"), c(3,1,2))
    }
    ####### (END) local functions 
    #
    X.spec=Spec(X, W=W, B.T=B.T, only.diag=FALSE, trace=FALSE,  subgrid=TRUE, subgrid.density=1, verbose=verbose, demean=demean, subgrid.density.relative.to.bandwidth=TRUE)
    Y.spec=Spec(Y, W=W, B.T=B.T, only.diag=FALSE, trace=FALSE,  subgrid=TRUE, subgrid.density=1, verbose=verbose, demean=demean, subgrid.density.relative.to.bandwidth=TRUE)
    #
    dim.spec<-dim(X.spec$spec)
    m=X.spec$m
    kappa.square=X.spec$kappa.square
    omega=X.spec$omega
    sigma.W=X.spec$sigma.W
    #
    pval.basis<-array(NA, dim.spec) ### the pointwise p-values on the basis coefficients
    ## compute pointwise p-values at each frequency
    for(freqi in 1:length(omega)){
        if(verbose)    cat("Frequency : ", freqi, "/", length(omega), "\n", sep="")
        if(freqi %in% c(1, length(omega))){
            is.pi.multiple=TRUE
        }else{
            is.pi.multiple=FALSE
        }
        pval.basis[freqi,,]=Marginal_basis_pval(X.spec$spec[freqi,,], Y.spec$spec[freqi,,], m=m, kappa.square=kappa.square, is.pi.multiple=is.pi.multiple)
    }
    ## adjust p-values within frequencies
    pval.basis.adj=P_array_upper_tri_adjust(pval.basis, method="fdr")
    ## apply Benjamini-Bogomolov criterion
    pval.tmp=pval.basis.adj
    min.adj.pval=apply(pval.tmp, 1, min)
    reject.freq=(stats::p.adjust(min.adj.pval, method='fdr') <= alpha)
    n.rej=sum(reject.freq)
    itmp=pval.tmp>n.rej*alpha/length(omega)
    pval.tmp[itmp]=accept
    pval.tmp[!itmp]=reject
    #
    sample.spec.diff <- list(values=pval.tmp, alpha=alpha, omega=omega, reject.freq=reject.freq, min.per.freq=min.adj.pval)
    class(sample.spec.diff) <- 'SampleSpecDiffFreqCurvelength'
    return(sample.spec.diff)
}

#' Printing method for class \code{SampleSpecDiffFreqCurvelength}
#' @param x Object of the class \code{SampleSpecDiffFreqCurvelength}
#' @param ... Additional arguments for print
#' @export
print.SampleSpecDiffFreqCurvelength <- function(x, ...){
    print("Comparison of Spectral Densities of Two Functional Time Series, Localization of Differences in Frequency and along curve length")
    print("--Localization in Frequencies--")
    print(data.frame(freq=x$omega, difference=x$reject.freq))
    print("--Localization along curve length, within Frequencies--")
    print(data.frame(freq=x$omega, difference=x$reject.freq, percentage.of.difference=format( 100*apply(x$values, 1, mean), digits=2)))
}

#' Plotting method for class \code{SampleSpecDiffFreqCurvelength}
#' @param x Object of the class \code{SampleSpecDiffFreqCurvelength}
#' @param ncolumns number of columns for the plots
#' @param ... additional parameters to be passed to \code{plot()}
#' @export
#' @importFrom grDevices colorRampPalette

plot.SampleSpecDiffFreqCurvelength <- function(x, ncolumns=3, ...){
    oop=graphics::par()
    cex.main=.8
    cex.axis=.8
    cex.lab.big1=1
    mai=c(0.0,.01,0.5,0.1)
    oma=c(1,.5, 1,0)
    pal.loc=grDevices::colorRampPalette(c("white", "black"), space="rgb")
    ncols=30
    diff.col="#00000015"
    {
        xlab=""
        xaxt='n'
        graphics::par(mfrow=c( ceiling(length(x$omega)/ncolumns) ,ncolumns), pty='s' )
        graphics::par(oma=oma)
        graphics::par(mai=mai)
        len.omega=length(x$omega)
        for(i in 1:len.omega){
            graphics::image(x$values[i,,], asp=1, cex.axis=cex.axis, xlim=c(0,1), xaxt=xaxt, yaxt='n', col=pal.loc(ncols), breaks=c(seq(0,.4, len=ncols), 2), cex.lab=cex.lab.big1)
            tmp.omega=round(x$omega[i], 2)
            graphics::mtext(bquote(omega == .(tmp.omega)), outer=FALSE, line=0, cex=cex.main)
            rm(tmp.omega)
        }
        #title("Main TITLE", outer=T)
    }
    graphics::par(oop)
}

## #############################################################################
## (END) Functions for localization of differences in frequency and along curve length
## #############################################################################

