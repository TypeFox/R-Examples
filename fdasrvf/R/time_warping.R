#' Group-wise function alignment
#'
#' This function aligns a collection of functions using the elastic square-root
#' slope (srsf) framework.
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param lambda controls the elasticity (default = 0)
#' @param method warp and caluclate to Karcher Mean or Median (options = "mean"
#' or "median", default = "mean")
#' @param showplot shows plots of functions (default = T)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param parallel enable parallel mode using \code{\link{foreach}} and
#'   \code{doParallel} pacakge
#' @param cores set number of cores to use with \code{doParallel} (default = 2)
#' @param omethod optimization method (DP,DP2,SIMUL,RBFGS)
#' @param MaxItr maximum number of iterations
#' @return Returns a list containing \item{f0}{original functions}
#' \item{fn}{aligned functions - matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples}
#' \item{qn}{aligned SRSFs - similar structure to fn}
#' \item{q0}{original SRSF - similar structure to fn}
#' \item{fmean}{function mean or median - vector of length \eqn{N}}
#' \item{mqn}{SRSF mean or median - vector of length \eqn{N}}
#' \item{gam}{warping functions - similar structure to fn}
#' \item{orig.var}{Original Variance of Functions}
#' \item{amp.var}{Amplitude Variance}
#' \item{phase.var}{Phase Variance}
#' \item{qun}{Cost Function Value}
#' @keywords srsf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' out = time_warping(simu_data$f,simu_data$time, MaxItr=1) # use more iterations for accuracy
time_warping <- function(f, time, lambda = 0, method = "mean",
                         showplot = TRUE, smooth_data = FALSE, sparam = 25,
                         parallel = FALSE, cores=2, omethod = "DP", MaxItr = 20){
    if (parallel){
        cl = makeCluster(cores)
        registerDoParallel(cl)
    } else
    {
        registerDoSEQ()
    }
    method <- pmatch(method, c("mean", "median")) # 1 - mean, 2 - median
    if (is.na(method))
        stop("invalid method selection")

    cat(sprintf("lambda = %5.1f \n",lambda))

    binsize = mean(diff(time))
    eps = .Machine$double.eps
    M = nrow(f)
    N = ncol(f)
    f0 = f
    w = 0.0;

    if (smooth_data){
        f = smooth.data(f,sparam)
    }

    if (showplot){
        matplot(time,f,type="l")
        title(main="Original data")
    }

    # Compute q-function of the functional data
    tmp = gradient.spline(f,binsize,smooth_data)
    f = tmp$f
    q = tmp$g/sqrt(abs(tmp$g)+eps)

    cat("\nInitializing...\n")
    mnq = rowMeans(q)
    dqq = sqrt(colSums((q - matrix(mnq,ncol=N,nrow=M))^2))
    min_ind = which.min(dqq)
    mq = q[,min_ind]
    mf = f[,min_ind]

    gam<-foreach(k = 1:N, .combine=cbind,.packages="fdasrvf") %dopar% {
        gam_tmp = optimum.reparam(mq,time,q[,k],time,lambda,omethod,w,mf[1],f[1,k])
    }

    gam = t(gam)
    gamI = SqrtMeanInverse(gam)
    gamI_dev = gradient(gamI, 1/(M-1))
    mf = approx(time,mf,xout=(time[length(time)]-time[1])*gamI + time[1])$y
    mq = f_to_srvf(mf,time)
    mq[is.nan(mq)] <- 0

    # Compute Mean
    if (method == 1)
        cat(sprintf("Computing Karcher mean of %d functions in SRSF space...\n",N))
    if (method == 2)
        cat(sprintf("Computing median of %d functions in SRSF space...\n",N))
    ds = rep(0,MaxItr+2)
    ds[1] = Inf
    qun = rep(0,MaxItr)
    tmp = matrix(0,M,MaxItr+2)
    tmp[,1] = mq
    mq = tmp
    tmp = matrix(0,M,MaxItr+2)
    tmp[,1] = mf
    mf = tmp
    tmp = array(0,dim=c(M,N,MaxItr+2))
    tmp[,,1] = f
    f = tmp
    tmp = array(0,dim=c(M,N,MaxItr+2))
    tmp[,,1] = q
    q = tmp
    qun = rep(0,MaxItr+1)
    for (r in 1:MaxItr){
        cat(sprintf("updating step: r=%d\n", r))
        if (r == MaxItr){
            cat("maximal number of iterations is reached. \n")
        }

        # Matching Step
        outfor<-foreach(k = 1:N, .combine=cbind,.packages='fdasrvf') %dopar% {
            gam = optimum.reparam(mq[,r],time,q[,k,1],time,lambda,omethod,w,mf[1,r],f[1,k,1])
            gam_dev = gradient(gam,1/(M-1))
            f_temp = approx(time,f[,k,1],xout=(time[length(time)]-time[1])*gam +
                time[1])$y
            q_temp = f_to_srvf(f_temp,time)
            list(gam,gam_dev,q_temp,f_temp)
        }
        gam = unlist(outfor[1,]);
        dim(gam)=c(M,N)
        gam = t(gam)
        gam_dev = unlist(outfor[2,]);
        dim(gam_dev)=c(M,N)
        gam_dev = t(gam_dev)
        q_temp = unlist(outfor[3,]);
        dim(q_temp)=c(M,N)
        f_temp = unlist(outfor[4,]);
        dim(f_temp)=c(M,N)
        q[,,r+1] = q_temp
        f[,,r+1] = f_temp
        tmp = (1-sqrt(gam_dev))^2

        if (method == 1){ # Mean
            ds_tmp  = sum(simpson(time,(matrix(mq[,r],M,N)-q[,,r+1])^2)) +
                lambda*sum(simpson(time, t(tmp)))
            if (is.complex(ds_tmp)){
                ds[r+1] = abs(ds_tmp)
            }
            else{
                ds[r+1] = ds_tmp
            }

            # Minimization Step
            # compute the mean of the matched function
            mq[,r+1] = rowMeans(q[,,r+1])
            mf[,r+1] = rowMeans(f[,,r+1])

            qun[r] = pvecnorm(mq[,r+1]-mq[,r],2)/pvecnorm(mq[,r],2)
        }

        if (method == 2){ # Median
            ds_tmp = sqrt(sum(simpson(time,(matrix(mq[,r],M,N)-q[,,r+1])^2))) +
                lambda*sum(simpson(time, t(tmp)))
            if (is.complex(ds_tmp)){
                ds[r+1] = abs(ds_tmp)
            }
            else{
                ds[r+1] = ds_tmp
            }

            # Minimization Step
            # compute the median of the matched function
            dist_iinv = (sum(1/ds[r+1]))^(-1)
            mq[,r+1] = (rowSums(q[,,r+1]/ds[r+1]))*dist_iinv
            mf[,r+1] = (rowSums(f[,,r+1]/ds[r+1]))*dist_iinv

            qun[r] = pvecnorm(mq[,r+1]-mq[,r],2)/pvecnorm(mq[,r],2)
        }
        if (qun[r] < 1e-4 || r >=MaxItr){
            break
        }
    }

    # One last run, centering of gam
    r = r+1
    outfor<-foreach(k = 1:N, .combine=cbind,.packages="fdasrvf") %dopar% {
        gam = optimum.reparam(mq[,r],time,q[,k,1],time,lambda,omethod,w,mf[1,r],f[1,k,1])
        gam_dev = gradient(gam,1/(M-1))
        list(gam,gam_dev)
    }
    gam = unlist(outfor[1,]);
    dim(gam)=c(M,N)
    gam = t(gam)
    gam_dev = unlist(outfor[2,]);
    dim(gam_dev)=c(M,N)
    gam_dev = t(gam_dev)

    gamI = SqrtMeanInverse(gam)
    gamI_dev = gradient(gamI, 1/(M-1))
    mq[,r+1] = approx(time,mq[,r],xout=(time[length(time)]-time[1])*gamI +
        time[1])$y*sqrt(gamI_dev)

    for (k in 1:N){
        q[,k,r+1] = approx(time,q[,k,r],xout=(time[length(time)]-time[1])*gamI +
            time[1])$y*sqrt(gamI_dev)
        f[,k,r+1] = approx(time,f[,k,r],xout=(time[length(time)]-time[1])*gamI +
            time[1])$y
        gam[k,] = approx(time,gam[k,],xout=(time[length(time)]-time[1])*gamI +
            time[1])$y
    }

    # Aligned data & stats
    fn = f[,,r+1]
    qn = q[,,r+1]
    q0 = q[,,1]
    mean_f0 = rowMeans(f[,,1]);
    std_f0 = apply(f[,,1], 1, sd)
    mean_fn = rowMeans(fn)
    std_fn = apply(fn, 1, sd)
    mqn = mq[,r+1]
    fmean = mean(f0[1,])+cumtrapz(time,mqn*abs(mqn));
    gam = t(gam)

    fgam = matrix(0,M,N)
    for (ii in 1:N){
        fgam[,ii] = approx(time,fmean,xout=(time[length(time)]-time[1])*gam[,ii] +
            time[1])$y
    }
    var_fgam = apply(fgam,1,var)

    orig.var = trapz(time,std_f0^2)
    amp.var = trapz(time,std_fn^2)
    phase.var = trapz(time,var_fgam)

    if (showplot){
        matplot((0:(M-1))/(M-1),gam,type="l",main="Warping functions",xlab="Time")

        matplot(time,fn,type="l",main=bquote(paste("Warped Data ",lambda ==
            .(lambda))))

        matplot(time,cbind(mean_f0,mean_f0+std_f0,mean_f0-std_f0),type="l",lty=1,
                        col=c("blue","red","green"),
                        ylab="",main=bquote(paste("Original Data: ", Mean %+-% STD)))
        legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),
                     col=c('blue','red','green'),lty=1)

        matplot(time,cbind(mean_fn,mean_fn+std_fn,mean_fn-std_fn),type="l",lty=1,
                        col=c("blue","red","green"),
                        ylab="",main=bquote(paste("Warped Data: ",lambda == .(lambda),": ",
                                                                            Mean %+-% STD)))
        legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),
                     col=c('blue','red','green'),lty=1)

        plot(time,fmean,type="l",col="green",main=bquote(paste(f[mean]," ",
                                                                                                                     lambda == .(lambda))))
    }

    if (parallel){
        stopCluster(cl)
    }


    return(list(f0=f[,,1],fn=fn,qn=qn,q0=q0,fmean=fmean,mqn=mqn,gam=gam,
                            orig.var=orig.var,amp.var=amp.var,phase.var=phase.var,qun=qun))

}
