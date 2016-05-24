#' Group-wise function alignment and PCA Extractions
#'
#' This function aligns a collection of functions while extracting pincipal
#' components.
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param num_comp number of principal components to extract (default = 3)
#' @param showplot shows plots of functions (default = T)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param parallel enable parallel mode using \code{\link{foreach}} and
#'   \code{doParallel} pacakge
#' @param cores set number of cores to use with \code{doParallel} (default = 2)
#' @param MaxItr maximum number of iterations
#' @return Returns a list containing \item{f0}{original functions}
#' \item{fn}{aligned functions - matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples}
#' \item{qn}{aligned srvfs - similar structure to fn}
#' \item{q0}{original srvfs - similar structure to fn}
#' \item{mqn}{srvf mean - vecotr of length \eqn{N}}
#' \item{gam}{warping functions - vecotr of length \eqn{N}}
#' \item{Dx}{cost function}
#' \item{vfpca}{list containing}
#' \item{q_pca}{srvf principal directions}
#' \item{f_pca}{f principal directions}
#' \item{latent}{latent values}
#' \item{coef}{coefficients}
#' \item{U}{eigenvectors}
#' @keywords srvf alignment, pca
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' out = align_fPCA(simu_data$f,simu_data$time,MaxItr = 1)  # use more iterations for accuracy
align_fPCA <- function(f, time, num_comp = 3, showplot = T, smooth_data = FALSE, sparam = 25,
                       parallel = FALSE, cores=8, MaxItr = 51){
    if (parallel){
        cl = makeCluster(cores)
        registerDoParallel(cl)
    } else
    {
        registerDoSEQ()
    }

    cat("\nInitializing...\n")
    binsize = mean(diff(time))
    eps = .Machine$double.eps
    M = nrow(f)
    N = ncol(f)
    f0 = f
    lambda = 0
    coef = -2:2
    NP = 1:num_comp  # number of principal components
    Nstd = length(coef)

    if (smooth_data){
        f = smooth.data(f,sparam)
    }

    if (showplot){
        matplot(time,f,type="l")
        title(main="Original data")
    }

    # Compute q-function of the plot
    tmp = gradient.spline(f,binsize)
    f = tmp$f
    q = tmp$g/sqrt(abs(tmp$g)+eps)

    # PCA Step
    mnq = rowMeans(q)
    dqq = sqrt(colSums((q - matrix(mnq,ncol=N,nrow=M))^2))
    min_ind = which.min(dqq)
    mq = q[,min_ind]
    qhat_cent = q-matrix(mq,M,N)
    K = 1/M * qhat_cent %*% t(qhat_cent)
    out = svd2(K)
    s = out$d
    U = out$u

    alpha_i = matrix(0,num_comp,N)
    for (ii in 1:num_comp){
        for (jj in 1:N){
            alpha_i[ii,jj] = simpson(time,qhat_cent[,jj]*U[,ii])
        }
    }

    tmp = U[,1:num_comp]%*%alpha_i
    if (smooth_data){
        mq = mq/pvecnorm(mq,2)
        for (ii in 1:N){
            if (sum(tmp[,ii]) != 0)
                tmp[,ii] = tmp[,ii]/pvecnorm(tmp[,ii],2)
        }
    }
    qhat = matrix(mq,M,N) + tmp

    # Matching Step
    gam_k<-foreach(k = 1:N, .combine=cbind,.packages="fdasrvf") %dopar% {
        gam0 = optimum.reparam(qhat[,k],time,q[,k],time)
    }

    cat(sprintf("Aligning %d functions in SRVF space to %d fPCA components...\n",N,num_comp))
    tmp = matrix(0,M,MaxItr+2)
    tmp[,1] = mq
    mq = tmp
    tmp = array(0,dim=c(M,N,MaxItr+2))
    tmp[,,1] = f
    f = tmp
    tmp = array(0,dim=c(M,N,MaxItr+2))
    tmp[,,1] = q
    q = tmp
    tmp = array(0,dim=c(M,N,MaxItr+2))
    tmp[,,1] = gam_k
    gam = tmp
    Dx = rep(0,MaxItr+1)
    for (r in 2:MaxItr){
        cat(sprintf("updating step: r=%d\n", r-1))
        if (r == MaxItr){
            cat("maximal number of iterations is reached. \n")
        }

        # PCA Step
        f_temp = matrix(0,M,N)
        q_temp = matrix(0,M,N)
        for (k in 1:N){
            gam_k = as.vector(gam[,k,r-1])
            f_temp[,k] = approx(time,f[,k,r-1],xout=(time[length(time)]-time[1])*gam_k + time[1])$y
            q_temp[,k] = f_to_srvf(f_temp[,k],time)
        }
        f[,,r] = f_temp
        q[,,r] = q_temp
        mq[,r] = rowMeans(q_temp)

        K = cov(t(q[,,r]))
        out = svd(K)
        s = out$d
        U = out$u

        qhat_cent = scale(t(q[,,r]),scale=F); qhat_cent = t(qhat_cent)
        alpha_i = matrix(0,num_comp,N)
        for (ii in 1:num_comp){
            for (jj in 1:N){
                alpha_i[ii,jj] = simpson(time,qhat_cent[,jj]*U[,ii])
            }
        }

        tmp = U[,1:num_comp]%*%alpha_i
        mq_c = mq[,r]
        if (smooth_data){
            mq_c = mq[,r]/pvecnorm(mq[,r],2)
            for (ii in 1:N){
                if (sum(tmp[,ii]) != 0)
                    tmp[,ii] = tmp[,ii]/pvecnorm(tmp[,ii],2)
            }
        }
        qhat = matrix(mq_c,M,N) + tmp

        # Matching Step
        gam_k<-foreach(k = 1:N, .combine=cbind,.packages="fdasrvf") %dopar% {
            gam0 = optimum.reparam(qhat[,k],time,q[,k,r],time)
        }
        gam[,,r] = gam_k

        psi =  sqrt(diff(gam_k)*(M-1))
        mu = rep(1,M-1)
        Dx1 = rep(0,N)
        for (ii in 1:N){
            Dx1[ii] = Re(acos(sum(mu*psi[,ii])/(M-1)))
        }
        Dx[r] = max(Dx1)

        if (abs(Dx[r]-Dx[r-1]) < 1e-4 || r >=MaxItr){
            break
        }
    }

    # Aligned data & stats
    fn = f[,,r]
    qn = q[,,r]
    q0 = q[,,1]
    mean_f0 = rowMeans(f[,,1]);
    std_f0 = apply(f[,,1], 1, sd)
    mqn = mq[,r]
    gamf = gam[,,1]
    for (k in 2:r){
        gam_k = gam[,,k]
        for (l in 1:N){
            gamf[,l] = approx(time,gamf[,l],xout=(time[length(time)]-time[1])*gam_k[,l] + time[1])$y
        }
    }

    # Center Mean
    gamI = SqrtMeanInverse(t(gamf))
    gamI_dev = gradient(gamI, 1/(M-1))
    mqn = approx(time,mqn,xout=(time[length(time)]-time[1])*gamI + time[1])$y*sqrt(gamI_dev)

    for (k in 1:N){
        qn[,k] = approx(time,qn[,k],xout=(time[length(time)]-time[1])*gamI + time[1])$y*sqrt(gamI_dev)
        fn[,k] = approx(time,fn[,k],xout=(time[length(time)]-time[1])*gamI + time[1])$y
        gamf[,k] = approx(time,gamf[,k],xout=(time[length(time)]-time[1])*gamI + time[1])$y
    }

    mean_fn = rowMeans(fn)
    std_fn = apply(fn, 1, sd)

    # Get Final PCA
    mq_new = rowMeans(qn)
    m_new = sign(fn[round(length(time)/2),])*sqrt(abs(fn[round(length(time)/2),]))  # scaled version
    mqn2 = c(mq_new,mean(m_new))
    K = cov(t(rbind(qn,m_new)))
    out = svd(K)
    s = out$d
    stdS = sqrt(s)
    U = out$u

    # compute the PCA in the q domain
    q_pca = array(0,dim=c((length(mq_new)+1),Nstd,num_comp))
    for (k in NP){
        for (i in 1:Nstd){
            q_pca[,i,k] = mqn2 + coef[i]*stdS[k]*U[,k]
        }
    }

    # compute the correspondence to the original function domain
    f_pca = array(0,dim=c((length(mqn)),Nstd,num_comp))
    for (k in NP){
        for (i in 1:Nstd){
            f_pca[,i,k] = cumtrapzmid(time,q_pca[1:(dim(q_pca)[1]-1),i,k]*abs(q_pca[1:(dim(q_pca)[1]-1),i,k]),sign(q_pca[dim(q_pca)[1],i,k])*(q_pca[dim(q_pca)[1],i,k]^2))
        }
    }

    N2 = dim(qn)[2]
    c = matrix(0,N2,num_comp)
    for (k in NP){
        for (i in 1:N2){
            c[i,k] = sum((c(qn[,i],m_new[i])-mqn2)*U[,k])
        }
    }

    vfpca = list()
    vfpca$q_pca = q_pca
    vfpca$f_pca = f_pca
    vfpca$latent = s
    vfpca$coef = c
    vfpca$U = U

    if (showplot){
        matplot((0:(M-1))/(M-1),gamf,type="l",main="Warping functions",xlab="Time")

        matplot(time,fn,type="l",main=bquote(paste("Warped Data ",lambda == .(lambda))))

        matplot(time,cbind(mean_f0,mean_f0+std_f0,mean_f0-std_f0),type="l",lty=1,col=c("blue","red","green"),
                        ylab="",main=bquote(paste("Original Data: ", Mean %+-% STD)))
        legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),col=c('blue','red','green'),lty=1)

        matplot(time,cbind(mean_fn,mean_fn+std_fn,mean_fn-std_fn),type="l",lty=1,col=c("blue","red","green"),
                        ylab="",main=bquote(paste("Warped Data: ",lambda == .(lambda),": ", Mean %+-% STD)))
        legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),col=c('blue','red','green'),lty=1)

        if (num_comp > 3)
            num_comp = 3
        layout(matrix(c(1:(2*num_comp)), 2, num_comp, byrow = TRUE))
        dims = dim(q_pca)
        for (ii in 1:num_comp){
            matplot(time,q_pca[1:(dim(q_pca)[1]-1),,ii],type="l")
            title(main=sprintf("q domain: PD %d", ii))
        }
        for (ii in 1:num_comp){
            matplot(time,f_pca[,,ii],type="l")
            title(main=sprintf("f domain: PD %d", ii))
        }

        layout(1)
        cumm_coef = 100*cumsum(s)/sum(s)
        plot(cumm_coef,type="l",col="blue",main="Coefficient Cumulative Percentage", ylab = "Percentage")
    }

    if (parallel){
        stopCluster(cl)
    }

    return(list(f0=f[,,1],fn=fn,qn=qn,q0=q0,mqn=mqn,gam=gamf,vfpca=vfpca,Dx=Dx))
}
