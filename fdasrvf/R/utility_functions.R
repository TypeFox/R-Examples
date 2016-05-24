ndims <- function(x){
    return(length(dim(x)))
}

meshgrid <- function(x, y = x) {
    if (!is.numeric(x) || !is.numeric(y))
        stop("Arguments 'x' and 'y' must be numeric vectors.")

    x <- c(x); y <- c(y)
    n <- length(x)
    m <- length(y)

    X <- matrix(rep(x, each = m),  nrow = m, ncol = n)
    Y <- matrix(rep(y, times = n), nrow = m, ncol = n)

    return(list(X = X, Y = Y))
}

gradient2 <- function(a,dx=1,dy=1){
    m = dim(a)[1]
    n = dim(a)[2]
    dxdu = matrix(0,m,n)
    dydv = matrix(0,m,n)

    for (i in 1:m) {
        dxdu[i,] = gradient(as.vector(a[i,]), dx)
    }

    for (i in 1:m) {
        dydv[,i] = gradient(as.vector(a[,i]), dy)
    }

    return(list(dxdu=dxdu,dydv=dydv))
}

cumtrapz <- function(x,y,dims=1){
    if ((dims-1)>0){
        perm = c(dims:max(ndims(y),dims), 1:(dims-1))
    } else {
        perm = c(dims:max(ndims(y),dims))
    }

    if (ndims(y) == 0){
        n = 1
        m = length(y)
    } else {
        if (length(x) != dim(y)[dims])
            stop('Dimension Mismatch')
        y = aperm(y, perm)
        m = nrow(y)
        n = ncol(y)
    }

    if (n==1){
        dt = diff(x)/2.0
        z = c(0, cumsum(dt*(y[1:(m-1)] + y[2:m])))
        dim(z) = c(m,1)
    } else {
        tmp = diff(x)
        dim(tmp) = c(m-1,1)
        dt = repmat(tmp/2.0,1,n)
        z = rbind(rep(0,n), apply(dt*(y[1:(m-1),] + y[2:m,]),2,cumsum))
        perm2 = rep(0, length(perm))
        perm2[perm] = 1:length(perm)
        z = aperm(z, perm2)
    }

    return(z)
}

trapz <- function(x,y,dims=1){
    if ((dims-1)>0){
        perm = c(dims:max(ndims(y),dims), 1:(dims-1))
    } else {
        perm = c(dims:max(ndims(y),dims))
    }

    if (ndims(y) == 0){
        m = 1
    } else {
        if (length(x) != dim(y)[dims])
            stop('Dimension Mismatch')
        y = aperm(y, perm)
        m = nrow(y)
    }

    if (m==1){
        M = length(y)
        out = sum(diff(x)*(y[-M]+y[-1])/2)
    } else {
        slice1 = y[as.vector(outer(1:(m-1), dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+')) ] 
        dim(slice1) = c(m-1, length(slice1)/(m-1)) 
        slice2 = y[as.vector(outer(2:m, dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+'))] 
        dim(slice2) = c(m-1, length(slice2)/(m-1)) 
        out = t(diff(x)) %*% (slice1+slice2)/2.
        siz = dim(y)
        siz[1] = 1
        out = array(out, siz)
        perm2 = rep(0, length(perm))
        perm2[perm] = 1:length(perm)
        out = aperm(out, perm2)
        ind = which(dim(out) != 1)
        out = array(out, dim(out)[ind])
    }

    return(out)
}

simpson <- function(x,y){
    M = nrow(y)
    if (is.null(M)){
        M = length(y)
        if (M < 3){
            out = trapz(x,y)
        }else{
            dx = diff(x)
            dx1 = dx[1:(length(dx)-1)]
            dx2 = dx[2:length(dx)]
            alpha = (dx1+dx2)/dx1/6
            a0 = alpha*(2*dx1-dx2)
            a1 = alpha*(dx1+dx2)^2/dx2
            a2 = alpha*dx1/dx2*(2*dx2-dx1)

            out = sum(a0[seq(1,length(a0),2)]*y[seq(1,M-2,2)] + a1[seq(1,length(a1),2)]*y[seq(2,M-1,2)]+a2[seq(1,length(a2),2)]*y[seq(3,M,2)])
            if (M %% 2 == 0){
                A = vandermonde.matrix(x[(length(x)-2):length(x)],3)
                C = solve(A[,3:1],y[(length(y)-2):length(y)])
                out = out + C[1]*(x[length(x)]^3-x[(length(x)-1)]^3)/3 + C[2]*(x[length(x)]^2-x[(length(x)-1)]^2)/2 + C[3]*dx[length(dx)]
            }
        }
    }else{
        M = nrow(y)
        N = ncol(y)

        # use  trapz if M < 3
        if (M < 3){
            out = trapz(x,y)
        }else{
            out = rep(0,N)
            dx = diff(x)
            dx1 = dx[1:(length(dx)-1)]
            dx2 = dx[2:length(dx)]
            alpha = (dx1+dx2)/dx1/6
            a0 = alpha*(2*dx1-dx2)
            a1 = alpha*(dx1+dx2)^2/dx2
            a2 = alpha*dx1/dx2*(2*dx2-dx1)
            for (i in 1:N){
                out[i] = sum(a0[seq(1,length(a0),2)]*y[seq(1,M-2,2),i] + a1[seq(1,length(a1),2)]*y[seq(2,M-1,2),i]+a2[seq(1,length(a2),2)]*y[seq(3,M,2),i])
                if (M %% 2 == 0){
                    A = vandermonde.matrix(x[(length(x)-2):length(x)],3)
                    C = solve(A[,3:1],y[(length(y)-2):length(y)])
                    out[i] = out[i] + C[1]*(x[length(x)]^3-x[(length(x)-1)]^3)/3 + C[2]*(x[length(x)]^3-x[(length(x)-1)]^2)/2 + C[3]*dx[length(dx)]
                }
            }
        }
    }

    return(out)
}

cumtraps <- function(x,y){
    M = length(y)
    dx = diff(x)
    dx1 = dx[1:(length(dx)-1)]
    dx2 = dx[2:length(dx)]
    alpha = (dx1+dx2)/dx1/6
    a0 = alpha*(2*dx1-dx2)
    a1 = alpha*(dx1+dx2)^2/dx2
    a2 = alpha*dx1/dx2*(2*dx2-dx1)

    A = vandermonde.matrix(x[1:3],3)
    C = solve(A[,3:1],y[1:3])
    z = rep(0,M)
    z[2] = C[1]*(x[2]^3-x[1]^3)/3 + C[2]*(x[2]^2-x[1]^2)/2 + C[3]*dx[1]
    z[seq(3,length(z),2)] = cumsum(a0[seq(1,length(a0),2)]*y[seq(1,M-2,2)] + a1[seq(1,length(a1),2)]*y[seq(2,M-1,2)] + a2[seq(1,length(a1),2)]*y[seq(3,M,2)])
    z[seq(4,length(z),2)] = cumsum(a0[seq(2,length(a0),2)]*y[seq(2,M-2,2)] + a1[seq(2,length(a1),2)]*y[seq(3,M-1,2)] + a2[seq(2,length(a1),2)]*y[seq(4,M,2)])+z[2]


    return(z)
}

pvecnorm <-function(v,p=2){
    sum(abs(v)^p)^(1/p)
}

pvecnorm2 <-function(dt,x){
    sqrt(sum(abs(x)*abs(x))*dt)
}

gradient.spline <- function(f,binsize,smooth_data=F){
    if (smooth_data==TRUE){
        n = nrow(f)
        if (is.null(n)){
            N = 1
            tmp.spline = smooth.spline(f)
            f.out = tmp.spline$y
            g = predict(tmp.spline,deriv=1)$y/binsize
        }else{
            N = ncol(f)
            f.out = matrix(0,nrow(f),ncol(f))
            g = matrix(0,nrow(f),ncol(f))
            for (jj in 1:N){
                tmp.spline = smooth.spline(f[,jj])
                f.out[,jj] = tmp.spline$y
                g[,jj] = predict(tmp.spline,deriv=1)$y/binsize
            }
        }
    }else{
        g = gradient(f,binsize)
        f.out = f
    }

    return(list(g=g, f=f.out))
}

SqrtMeanInverse <- function(gam){
    n = nrow(gam)
    T1 = ncol(gam)
    dt = 1/(T1-1)
    psi = matrix(0,n,T1-1)
    for (i in 1:n){
        psi[i,] = sqrt(diff(gam[i,])/dt+.Machine$double.eps)
    }

    # Find direction
    mnpsi = colMeans(psi)
    dqq = sqrt(colSums((t(psi) - matrix(mnpsi,ncol=n,nrow=T1-1))^2))
    min_ind = which.min(dqq)
    mu = psi[min_ind,]
    tt = 1
    maxiter = 20
    eps = .Machine$double.eps
    lvm = rep(0,1,maxiter)
    vec = matrix(0,n,T1-1)
    for (iter in 1:maxiter){
        for (i in 1:n){
            v = psi[i,] - mu
            dot<- simpson(seq(0,1,length.out=T1-1),mu*psi[i,])
            dot.limited<- ifelse(dot>1, 1, ifelse(dot<(-1), -1, dot))
            len = acos(dot.limited)
            if (len > 0.0001){
                vec[i,] = (len/sin(len))*(psi[i,] - cos(len)*mu)
            }else{
                vec[i,] = rep(0,T1-1)
            }
        }
        vm = colMeans(vec)
        lvm[iter] = sqrt(sum(vm*vm)*dt)
        if (sum(vm) == 0){ # we had a problem, pick id (i.e., they are all id)
            mu = rep(1,T1-1)
            break
        }else{
            mu = cos(tt*lvm[iter])*mu + (sin(tt*lvm[iter])/lvm[iter])*vm
            if (lvm[iter] < 1e-6 || iter >=maxiter){
                break
            }
        }

    }
    gam_mu = c(0,cumsum(mu*mu))/T1
    gam_mu = (gam_mu - min(gam_mu))/(max(gam_mu)-min(gam_mu))
    gamI = invertGamma(gam_mu)
    return(gamI)
}

invertGamma <- function(gam){
    N = length(gam)
    x = (0:(N-1))/(N-1)
    gamI = approx(gam,x,xout=x)$y
    gamI[N] = 1
    gamI = gamI/gamI[N]
    return(gamI)
}

randomGamma <- function(gam,num){
    out = SqrtMean(gam)
    mu = out$mu
    psi = out$psi
    vec = out$vec

    K = cov(t(vec))
    out = svd(K)
    s = out$d
    U = out$u
    n = 5
    TT = nrow(vec) + 1
    vm = rowMeans(vec)

    rgam = matrix(0,num,TT)
    for (k in 1:num){
        a = rnorm(n)
        v = rep(0,length(vm))
        for (i in 1:n){
            v = v + a[i]*sqrt(s[i])*U[,i]
        }
        vn = pvecnorm(v,2)/sqrt(TT)
        psi = cos(vn)*mu + sin(vn)*v/vn
        tmp = rep(0,TT)
        tmp[2:TT] = cumsum(psi*psi)
        rgam[k,] = tmp/TT
    }
    return(rgam)
}

f_K_fold <- function(Nobs,K=5){
    rs <- runif(Nobs)
    id <- seq(Nobs)[order(rs)]
    k <- as.integer(Nobs*seq(1,K-1)/K)
    k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
    k[,1] <- k[,1]+1
    l <- lapply(seq.int(K),function(x,k,d)
        list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
                 test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
    return(l)
}

cov_samp <- function(x,y=NULL){
    x = scale(x,scale=F)
    N = dim(x)[1]
    if (length(y) == 0){
        sigma = 1/N * t(x) %*% x
    }else{
        y = scale(y,scale=F)
        sigma = 1/N * t(x) %*% y
    }

    return(sigma)
}

diffop <- function(n,binsize = 1){
    m = matrix(0,nrow=n,ncol=n)
    diag(m[-1,]) <- 1
    diag(m) <- -2
    diag(m[,-1]) <- 1
    m = t(m) %*% m
    m[1,1] = 6
    m[n,n] = 6
    m = m/(binsize^4)
    return(m)
}

geigen <- function (Amat, Bmat, Cmat)
{
    Bdim <- dim(Bmat)
    Cdim <- dim(Cmat)
    if (Bdim[1] != Bdim[2])
        stop("BMAT is not square")
    if (Cdim[1] != Cdim[2])
        stop("CMAT is not square")
    p <- Bdim[1]
    q <- Cdim[1]
    s <- min(c(p, q))
    if (max(abs(Bmat - t(Bmat)))/max(abs(Bmat)) > 1e-10)
        stop("BMAT not symmetric.")
    if (max(abs(Cmat - t(Cmat)))/max(abs(Cmat)) > 1e-10)
        stop("CMAT not symmetric.")
    Bmat <- (Bmat + t(Bmat))/2
    Cmat <- (Cmat + t(Cmat))/2
    Bfac <- chol(Bmat)
    Cfac <- chol(Cmat)
    Bfacinv <- solve(Bfac)
    Cfacinv <- solve(Cfac)
    Dmat <- t(Bfacinv) %*% Amat %*% Cfacinv
    if (p >= q) {
        result <- svd2(Dmat)
        values <- result$d
        Lmat <- Bfacinv %*% result$u
        Mmat <- Cfacinv %*% result$v
    }
    else {
        result <- svd2(t(Dmat))
        values <- result$d
        Lmat <- Bfacinv %*% result$v
        Mmat <- Cfacinv %*% result$u
    }
    geigenlist <- list(values, Lmat, Mmat)
    names(geigenlist) <- c("values", "Lmat", "Mmat")
    return(geigenlist)
}

svd2 <- function (x, nu = min(n, p), nv = min(n, p), LINPACK = FALSE)
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    svd.x <- try(svd(x, nu, nv, LINPACK))
    if (class(svd.x) == "try-error") {
        nNA <- sum(is.na(x))
        nInf <- sum(abs(x) == Inf)
        if ((nNA > 0) || (nInf > 0)) {
            msg <- paste("sum(is.na(x)) = ", nNA, "; sum(abs(x)==Inf) = ",
                                     nInf, ".  'x stored in .svd.x.NA.Inf'", sep = "")
            stop(msg)
        }
        attr(x, "n") <- n
        attr(x, "p") <- p
        attr(x, "LINPACK") <- LINPACK
        .x2 <- c(".svd.LAPACK.error.matrix", ".svd.LINPACK.error.matrix")
        .x <- .x2[1 + LINPACK]
        msg <- paste("svd failed using LINPACK = ", LINPACK,
                                 " with n = ", n, " and p = ", p, ";", sep = "")
        warning(msg)
        svd.x <- try(svd(x, nu, nv, !LINPACK))
        if (class(svd.x) == "try-error") {
            .xc <- .x2[1 + (!LINPACK)]
            stop("svd also failed using LINPACK = ", !LINPACK)
        }
    }
    svd.x
}

cumtrapzmid <- function(x,y,c){
    a = length(x)
    mid = round(a/2)

    # case < mid
    fn = rep(0,a)
    tmpx = x[seq(mid-1,1,-1)]
    tmpy = y[seq(mid-1,1,-1)]
    tmp = c + cumtrapz(tmpx,tmpy)
    fn[1:mid-1] = rev(tmp)

    # case >= mid
    fn[mid:a] = c + cumtrapz(x[mid:a],y[mid:a])

    return(fn)
}

warp_q_gamma <- function(time, q, gam){
  M = length(gam)
  gam_dev = gradient(gam, 1/(M-1))
  q_tmp = approx(time,q,xout=(time[length(time)]-time[1])*gam +
               time[1])$y*sqrt(gam_dev)
  return(q_tmp)
}

zero_crossing <- function(Y, q, bt, time, y_max, y_min, gmax, gmin){
  # finds zero-crossing of optimal gamma, gam = s*gmax + (1-s)*gmin
  # from elastic regression model
  max_itr = 100
  a = rep(0, max_itr)
  a[1] = 1
  f = rep(0, max_itr)
  f[1] = y_max - Y
  f[2] = y_min - Y
  mrp = f[1]
  mrn = f[2]
  mrp_ind = 1  # most recent positive index
  mrn_ind = 2  # most recent negative index

  for (ii in 3:max_itr){
    x1 = a[mrp_ind]
    x2 = a[mrn_ind]
    y1 = mrp
    y2 = mrn
    a[ii] = (x1 *y2 - x2 * y1) / (y2-y1)
    gam_m = a[ii] * gmax + (1-a[ii]) * gmin
    qtmp = warp_q_gamma(time, q, gam_m)
    f[ii] = trapz(time, qtmp*bt) - Y

    if (abs(f[ii]) < 1e-5){
      break
    } else if (f[ii]> 0){
      mrp = f[ii]
      mrp_ind = ii
    } else{
      mrn = f[ii]
      mrn_ind = ii
    }
  }

  gamma = a[ii] * gmax + (1-a[ii]) * gmin

  return(gamma)
}

repmat <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  if (is.null(mx)){
    mx = 1
    nx = length(X)
    mat = matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }else {
    nx = dim(X)[2]
    mat = matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }

  return(mat)
}

simul_align <- function(f1, f2){
    # parameterize by arc-length
    s1 = arclength(f1)
    s2 = arclength(f2)

    len1 = max(s1)
    len2 = max(s2)

    f1 = f1/len1
    s1 = s1/len1
    f2 = f2/len2
    s2 = s2/len2

    # get srvf (should be +/-1)
    q1 = diff(f1)/diff(s1)
    q1[diff(s1)==0] = 0;
    q2 = diff(f2)/diff(s2);
    q2[diff(s2)==0] = 0;

    # get extreme points
    out = extrema_1s(s1, q1);
    ext1 = out$ext2
    d1 = out$d
    out = extrema_1s(s2,q2);
    ext2 = out$ext2
    d2 = out$d

    out = match_ext(s1,ext1,d1,s2,ext2,d2)
    D = out$D
    P = out$P
    mpath = out$mpath

    te1 = s1[ext1]
    te2 = s2[ext2]

    fe1 = f1[ext1]
    fe2 = f2[ext2]

    out = simul_reparam(te1, te2, mpath)

    g1 = out$g1
    g2 = out$g2

    return(list(s1=s1,s2=s2,g1=g1,g2=g2,ext1=ext1,ext2=ext2,mpath=mpath))
}

arclength <- function(f){
    t1 = rep(0,length(f))
    t1[2:length(f)] = abs(diff(f))
    t1 = cumsum(t1)

    return(t1)
}

extrema_1s <- function(t, q){
    q = round(q)

    if (q[1] != 0){
        d = -q[1]
    } else{
        d = q[q!=0]
        d = d[1]
    }

    ext = which(diff(q) != 0) + 1
    ext2 = rep(0, length(ext)+2)
    ext2[1] = 1
    ext2[2:(length(ext2)-1)] = round(ext)
    ext2[length(ext2)] = length(t)

    return(list(ext2=ext2,d=d))
}

match_ext <- function(t1, ext1, d1, t2, ext2, d2){
    te1 = t1[ext1];
    te2 = t2[ext2];

    # We'll pad each sequence to start on a 'peak' and end on a 'valley'
    pad1 = rep(0, 2)
    pad2 = rep(0, 2)
    if (d1 == -1){
        te1a = rep(0, length(te1)+1)
        te1a[2:length(te1a)] = te1
        te1a[1] = te1a[2]
        te1 = te1a
        pad1[1] = 1
    }

    if ((length(te1)%%2)==1){
        te1a = rep(0,length(te1)+1)
        te1a[1:length(te1a)-1] = te1
        te1a[length(te1a)] = te1[length(te1)]
        te1 = te1a
        pad1[2] = 1
    }

    if (d2 == -1){
        te2a = rep(0, length(te2)+1)
        te2a[2:length(te2a)] = te2
        te2a[1] = te2a[2]
        te2 = te2a
        pad2[1] = 1
    }

    if ((length(te2)%%2)==1){
        te2a = rep(0,length(te2)+1)
        te2a[1:length(te2a)-1] = te2
        te2a[length(te2a)] = te2[length(te2)]
        te2 = te2a
        pad2[2] = 1
    }

    n1 = length(te1)
    n2 = length(te2)

    # initialize weight and path matrices
    D = matrix(0,n1,n2)
    P = array(0,c(n1,n2,2))

    for (i in 1:n1){
        for (j in 1:n2){
            if (((i+j)%%2)==0){
                if ((i-1)>=(1+(i%%2))){
                    for (ib in seq(i-1,1+(i%%2),by=-2)){
                        if ((j-1)>=(1+(j%%2))){
                            for (jb in seq(j-1,1+(j%%2),by=-2)){
                                icurr = seq(ib+1,i,2)
                                jcurr = seq(jb+1,j,2)
                                W = sqrt(sum(te1[icurr]-te1[icurr-1])) *
                                    sqrt(sum(te2[jcurr]-te2[jcurr-1]))
                                Dcurr = D[ib,jb] + W
                                if (Dcurr > D[i,j]){
                                    D[i,j] = Dcurr
                                    P[i,j,] = c(ib,jb)
                                }
                            }

                        }
                  }

                }
            }
        }
    }

    D = D[(1+pad1[1]):(n1-pad1[2]), (1+pad2[1]):(n2-pad2[2])]
    P = P[(1+pad1[1]):(n1-pad1[2]), (1+pad2[1]):(n2-pad2[2]),]
    P[,,1] = P[,,1]-pad1[1]
    P[,,2] = P[,,2]-pad2[1]

    # Retrieve Best Path
    if (pad1[2] == pad2[2]){
        mpath = dim(D)
    } else if (D[nrow(D)-1,ncol(D)] > D[nrow(D),ncol(D)-1]){
        mpath = dim(D) - c(1,0)
    } else {
        mpath = dim(D) - c(0,1)
    }
    mpath = round(mpath)
    P = round(P)
    prev_vert = P[mpath[1],mpath[2],]

    while (any(prev_vert>0)){
        mpath = rbind(prev_vert,mpath,deparse.level=0)
        prev_vert = P[mpath[1,1],mpath[1,2],]
    }

    return(list(D=D,P=P,mpath=mpath))
}

simul_reparam <- function(te1, te2, mpath){
    g1 = 0
    g2 = 0

    if (mpath[1,2] == 2){
        g1 = c(g1,0)
        g2 = c(g2,te2[2])
    }else if (mpath[1,1] == 2){
        g1 = c(g1, te1[2])
        g2 = c(g2,0)
    }

    m = nrow(mpath)
    for (ii in 1:(m-1)){
        out = simul_reparam_segment(mpath[ii,], mpath[ii+1,],te1,te2)

        g1 = c(g1,out$gg1)
        g2 = c(g2,out$gg2)
    }

    n1 = length(te1)
    n2 = length(te2)

    if ((mpath[nrow(mpath),1]==(n1-1)) || (mpath[nrow(mpath),2] == (n2-1))){
        g1 = c(g1,1)
        g2 = c(g2,1)
    }

    return(list(g1=g1,g2=g2))
}

simul_reparam_segment <- function(src, tgt, te1, te2){
    i1 = seq(src[1]+1,tgt[1],2)
    i2 = seq(src[2]+1,tgt[2],2)

    v1 = sum(te1[i1]-te1[i1-1])
    v2 = sum(te2[i2]-te2[i2-1])
    R = v2/v1

    a1 = src[1]
    a2 = src[2]
    t1 = te1[a1]
    t2 = te2[a2]
    u1 = 0.
    u2 = 0.

    gg1 = vector()
    gg2 = vector()

    while ((a1<tgt[1]) && (a2<tgt[2])){
        if ((a1==(tgt[1]-1)) && (a2==(tgt[2]-1))){
            a1 = tgt[1]
            a2 = tgt[2]
            gg1 = c(gg1, te1[a1])
            gg2 = c(gg2, te2[a2])
        } else {
            p1 = (u1+te1[a1+1]-t1)/v1
            p2 = (u2+te2[a2+1]-t2)/v2
            if (p1<p2){
                lam = t2+R*(te1[a1+1]-t1)
                gg1 = c(gg1,te1[a1+1],te1[a1+2])
                gg2 = c(gg2, lam, lam)
                u1 = u1 + te1[a1+1] - t1
                u2 = u2 + lam - t2
                t1 = te1[a1+2]
                t2 = lam
                a1 = a1 + 2
            } else {
                lam = t1 + (1./R)*(te2[a2+1]-t2)
                gg1 = c(gg1,lam,lam)
                gg2 = c(gg2, te2[a2+1], te2[a2+2])
                u1 = u1 + lam - t1
                u2 = u2 + te2[a2+1] - t2
                t1 = lam
                t2 = te2[a2+2]
                a2 = a2 + 2
            }
        }
    }

    return(list(gg1=gg1,gg2=gg2))
}

simul_gam <- function(u,g1,g2,t,s1,s2,tt){
    gs1 = interp1_flat(u,g1,tt)
    gs2 = interp1_flat(u,g2,tt)

    gt1 = interp1_flat(s1,t,gs1)
    gt2 = interp1_flat(s2,t,gs2)

    gam = interp1_flat(gt1,gt2,tt)

    return(gam)
}

interp1_flat <- function(x,y,xx){
    flat = which(diff(x) <=0)
    n = length(flat)

    if (n==0){
        yy = approx(x,y,xx)$y
    } else {
        yy = rep(0,length(xx))
        i1 = 1
        if (flat[1] == 1){
            i2 = 1
            j = xx==x[i2]
            yy[j] = min(y[i2:i2+1])
        } else {
            i2 = flat[1]
            j = (xx>=x[i1]) & (xx<=x[i2])
            yy[j] = approx(x[i1:i2],y[i1:i2],xx[j])$y
            i1 = i2
        }
        if (n>1){
          for (k in 2:n){
            i2 = flat[k]
            if (i2 > (i1+1)){
              j = (xx>=x[i1]) & (xx<=x[i2])
              yy[j] = approx(x[i1+1:i2],y[i1+1:i2],xx[j])$y
            }
            j = xx == x[i2]
            yy[j] = min(y[i2:i2+1])
            i1 = i2
          }
        }

        i2 = length(x)
        j = (xx>=x[i1]) & (xx<=x[i2])
        if ((i1+1) == i2){
            yy[j] = y[i2]
        } else {
            yy[j] = approx(x[i1+1:i2],y[i1+1:i2],xx[j])$y
        }
    }

    return(yy)
}


circshift <- function(a, sz) {
    if (is.null(a)) return(a)

    if (is.vector(a) && length(sz) == 1) {
        n <- length(a)
        s <- sz %% n
        a <- a[(1:n-s-1) %% n + 1]

    } else if (is.matrix(a) && length(sz) == 2) {
        n <- nrow(a); m <- ncol(a)
        s1 <- sz[1] %% n
        s2 <- sz[2] %% m
        a <- a[(1:n-s1-1) %% n + 1, (1:m-s2-1) %% m + 1]
    } else
        stop("Length of 'sz' must be equal to the no. of dimensions of 'a'.")

    return(a)
}
