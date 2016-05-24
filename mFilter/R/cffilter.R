### Christiano-Fitzgerald filter
cffilter <- function(x,pl=NULL,pu=NULL,root=FALSE,drift=FALSE,
                     type=c("asymmetric","symmetric","fixed","baxter-king",
                     "trigonometric"),nfix=NULL,theta=1)
{

    type = match.arg(type)
    if(is.null(root)) root <- FALSE
    if(is.null(drift)) drift <- FALSE
    if(is.null(theta)) theta <- 1
    if(is.null(type)) type <- "asymmetric"

    if(is.ts(x))
        freq=frequency(x)
    else
        freq=1

    if(is.null(pl))
    {
        if(freq > 1)
            pl=trunc(freq*1.5)
        else
            pl=2
    }

    if(is.null(pu))
        pu=trunc(freq*8)

    if(is.null(nfix))
        nfix = freq*3

    nq=length(theta)-1;
    b=2*pi/pl;
    a=2*pi/pu;

    xname=deparse(substitute(x))
    xo = x
    x = as.matrix(x)
    n = nrow(x)
    nvars = ncol(x)

    if(n < 5)
        warning("# of observations < 5")

    if(n < (2*nq+1))
        stop("# of observations must be at least 2*q+1")

    if(pu <= pl)
        stop("pu must be larger than pl")

    if(pl < 2)
    {
        warning("pl less than 2 , reset to 2")
        pl = 2
    }

    if(root != 0 && root != 1)
        stop("root must be 0 or 1")

    if(drift<0 || drift > 1)
        stop("drift must be 0 or 1")

    if((type == "fixed" || type == "baxter-king") && nfix < 1)
        stop("fixed lag length must be >= 1")

    if(type == "fixed" & nfix < nq)
        stop("fixed lag length must be >= q")

    if((type == "fixed" || type == "baxter-king") && nfix >= n/2)
        stop("fixed lag length must be < n/2")

    if(type == "trigonometric" && (n-2*floor(n/2)) != 0)
        stop("trigonometric regressions only available for even n")

    theta = as.matrix(theta)
    m1 = nrow(theta)
    m2 = ncol(theta)
    if(m1 > m2)
        th=theta
    else
        th=t(theta)

    ##   compute g(theta)
    ##   [g(1) g(2) .... g(2*nq+1)] correspond to [c(q),c(q-1),...,c(1),
    ##                                        c(0),c(1),...,c(q-1),c(q)]
    ##   cc = [c(0),c(1),...,c(q)]
    ## ?? thp=flipud(th)
    g=convolve(th,th,type="open")
    cc = g[(nq+1):(2*nq+1)]
    ##   compute "ideal" Bs
    j = 1:(2*n)
    B = as.matrix(c((b-a)/pi, (sin(j*b)-sin(j*a))/(j*pi)))
    ##    compute R using closed form integral solution
    R = matrix(0,n,1)
    if(nq > 0)
    {
        R0 = B[1]*cc[1] + 2*t(B[2:(nq+1)])*cc[2:(nq+1)]
        R[1] = pi*R0
        for(i in 2:n)
        {
            dj = Bge(i-2,nq,B,cc)
            R[i] = R[i-1] - dj
        }
    }
    else
    {
        R0 = B[1]*cc[1]
        R[1] = pi*R0;
        for(j in 2:n)
        {
            dj = 2*pi*B[j-1]*cc[1];
            R[j] = R[j-1] - dj;
        }
    }

    AA = matrix(0,n,n)

###  asymmetric filter

    if(type == "asymmetric")
    {
        if(nq==0)
        {
            for(i in 1:n)
            {
                AA[i,i:n] = t(B[1:(n-i+1)])
                if(root)
                    AA[i,n] = R[n+1-i]/(2*pi)
            }
            AA[1,1] = AA[n,n]
            ##  Use symmetry to construct bottom 'half' of AA
            AAu = AA
            AAu[!upper.tri(AAu)] <- 0
            AA = AA + flipud(fliplr(AAu))
        }
        else
        {
            ## CONSTRUCT THE A MATRIX size n x n
            A = Abuild(n,nq,g,root)
            Ainv = solve(A)
            ## CONSTRUCT THE d MATRIX size n x 1
            for(np in 0:ceiling(n/2-1))
            {
                d = matrix(0,n,1)
                ii = 0
        	for(jj in (np-root):(np+1+root-n))
                {
                    ii = ii+1
                    d[ii] = Bge(jj,nq,B,cc)
                }
                if (root == 1)
                    d[n-1] = R[n-np]
                ##  COMPUTE Bhat = inv(A)*d
                Bhat = Ainv%*%d
                AA[np+1,] = t(Bhat)
            }
            ##  Use symmetry to construct bottom 'half' of AA
            AA[(ceiling(n/2)+1):n,] = flipud(fliplr(AA[1:floor(n/2),]))
        }
    }


###  symmetric filter

    if (type == "symmetric")
    {
        if(nq==0)
        {
            for(i in 2:ceiling(n/2))
            {
                np = i-1
                AA[i,i:(i+np)] = t(B[1:(1+np)])
                if(root)
                    AA[i,i+np] = R[np+1]/(2*pi);
                AA[i,(i-1):(i-np)] = AA[i,(i+1):(i+np)];
            }
            ##  Use symmetry to construct bottom 'half' of AA
            AA[(ceiling(n/2)+1):n,] = flipud(fliplr(AA[1:floor(n/2),]))
        }
        else
        {
            for(np in nq:ceiling(n/2-1))
            {
                nf = np
                nn = 2*np+1
                ## CONSTRUCT THE A MATRIX size nn x nn
                A = Abuild(nn,nq,g,root)
                Ainv = solve(A)
                ## CONSTRUCT THE d MATRIX size nn x 1
                d = matrix(0,nn,1)
                ii = 0
                for(jj in (np-root):(-nf+root))
                {
                    ii = ii+1
                    d[ii] = Bge(jj,nq,B,cc)
                }
                if(root)
                    d[nn-1] = R[nf+1]
                ##  COMPUTE Bhat = inv(A)*d
                Bhat = Ainv%*%d
                AA[np+1,1:(2*np+1)] = t(Bhat)
            }
            ##  Use symmetry to construct bottom 'half' of AA
            AA[(ceiling(n/2)+1):n,] = flipud(fliplr(AA[1:floor(n/2),]))
        }
    }


###  fixed length symmetric filter

    if (type == "fixed")
    {
        if(nq==0)
        {
            bb = matrix(0,2*nfix+1,1)
            bb[(nfix+1):(2*nfix+1)] = B[1:(nfix+1)]
            bb[nfix:1] = B[2:(nfix+1)]
            if(root)
            {
                bb[2*nfix+1] = R[nfix+1]/(2*pi)
                bb[1] = R[nfix+1]/(2*pi)
            }
            for(i in (nfix+1):(n-nfix))
                AA[i,(i-nfix):(i+nfix)] = t(bb)
        }
        else
        {
            nn = 2*nfix+1
            ## CONSTRUCT THE A MATRIX size nn x nn
            A = Abuild(nn,nq,g,root)
            Ainv = solve(A)
            ## CONSTRUCT THE d MATRIX size nn x 1
            d = matrix(0,nn,1)
            ii = 0
            for(jj in (nfix-root):(-nfix+root))
            {
                ii = ii+1
                d[ii] = Bge(jj,nq,B,cc)
            }
            if(root)
                d[nn-1] = R[nn-nfix]
            ##  COMPUTE Bhat = inv(A)*d
            Bhat = Ainv%*%d
            for(ii in (nfix+1):(n-nfix))
                AA[ii,(ii-nfix):(ii+nfix)] = t(Bhat)
        }
    }

###  Baxter-King filter

    if (type == "baxter-king")
    {
        bb = matrix(0,2*nfix+1,1)
        bb[(nfix+1):(2*nfix+1)] = B[1:(nfix+1)]
        bb[nfix:1] = B[2:(nfix+1)]
        bb = bb - sum(bb)/(2*nfix+1)
        for(i in (nfix+1):(n-nfix))
            AA[i,(i-nfix):(i+nfix)] = t(bb)
    }

###  Trigonometric Regression filter

    if(type == "trigonometric")
    {
        jj = 1:(n/2)
	## find frequencies in desired band omitting n/2;
	jj = jj[((n/pu)<=jj & jj<=(n/pl) & jj<(n/2))]
        if(!any(jj))
	   stop("frequency band is empty in trigonometric regression")
        om = 2*pi*jj/n
	if(pl > 2)
        {
            for(t in 1:n)
            {
   		for(k in n:1)
                {
                    l = t-k
                    tmp = sum(cos(om*l))
                    AA[t,k] = tmp
                }
            }
        }
	else
        {
            for(t in 1:n)
            {
   		for(k in n:1)
                {
                    l = t-k
                    tmp = sum(cos(om*l))
                    tmp2 = (cos(pi*(t-l))*cos(pi*t))/2
                    AA[t,k] = tmp + tmp2
                }
            }
        }
	AA = AA*2/n
    }


###  check that sum of all filters equal 0 if assuming unit root

    if(root)
    {
        tst = max(abs(c(apply(AA,1,sum))))
        if((tst > 1.e-09) && root)
        {
            warning("Bhat does not sum to 0 ")
            cat("test =",tst,"\n")
        }
    }

###	compute filtered time series using selected filter matrix AA

    if(drift)
        x = undrift(x)

    x.cycle = AA%*%as.matrix(x)

    if(type=="fixed" || type=="symmetric" || type=="baxter-king")
    {
        if(nfix>0)
            x.cycle[c(1:nfix,(n-nfix+1):n)] = NA
    }
    x.trend = x-x.cycle
    if(is.ts(xo))
    {
        tsp.x = tsp(xo)
        x.cycle=ts(x.cycle,star=tsp.x[1],frequency=tsp.x[3])
        x.trend=ts(x.trend,star=tsp.x[1],frequency=tsp.x[3])
        x=ts(x,star=tsp.x[1],frequency=tsp.x[3])
    }

    if(type=="asymmetric")
        title = "Chiristiano-Fitzgerald Asymmetric Filter"
    if(type=="symmetric")
        title = "Chiristiano-Fitzgerald Symmetric Filter"
    if(type=="fixed")
        title = "Chiristiano-Fitzgerald Fixed Length Filter"
    if(type=="baxter-king")
        title = "Baxter-King Fixed Length Filter"
    if(type=="trigonometric")
        title = "Trigonometric Regression Filter"

    res <- list(cycle=x.cycle,trend=x.trend,fmatrix=AA,title=title,
                xname=xname,call=as.call(match.call()),
                type=type,pl=pl,pu=pu,nfix=nfix,root=root,drift=drift,
                theta=theta,method="cffilter",x=x)

    return(structure(res,class="mFilter"))
}


###======================================================================
###				Functions
###======================================================================

Bge <- function(jj,nq,B,cc)
{
###
###  closed form solution for integral of B(z)g(z)(1/z)^j  (eqn 16)
###     nq > 0, jj >= 0
###     if nq = 0, y = 2*pi*B(absj+1)*cc(1);
###
    absj =abs(jj)
    if(absj >= nq)
    {
        dj = B[absj+1]*cc[1] + t(B[(absj+2):(absj+nq+1)])%*%cc[2:(nq+1)]
        dj = dj + t(flipud(B[(absj-nq+1):absj]))%*%cc[2:(nq+1)]
    }
    else if(absj >= 1)
    {
        dj = B[absj+1]*cc[1] + t(B[(absj+2):(absj+nq+1)])%*%cc[2:(nq+1)]
        dj = dj + t(flipud(B[1:absj]))%*%cc[2:(absj+1)]
        dj = dj + t(B[2:(nq-absj+1)])%*%cc[(absj+2):(nq+1)]
    }
    else
        dj = B[absj+1]*cc[1] + 2*t(B[2:(nq+1)])%*%cc[2:(nq+1)]

    y = 2*pi*dj
    return(y)
}

###
###-----------------------------------------------------------------------
###

Abuild <- function(nn,nq,g,root)
{
###
###  builds the nn x nn A matrix in A.12
###   if root == 1 (unit root)
###      Abig is used to construct all but the last 2 rows of the A matrix
###   elseif root == 0 (no unit root)
###      Abig is used to construct the entire A matrix
###
    if(root)
    {
	Abig=matrix(0,nn,nn+2*(nq-1))
        for(j in 1:(nn-2))
	   Abig[j,j:(j+2*nq)] = t(g)
	A = Abig[,nq:(nn+nq-1)]
	##   construct A(-f)
	Q = -matrix(1,nn-1,nn)
	##Q = tril(Q);
        Q[upper.tri(Q)] <- 0
	F = matrix(0,1,nn-1)
	F[(nn-1-nq):(nn-1)] = g[1:(nq+1)]
	A[(nn-1),] = F%*%Q
	##    construct last row of A
        A[nn,] = matrix(1,1,nn)
    }
    else
    {
	Abig=matrix(0,nn,nn+2*(nq-0))
	for(j in 1:nn)
        {
            Abig[j,j:(j+2*nq)] = c(g)
        }
        A = Abig[,(nq+1):(nn+nq)]
    }
    ##    multiply A by 2*pi
    A = 2*pi*A
    return(A)
}
###
###-----------------------------------------------------------------------
###

undrift <- function(x)
{
###
###  This function removes the drift or a linear time trend from a time series using the formula
###			drift = (x(n) - x(1)) / (n-1).
###
###  Input:  x - data matrix x where columns represent different variables, x is (n x # variables).
###  Output: xun - data matrix same size as x with a different drift/trend removed from each variable.
###
    x = as.matrix(x)
    nv = dim(x)
    n = nv[1]
    nvars = nv[2]
    xun = matrix(0,n,nvars)
    dd = as.matrix(0:(n-1))
    for(ivar in 1:nvars)
    {
        drift = (x[n,ivar]-x[1,ivar]) / (n-1)
        xun[,ivar] = x[,ivar] - dd*drift
    }
    if(nvars==1) xun = c(xun)

    return(xun)
}

###
###-----------------------------------------------------------------------
###

### function that reverses the columns of a matrix (matlab equivalent)
flipud <- function(x) {apply(as.matrix(x),2,rev)}

### function that reverses the rows of a matrix (matlab equivalent)
fliplr <- function(x) {t(apply(as.matrix(x),1,rev))}
