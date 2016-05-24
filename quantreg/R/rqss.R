"rqss.fit" <-
function (x, y, tau = 0.5, method = "sfn", rhs = NULL, control, ...)
{
    if(is.null(rhs)) rhs <- (1 - tau) * t(x) %*% rep(1,length(y))
    else tau <- 0.5
    fit <- switch(method,
	sfn = rq.fit.sfn(x, y, tau = tau,  rhs = rhs, control = control, ...),
        sfnc = rq.fit.sfnc(x, y, tau = tau,  rhs = rhs, control = control, ...), {
            what <- paste("rq.fit.", method, sep = "")
            if (exists(what, mode = "function"))
                (get(what, mode = "function"))(x, y, ...)
            else stop(paste("unimplemented method:", method))
        })
    fit$contrasts <- attr(x, "contrasts")
    fit$resid <- c(y - x %*% fit$coef)

    fit
}
"untangle.specials" <-
function (tt, special, order = 1)
{
    spc <- attr(tt, "specials")[[special]]
    if (length(spc) == 0)
        return(list(vars = character(0), terms = numeric(0)))
    facs <- attr(tt, "factor")
    fname <- dimnames(facs)
    ff <- apply(facs[spc, , drop = FALSE], 2, sum)
    list(vars = (fname[[1]])[spc], terms = seq(ff)[ff & match(attr(tt,
        "order"), order, nomatch = 0)])
}

"qss" <-
function (x, constraint = "N", lambda = 1, ndum = 0, dummies = NULL, 
    w = rep(1, length(x))) 
{
    if (is.matrix(x)) {
        if (ncol(x) == 2) 
            qss <- qss2(x, constraint = constraint, dummies = dummies, 
                lambda = lambda, ndum = ndum, w = w)
        else if (ncol(x) == 1) 
            x <- as.vector(x)
        else stop("qss objects must have dimension 1 or 2")
    }
    if (is.vector(x)) 
        qss <- qss1(x, constraint = constraint, lambda = lambda, 
            dummies = dummies, ndum = ndum, w = w)
    qss
}

"qss2" <-
function(x, y, constraint = "N", lambda = 1, ndum= 0, dummies = NULL, w=rep(1,length(x))){
#
# Sparse Additive Quantile Smoothing Spline Models - Bivariate (Triogram) Module
#
# This function returns a structure intended to make model.matrix for a bivariate
# nonparametric component of a model formula specified by a call to rqss().  A sparse form
# of the Frisch Newton algorithm is eventually called to compute the estimator.
# An optional convexity/concavity constraint can be specified.  If
# the formula consists of a single qss component then the estimator solves the
# following variational problem:
#
#       min sum rho_tau (y_i - g(x_i)) + lambda V(grad(g))
#
# where V(f) denotes the total variation of the function f.  The solution is a piecewise
# linear function on the Delaunay triangulation formed by the observed (x_i,y_i) points.
# Additive models can consist
# of several components of this form plus partial linear and univariate qss components.
# To resolve the identifiability problem we delete the first column of the qss design
# components.  On return F contains the fidelity portion of the design, A the penalty
# contribution of the design, R the constraint portion, and r the rhs of the constraints.
#
# Constraints are specified by the constraint argument:
#
#	N	none
#	V	convex
#	C	concave
#
# Author:  Roger Koenker   April 2, 2003
#
# For a prototype see triogram in ~roger/projects/tv/cobar/.RData on ysidro.
#
# Warning:   Under development...todo:
#
#       o  weights
#       o  dummy x's
#       o  tau's
#       o  lambda's
#       o  ...
#
   stopifnot(requireNamespace("tripack"))
#
    y <- x[,2]
    x <- x[,1]
    n <- length(x)
    if (n != length(y))
        stop("xy lengths do not match")
    f <- triogram.fidelity(x, y, ndum = ndum, dummies = dummies)
    F <- f$F
    A <- triogram.penalty(f$x, f$y)
    switch(constraint, V = {
        R <- A
        r <- rep(0, nrow(R))
    }, C = {
        R <- -A
        r <- rep(0, nrow(R))
    }, N = {
        R = NULL
        r = NULL
    })
    list(x = list(x = f$x, y = f$y), F = F[, -1], dummies = f$dummies,
        lambda = lambda, A = A[, -1], R = R[, -1], r = r)
}


"qss1" <-
function (x, constraint = "N", lambda = 1, dummies = dummies,
        ndum = 0, w = rep(1, length(x))){
#
# Sparse Additive Quantile Smoothing Spline Models - Univariate Module
#
# This function returns a structure intended to make model.matrix for a univariate
# nonparametric component of a model formula specified by a call to rq().  A sparse form
# of the Frisch Newton algorithm is eventually called to compute the estimator.
# Optional monotonicity and/or convexity/concavity constraints can be specified.  If
# the formula consists of a single qss component then the estimator solves the
# following variational problem:
#
#       min sum rho_tau (y_i - g(x_i)) + lambda V(g')
#
# where V(f) denotes the total variation of the function f.  The solution is a piecewise
# linear function with "knots" at the observed x_i points.  Additive models can consist
# of several components of this form plus partial linear and triogram components.
# To resolve the identifiability problem we delete the first column of the qss design
# components.  On return F contains the fidelity portion of the design, A the penalty
# contribution of the design, R the constraint portion, and r the rhs of the constraints.
#
# Constraints are specified by the constraint argument:
#
#	N	none
#	I	monotone increasing
#	D	monotone decreasing
#	V	convex
#	C	concave
#	CI	concave and monotone increasing
#	...	etc
#
# Author:  Roger Koenker   February 27, 2003
#
# Warning:   Under development...todo:
#
#       o  weights
#       o  dummy x's
#       o  tau's
#       o  lambda's
#       o  ...
#
#
    xun <- unique(x[order(x)])
    h <- diff(xun)
    nh <- length(h)
    nx <- length(x)
    p <- nh + 1
    B <- new("matrix.csr", ra = c(rbind(-1/h, 1/h)), ja = as.integer(c(rbind(1:nh,
        2:(nh + 1)))), ia = as.integer(2 * (1:(nh + 1)) - 1),
        dimension = as.integer(c(nh, nh + 1)))
    makeD <- function(p) {
        new("matrix.csr", ra = c(rbind(rep(-1, (p - 1)), rep(1,
            (p - 1)))), ja = as.integer(c(rbind(1:(p - 1), 2:p))),
            ia = as.integer(2 * (1:p) - 1), dimension = as.integer(c(p -
                1, p)))
    }
    D <- makeD(nh)
    A <- D %*% B
    if (length(xun) == length(x)){
        F <- new("matrix.csr", ra = rep(1, nx), ja = as.integer(rank(x)),
            ia = 1:(nx + 1), dimension = as.integer(c(nx, nx)))
	}
    else {
        F <- new("matrix.csr", ra = rep(1, nx), ja = as.integer(factor(x)),
            ia = 1:(nx + 1), dimension = as.integer(c(nx, length(xun))))
	}

   switch(constraint,
        V = {   R <- A;
                r <- rep(0,nrow(R))
                },
        C = {   R <- -A;
                r <- rep(0,nrow(R))
                },
        I = {   R <- makeD(p)
                r <- rep(0,p-1)
                },
        D = {   R <- -makeD(p)
                r <- rep(0,p-1)
                },
        VI = {  R <- makeD(p)
		R <- rbind(R,A)
                r <- rep(0,nrow(R))
		},
        VD = {  R <- -makeD(p)
		R <- rbind(R,A)
                r <- rep(0,nrow(R))
		},
        CI = {  R <- makeD(p)
		R <- rbind(R,-A)
                r <- rep(0,nrow(R))
		},
        CD = {  R <- -makeD(p)
		R <- rbind(R,-A)
                r <- rep(0,nrow(R))
		},
	N = { R=NULL; r=NULL}
	)
   list(x = list(x=xun), F=F[,-1], lambda = lambda, A=A[,-1], R=R[,-1], r=r)
}
"plot.qss1" <-
function(x, rug = TRUE, jit = TRUE, add = FALSE, ...)
{
if(!add) plot(x[,1],x[,2],type = "n", ...)
lines(x[,1],x[,2], ...)
if(rug) {
    if(jit)
       rug(jitter(x[,1]))
    else rug(x[,1])
    }
}
"plot.qss2" <-
function (x, render = "contour", ncol = 100, zcol = NULL, ...)
{
    stopifnot(requireNamespace("tripack"))
    y <- x[, 2]
    z <- x[, 3]
    x <- x[, 1]
    tri <- tripack::tri.mesh(x, y)
    if (render == "rgl") {
        if(!requireNamespace("rgl",quietly=TRUE))
		stop("The package rgl is required")
        collut <- terrain.colors(ncol)
        if (!length(zcol))
            zcol <- z
        if (max(z) > max(zcol) || min(z) < min(zcol))
            warning("fitted z values out of range of zcol vector")
        zlim <- range(zcol)
        colz <- ncol * (z - zlim[1])/(zlim[2] - zlim[1]) + 1
        colz <- collut[colz]
        s <- c(t(tripack::triangles(tri)[, 1:3]))
        rgl::rgl.triangles(x[s], y[s], z[s], col = colz[s])
    }
    else {
        stopifnot(requireNamespace("akima"))
        if(render == "contour"){
                plot(x, y, type = "n", ...)
                contour(akima::interp(x, y, z), add = TRUE, frame.plot = TRUE, ...)
                tripack::convex.hull(tri, plot.it = TRUE, add = TRUE)
                }
        else if(render == "persp")
                persp(akima::interp(x, y, z, ), theta = -40, phi = 20, xlab = "x",
                        ylab = "y", zlab = "z", ...)
        else stop(paste("Unable to render: ",render))
    }
}
plot.rqss <-
function (x, rug = TRUE, jit = TRUE, bands = NULL, coverage = 0.95, 
    add = FALSE, shade = TRUE, select = NULL, pages = 0, titles = NULL, ...) 
{
    SetLayout <- function(m, p) {
            # Shamelessly cribbed from mgcv
            # m is the number of plots
            # p is the number of pages
        if (p > m) 
            p <- m
        if (p < 0) 
            p <- 0
        if (p != 0) {
            q <- m%/%p
            if ((m%%p) != 0) {
                q <- q + 1
                while (q * (p - 1) >= m) p <- p - 1
            }
            c <- trunc(sqrt(q))
            if (c < 1) 
                c <- 1
            r <- q%/%c
            if (r < 1) 
                r <- 1
            while (r * c < q) r <- r + 1
            while (r * c - q > c && r > 1) r <- r - 1
            while (r * c - q > r && c > 1) c <- c - 1
            oldpar <- par(mfrow = c(r, c))
        }
        else oldpar <- par()
        return(oldpar)
    }
    m <- length(x$qss)
    if (m == 0) 
        stop("No qss object to plot")
    if(length(select)) {
        if(all(select %in% 1:m))
            oldpar <- SetLayout(length(select), pages)
        else
            stop(paste("select must be in 1:",m,sep="")) 
        }
    else
        oldpar <- SetLayout(m, pages)
    if ((pages == 0 && prod(par("mfrow")) < m && dev.interactive()) || 
        pages > 1 && dev.interactive()) 
        ask <- TRUE
    else ask <- FALSE
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    qssnames <- names(x$qss)
    if(length(titles)){
        if(length(titles) != length(qssnames))
            stop("Length of titles doesn't match length of qssnames")
	}
    else{
        titles <- paste("Effect of ", qssnames)
	}
    if (length(bands)) {
        rdf <- x$n - x$edf
        if (any(unlist(lapply(x$qss, function(x) ncol(x$xyz) == 3)))) 
            warning("Can't plot confidence bands in 3D (yet)")
        band <- as.list(rep(NA, m))
        V <- summary(x, cov = TRUE, ...)$Vqss
        "summary.qss1" <- function(object, V, ngrid = 400, ...) {
            x <- object$xyz[, 1]
            eps <- 0.01
            newd <- data.frame(x = seq(min(x) + eps, max(x) - 
                eps, length = ngrid))
            G <- predict(object, newd, ...)
            ones <- as.matrix.csr(matrix(1, nrow(G$D), 1))
            D <- cbind(ones, G$D)
            S <- as.matrix(D %*% V %*% t(D))
            se <- sqrt(diag(S))
            cv <- qt(1 - (1 - coverage)/2, rdf)
            if (bands %in% c("uniform","both")) {
                E <- eigen(as.matrix(V))
                B <- E$vectors %*% diag(sqrt(pmax(0,E$values))) %*% t(E$vectors)
                D <- as.matrix(D)
                BX1 <- B %*% t(D[-1, ])
                BX1 <- BX1/sqrt(apply(BX1^2, 2, sum))
                BX0 <- B %*% t(D[-nrow(D), ])
                BX0 <- BX0/sqrt(apply(BX0^2, 2, sum))
                kappa <- sum(sqrt(apply((BX1 - BX0)^2, 2, sum)))
                cvu <- critval(kappa, alpha = 1 - coverage, rdf = rdf)
                if(bands == "both") cv <- c(cvu,cv)
                else cv <- cvu
                }
            list(pred = data.frame(x = G$x, y = G$y, se = se), 
                cv = cv)
        }
    }
    if(!length(select)) 
        select <- 1:m
    for (i in select) {
        qss <- x$qss[[i]]$xyz
        if (ncol(qss) == 3) {
            qss[, 3] <- x$coef[1] + qss[, 3]
            plot.qss2(qss, ...)
        }
        else if (ncol(qss) == 2) {
            if (length(bands)) {
                if (is.na(x$coef["(Intercept)"])) 
                  stop("rqss confidence bands require an intercept parameter")
                B <- summary(x$qss[[i]], V[[i]], ...)
                cv <- B$cv
                B <- B$pred
                B$y <- B$y + x$coef["(Intercept)"]
                bandcol <- c("grey85","grey65")
		for(k in 1:length(cv)){
                 if (add || k > 1) 
                  if(shade){
                     polygon(c(B$x,rev(B$x)),
                        c(B$y - cv[k] * B$se,rev(B$y + cv[k] * B$se)),
			col = bandcol[k], border = FALSE)
			}
                  else
                    matlines(B$x, cbind(B$y, B$y + cv[k] * cbind(-B$se, 
                        B$se)), lty = c(1, 2, 2), col = c("black", "blue", "blue"))
                else {
                  matplot(B$x, B$y + cv[k] * cbind(-B$se, B$se), 
                    xlab = paste(qssnames[i]), ylab = "Effect", 
                    type = "n", ...)
                  if(shade){
                     polygon(c(B$x,rev(B$x)),
                        c(B$y - cv[k] * B$se,rev(B$y + cv[k] * B$se)),
			col = bandcol[k], border = FALSE)
			}
                  else{
                     lines(B$x, B$y + cv * B$se, lty = 2, ...)
                     lines(B$x, B$y - cv * B$se, lty = 2, ...)
		     }
	           }
                  lines(B$x, B$y, ...)
                }
                band[[1]] <- list(x = B$x, blo = B$y - B$se %o% cv, 
                  bhi = B$y + B$se %o% cv)
                if (rug) {
                  if (jit) 
                    rug(jitter(qss[, 1]))
                  else rug(qss[, 1])
                }
            }
            else {
                qss[, 2] <- x$coef[1] + qss[, 2]
                plot.qss1(qss, xlab = paste(qssnames[i]), ylab = "Effect", 
                  rug = rug, jit = jit, add = add, ...)
            }
            title(titles[i])
        }
        else stop("invalid fitted qss object")
    }
    if (pages > 0) 
        par(oldpar)
    if (length(bands)) 
        class(band) <- "rqssband"
    else
        band <- NULL
    invisible(band)
}

"triogram.fidelity" <- function (x, y, ndum=0, dummies = NULL)
{
#Make fidelity block of the triogram design in sparse matrix.csr form
#The rather esoteric match call identifies and handles duplicated xy points
n <- length(x)
A <- as.data.frame(cbind(x,y))
dupA <- duplicated(A)
if(any(dupA)){
  x <- x[!dupA]
  y <- y[!dupA]
  J <- match(do.call("paste",c(A,"\r")),do.call("paste",c(A[!dupA,],"\r")))
  z <- new("matrix.csr",ra=rep(1,n), ja=J, ia=1:(n+1),dimension=as.integer(c(n,max(J))))
  }
else{
  z <- as(n,"matrix.diag.csr")
  z <- as(z,"matrix.csr")
 }
#Augment with dummy vertices, if any...
if(length(dummies)){
	if (is.list(dummies)){
        	if (all(!is.na(match(c("x", "y"), names(dummies))))){
        		ndum <- length(dummies$x)
			if(length(dummies$y) == ndum){
        			x <- c(x,dummies$x)
        			y <- c(y,dummies$y)
        			zdum <- as.matrix.csr(0,n,ndum)
        			z <- cbind(z,zdum)
				}
        		else stop("dummies x and y components differ in length")
			}
        	else stop("dummies list lacking x and y elements")
		}
    	else stop("dummies argument invalid (not a list) in triogram.fidelity")
	}
else if(ndum > 0){
        u <- runif(ndum); v <- runif(ndum)
        xd <- min(x) + u * (max(x)-min(x))
        yd <- min(y) + v * (max(y)-min(y))
        T <- tripack::tri.mesh(x,y)
        s <- tripack::in.convex.hull(T,xd,yd)
        x <- c(x,xd[s])
        y <- c(y,yd[s])
        ndum <- sum(s)
        zdum <- as.matrix.csr(0,n,ndum)
        z <- cbind(z,zdum)
	dummies <- list(x = xd[s],y = yd[s])
        }
list(x=x,y=y,F=z, dummies = dummies)
}
"triogram.penalty" <- function (x, y, eps = .Machine$double.eps)
{
    n <- length(x)
    tri <- tripack::tri.mesh(x, y)
    bnd <- tripack::on.convex.hull(tri,x,y)
    q <- length(tri$tlist)
    m <- 13 * n
    z <- .Fortran("penalty", as.integer(n), as.integer(m), as.integer(q),
        as.double(x), as.double(y), as.integer(bnd),as.integer(tri$tlist),
        as.integer(tri$tlptr), as.integer(tri$tlend), rax = double(m),
	jax = integer(m), ned = integer(1), as.double(eps), ierr = integer(1),
	PACKAGE = "quantreg")[c("rax", "jax", "iax", "ned", "ierr")]
    if (z$ierr == 1)
        stop("collinearity in ggap")
    nnz <- 4 * z$ned
    ra <- z$rax[1:nnz]
    ja <- z$jax[1:nnz]
    ia <- as.integer(1 + 4 * (0:z$ned))
    dim <- as.integer(c(z$ned, n))
    new("matrix.csr",ra=ra,ja=ja,ia=ia,dimension=dim)
}

predict.rqss <-
function (object, newdata, interval = "none",  level = 0.95, ...) 
{
    ff <- object$fake.formula
    Terms <- delete.response(terms(object$formula, "qss"))
    Names <- all.vars(parse(text = ff))
    if (any(!(Names %in% names(newdata)))) 
        stop("newdata doesn't include some model variables")
    ff <- reformulate(ff)
    nd <- eval(model.frame(ff, data = newdata), parent.frame())
    qssterms <- attr(Terms, "specials")$qss
    if (length(qssterms)) {
        tmp <- untangle.specials(Terms, "qss")
        dropv <- tmp$terms
        m <- length(dropv)
        if (length(dropv)) 
            PLTerms <- Terms[-dropv]
        attr(PLTerms, "specials") <- tmp$vars
    }
    else {
        PLTerms <- Terms
        m <- 0
    }
    if(requireNamespace("MatrixModels") && requireNamespace("Matrix"))
        X <- as(MatrixModels::model.Matrix(PLTerms, data = nd, sparse = TRUE),"matrix.csr")
    else
        X <- model.matrix(PLTerms, data = nd)
    p <- ncol(X)
    y <- X %*% object$coef[1:p]
    X <- as.matrix.csr(X)
    if (m > 0) {
        for (i in 1:m) {
            qss <- object$qss[[i]]
            names <- all.vars(Terms[dropv[i]])
            names <- names[names %in% Names]
            dimnames(qss$xyz)[[2]] <- c(names, "zfit")
            newd <- nd[names]
            if (ncol(qss$xyz) == 3) {
                g <- predict.qss2(qss$xyz, newdata = newd, ...)
                y <- y + g$z
                if(interval == "confidence")
                    X <- cbind(X,g$D)
            }
            else if (ncol(qss$xyz) == 2) {
                g <- predict(qss, newdata = newd, ...)
                y <- y + g$y
                if(interval == "confidence")
                    X <- cbind(X,g$D)
            }
            else stop("invalid fitted qss object")
        }
    }
    if(interval == "confidence"){
        v <- sqrt(diag(X %*% summary(object, cov = TRUE)$V %*% t(X)))
        calpha <- qnorm(1 - (1-level)/2)
        y <- cbind(y,y - v*calpha,y + v*calpha)
        dimnames(y)[[2]] <- c("yhat","ylower","yupper")
	}
    y
}

"predict.qss1" <-
function (object, newdata, ...)
{
    x <- object$xyz[, 1]
    y <- object$xyz[, 2] 
    if(ncol(newdata)==1)
        newdata <- newdata[,1]
    else
        stop("newdata should have only one column for predict.qss1")
    if (any(diff(x) < 0))
        stop("x coordinates in qss1 object not monotone")
    if (max(newdata) > max(x) || min(newdata) < min(x))
        stop("no extrapolation allowed in predict.qss")
    bin <- cut(newdata, unique(x), label = FALSE, include.lowest = TRUE)
    p <- length(x)
    m <- length(newdata)
    V <- cbind(bin, bin + 1)
    B <- cbind(x[bin + 1] - newdata, newdata - x[bin])/(x[bin +
        1] - x[bin])
    ra <- c(t(B))
    ja <- as.integer(c(t(V)))
    ia <- as.integer(c(2 * (1:m) - 1, 2 * m + 1))
    dim <- c(m, p)
    D <- new("matrix.csr", ra = ra, ja = ja, ia = ia, dimension = dim)
    list(x = newdata, y = D %*% y, D = D[, -1])
}

predict.qss2 <- function (object, newdata, ...) 
{
    x <- object[, 1]
    y <- object[, 2]
    z <- object[, 3]
    tri.area <- function(v) {
        0.5 * ((v[2, 1] - v[1, 1]) * (v[3, 2] - v[1, 2]) - (v[3, 
            1] - v[1, 1]) * (v[2, 2] - v[1, 2]))
    }
    barycentric <- function(v) {
        b <- rep(0, 3)
        Area <- tri.area(v[1:3, ])
        b[1] <- tri.area(v[c(4, 2, 3), ])/Area
        b[2] <- tri.area(v[c(1, 4, 3), ])/Area
        b[3] <- tri.area(v[c(1, 2, 4), ])/Area
        if (any(b < 0 || b > 1)) 
            stop("barycentric snafu")
        b
    }
    if (is.list(newdata)) {
	fnames <- (dimnames(object)[[2]])[1:2]
        if (all(!is.na(match(fnames, names(newdata))))) {
            newx <- newdata[[fnames[1]]]
            newy <- newdata[[fnames[2]]]
	   }
        else (stop("qss object and newdata frame names conflict"))
        }
    else if (is.matrix(newdata)) 
        if (ncol(newdata) == 2) {
            newx <- newdata[, 1]
            newy <- newdata[, 2]
        }
        else (stop("newdata matrix must have 2 columns"))
    tri <- tripack::tri.mesh(x, y)
    if (!all(tripack::in.convex.hull(tri, newx, newy))) 
        stop("some newdata points outside convex hull")
    p <- length(x)
    m <- length(newx)
    V <- matrix(0, m, 3)
    B <- matrix(0, m, 3)
    for (i in 1:m) {
        V[i, ] <- unlist(tripack::tri.find(tri, newx[i], newy[i]))
        v <- rbind(cbind(x[V[i, ]], y[V[i, ]]), c(newx[i], newy[i]))
        B[i, ] <- barycentric(v)
    }
    ra <- c(t(B))
    ja <- as.integer(c(t(V)))
    ia <- as.integer(3 * (0:m) + 1)
    D <- new("matrix.csr", ra = ra, ja = ja, ia = ia, dimension = c(m, 
        p))
    list(x = newx, y = newy, z = c(D %*% z), D = D[,-1])
}

fitted.rqss  <-  function(object, ...) (object$X %*% object$coef)[1:object$n]
resid.rqss <- function(object, ...) object$resid[1:object$n]

"rqss" <- function (formula, tau = 0.5, data = parent.frame(), weights, 
    na.action, method = "sfn", lambda = NULL, contrasts = NULL, 
    ztol = 1e-05,  control = sfn.control(), ...) 
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "weights", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    m[[1]] <- as.name("model.frame")
    special <- "qss"
    Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
    qssterms <- attr(Terms, "specials")$qss
    dropx <- NULL
    if(length(tau) > 1){
        tau <- tau[1]
        warning("multiple taus not supported, using first element")
        }
    if (length(qssterms)) {
        tmpc <- untangle.specials(Terms, "qss")
        ord <- attr(Terms, "order")[tmpc$terms]
        if (any(ord > 1)) 
            stop("qss can not be used in an interaction")
        dropx <- tmpc$terms
        if (length(dropx)) 
            Terms <- Terms[-dropx]
        attr(Terms, "specials") <- tmpc$vars
        fnames <- function(x){
           fy <- all.names(x[[2]])
           if(fy[1] == "cbind") fy <- fy[-1]
           fy
           }
        fqssnames <- unlist(lapply(parse(text = tmpc$vars), fnames))
        qssnames <- unlist(lapply(parse(text = tmpc$vars), function(x) deparse(x[[2]])))
    }
    m$formula <- Terms
    ff <- delete.response(terms(formula(Terms)))
    if(exists("fqssnames")){
       ffqss <- paste(fqssnames, collapse = "+")
       ff <- paste(deparse(formula(ff)), "+", ffqss)
    }
    m <- eval(m, parent.frame())
    weights <- model.extract(m, weights)
    process <- (tau < 0 || tau > 1)
    Y <- model.extract(m, "response")
    if(requireNamespace("MatrixModels") && requireNamespace("Matrix")){
         X <- MatrixModels::model.Matrix(Terms, m, contrasts, sparse = TRUE)
         vnames <- dimnames(X)[[2]]
         X <- as(X ,"matrix.csr")
         }
    else{
         X <- model.matrix(Terms, m, contrasts)
         vnames <- dimnames(X)[[2]]
         }
    p <- ncol(X)
    pf <- environment(formula)
    nrL <- 0
    if (method == "lasso") {
        if (!length(lambda)) 
            stop("No lambda specified for lasso constraint")
        if (length(lambda) == 1) 
            lambda <- c(0, rep(lambda, p - 1))
        if (length(lambda) != p) 
            stop("lambda must be either of length p, or length 1")
        if (any(lambda < 0)) 
            stop("negative lambdas disallowed")
        L <- diag(lambda, nrow = length(lambda))
        L <- L[which(lambda != 0), , drop = FALSE]
        L <- as.matrix.csr(L)
        nrL <- nrow(L)
        ncL <- ncol(L)
    }
    if (length(qssterms) > 0) {
        F <- as.matrix.csr(X)
        qss <- lapply(tmpc$vars, function(u) eval(parse(text = u), 
            data, enclos = pf))
        mqss <- length(qss)
        ncA <- rep(0, mqss + 1)
        nrA <- rep(0, mqss + 1)
        nrR <- rep(0, mqss + 1)
        for (i in 1:mqss) {
            F <- cbind(F, qss[[i]]$F)
            ncA[i + 1] <- ncol(qss[[i]]$A)
            nrA[i + 1] <- nrow(qss[[i]]$A)
            nrR[i + 1] <- ifelse(is.null(nrow(qss[[i]]$R)), 0, 
                nrow(qss[[i]]$R))
            vnames <- c(vnames, paste(qssnames[i], 1:ncA[i + 
                1], sep = ""))
        }
        A <- as.matrix.csr(0, sum(nrA), sum(ncA))
        if (sum(nrR) > 0) {
            R <- as.matrix.csr(0, sum(nrR), sum(ncA))
            nrR <- cumsum(nrR)
        }
        ncA <- cumsum(ncA)
        nrA <- cumsum(nrA)
        lambdas <- rep(0, mqss)
        for (i in 1:mqss) {
            lambdas[i] <- qss[[i]]$lambda
            Arows <- (1 + nrA[i]):nrA[i + 1]
            Acols <- (1 + ncA[i]):ncA[i + 1]
            A[Arows, Acols] <- qss[[i]]$lambda * qss[[i]]$A
            if (nrR[i] < nrR[i + 1]) 
                R[(1 + nrR[i]):nrR[i + 1], (1 + ncA[i]):ncA[i + 
                  1]] <- qss[[i]]$R
        }
        A <- cbind(as.matrix.csr(0, nrA[mqss + 1], p), A)
        if (nrR[mqss + 1] > 0) {
            R <- cbind(as.matrix.csr(0, nrR[mqss + 1], p), R)
            r <- rep(0, nrR[mqss + 1])
        }
        else {
            R <- NULL
            r <- NULL
        }
        if (method == "lasso") {
            A <- rbind(cbind(L, as.matrix.csr(0, nrL, ncol(F) - ncL)), A)
        }
        X <- rbind(F, A)
        Y <- c(Y, rep(0, nrow(A)))
        rhs <- t(rbind((1 - tau) * F, 0.5 * A)) %*% rep(1, nrow(X))
        XpX <- t(X) %*% X
        nnzdmax <- XpX@ia[length(XpX@ia)] - 1
	if(is.null(control[["nsubmax"]]))
           control[["nsubmax"]] <- max(nnzdmax, floor(1000 + exp(-1.6) * nnzdmax^1.2))
        if(is.null(control[["nnzlmax"]]))
           control[["nnzlmax"]] <- floor(2e+05 - 2.8 * nnzdmax + 7e-04 * nnzdmax^2)
        if(is.null(control[["tmpmax"]]))
           control[["tmpmax"]] <- floor(1e+05 + exp(-12.1) * nnzdmax^2.35)
        fit <- if (length(r) > 0) 
            rqss.fit(X, Y, tau = tau, rhs = rhs, method = "sfnc", 
                R = R, r = r, control = control, ...)
        else rqss.fit(X, Y, tau = tau, rhs = rhs, method = "sfn", control = control, ...)
        for (i in 1:mqss) {
            ML <- p + 1 + ncA[i]
            MU <- p + ncA[i + 1]
            qss[[i]] <- list(xyz = cbind(qss[[i]]$x$x, qss[[i]]$x$y, 
                c(0, fit$coef[ML:MU])), dummies = qss[[i]]$dummies)
            if (ncol(qss[[i]]$xyz) == 2) 
                class(qss[[i]]) <- "qss1"
            else class(qss[[i]]) <- "qss2"
        }
        names(qss) <- qssnames
        fit$qss <- qss
    }
    else {
        X <- as.matrix.csr(X)
        nrA <- 0
        if (method == "lasso") {
            rhs <- t(rbind((1 - tau) * X, 0.5 * L)) %*% rep(1, nrow(X) + nrow(L))
            X <- rbind(X, L)
            Y <- c(Y, rep(0, nrL))
            nrA <- c(nrA, nrL)
        }
        else
            rhs <- NULL
        if (length(weights)) {
            if (any(weights < 0)) 
                stop("negative weights not allowed")
            X <- X * weights
            Y <- Y * weights
        }
        fit <- rqss.fit(X, Y, tau = tau,  rhs = rhs, control = control, ...)
        fit$nrA <- nrA
    }
    names(fit$coef) <- vnames
    n <- length(fit$resid) - nrL - nrA[length(nrA)]
    uhat <- fit$resid[1:n]
    Rho <- function(u, tau) sum(u * (tau - (u < 0)))
    fit$fidelity <- Rho(uhat, tau)
    fit$edf <- sum(abs(uhat) < ztol)
    fit$X <- X
    fit$n <- n
    fit$nrL <- nrL
    fit$terms <- Terms
    fit$fake.formula <- ff
    fit$formula <- formula
    fit$method <- method
    fit$call <- call
    fit$tau <- tau
    if (length(qssterms)) {
        fit$lambdas <- lambdas
        fit$qssnames <- qssnames
        fit$nrA <- nrA
        fit$ncA <- cumsum(c(p, diff(ncA)))
    }
    else fit$ncA <- p
    attr(fit, "na.message") <- attr(m, "na.message")
    class(fit) <- "rqss"
    fit
}


"summary.rqss" <- function(object, cov = FALSE, ztol = 1e-5, ...){
    resid <- object$resid
    coef <- object$coef
    lambdas <- object$lambdas
    formula <- object$formula
    fidelity <- object$fidelity
    edf <- object$edf
    nrA <- object$nrA
    ncA <- object$ncA
    nrL <- object$nrL
    tau <- object$tau
    X <- object$X
    p <- ncol(object$X)
    m <- length(ncA)
    n <- length(resid) - nrL - nrA[m]
    uhat <- resid[1:n]
    h <- bandwidth.rq(tau, n, hs = TRUE)
    if(tau + h > 1)
       stop("tau + h > 1:  error in summary.rq")
    if(tau - h < 0)
       stop("tau - h < 0:  error in summary.rq")
    h <- (qnorm(tau + h) - qnorm(tau - h)) * max(ztol,min(sqrt(var(uhat)),
        (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34))
    f <- c(dnorm(uhat/h)/h,rep(1, length(resid) - n))
    D <- t(X) %*% (f * X)
    D <- chol(.5 * (D + t(D)), ...)
    D <- backsolve(D,diag(p))
    D0 <- tau * (1 - tau) * t(X) %*% X
    V <- D %*% D0 %*% D
    scale <- mean(f)
    serr <- sqrt(diag(V))
    ptab <- array(coef[1:ncA[1]], c(ncA[1],4))
    pnames <- names(coef[1:ncA[1]])
    dimnames(ptab) <- list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ptab[, 2] <- serr[1:ncA[1]]
    ptab[, 3] <- ptab[, 1]/ptab[, 2]
    ptab[, 4] <- 2 * (1 - pt(abs(ptab[, 3]), n - edf))
    if(cov) {
        Vcov <- V[1:ncA[1],1:ncA[1]]
        Vqss <- as.list(1:(m-1))
        }
    if(m > 1) {
        penalty <- rep(NA, m - 1)
        qssedfs <- rep(NA, m - 1)
        ntab <- matrix(0, m - 1, 5)
        vnames <- sub(".$","",names(coef[ncA[-length(ncA)]+1]))
        dimnames(ntab) <- list(vnames, c("EDF", "Lambda", "Penalty", "F value", "Pr(>F)"))
        for(i in 2:m) {
           ntab[i - 1, 3] <- sum(abs(resid[n + nrL + ((nrA[i - 1] + 1):nrA[i])]))
           ntab[i - 1, 1] <- sum(abs(resid[n + nrL + ((nrA[i - 1] + 1):nrA[i])]) > ztol)
           v <- V[c(1, (ncA[i-1] + 1):ncA[i]), c(1, (ncA[i-1] + 1):ncA[i])]
           v <- .5 * (v + t(v))
           if(cov) Vqss[[i-1]] <- v
           b <- coef[(ncA[i-1] + 1):ncA[i]]
           ntab[i-1, 4]  <- t(b) %*% solve(v[-1,-1],b)/ntab[i-1,1]
           ntab[i-1, 5]  <- 1 - pf(ntab[i-1,4],ntab[i-1,1], n - edf)
           }
        ntab[,2] <- lambdas
        ntab[,3] <- ntab[,3]/lambdas
        if(cov)
           z <- list(coef = ptab, qsstab = ntab, fidelity = fidelity, tau = tau,
		formula = formula, edf = edf, n = n, Vcov = Vcov, Vqss = Vqss, V = V)
        else
           z <- list(coef = ptab, qsstab = ntab, fidelity = fidelity, tau = tau,
		formula = formula, edf = edf, n = n)
        }
    else
        if(cov)
           z <- list(coef = ptab, fidelity = fidelity, formula = formula, tau = tau,
		edf = edf, n = n, Vcov = Vcov, V = V)
        else
           z <- list(coef = ptab, fidelity = fidelity, formula = formula, tau = tau,
		edf = edf, n = n)

    class(z) <- "summary.rqss"
    return(z)
}
"plot.summary.rqss" <- function(x, ...)
    warning("No plot method for summary.rqss objects:   plot the rqss object instead")

"print.rqss" <- function(x, ...) {
    cat("Formula:\n")
    print(x$formula)
    sx <- summary(x)
    cat("Quantile fidelity at tau = ", x$tau, "is", sx$fidelity, "\n")
    cat("Estimated Model Dimension is", sx$edf, "\n")
}
 
"logLik.rqss" <- function(object, ...){
	n <- object$n
	tau <- object$tau
	val <- n * (log(tau * (1-tau)) - 1 - log(object$fidelity/n))
	attr(val,"n") <- n
	attr(val,"df") <- object$edf
	class(val) <- "logLik"
	val
	}
"AIC.rqss" <- function(object, ... , k = 2){
	v <- logLik(object)
	if(k < 0) 
		k <- log(attr(v,"n"))
	val <- AIC(logLik(object), k = k)
	attr(val,"edf") <- attr(v,"df")
	val
	}
"dither"  <- function(x, type = "symmetric", value = NULL) {
	if(length(x) == 0)
		return(x)
	if(!is.numeric(x))
		stop("'x' must be numeric")
	if(!length(value))
		value <- min(diff(sort(unique(x))))
	if(type == "symmetric")
		v <- x + runif(length(x), -value/2, value/2)
	else if(type == "right")
		v <- x + runif(length(x), 0, value)
	else
		stop("invalid type")
	v
}
	
print.summary.rqss <-
function (x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), 
    ...) 
{
    cat("Formula:\n")
    print(x$formula)
    if (length(x$coef) > 0) {
        cat("\nParametric coefficients:\n")
        printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    cat("\n")
    if (length(x$qsstab) > 0) {
        cat("Approximate significance of qss terms:\n")
        printCoefmat(x$qsstab, digits = digits, signif.stars = signif.stars, 
            has.Pvalue = TRUE, na.print = "NA", cs.ind = 1, ...)
    }
    cat("\n")
    if (length(x$fidelity) > 0) 
        cat("  Quantile Fidelity at tau = ", x$tau, "  is  ", formatC(x$fidelity, 
		digits = 6, width = 11), "\n", sep = "")
    cat("  Effective Degrees of Freedom = ", formatC(x$edf, digits = 5, width = 8, 
        flag = "-"), "  Sample Size = ", x$n, "\n", sep = "")
    invisible(x)
}
critval <- function(kappa, alpha = 0.05, rdf = 0){
# Hotelling tube critical value for uniform confidence bands
# This should agree to about 6 digits with the following call to locfit:
# crit(const = c(kappa,1), cov = 1-alpha,rdf = rdf)$crit.val
   tube <- function(x,alpha,rdf){
      if(rdf <= 0)
         kappa * exp(-x^2/2)/pi + 2 * (1 - pnorm(x)) -  alpha
      else
         kappa * (1+x^2/rdf)^(-rdf/2)/pi + 2 * (1 - pt(x,rdf)) -  alpha
      }
   uniroot(tube,c(1,5),alpha = alpha, rdf = rdf)$root
   }

