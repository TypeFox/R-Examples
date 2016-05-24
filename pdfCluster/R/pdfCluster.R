#package
#pdfCluster +pdfClusterPairs

#dyn.load("C:/Documents and Settings/menardi/Documenti/assegnoAA/pacchettoNew/pdfClusterNew/src/macon.dll")
#dyn.load("C:/Documents and Settings/menardi/Documenti/assegnoAA/pacchettoNew/pdfClusterNew/src/kepdf.dll")
#dyn.load("C:/Documents and Settings/menardi/Documenti/assegnoAA/pacchettoNew/pdfClusterNew/src/nc.dll")
#####
# R code for kernel density estimation
# R interfaced with C
		  
kepdf <- function (x, eval.points = x, kernel = "gaussian", bwtype = "fixed", h = h.norm(x), hx = NULL, alpha = 1/2) 
{
    bwtype <- match.arg(bwtype, choices=c("fixed", "adaptive"))
	if (!bwtype%in% c("fixed", "adaptive")) stop("bwtype should be one of 'fixed' or 'adaptive'")
    kernel <- match.arg(kernel, choices=c("gaussian", "t7"))
    if (!kernel%in% c("gaussian", "t7")) stop("kernel should be one of 'gaussian', 't7'")
    x <- data.matrix(x)
    eval.points <- data.matrix(eval.points)
    nx <- as.integer(nrow(x))
    ndim <- as.integer(ncol(x))
    #x <- matrix(as.double(x), nx, ndim)
    neval <- as.integer(nrow(eval.points))
    eval.points <- matrix(as.double(eval.points), neval, ndim)
    f <- as.double(rep(0, neval))
    if (bwtype == "fixed"){ 
        hm <- matrix(as.double(h), nx, ndim, byrow = TRUE)
        par <- list(h = h, hx = hx )
        }
    else if (bwtype == "adaptive") {
           hm <- hx
           if (is.null(hx)|| any(dim(hx)!=c(nx, ndim)))   
           hx <- hm <- matrix(as.double(hprop2f(x = x, h = h, alpha, kernel = kernel)), nx, ndim)
           par <- list(h = h, hx = hx, alpha= alpha )
           }
    rec.hm <- 1/(hm)
    prodhm <-apply(rec.hm, 1, prod) 
    pdf <- switch(kernel,  
        gaussian = .C("c_kepdfN", as.double(x), as.double(rec.hm/sqrt(2)), as.double(eval.points), 
              as.integer(nx), as.integer(ndim), as.integer(neval), 
              double(neval), as.double(prodhm), PACKAGE="pdfCluster"),
        t7 = .C("c_kepdft7", as.double(x), as.double(rec.hm), as.double(eval.points), 
              as.integer(nx), as.integer(ndim), as.integer(neval), 
              double(neval), as.double(prodhm), PACKAGE="pdfCluster")
	)
	new("kepdf", call = match.call(), x = x, eval.points = eval.points, estimate = pdf[[7]], kernel= kernel, bwtype = bwtype, par = par)
	}	

		
 
setClass("kepdf", representation(call="call", x="matrix", eval.points="matrix", 
estimate="numeric", kernel="character", bwtype="character", par="list"))

#function to determine bandwiths

# normal reference rule
h.norm <- function(x){
   x <- as.matrix(x)
   nd <- ncol(x)
   sd <- sqrt(diag(var(x)))
   sd*(4/((nd+2)*nrow(x)))^(1/(nd+4))
 }
 
#smoothing matrix of dimension n x d
		
 hprop2f <- function(x, h=h.norm(x), alpha=1/2, kernel="gaussian" ){
    x<-data.matrix(x)
    nx <- as.integer(nrow(x))
    ndim <- as.integer(ncol(x))
    hm <- matrix(as.double(h), nx, ndim, byrow = TRUE)
    rec.hm <- 1/(hm)
    prodhm <-apply(rec.hm, 1, prod) 
    pdfpil <- switch(kernel,  
        gaussian = .C("c_kepdfN", as.double(x), as.double(rec.hm/sqrt(2)), as.double(x), 
              as.integer(nx), as.integer(ndim), as.integer(nx), 
              double(nx), as.double(prodhm), PACKAGE="pdfCluster"),
        t7= .C("c_kepdft7", as.double(x), as.double(rec.hm), as.double(x), 
              as.integer(nx), as.integer(ndim), as.integer(nx), 
              double(nx), as.double(prodhm), PACKAGE="pdfCluster") 			  
        )
        #cons <- switch (kernel, gaussian = 1/((2 * pi)^(ndim/2) * nx), biweight = (15/16)/nx, triweight = (35/32)/nx)
        #g <- exp(mean(log(pdfpil[[7]]*cons)))
    #lambda <- ((pdfpil[[7]]*cons)/g)^(-alpha)
    g <- exp(mean(log(pdfpil[[7]])))
    lambda <- ((pdfpil[[7]])/g)^(-alpha)
    hm <- outer(lambda, h)
    hm <- matrix(as.double(hm), nx, ndim)
    hm
 }
		

summary.kepdf <- function (object, ..., props=c(75,50,25)) {
	#compute
	if(dim(object@eval.points)[1]==dim(object@x)[1] && all(object@x==object@eval.points)) pdf.obs <- object@estimate else pdf.obs <- kepdf(object@x, bwtype=object@bwtype, kernel=object$kernel, h=object@par$h, hx = object@par$hx, eval.points=object@x)@estimate 
		hts <- quantile(pdf.obs, prob = (100 - props)/100)
	cat("An S4 object of class \"", class(object),"\"","\n","\n",sep="")
	modepos = which.max(pdf.obs)				
	cat("The highest density data point has position", modepos, "in the sample data", "\n","\n",sep=" ")
	indices <- list()
	for(i in 1:length(props)){
		indices[[i]] <- (which(pdf.obs>hts[i]))
		cat("Rows of ", props[i], "% top density data points:", indices[[i]],"\n","\n",sep=" ")
	}				
	#names(indices) <- paste(as.character(props),"%",sep="")
	#new("summary.kepdf", obj.class=class(object)[1], mode = modepos, props = props, indices = indices)
    #show 
    invisible(object)
    }

#setClass("summary.kepdf", representation(obj.class = "character", mode = "numeric", props = "numeric", indices = "list"))

setGeneric("summary", function(object, ...)standardGeneric("summary"))
setMethod("summary", "kepdf", summary.kepdf)


#show methods
setMethod("show","kepdf", function(object){
	cat("An S4 object of class \"",class(object),"\"","\n","\n",sep="")
	cat("Call: ")
	print(object@call)
	cat("\n")
	cat("Kernel: ", "\n")
	print(object@kernel)
	cat("\n")
	if(object@bwtype == "fixed") {
		cat("Estimator type: fixed bandwidth", "\n")
		cat("\n")
		if(ncol(object@x)==1) 
		cat("Smoothing parameter: ", object@par$h, "\n") else 
		cat("Diagonal elements of the smoothing matrix: ", object@par$h, "\n")
		} else {
		cat("Estimator type: adaptive bandwidth", "\n")
		cat("\n")
		cat("Smoothing parameters of observations: ", "\n")
		print(object@par$hx)	
		}
	cat("\n")
	cat("Density estimate at evaluation points: ", "\n")
	print(object@estimate)
	cat("\n")
	})



plot.kepdf1d <- function(x, eval.points=NULL, n.grid=NULL, data=NULL, add=FALSE,
                         main	= NULL, xlab=NULL, ylab=NULL, col=NULL, col.data=2, type="l",...){
  ##########################################
  ##########################################
  # comments by GM 03/02/2011
  # plot of univariate kernel density estimate 
  # DEFAULT: density is evaluated on a grid defined on the range 
  # of evaluation points of the density
  ##########################################
  ##########################################
  if(is.null(n.grid)) n.grid <- 30
  if(is.null(main))	main="Kernel density estimate of data"
  if(is.null(xlab))	xlab="eval.points"
  if(is.null(ylab))	ylab="density function"
  if(is.null(eval.points)) {
    datarange <- apply(x@eval.points, 2, extendrange)
    eval.points <- seq(datarange[1, 1], datarange[2, 1], length=n.grid)
  } 
  pdf <- kepdf(x@x, eval.points=eval.points, kernel = x@kernel, 
               bwtype = x@bwtype, h=x@par$h, hx=x@par$hx)
  pdf.est.sort <- pdf@estimate[order(eval.points)]
  eval.points.sort=sort(eval.points)
  if(is.null(col)) col=1
  if (add){lines(eval.points.sort, pdf.est.sort, col=col,...)} else {
    plot(eval.points.sort, pdf.est.sort, type=type, main=main, xlab=xlab, 
         ylab=ylab, col=col, ...)}
  if (!is.null(data)) rug(data, col=col.data)
  invisible(list(eval.points=as.vector(eval.points), estimate=pdf@estimate))
}

plot.kepdf2d <- function (x, eval.points=NULL, n.grid=NULL, data=NULL,
                          add=FALSE, main=NULL, xlab=NULL, ylab=NULL, zlab=NULL, col=NULL, 
                          col.data=2, props=c(75,50,25), method="contour", ticktype = "detailed",
                          ...) {
  if(is.null(n.grid)) n.grid <-c(30,30)
  if(is.null(main)) main="Kernel density estimate of data"
  if(is.null(eval.points)) {  
    datarange <- apply(x@eval.points, 2, extendrange, f= 0.15)
    xgrid <- seq(datarange[1, 1], datarange[2, 1], length=n.grid[1]) 
    # sort(x@eval.points[,1])
    ygrid <- seq(datarange[1, 2], datarange[2, 2], length=n.grid[2]) 
    # sort(x@eval.points[,2])
  } else { 
    xgrid=sort(eval.points[, 1])
    ygrid=sort(eval.points[, 2])
  }
  method <- match.arg(method,choices=c("contour","perspective","image"))
  nx <- length(xgrid)
  ny <- length(ygrid)
  grid.points <- as.matrix(expand.grid(xgrid,ygrid))
  pdf <- kepdf(x@x, eval.points=grid.points, kernel = x@kernel, 
               bwtype=x@bwtype, h=x@par$h, hx=x@par$hx)
  lab <- colnames(x@x) 
  if(is.null(xlab)) if(is.null(lab)) xlab <- "Var 1" else xlab <- lab[1]
  if(is.null(ylab)) if(is.null(lab)) ylab <- "Var 2" else ylab <- lab[2]
  if(method=="image"){
    if(is.null(col)) col=heat.colors(12)		
    image(xgrid, ygrid, matrix(pdf@estimate, nx,ny), add = add, main=main,
          col=col, xlab=xlab, ylab=ylab, ...)
  } else if(method=="perspective"){
    if (is.null(zlab)) zlab="density function"
    if(is.null(col)) col="lightblue"		
    persp(xgrid,ygrid, matrix(pdf@estimate, nx,ny), main=main, xlab = xlab,
          ylab = ylab, zlab=zlab, col=col, ticktype = "detailed", ...)
  } else {
    if(!add) plot(grid.points, type = "n", main=main, xlab = xlab, 
                  ylab = ylab, ...)
    if(dim(x@eval.points)[1]==dim(x@x)[1] && all(x@x==x@eval.points)) pdf.obs <- x@estimate else pdf.obs <- kepdf(x@x, eval.points=x@x, kernel = x@kernel, bwtype=x@bwtype, h=x@par$h, hx=x@par$hx)@estimate
    hts <- quantile(pdf.obs, prob = (100 - props)/100)
    if(is.null(col)) col=1		
    for (i in 1:length(props)) {
      scale <- props[i]/hts[i]
      contour(xgrid, ygrid, matrix(pdf@estimate, nx,ny) * scale, 
              level = hts[i] *scale, add = TRUE, main=main, col=col, ...)
    }
    if(!is.null(data)) points(data, col=col.data, ...)#aggiunge i punti	
  }
  invisible(list(eval.points = pdf@eval.points, estimate = pdf@estimate))
}		


plot.kepdfdd <- function (x, eval.points= NULL, n.grid=NULL, data=NULL,
                          main=NULL, method="contour", indcol=NULL,
                          text.diag.panel = NULL, gap = 0.5, col=NULL, col.data=2, ...) {
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) 
      text(x, y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, xpd, bg, main,  oma, ...) {
      if (side%%2 == 1) Axis(x, side = side, xpd = NA, ...) else 
      Axis(y, side = side, xpd = NA, ...)
      }
  localPlot <- function(..., oma, font.main, cex.main) plot.kepdf2d(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is(x, "kepdf")) stop("argument must be of class \"kepdf\"")
  nc <- ncol(x@x)
  if (is.null(text.diag.panel)) if(!is.null(lab <- colnames(x@x))) text.diag.panel <- lab else text.diag.panel <- paste(rep("V", nc),1:nc,sep="")
  oma <- if ("oma" %in% nmdots) dots$oma else NULL
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) 
      oma[3L] <- 6
  }    
  #if (!is.null(main)) main="Pairwise marginal kernel density estimates of data"
  if(is.null(indcol)) indcol <- 1L:nc
  opar <- par(mfrow = c(length(indcol), length(indcol)), mar = rep(c(gap,gap/2),each=2),
 oma=oma)
  on.exit(par(opar))
  out <- list()
  count <- 1
  method <- match.arg(method,choices=c("contour","perspective","image"))
  for (i in indcol){ 
    for (j in indcol) 
      if(i==j) {
        plot(1,type="n",xlab="",ylab="", axes=FALSE)
        text(1,1,text.diag.panel[i],cex=2)
        box()} else {
          if (!is.null(data)) dataji <- data[,c(j,i)] else dataji<-data
          if (!is.null(eval.points)) eval.pointsji <- eval.points[,c(j,i)] else eval.pointsji<-eval.points
          if (!is.null(n.grid)) n.gridji <- n.grid[c(j,i)] else n.gridji<-n.grid
          if (!is.null(x@par$hx)) hxji <- matrix(x@par$hx, ncol=nc)[,c(j,i)] else hxji <-NULL 
          marg <- kepdf(x@x[,c(j,i)], kernel = x@kernel, bwtype=x@bwtype, h = x@par$h[c(j,i)], hx = hxji)
          out[[count]] <- localPlot(x=marg, eval.points = eval.pointsji,
                                    n.grid=n.gridji, data=dataji, xlab = "",
                                    ylab = "", main="", yaxt="n", xaxt="n",
                                    method=method, col=col, col.data=col.data,...)
          count <- count + 1     
          if (method!="perspective"){
            if(i==indcol[1]) {axis(3); if(j==length(indcol)) axis(4)}
            if(j==indcol[1]) {axis(2); if(i==length(indcol)) axis(1)}
            box()} #else if (i == j) {j=j+1}
        }
    count <- count + 1}			
  par(new = FALSE)
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) dots$font.main
    else par("font.main") 
    cex.main <- if ("cex.main" %in% nmdots) dots$cex.main
    else par("cex.main") 
    mtext(main, side=3, TRUE, line = 3, outer = TRUE, at = NA, cex = cex.main, font = font.main, adj= 0.5)
  }	
  invisible(out)
}

plot.kepdf <- function(x, y, eval.points = NULL, n.grid = NULL, data = NULL,
                       add = FALSE, main = NULL, xlab = NULL, ylab = NULL, zlab = NULL, col=NULL,
                       col.data=2, type="l", props = c(75,50,25), method="contour", 
                       ticktype = "detailed", indcol = NULL, text.diag.panel = NULL, 
                       gap = 0.5, ...){
  d <- NCOL(x@x)
  if (d==1) plot.kepdf1d(x, eval.points = eval.points, n.grid = n.grid,
                         data=data, add=add, main = main, xlab = xlab, ylab = ylab,
                         col= col, col.data = col.data, type=type, ...) else 
  if (d==2) plot.kepdf2d(x, eval.points = eval.points, n.grid = n.grid,
                         data=data, add=add, main=main, xlab=xlab, ylab=ylab, 
                         zlab=zlab, col=col, col.data = col.data, props = props, 
                         method=method, ticktype = ticktype,	...) else 
  plot.kepdfdd(x, eval.points = eval.points, n.grid = n.grid, data=data,
               main = main, props= props, method = method, indcol = indcol,
               text.diag.panel = text.diag.panel, gap = gap, col=col, col.data=col.data,...) 
}

setMethod("plot", signature(x="kepdf", y="missing"), plot.kepdf)		

	
##########################################################
##### functions, classes and methods related to clustering
##########################################################

setOldClass("dendrogram")  

setClass("pdfCluster", representation(call="call", x="matrix", pdf="list", 
	nc="list",  graph= "list", cluster.cores="ANY", tree="dendrogram", noc="numeric", stages="ANY", clusters="ANY"
	))


#library(geometry)

pdfCluster.data <- function (x, graphtype, hmult, Q="QJ", lambda = 0.10, grid.pairs = 10, n.grid=min(round((5 + sqrt(NROW(x)))*4), NROW(x)), ... ) 
{
 call<-match.call()
 x<-data.matrix(x)
 if(any(sapply(x[1,],is.factor))) stop("All variables must be numeric")
 if(any(apply(x,2,check.discrete)<10)) warning("One or more variables look to be discrete: pdfCluster is designed for continuous data")
 if(any(is.na(x))) {warning(cat("NA in object",deparse(substitute(x)), "have been omitted", "\n"));  x<- na.omit(x)}
  #manage dots
	dots <- list(...)
	ndots <- names(dots)
	args.kepdf <- ndots%in%(names(formals(kepdf)))
        dots.kepdf <- dots[args.kepdf]
        args.pdfClassification <- ndots%in%(names(formals(pdfClassification)))
        dots.pdfClassification <- dots[args.pdfClassification]
	#x <- data.matrix(x)  
	if (missing(hmult)) if (ncol(x)>6) hmult=1 else hmult=0.75          
	if ("h" %in% names(dots.kepdf)) dots.kepdf$h <-dots.kepdf$h*hmult else dots.kepdf$h <- hmult*h.norm(x)	
	pdf <- do.call(kepdf, c(list(x = x, eval.points = x), dots.kepdf))
	estimate <- pdf@estimate
        #check given arguments 
	hasQ<-hasArg(Q)
	haslambda<-hasArg(lambda)
	hasgrid.pairs<-hasArg(grid.pairs)
	n <- nrow(x)
	d <-ncol(x)
    if (n.grid > n) {
        warning("n.grid too large, set equal to n")
        n.grid <- min(n.grid, nrow(x))
    	}
	# build graph
	if (missing(graphtype)) if(d==1) graphtype ="unidimensional" else if (d <= 6) graphtype ="delaunay" else graphtype ="pairs"
	graph.par <- list()
	graph.par$type <- graphtype
	if (NCOL(x)==1) {
		if(hasQ) message("Unused argument 'Q' when graphtype= 'unidimensional'")
		if(haslambda) message("Unused argument 'lambda' when graphtype= 'unidimensional'")
		if(hasgrid.pairs) message("Unused argument 'grid.pairs' when graphtype= 'unidimensional'")
		graph.nb <- graph.uni(x) 
		} else  
	if (graphtype == "delaunay") {
		if(haslambda) message("Unused argument 'lambda' when graphtype= 'delaunay'")
		if(hasgrid.pairs) message("Unused argument 'grid.pairs' when graphtype= 'delaunay'")
      if(NCOL(x)> 6) warning("Running pdfCluster with graphtype='delaunay' may require a long time when the number of variables is so large")		
		graph.nb <- graph.del(x, Q = Q) 
		} else
	if (graphtype == "pairs") {
		if(hasQ) message("Unused argument 'Q' when graphtype= 'pairs'")
		graph.nb <- graph.pairs(x, pdf, lambda = lambda, grid.pairs=grid.pairs)
		graph.par$lambda = lambda; graph.par$comp.area <- graph.nb$comp.areas
		} else 
	stop("graphtype should be one of 'unidimensional', 'delaunay', or 'pairs'")
	nc <- num.con(x, estimate, graph.nb$graph, profile.grid = n.grid-2, correct=TRUE)  
  	    #determine connected components
        struct <- con2tree(nc, estimate)
    if (struct$bad) {
        message("No output given")
		      if (graphtype == "pairs") message ("The grid is too coarse: re-run with a larger value of one among 'n.grid', 'hmult' or 'lambda'. See the 'Warning' Section in 'help(pdfCluster)' for further details.") else message("The grid is too coarse: re-run with larger value of 'n.grid' or hmult. See the 'Warning' Section in 'help(pdfCluster)' for further details.")
      }
	    else {
    	g <- struct$g
        g[struct$g == 0] <- NA
		pdf.comp <-list()
		pdf.comp$kernel <- pdf@kernel
		pdf.comp$bwtype <- pdf@bwtype
		pdf.comp$par <- pdf@par
		pdf.comp$estimate <- estimate
		out <- new("pdfCluster", call = call, x = data.matrix(x), 
            pdf = pdf.comp,  nc = nc, graph = graph.par, cluster.cores = g, tree = struct$tree, noc = struct$noc, 
            stages = NULL, clusters = NULL)
        if ((!"n.stage" %in% ndots)||("n.stage" %in% ndots && dots.pdfClassification$n.stage > 0)) {
           out <- do.call(pdfClassification, c(list(obj = out), dots.pdfClassification))
}
        out
    }
}

# function to find clusters starting from a graph of data

#pdfCluster.graphpairs <- function (x, lambda = 0.10, n.grid=min(50, nrow(as.matrix(x@x))), ... ) 
pdfCluster.graphpairs <- function (x, graphtype, hmult, Q, lambda=0.10, grid.pairs, n.grid=min(round((5 + sqrt(NROW(x@x)))*4), NROW(x@x)), ...)
{
    #manage dots
	dots <- list(...)
	ndots <- names(dots)
        args.pdfClassification <- ndots%in%(names(formals(pdfClassification)))
        dots.pdfClassification <- dots[args.pdfClassification]
	#x is an object of class pdfCluster 
	#containing a graph built with graph.pairs
	if (n.grid > nrow(x@x)) {
        warning("n.grid too large, set equal to n")
        n.grid <- min(n.grid, nrow(x@x))
    }  
	# build graph
	if (x@graph$type!= "pairs") stop (cat(deparse(substitute(x)), "should contain a graph of type 'pairs'"))
	graph.nb <- area.to.graph(x@graph$comp.area, lambda = lambda)	
	graph.par <- list()
	graph.par$type <- graph.nb$type
	graph.par$comp.area <- graph.nb$comp.areas
	graph.par$lambda=lambda
	nc <- num.con(x@x, x@pdf$estimate, graph.nb$graph, profile.grid = n.grid-2, correct=TRUE)  
  	    #determine connected components
        struct <- con2tree(nc,x@pdf$estimate)
    if (struct$bad) {
        message("No output given")
		message ("The grid is too coarse: re-run with a larger value of one among 'n.grid', 'hmult' or 'lambda'. See the 'Warning' Section in 'help(pdfCluster)' for further details.")
    }
	    else {
    	g <- struct$g
        g[struct$g == 0] <- NA
		pdf.comp <-x@pdf
		out <- new("pdfCluster", call = match.call(), x = data.matrix(x@x), 
            pdf = pdf.comp,  nc = nc, graph = graph.par, cluster.cores = g, tree = struct$tree, noc = struct$noc, 
            stages = NULL, clusters = NULL)
        if ((!"n.stage" %in% ndots)||("n.stage" %in% ndots && dots.pdfClassification$n.stage > 0)) {
           out <- do.call(pdfClassification, c(list(obj = out), dots.pdfClassification))
}
        out
    }
}

#setGeneric("pdfCluster", function(x, graphtype, Q="QJ",   ...) standardGeneric("pdfCluster"))
#setMethod("pdfCluster", signature(x="matrix", graphtype = "ANY", Q = "ANY"), pdfCluster.data)
#setMethod("pdfCluster", signature(x="data.frame", graphtype = "ANY", Q = "ANY"), pdfCluster.data)
#setMethod("pdfCluster", signature(x="numeric", graphtype = "ANY", Q = "ANY"), pdfCluster.data)
#setMethod("pdfCluster", signature(x="pdfCluster", graphtype = "missing", Q = "missing"), pdfCluster.graphpairs)

##versione 1:0-0
#setGeneric("pdfCluster", function(x, graphtype, hmult, Q="QJ", lambda = 0.10, grid.pairs = 10, n.grid = min(round((6 + sqrt(NROW(x)))*4), NROW(x)),...) standardGeneric("pdfCluster"))
#setMethod("pdfCluster", signature(x="matrix", graphtype="ANY", hmult= "ANY", Q="ANY", lambda = "ANY", grid.pairs = "ANY", n.grid="ANY"), pdfCluster.data)
#setMethod("pdfCluster", signature(x="data.frame", graphtype="ANY", hmult= "ANY", Q="ANY", lambda = "ANY", grid.pairs = "ANY", n.grid="ANY"), pdfCluster.data)
#setMethod("pdfCluster", signature(x="numeric", graphtype="ANY", hmult= "ANY", Q="ANY", lambda = "ANY", grid.pairs = "ANY", n.grid="ANY"), pdfCluster.data)
#setMethod("pdfCluster", signature(x="pdfCluster", graphtype="missing", hmult= "ANY", Q="missing",lambda="ANY", grid.pairs="missing", n.grid="ANY"), pdfCluster.graphpairs)

##versione 1:0-1
#setGeneric("pdfCluster",function(x, graphtype, hmult, Q="QJ", lambda = 0.10, grid.pairs = 10, n.grid=min(round((6 + sqrt(NROW(x)))*4), NROW(x)),...) standardGeneric("pdfCluster"))
#setMethod("pdfCluster", signature(x="matrix", graphtype="character", hmult="numeric", Q="character", lambda = "numeric", grid.pairs = "numeric", n.grid="numeric"), pdfCluster.data)
#setMethod("pdfCluster", signature(x="data.frame", graphtype="character", hmult="numeric", Q="character", lambda = "numeric", grid.pairs = "numeric", n.grid="numeric"), pdfCluster.data)
#setMethod("pdfCluster", signature(x="numeric", graphtype="character", hmult="numeric", Q="character", lambda = "numeric", grid.pairs = "numeric", n.grid="numeric"), pdfCluster.data)
#setMethod("pdfCluster", signature(x="pdfCluster", graphtype="missing", hmult="missing", Q="missing",lambda="numeric", grid.pairs="missing", n.grid="numeric"), pdfCluster.graphpairs)

#versione 1:0-1
setGeneric("pdfCluster", function(x, graphtype, hmult, Q="QJ", lambda = 0.10, grid.pairs = 10, n.grid = min(round((5 + sqrt(NROW(x)))*4), NROW(x)),...) standardGeneric("pdfCluster"))
setMethod("pdfCluster", signature(x="matrix"), pdfCluster.data)
setMethod("pdfCluster", signature(x="data.frame"), pdfCluster.data)
setMethod("pdfCluster", signature(x="numeric"), pdfCluster.data)
setMethod("pdfCluster", signature(x="pdfCluster"), pdfCluster.graphpairs)

#setClass("summary.pdfCluster", representation(obj.class="character", cluster.cores="numeric", clusters="ANY", tree="dendrogram"))

#summary.pdfCluster <- function(object, ...){
summary.pdfCluster <- function(object,...){
    noc <-object@noc
    t.cluster.cores<-table(groups(object, stage=0),exclude=NULL)
    t.clusters<-table(groups(object))
    maxdigcc<-nchar(max(t.cluster.cores))	
    maxdigc<-nchar(max(t.clusters))	
    cat("An S4 object of class \"", class(object), "\"","\n","\n",sep="")
    cat("Call: ")
    print(object@call)
    cat("\n")
    cat("Initial groupings:","\n")
    cat(" label ", format(c(1:noc, NA), width=maxdigcc),"\n","count ", format(t.cluster.cores, width=maxdigcc),"\n")
    cat("\n")
	cat("Final groupings:","\n");
    cat(" label ", format(1:noc, width=maxdigc),"\n","count ", format(t.clusters, width=maxdigc),"\n")
	cat("\n")
	cat("Groups tree (here 'h' denotes 'height'):\n")
	str(object@tree)
	invisible(object)
}

setMethod("summary",signature("pdfCluster"), summary.pdfCluster)


setMethod("show",signature("pdfCluster"), function(object){
	cat("Clustering via nonparametric density estimation", "\n","\n")
	cat("Call: ")
    print(object@call)
    cat("\n")
	#h  <-  object@hmult*object@h
	#if(ncol(object@x)==1) 
	#cat("Smoothing parameter: ", h, "\n", "\n")
	#else 
	#cat("Diagonal elements of the smoothing matrix: ", h, "\n","\n")
	#cat("Density estimate at data points: ", "\n")
	#print(object@estimate)
	#cat("\n")
    cat("Kernel estimator type:\n")
		print(object@pdf$bwtype)
		cat("\n")
	cat("Graph type:\n")
	    print(object@graph$type)
		cat("\n")
	cat("Groups tree (here 'h' denotes 'height'):\n")
	str(object@tree)
	cat("\n")
  	cat("Initial groupings:", "\n")
      print(groups(object, stage=0))
    cat("\n")
	if(!is.null(groups(object))){
    if (length(object@stages)>1) for(i in seq(1:(length(object@stages)-1))){
	cat("Stage", i, "groupings:","\n")
    print(object@stages[[i]])
    cat("\n")}
	cat("Final groupings:", "\n")
      print(groups(object))
    cat("\n")
		}
   })

groups <- function(obj, stage=length(obj@stages)) {
 if(class(obj)!="pdfCluster") stop("Function 'groups' extracts groups from objects of class 'pdfCluster'")
 if (stage==0) out <- obj@cluster.cores else {if (is.null(obj@stages)) out <- NULL else out <- obj@stages[[stage]]}
 out
	}

pdfClassification <- function (obj, n.stage = 5, se=TRUE, hcores=FALSE) 
{
    x <- obj@x
    lista <- list()
    g <- obj@cluster.cores
    n <- length(g)
    g[is.na(g)] <- 0
    M <- obj@noc
    kernel <- obj@pdf$kernel 
    if (obj@noc == 1) {
        obj@cluster.cores[g == 0] <- 1
        for (stage in 1:n.stage) lista[[stage]] <- rep(1, nrow(x))
    } else {
        stage <- 1
        while (stage <= n.stage) {
            unallocated <- which(g == 0)
			alp<-as.vector(table(g[g!=0]))/sum(g!=0)
            f <- matrix(NA, nrow = sum(g == 0), ncol = M)
			se2.logf <- matrix(NA, nrow = sum(g == 0), ncol = M)
            if (obj@pdf$bwtype == "fixed")
			{ 
                h0 <- obj@pdf$par$h
				if (hcores){
					for (m in 1:M) {
						indm <- which(g==m)
						h <- h0 * (1 + (stage)/n.stage)
						f[, m] <- kepdf(x = x[indm, ], kernel=kernel, bwtype="fixed", h = h, eval.points = matrix(x[g == 0, ], ncol=ncol(x)))@estimate
						se2.logf[, m] <- 1/(f[, m]*length(indm)*prod(h))	
					}
				} else
					for (m in 1:M) {
					indm <- which(g==m)
					h1 <- max(h.norm(x[g==m,]), 0.0001)
					h <-exp((1-alp[m])*log(h0)+(alp[m])*log(h1))
					f[, m] <- kepdf(x = x[g == m, ], kernel=kernel, bwtype="fixed", h = h, eval.points = matrix(x[g == 0, ], ncol=ncol(x)))@estimate
					se2.logf[, m] <- 1/(f[, m]*length(indm)*prod(h))	
					}				
			} else
			if (obj@pdf$bwtype == "adaptive") 
			{
				h0 <- obj@pdf$par$h 
				hx <- obj@pdf$par$hx
				alpha=obj@pdf$par$alpha
				if (hcores) {
					for (m in 1:M) {
						indm <- which(g==m)
						hx <- hx* (1 + (stage)/n.stage)
						f[, m] <- kepdf(x = x[indm, ], kernel=kernel, bwtype="adaptive", hx = hx[indm, ], eval.points = matrix(x[g == 0, ], ncol=ncol(x)))@estimate
						se2.logf[, m] <- 1/(f[, m]*length(indm)*prod(h0))	
					}
				} else 
					for (m in 1:M) {
						indm <- which(g==m)
						h0x<-hx[indm, ]
						h1x <- hprop2f(x=x[indm,], h=h.norm(x), alpha)
						hxmed <-exp((1-alp[m])*log(h0x)+alp[m]*log(h1x))
						f[, m] <- kepdf(x = x[indm, ], kernel=kernel, bwtype="adaptive", hx = hxmed, eval.points = matrix(x[g == 0, ], ncol=ncol(x)))@estimate
						se2.logf[, m] <- 1/(f[, m]*length(indm)*prod(h0))	
					}
			} 
            unalloc.outliers <- (apply(f, 1, sum) == 0)
            unallocated <- unallocated[!unalloc.outliers]
            f <- matrix(f[!unalloc.outliers, ], ncol = M)
            ind.max <- t(apply(f, 1, order)[rev(1:M), ])
            i1 <- ind.max[, 1]
            i2 <- ind.max[, 2]
            f1 <- diag(as.matrix(f[, i1]))
			f2 <- diag(as.matrix(f[, i2]))
            log.ratio <- log(f1) - log(f2)
            se2.logf1 <- diag(as.matrix(se2.logf[, i1])) 
            se2.logf2 <- diag(as.matrix(se2.logf[, i2]))
			se.log.ratio <- sqrt(se2.logf1 + se2.logf2)
			if (se) log.ratio.sc <- (log(f1) - log(f2))/se.log.ratio else log.ratio.sc <- log.ratio
			if (stage < n.stage) {
                alti <- (log.ratio.sc >= quantile(log.ratio.sc, (n.stage - 
                  stage)/n.stage,na.rm = TRUE))
            } else {
                alti <- rep(TRUE, (sum(g == 0) - sum(unalloc.outliers)))
            }
            g.alti <- ind.max[alti, 1]
            for (m in 1:M) {
                nuovi <- unallocated[alti][g.alti == m]
                g[nuovi] <- m
            }
            lista[[stage]] <- g
            if (sum(g == 0) == 0 & stage < n.stage) {
                message(paste("Classification accomplished in", 
                  stage, "stages"))
                n.stage <- stage
            }
            stage <- stage + 1
        }
        if (any(unalloc.outliers)) 
            warning(cat(sum(unalloc.outliers), "outliers have been found and not classified", 
                "\n"))
    }
    unlista <- unlist(lista)
    unlista[unlista == 0] <- NA
    obj@stages <- split(unlista, rep(1:n.stage, each = nrow(x)))
    obj@clusters <- lista[[n.stage]]
    obj
}

con2tree <- function(object, f){
  # object is the output of num.con()
  # f density estimate
  # 
  ow <- options("warn")
  nc <- object$nc
  p <- object$p
  index <- which(diff(nc) != 0)         # posizione punti di salto
  K <- length(index)  					# numero di punti di salto
  ps <- p[index]      					# frazione di dati inclusi ai vari punti di salto
  lista <- list()     					
  if (K == 1){
  gruppi <- as.vector(object$id[,ncol(object$id)])
  M <-1 } else {
  for(j in 1:K) lista[[j]] <-  as.vector(object$id[,index[j]])  # list of connected sets at the jumps
   gruppi <- rep(0,length(lista[[1]]))    
  M <- 0                                 
  insiemi <- list()
  # step
  k <- 1
  allocated <- (gruppi>0)
  insiemi[[k]] <- setdiff(unique(gruppi),0)
  # loop k=2,...,K
  while(k < K) {
    k <- k+1
    sets  <- lista[[k]]                #elenco connessi al salto k 
    insieme <- list()                   
    for(m in 1:max(sets)) {
      set <- which(sets==m)            
	#which objects are in connected component m?         
      new <- setdiff(set, which(allocated))    
	#which objects in connected component m are not allocated yet?
      if(length(new)>0){               
        #are there new objects in the connected component?  
        g.new <- unique(gruppi[intersect(set, which(allocated))]) 
	# which is the label group of objects in connected component m and has been already allocated?
        if(length(g.new) == 0)  gruppi[set] <- M <- M+1           
        if(length(g.new) == 1)  gruppi[set] <- g.new 
        allocated <- (gruppi>0)
        }
      gg <- sort(setdiff(unique(gruppi[set]),0))
      if(length(gg)>0) insieme[[length(insieme)+1]] <- gg 
    }
    insiemi[[k]] <- insieme
    }
  if(!missing(f)){
    u <- unique(gruppi[rev(order(f))])
    g <- rep(0,length(f))
    u0 <- u[u>0]
    for(i in 1:max(gruppi)) g[gruppi==u0[i]] <- i
    gruppi <- g
    } 
  }
  salti <- diff(nc)
  salta.giu <- rev(which(diff(nc)[index]<0))
  altezza <- numeric(M)
  m <- 0
  salti.su <- salti[salti>0]
       options(warn = -1)
  while(m <M){
    m <- m+1
    r <- min(which(cumsum(salti.su) >= m))
    altezza[m] <- p[salti>0][r]
  }
  bad.grid <- any(is.na(altezza))
  sotto.albero <- function(tree, k, set){
      insieme <- insiemi[[salta.giu[k]]]      
      r <- 0
      branch <- list()
      for(item0 in insieme){
        item <- intersect(set, unlist(item0))
        if(length(item)==1){
          r <- r+1
          u <- item
          attr(u,"members") <- 1
          attr(u,"height") <- altezza[item]
          attr(u,"label")  <- paste(as.character(item), " ", sep="")
          attr(u,"leaf")   <- TRUE
          branch[[r]] <- u
        }     
      if(length(item) > 1) {
        r <- r+1
        u <- sotto.albero(list(), k+1, item)
        attr(u,"members") <- length(unlist(item))
        attr(u,"height") <- max(ps[salta.giu[k+1]])
        attr(u,"label")  <- paste("{", paste(item, collapse=","),"}", sep="")
        attr(u,"leaf")   <- FALSE   
        branch[[r]] <- u
        # browser()
      }
      }
      branch
    }        
  if(M > 1) {
	tree <- sotto.albero(list(), 1, 1:M)
	attr(tree, "members") <- M
	attr(tree, "height") <-  max(ps[salta.giu[1]])
	attr(tree, "label") <- paste("{", paste(1:M, collapse=","),"}", sep="")
	noc <- M } else {
	#bassi.connect<- object$id[ ,(ncol(object$id)-1)] 
	# when there is one group only, cluester cores are all the allocated objects at the lower considered threshold of density
	#gruppi[bassi.connect>0] <- bassi.connect[bassi.connect>0] 
	tree <- list()
	tree[[1]]<-1
	attr(tree[[1]], "members") <- 1
	attr(tree[[1]], "height") <- 0
	attr(tree[[1]], "label") <- "1"
	attr(tree[[1]], "leaf") <- TRUE
    attr(tree,"members")<-1
	attr(tree,"height")<-1
	attr(tree,"label")<-paste("{", paste(1, collapse=","),"}", sep="")
	noc <- 1
	}
	tree <- list(tree)
	attr(tree, "members") <- M
	attr(tree, "height") <- 1
    attr(tree, "class") <- "dendrogram"
  options(warn = ow$warn)
  invisible(list(g=gruppi, tree=tree, bad=bad.grid, noc=noc))
 }
 
#####
plot.pdfCluster <- function(x, y, which=1:4, stage=Inf, pch=NULL, col=NULL, ...){
  if (is.element(4, which) & x@noc==1) which <- setdiff(which, 4)
  w12 <- sort(which)[1:min(2,length(which))]
  if(setequal(1:2,w12)) par(mfrow=c(1,2))
  if(is.element(1, which)){
     nc <- x@nc$nc
     p <- x@nc$p
     plot(c(0,p,1,1), c(0,nc,nc[length(nc)],0), 
           ylim=c(0,max(nc)), yaxt="n", type="S", cex=0.7,
           xlab="fraction of data points",
           ylab="mode function", ...)
     axis(2,at=c(0:max(nc)))
  }
  if(is.element(2, which)) {
    plot(x@tree, center=TRUE, ...)
    title(sub="groups tree")
    }
  if(setequal(1:2,w12)) par(mfrow=c(1,1))
 # if(length(setdiff(which,w12)) > 0) {
  #if(length(setdiff(w12,c(3,4))) > 0 ) {
 #    cat("press <enter> to continue...")
 #    readline()
 #  }
  if(is.element(3, which))    {
    #browser()  
	if(sum(which<3)) {
     cat("press <enter> to continue...")
     readline()
   }
    stage <- min(stage, length(x@stages))
    if(stage==0) {g <- x@cluster.cores; g[is.na(x@cluster.cores)] <- 0}
    else {g <- x@stages[[stage]]; g[is.na(x@stages[[stage]])] <- 0} 
	dat <- x@x
    M <-  length(table(g))
    if (is.null(pch)) #|(length(pch) < M & length(pch) > 1)) {
	pch <- as.character(g)
	if (length(pch) < M & length(pch) > 1) warning("pch argument should have length 1 or equal to the number of clusters")
	if(length(pch) > 1){ 
	newpch <- as.factor(g)
	levels(newpch) <- rep(levels(as.factor(pch)), length(levels(as.factor(pch))))[1:M]
	if (is.numeric(pch)) pch <- as.numeric(as.character(newpch)) else pch <- as.character(newpch)
	} 
	if (is.null(col)){#|(length(col) < M & length(col) > 1)) {
	col=g 
	if (0 %in% g) {
	g0<- g
	g0[g0==0]<-x@noc+1	
	col=g0
	}	
	} 
	if (length(col) < M & length(col) > 1) warning("col argument should have length 1 or equal to the number of clusters")
	if(length(col) > 1){ 
	newc <- col
	newc <- as.factor(newc)
	levels(newc) <- rep(levels(as.factor(col)), length(levels(as.factor(col))))[1:M]
	if (is.numeric(col)) col <- as.numeric(as.character(newc)) else col <- as.character(newc)
	}
	if(ncol(dat)==1){
	if(is.null(list(...)$add))  plot(dat[,1], pch=pch, col=col, cex=.75, ylab="x", ...)
#	if(is.null(list(...)$add))  plot(dat[,1], type="n", ...)
#    for(m in 0:M) {
#      gr <- (g==m)
#      if(m==0) points(seq(along=x)[gr], dat[gr,],  pch=0, cex=0.2)
#      else points(seq(along=x)[gr], dat[gr,], pch=m, col=m+1, cex=0.75)
#      } 
	} else if (ncol(dat)==2){
if(is.null(lab <- colnames(dat))) {lab <- paste(rep("V", ncol(dat)), 1:ncol(dat), sep="")}	
if(is.null(list(...)$add))  plot(dat, pch=pch, col=col, xlab=lab[1], ylab=lab[2], cex=.75, ...)
#	if(is.null(list(...)$add))  plot(dat[,coord], type="n", ...)
#    for(m in 0:M) {
#      gr <- (g==m)
#      if(m==0) {pch=pch[1]; points(dat[gr,coord],  pch=pch , cex=0.2)}
#      else points(dat[gr,coord], pch=pch[m], col=m, cex=0.75)
#      }
	  } else pairs(dat, pch=pch, col=col, ...)
  }
    if(is.element(4, which))    {
if(sum(which<4)) {
     cat("press <enter> to continue...")
     readline()
   }
   plot(dbs(x))
   }
  #browser()
  invisible()
}

setMethod("plot",signature(x="pdfCluster", y="missing"), plot.pdfCluster)

###############################################################################
##### functions, classes and methods related to graph building and CC detection
###############################################################################
  
graph.del <- function(x, Q="QJ") {
   # ----------------------------------------------------------------------
   # Crea la triangolazione e elenca i VICINI PER OGNI PUNTO
 
   # REQUIRES geometry
   # Options implemented for qdelaunay.exe:
   # Qt  - all regions are simplicial (in 2D triangles)
   # Qj  - joggle input all regions are simplicial (in 2D triangles)
   # Fv  - output is a list of points of each triangle
 
   # OUTPUT - for each point, list of connected points
   # ----------------------------------------------------------------------
   x <- as.matrix(x)
   x0 <- x
   sd <- sqrt(diag(var(x)))
   m  <- apply(x,2,mean)
   x.scaled <- t((t(x)-m)/sd)
   ncx <- dim(x.scaled)[2]
   nrx <- dim(x.scaled)[1]
   if (Q=="QT") xt <- delaunayn(x.scaled, options="Qt")
   if (Q=="QJ") xt <- delaunayn(x.scaled, options="QJ")
   on.exit(unlink("qhull_out.txt"))
   x.nb <- vector(mode="list", length=nrx)
   for (i in 1:nrx) {
		x.nb[[i]] <- as.integer(sort(unique(as.vector(xt[(xt==i)%*%rep(1,ncx+1)>0,]), FALSE))) 
	}
   rownames <- as.character(1:length(x.nb))
   list(graph = x.nb, type="delaunay") 
	}
  
#GM 15/03: forse serve ordinare la lista?
graph.uni <- function(x) {
	n <- NROW(x)
	x.nb <- cbind(1:n,1:n,1:n)
	x.nb[order(x)[1:(n-1)],2]<-order(x)[2:n]
	x.nb[order(x)[2:n],3]<-order(x)[(1:(n-1))]
	x.nb <-apply(x.nb,1,unique)
	list(graph = x.nb, type="unidimensional") 
	}

graph.pairs <- function(x, xd, lambda = 0.10, grid.pairs = 10) {
    graph.build <-function(val, ind, pos) sort(c(pos, val[which(ind==pos)], ind[which(val==pos)]))
	comp.areas<-matrix.connected.C(x, xd, grid.pairs = grid.pairs)
	ind.novalley <- which(comp.areas$area <= lambda) 
	x.nb <-lapply(c(1:nrow(x)), graph.build, val=comp.areas$pairs.ord[2,][ind.novalley],ind=comp.areas$pairs.ord[1,][ind.novalley])
	list(graph = x.nb, type="pairs", comp.areas = comp.areas) 
	}	

area.to.graph <- function(comp.areas, lambda = 0.10)	{
	graph.build <-function(val, ind, pos) sort(c(pos, val[which(ind==pos)], ind[which(val==pos)]))
	ind.novalley <- which(comp.areas$area <= lambda)  
	x.nb <-lapply(c(1:max(comp.areas$pairs.ord)), graph.build, val=comp.areas$pairs.ord[2,][ind.novalley],ind=comp.areas$pairs.ord[1,][ind.novalley])
	list(graph = x.nb, type="pairs", comp.areas = comp.areas) 
	}

	
valley.measureC <- function(F, output.complete=FALSE){
	RES <- .C("valley_measure", area=double(1), F0valley=as.double(F), npoints=as.integer(length(F)), PACKAGE="pdfCluster")
	if(output.complete) out <- list(area = RES$area, F0valley = RES$F0valley, F = F) else out=RES$area 
	out
}

apply_valley_measureC <- function(F){
	N <- nrow(F)
	n <- ncol(F)
	.C("apply_valley_measure",	aree=double(N), as.double(t(F)), as.integer(n), as.integer(N), PACKAGE="pdfCluster")[[1]]
}

matrix.connected.C <- function(x, pdf, grid.pairs=10, output.complete=FALSE ){
	n<-nrow(x)
	d<-ncol(x)
	pairs.ord <- combn(1:n,2)
	alpha<-seq(0,1,length=grid.pairs)
	#set a matrix as follows: nrow = n*grid.pairs ncol=d
	#each k-ple of grid.pairs rows is a set of points along the segment joining a pair
	#pairs indicated the order--> 12 13 14 15...23 24 25...34 35...
	new <- as.vector(x[t(pairs.ord)[,1],])%*%t(alpha)+as.vector(x[t(pairs.ord)[,2],])%*%t(1-alpha)
	new <- matrix(as.vector(t(new)),ncol=d)
	F <- kepdf(pdf@x, bwtype=pdf@bwtype, kernel=pdf@kernel,eval.points=new, h=pdf@par$h, hx = pdf@par$hx)@estimate
	F <- matrix(F,ncol=grid.pairs,byrow=T)
	res <- apply_valley_measureC(F)
	list( area = res, pairs.ord = pairs.ord)
}	

num.con<-function (x, estimate, x.nb, pn = 0.9, profile = TRUE, profile.grid = 25, correct=FALSE) {
    n <- length(estimate)
    ngrid <- profile.grid
    if (profile) pn <- seq(0, 1, length = ngrid + 2)[-c(1, ngrid + 2)]
    qn <- as.numeric(quantile(estimate, 1 - pn))
    xd.id <- (matrix(estimate, ncol = 1) %*% (1/qn)) < 1 # TRUE if f < quantile 
    nc.nc <- rep(0, length(qn))
    names(nc.nc) <- format(pn, digit = 3)
    nc.id <- matrix(0, nrow = n, ncol = length(qn))   #14/03 e quindi anche nc.id e nc.id1 al momento
    colnames(nc.id) <- format(pn, digit = 3)
    for (i in 1:length(qn)) {
        ni <- sum(xd.id[, i])
        if (ni < n) {
            x.nbc <- x.nb
            ind <- which(!xd.id[, i])
            for (j in ind) x.nbc[[j]] <- intersect(x.nb[[j]], 
                ind)
            x.nbc[xd.id[, i]] <- as.integer(0)
            nc <- find.nc(x.nbc)
            nc.nc[i] <- nc$nc - ni
            nc.id[xd.id[, i], i] <- -1
            nc.id[!xd.id[, i], i] <- unclass(factor(nc$comp.id[!xd.id[, i]]))
        }
        if (ni == n) {
            nc.nc[i] <- 0
            nc.id[, i] <- -1
        }
        if (correct) {
            tt <- table(nc.id[nc.id[, i] != -1, i])
            ttc <- as.integer(names(tt)[tt <= 1])
            ttp <- as.integer(names(tt)[tt > 1])
            if (length(ttc) > 0) 
                for (l in 1:length(ttc)) nc.id[nc.id[, i] == 
                  ttc[l], i] <- -1
            nc.nc[i] <- length(ttp)
            nc.nc[i] <- length(tt)
        }
        nc.nc[i] <- length(unique(nc.id[, i])) - 1
    }
	#correction 14/06/2011
	if (profile) {
		newnames <- c("0",names(nc.nc),"1")
		nc.nc <- c(0,nc.nc,1)
		names(nc.nc) <- newnames
		nc.id <- cbind(rep(-1,nrow(nc.id)),nc.id,rep(1,nrow(nc.id)))
		pn <- c(0,pn,1)
		qn <- c(max(estimate),qn,min(estimate))
		}
    list(nc = nc.nc, p = pn, id = nc.id, pq = cbind(p = pn, q = qn))
}

 find.nc <- function (nb.obj) 
{
    comp <- integer(length(nb.obj))
    #comp <- .Call("g_components", nb.obj, as.integer(comp), PACKAGE="pdfCluster")
    comp <- .Call("g_components", nb.obj, as.integer(comp), PACKAGE="pdfCluster")
    answ <- list(nc = length(unique(comp)), comp.id = comp)
    answ
}


###################################################
##### functions, classes and methods related to dbs
###################################################


setClass("dbs", representation(call="call", x="matrix", prior="numeric", dbs="numeric", clusters="numeric", noc="numeric", stage="ANY"))


dbs.cluster <- function(x, clusters, h.funct="h.norm", hmult = 1, prior, ...){
	x <- as.matrix(x)
	h.fn <- get(h.funct, inherits = TRUE)
	M <- max(clusters)
	if(M==1) stop("dbs can be computed for 2 groups at least")
	lik.fine <- matrix(NA,nrow=nrow(x),ncol=M)
	if(missing(prior)) prior <- rep(1/M,M) 
	for(m in 1:M){
		h=(h.fn(x=x[clusters==m,]))
		if (sum(h==0)>0) h[h==0]<-1e-10 
		lik.fine[,m] <- kepdf(x=x[clusters==m,],h=h*hmult, eval.points=x, ...)@estimate
		}
	ind <- t(t(lik.fine)*prior)
	den <- apply(ind,1,sum)
	tau_m <- ind/den
	ordered <- t(apply(tau_m,1,order,decreasing = T))
	ind.ug <- which(clusters==ordered[,1])
	dbs.ind <- diag(log(tau_m[,clusters]/tau_m[,ordered[,1]]))
	dbs.ind[ind.ug] <- diag(log(tau_m[ind.ug,ordered[ind.ug,1]]/ind[ind.ug,ordered[ind.ug,2]]))
	dbs.ind[is.infinite(dbs.ind)]<-max(dbs.ind[is.finite(dbs.ind)])+1
	new("dbs", call=match.call(), x=x, prior=prior, dbs=dbs.ind/max(dbs.ind), clusters=clusters, noc=M, stage=NULL)
	}

dbs.pdfCluster <- function(x, h.funct="h.norm", hmult=1, prior=as.vector(table(x@cluster.cores)/sum(table(x@cluster.cores))), stage=NULL, ...){
	if(is.null(stage)) stage <- length(x@stages)
	if (stage == 0) {
		ind.all <- which(!is.na(x@cluster.cores))
		clusters <- x@cluster.cores[ind.all] 
		} else {
		ind.all <- which(!is.na(x@stages[[stage]]))
		clusters <- x@stages[[stage]][ind.all] 
		}
	y <- as.matrix((as.matrix(x@x))[ind.all,])
	h.fn <- get(h.funct, inherits = TRUE)
	M <- max(clusters)
	if(M==1) stop("dbs can be computed for 2 groups at least")
	lik.fine <- matrix(NA,nrow=nrow(y),ncol=M)
	if(missing(prior)) prior <- rep(1/M,M) 
	for(m in 1:M){
		h=(h.fn(x=y[clusters==m,]))
		if (sum(h==0)>0) h[h==0] <- 1e-10 
		lik.fine[, m] <- kepdf(x=y[clusters==m, ],h=h*hmult, eval.points=y, ...)@estimate
		}
	ind <- t(t(lik.fine)*prior)
	den <- apply(ind,1,sum)
	tau_m <- ind/den
	ordered <- t(apply(tau_m,1,order,decreasing = T))
	ind.ug <- which(clusters==ordered[,1])
	dbs.ind <- diag(log(tau_m[,clusters]/tau_m[,ordered[,1]]))
	dbs.ind[ind.ug] <- diag(log(tau_m[ind.ug,ordered[ind.ug,1]]/ind[ind.ug,ordered[ind.ug,2]]))
	dbs.ind[is.infinite(dbs.ind)]<-max(dbs.ind[is.finite(dbs.ind)])+1
	dbs.all <- rep(NA, length(x@cluster.cores))
	dbs.all[ind.all] <- dbs.ind/max(dbs.ind)
	cl <- rep(NA, length(x@cluster.cores))
	cl[ind.all] <- clusters
	new("dbs", call=match.call(), x=x@x, prior=prior, dbs=dbs.all, clusters=cl, noc=M, stage=stage)
	}	

setGeneric("dbs", function(x, clusters, h.funct="h.norm", hmult=1, prior, ...) standardGeneric("dbs"))
setMethod("dbs", signature(x="matrix", clusters="numeric"), dbs.cluster)
setMethod("dbs", signature(x="pdfCluster",clusters="missing"), dbs.pdfCluster)
	
plot.dbs <- function(x, y, xlab="", ylab="", col= NULL, lwd =3, cex=.9, cex.axis=.5, main=NULL, labels=FALSE,...){
	dbs=x@dbs
	noc=x@noc
	if (is.null(col)) cl <- 1:noc else cl=rep(col,noc)
	if (is.null(main)) main="dbs plot" 
	plot(0,axes=F, xlim=c(min(0, min(dbs, na.rm=T)),1), ylim=c(0,sum(!is.na(dbs))),type="n",ylab=ylab, xlab=xlab, main=main, ...)
	ini <- 0
	#lb <- 0
	pos <- numeric(noc)
	lab <- numeric(0)
	ind <- numeric(0)
	for(m in 1:noc){
		lab.gr <- which(x@clusters==m)
		ind.gr <- dbs[lab.gr]
		lab <- c(lab,lab.gr[order(ind.gr)])
		ind.gr <- c(sort(ind.gr))
		ind <- c(ind,ind.gr)
		median.gr <- median(ind.gr)
		fin <- length(lab)
		ini <- (fin-(length(lab.gr)))+1
		segments(rep(0,(length(lab.gr))), ini:fin,ind.gr,ini:fin,lwd=lwd,col=cl[m])
		text(0.65,(ini+fin)/2, paste("cluster median dbs = ", round(median.gr,2)), cex=cex,pos=4,xpd=T)
		#segments(0,ini,0,fin,lwd=4,col=cl[m])
		}
	axis(1, seq(round(min(0, min(dbs, na.rm=T)), 1), 1, by=.1), seq(round(min(0, min(dbs, na.rm=T)), 1), 1, by=.1), cex.axis=cex.axis)
	if (labels) axis(2, at=(1:sum(!is.na(dbs))), lab, cex.axis=cex.axis)
	}

setMethod("plot", signature(x="dbs",y="missing"),  plot.dbs)

summary.dbs <- function(object, ...){
	cat("An S4 object of class", class(object), "\n", "\n")
	cat("Density based silhouette summary of cluster", "\n")
	dbs.groups <- list()
	for(m in 1:max(object@clusters,na.rm=T)){
		ind <- which(object@clusters==m)
		dbs.groups[[m]] <- object@dbs[ind]
		names(dbs.groups[[m]]) <- ind
		cat(m, ":", "\n")
		print(summary(dbs.groups[[m]]))
		cat("\n")
		}
	cat("Density based silhouette summary of data", "\n")
	print(summary(object@dbs))
	invisible(object)
	}

setMethod("summary", "dbs", summary.dbs)



	
setMethod("show", "dbs", function(object){
	cat("Density-based silhouette ", "\n", "\n")
	cat("dbs index: ", "\n")
	print(round(object@dbs,digits=4))
	cat("\n")
	cat("clusters: ", "\n")
	print(object@clusters)
	})






###############################
#other functions


adj.rand.index <- function (cl1, cl2){
	tab <- table(cl1, cl2)
	f2 <- function(n) n*(n-1)/2
	sum.f2 <- function(v) sum(f2(v))
	marg.1 <- apply(tab,1,sum)
	marg.2 <- apply(tab,2,sum)  
	n	<- sum(tab)
	prod <- sum.f2(marg.1)*sum.f2(marg.2)/f2(n)
	num <- (sum.f2(as.vector(tab))- prod)
	den <- 0.5*(sum.f2(marg.1)+sum.f2(marg.2))-prod
	num/den
}


check.discrete <-function (vec) length(unique(vec))





