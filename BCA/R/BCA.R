## Created 26Feb11 by Dan Putler
## Last modified on 02Sep14 by Dan Putler

## Data related items
create.samples <- function(x, est=0.34, val=0.33, rand.seed=1) {
    if((est + val) > 1) stop("The estimation and validation samples exceed 100%.")
    if(est < 0 | val < 0) stop("A negative sample size was provided.")
    nEst <- round(est*nrow(x))
    nVal <- round(val*nrow(x))
    if((nEst + nVal) < nrow(x)) {
        assignmnts <- c(rep("Estimation", nEst), rep("Validation", nVal),
          rep("Holdout", (nrow(x) - nEst - nVal)))
        }
    else {
        assignmnts <- c(rep("Estimation", nEst), rep("Validation", nVal))
        }
    set.seed(rand.seed)
    randVar <- runif(nrow(x))
    assignmnts <- assignmnts[order(randVar)]
    return(assignmnts)
    }

variable.summary <- function(dframe) {
    if(!is.data.frame(dframe)) stop("The object is not a data frame.")
    Class <- unlist(lapply(1:ncol(dframe), function(i) class(dframe[,i])))
    NAs <- 100*(apply(dframe, 2, function(x) sum(as.numeric(is.na(x))))/nrow(dframe))
    varNum <- 1:length(names(dframe))
    Levels <- rep(NA, length(varNum))
    Min.Level <- rep(NA, length(varNum))
    Mean <- rep(NA, length(varNum))
    SD <- rep(NA, length(varNum))
    if(any(Class == "factor")) {
        varNumFac <- varNum[Class=="factor"]
        totLevel <- unlist(lapply(varNumFac, 
          function(i) length(summary(dframe[,i], maxsum=nrow(dframe)))))
        minLevel <- unlist(lapply(varNumFac, 
          function(i) min(summary(dframe[,i], maxsum=nrow(dframe)))))
        Levels <- rep(NA, length(varNum))
        Levels[varNumFac] <- totLevel
        Min.Level[varNumFac] <- minLevel

	}
    if(any(Class == "numeric") | any(Class == "integer")) {
        varNumNum <- varNum[Class=="numeric" | Class=="integer"]
        Means <- unlist(lapply(varNumNum,
          function(i) mean(dframe[,i], na.rm=TRUE)))
        Mean <- rep(NA, length(varNum))
        Mean[varNumNum] <- Means
        SDs <- unlist(lapply(varNumNum,
          function(i) sd(dframe[,i], na.rm=TRUE)))
        SD[varNumNum] <- SDs
	}
    outInfo <- data.frame(Class=Class, NAs=NAs, Levels=Levels,
      Min.Level=Min.Level, Mean=Mean, SD=SD)
    outInfo <- outInfo[order(outInfo$NAs),]
    names(outInfo) <- c("Class", "%.NA", "Levels", "Min.Level.Size", "Mean",
      "SD")
    return(outInfo)
    }

relabel.factor <- function(x, new.labels, old.labels=levels(x)) {
    if(!is.factor(x)) stop("The vector given is not a factor")
    if(length(new.labels) != length(old.labels)) { stop(
    "The number of new level labels must equal the original number of labels.")
    }
    xx <- rep(NA, length(x))
    for(i in 1:length(old.labels)) {
        xx[x == old.labels[i]] <- new.labels[i]
        }
    return(as.factor(xx))
    }


## Clustering
# 2D biplot
bpCent <- function(pc, clsAsgn, data.pts = TRUE, centroids = TRUE, 
  choices = 1:2, scale = 1, pc.biplot=FALSE, var.axes = TRUE, col=palette()[1:2],
  cex = rep(par("cex"), 2), xlabs = NULL, ylabs = NULL, expand=1, xlim = NULL,
  ylim = NULL, arrow.len = 0.1, main = NULL, sub = NULL, xlab = NULL,
  ylab = NULL, ...)
{
    if(length(choices) != 2) stop("length of choices must be 2")
    if(!length(scores <- pc$x))
	stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
             domain = NA)
    if(is.complex(scores))
        stop("biplots are not defined for complex PCA")
    lam <- pc$sdev[choices]
    n <- NROW(scores)
    lam <- lam * sqrt(n)
    if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
    if(scale != 0) lam <- lam^scale else lam <- 1
    if(pc.biplot) lam <- lam / sqrt(n)
    cntrs <- apply(data.frame(scores[, choices]), MARGIN=2,
      function(x) tapply(x,as.factor(clsAsgn),mean))
    x <- t(t(scores[, choices]) / lam)
    cntrs <- t(t(cntrs) / lam)
    y <- t(t(pc$rotation[, choices]) * lam)
    n <- nrow(pc)
    ncls <- nrow(cntrs)
    p <- nrow(y)
    if(missing(xlabs)) {
	xlabs <- dimnames(x)[[1]]
	if(is.null(xlabs)) xlabs <- 1:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2]])
    if(missing(ylabs)) {
	ylabs <- dimnames(y)[[1]]
	if(is.null(ylabs)) ylabs <- paste("Var", 1:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2]])

    if(length(cex) == 1) cex <- c(cex, cex)
    if(missing(col)) {
	col <- par("col")
	if (!is.numeric(col)) col <- match(col, palette(), nomatch=1)
	col <- c(col, col + 1)
    }
    else if(length(col) == 1) col <- c(col, col)
    else col <- col

    unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
    rangx1 <- unsigned.range(x[, 1])
    rangx2 <- unsigned.range(x[, 2])
    rangy1 <- unsigned.range(y[, 1])
    rangy2 <- unsigned.range(y[, 2])

    if(missing(xlim) && missing(ylim))
	xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if(missing(xlim)) xlim <- rangx1
    else if(missing(ylim)) ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")
    if(!is.null(main))
        op <- c(op, par(mar = par("mar")+c(0,0,1,0)))
    plot(x, type = "n", xlim = xlim, ylim = ylim, col = col[1],
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    if(data.pts) text(x, xlabs, cex = cex[1], col = col[1], ...)
    if(centroids) text(cntrs, as.character(1:ncls), col="blue", cex=1.5)
    
    par(new = TRUE)
    plot(y, axes = FALSE, type = "n", xlim = xlim*ratio, ylim = ylim*ratio,
	 xlab = "", ylab = "", col = col[1], ...)
    axis(3, col = col[2])
    axis(4, col = col[2])
    box(col = col[1])
    text(y, labels=ylabs, cex = cex[2], col = "red", ...)
    if(var.axes)
	arrows(0, 0, y[,1] * 0.8, y[,2] * 0.8, col = "red", length=arrow.len)
    invisible()
}

# 3D biplot
bpCent3d <- function(pc, clsAsgn, data.pts = TRUE, centroids = TRUE, 
  choices = 1:3, scale = 1, pc.biplot=FALSE, var.axes = TRUE, 
  col.score = "black", col.cntrs = "blue", col.load = "red",
  xlabs = NULL, ylabs = NULL, xlim = NULL, ylim = NULL, zlim = NULL,
  xlab = NULL,  ylab = NULL, dim.lab = NULL, fov = 60)
{
	if (requireNamespace("rgl", quietly = TRUE)) {
		if(length(choices) != 3) stop("length of choices must be 2")
		if(!length(scores <- pc$x))
			stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
             domain = NA)
		if(is.complex(scores))
			stop("biplots are not defined for complex PCA")
		lam <- pc$sdev[choices]
		n <- NROW(scores)
		lam <- lam * sqrt(n)
		if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
		if(scale != 0) lam <- lam^scale else lam <- 1
		if(pc.biplot) lam <- lam / sqrt(n)
		dim.lab <- paste("PC", choices, sep="")
		use.rgl <- options("Rcmdr")[[1]]$use.rgl
		if(length(use.rgl) == 0 || use.rgl) require(rgl)
			cntrs <- apply(data.frame(scores[, choices]), MARGIN=2,
		function(x) tapply(x,as.factor(clsAsgn),mean))
		x <- t(t(scores[, choices]) / lam)
		cntrs <- t(t(cntrs) / lam)
		y <- t(t(pc$rotation[, choices]) * lam)
		n <- nrow(pc)
		ncls <- nrow(cntrs)
		p <- nrow(y)
		if(missing(xlabs)) {
			xlabs <- dimnames(x)[[1]]
			if(is.null(xlabs)) xlabs <- 1:n
		}
		xlabs <- as.character(xlabs)
		dimnames(x) <- list(xlabs, dimnames(x)[[2]])
		if(missing(ylabs)) {
			ylabs <- dimnames(y)[[1]]
			if(is.null(ylabs)) ylabs <- paste("Var", 1:p)
		}
		ylabs <- as.character(ylabs)
		dimnames(y) <- list(ylabs, dimnames(y)[[2]])

		unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
		rangx1 <- unsigned.range(x[, 1])
		rangx2 <- unsigned.range(x[, 2])
		rangx3 <- unsigned.range(x[, 3])
		rangy1 <- unsigned.range(y[, 1])
		rangy2 <- unsigned.range(y[, 2])
		rangy3 <- unsigned.range(y[, 3])
		ratio <- max(rangy1/rangx1, rangy2/rangx2, rangy3/rangx2)
		rangy1 <- rangy1/ratio
		rangy2 <- rangy2/ratio
		rangy3 <- rangy3/ratio
		rangx <- c(min(rangx1[1], rangy1[1]), max(rangx1[2], rangy1[2]))
		rangy <- c(min(rangx2[1], rangy2[1]), max(rangx2[2], rangy2[2]))
		rangz <- c(min(rangx3[1], rangy3[1]), max(rangx3[2], rangy3[2]))
		if (missing(xlim) && missing(ylim) && missing(zlim))
			xlim <- ylim <- zlim <- range(rangx, rangy, rangz)
		else if (missing(xlim) && missing(ylim))
			xlim <- ylim <- range(rangx, rangy)
		else if (missing(xlim) && missing(zlim))
			xlim <- zlim <- range(rangx, rangz)
		else if (missing(ylim) && missing(zlim))
			ylim <- zlim <- range(rangy, rangz)
		else if (missing(xlim))
			xlim <- rangx
		else if (missing(ylim))
			ylim <- rangy
		else if (missing(zlim))
			zlim <- rangz
		cMin <- min(xlim, ylim, zlim)
		rgl::rgl.clear()
		rgl::rgl.viewpoint(fov=fov)
		rgl::rgl.bg(color="white")
		rgl::rgl.lines(c(cMin,rangx[2]), c(cMin,cMin), c(cMin,cMin), color="black")
		rgl::rgl.lines(c(cMin,cMin), c(cMin,rangy[2]), c(cMin,cMin), color="black")
		rgl::rgl.lines(c(cMin,cMin), c(cMin,cMin), c(cMin,rangz[2]), color="black")
		if(!is.null(dim.lab)) {
			rgl::rgl.texts(rangx[2], cMin, cMin, dim.lab[1], adj = 1, color = "black")
			rgl::rgl.texts(cMin, rangy[2], cMin, dim.lab[2], adj = 1, color = "black")
			rgl::rgl.texts(cMin, cMin, rangz[2], dim.lab[3], adj = 1, color = "black")
		}
		if(data.pts) rgl::texts3d(x[,1], x[,2], x[,3], texts = xlabs, color= col.score)
		if(centroids) rgl::texts3d(cntrs[,1], cntrs[,2], cntrs[,3], 
		  texts = as.character(1:ncls), color = col.cntrs)
		if (var.axes) {
			zrs <- rep(0, nrow(y))
			for(i in 1:nrow(y)) {
				rgl::lines3d(c(0,y[i,1]/ratio), c(0,y[i,2]/ratio), c(0,y[i,3]/ratio),
				color = "red")
	        }
	    }
		rgl::texts3d(1.05*y[,1]/ratio, 1.05*y[,2]/ratio, 1.05*y[,3]/ratio, texts = ylabs,
		color = "red")
		return(invisible())
	} else {
		stop("The needed rgl package is not available.")
	}
}

# The SD Index plot functions for k-means
SD.clv <- function(x, clus, alpha) {
    if(!is.data.frame(x)) x <- as.data.frame(x)
    scatt <- clv.Scatt(x, clus)
    dis <- clv.Dis(scatt$cluster.center)
    SD <- clv.SD(scatt$Scatt, dis, alfa=alpha)
    return(SD)
    }

SDIndex <- function(x, minClust, maxClust, iter.max=10, num.seeds=10) {
    if(minClust < 2) stop("The minimum number of clusters must be two or more")
    cluster.vec <- as.integer(seq(minClust, maxClust, by=1))
    maxClus <- KMeans(x, centers=maxClust, iter.max=iter.max,
      num.seeds=num.seeds)$cluster
    maxScatt <- clv.Scatt(x, maxClus)
    maxDis <- clv.Dis(maxScatt$cluster.center)
    SD_Index <- rep(NA,length(cluster.vec))
    for(i in 1:length(cluster.vec)) {
        kmeansClus <- KMeans(x, centers=cluster.vec[i], iter.max=iter.max,
          num.seeds=num.seeds)$cluster
        SD_Index[i] <- SD.clv(x, kmeansClus, maxDis)
        }
    plot(cluster.vec, SD_Index, type="l", col=palette()[1], lwd=2,
      main= "K-Means SD Index for Different Numbers of Clusters" ,
      xlab="Number of Clusters",
      ylab="SD Index Value")
    invisible()
    }

# The cluster diagnostic functions for K-Centroids clustering

# This function calculates the Calinski-Harabas index for all solutions
# contained in a bootFlexclust object. It does it for each Rand test paired
# comparison, so 100 bootstrap replicates of the Rand index will result in 200
# Calinski-Harabase index values. This code borrows from the function index.G1
# of Marek Walesiak and Andrzej Dudek's clusterSim package in terms of
# implementing the sum of squares components
bootCH <- function(xdat, k_vals, clstr1, clstr2, cntrs1, cntrs2, 
  method = c("kmn", "kmd", "neuralgas")) {
  method = match.arg(method)
  if(method == "kmd") all_centers <- apply(xdat, 2, median)
  else all_centers <- apply(xdat, 2, mean)
  all_dif <- sweep(xdat, 2, all_centers,"-")
  tss <- sum(all_dif^2)
  n_solu <- dim(clstr1)[2] # total number of cluster solutions considered
  nboot <- dim(clstr1)[3]
  n_obs <- dim(clstr1)[1]
  ch_mat <- matrix(NA, nrow=2*nboot, ncol=n_solu)
  for(k_ind in 1:n_solu) {
    cent_array1 <- cntrs1[[k_ind]]
    cent_array2 <- cntrs2[[k_ind]]
    cls_asgn_m1 <- clstr1[,k_ind,]
    cls_asgn_m2 <- clstr2[,k_ind,]
    k <- k_vals[k_ind]
    ch_reps1 <- rep(NA, nboot)
    ch_reps2 <- rep(NA, nboot)
    for(b_ind in 1:nboot) {
      clus1 <- cls_asgn_m1[,b_ind]
      clus2 <- cls_asgn_m2[,b_ind]
      centrds1 <- cent_array1[,,b_ind]
      centrds2 <- cent_array2[,,b_ind]
      wss1a <- (xdat - centrds1[clus1,])^2
      wss2a <- (xdat - centrds2[clus2,])^2
      wss1 <- 0
      wss2 <- 0
      for(m in 1:k) {
        wss1 <- wss1 + sum(wss1a[clus1 == m,])
        wss2 <- wss2 + sum(wss2a[clus2 == m,])
      }
      bss1 <- tss - wss1
      bss2 <- tss - wss2
      ch_reps1[b_ind] <- ((n_obs - k)/(k - 1))*(bss1/wss1)
      ch_reps2[b_ind] <- ((n_obs - k)/(k - 1))*(bss2/wss2)
    }
    ch_mat[,k_ind] <- c(ch_reps1, ch_reps2)
    dimnames(ch_mat)[[2]] <- as.character(k_vals)
  }
  return(ch_mat)
}

# A function to produce the diagnostic plots for a bootFlexclust object
bootPlot <- function(fc, ch, col1="blue", col2="green") {
  if(class(fc) != "bootFlexclust") {
    stop("An object of class bootFlexclust was not provided")
  }
  if(class(ch) != "ch_index") {
    stop("An object of class ch_index was not provided")
  }
  omfrow <- par()$mfrow
  par(mfrow = c(2, 1))
  boxplot(fc@rand, main=paste("Adj Rand Index for", ch$data, "using",
    ch$method), xlab="Number of Clusters", ylab="Adjusted Rand", col=col1)
  boxplot(ch, main=paste("C-H Index for", ch$data, "using", ch$method),
    xlab="Number of Clusters", ylab="Calinski-Harabas", col=col2)
  par(mfrow = omfrow)
  return(invisible())
}

# A function needed as a workaround for Rcmdr namespace issues, really a
# wrapper for bootFlexclust, bootCH, and bootPlot that enables a single call
# to be made in RcmdrPlugin.BCA
bootCVD <- function(x, k, nboot=100, nrep=1, method = c("kmn", "kmd", "neuralgas"),
  col1, col2, dsname) {
  print(class(x))
  method = match.arg(method)
  if(method == "kmn") {
    bfc <- bootFlexclust(x=x, k=k, nboot=nboot, nrep=nrep, FUN = cclust,
      dist = "euclidean", method = "kmeans")
    the_method <- "K-Means"
  }
  else if(method == "kmd") {
    bfc <- bootFlexclust(x=x, k=k, nboot=nboot, nrep=nrep, FUN = kcca,
      family = kccaFamily("kmedians"))
    the_method <- "K-Medians"
  }
  else {
    bfc <- bootFlexclust(x=x, k=k, nboot=nboot, nrep=nrep, FUN = cclust,
      dist = "euclidean", method = "neuralgas")
    the_method <- "Neural Gas"
  } 
  cat("\nSummary of Rand Indices:\n")
  print(summary(bfc@rand))
  ch <- bootCH(x, k, bfc@cluster1, bfc@cluster2, bfc@centers1, bfc@centers2, method)
  cat("\nSummary of Calinski-Harabas Indices:\n")
  print(summary(ch))
  omfrow <- par()$mfrow
  par(mfrow = c(2, 1))
  boxplot(bfc@rand, main=paste("Adj Rand Index for", dsname, "using", the_method),
    xlab="Number of Clusters", ylab="Adjusted Rand", col=col1)
  boxplot(ch, main=paste("C-H Index for", dsname, "using", the_method),
    xlab="Number of Clusters", ylab="Calinski-Harabas", col=col2)
  par(mfrow = omfrow)
}

## Lift charts and scoring
lift.chart <- function(modelList, data, targLevel, trueResp, type="cumulative",
  sub="") {
    if(type != "cumulative" & type != "incremental") {
        stop("An improper lift chart type is specified.")
        }
    set.seed(1)
    data <- data[order(runif(nrow(data))), ]
    yvar1 <- rep(NA, length(modelList))
    modAvail <- rep(NA, length(modelList))
    probVar <- NULL
    for(i in 1:length(modelList)) {
        mod <- eval(parse(text=modelList[i]))
        modtype <- class(mod)[1]
        if(modtype != "glm" & modtype != "rpart" & modtype != "nnet.formula") {
      stop("Models can only be estimated using glm, rpart, or nnet.")
        }
        yvar1[i] <- as.character(mod$call$formula)[2]
        xvars <- unlist(strsplit(as.character(mod$call$formula)[3]," + ",
          fixed=TRUE))
        if(!all(xvars %in% names(data))) {
            probVar <- c(probVar, xvars[!(xvars %in% names(data))])
            modAvail[i] <- FALSE
            }
        else {
            modAvail[i] <- TRUE
            }
        }
    if(any(yvar1!=yvar1[1])) {
          stop("Not all the models have the same dependent variable")
          }
    yvar2 <- data[[yvar1[1]]]
    if(!is.factor(yvar2)) {
        stop("The y variable must be a two-level factor.")
        }
    if(length(levels(yvar2)) != 2) {
        stop("The y variable must be a two-level factor.")
        }
    if(!any(as.character(yvar2) == targLevel)) {
        stop(paste('None of the levels of the response variable is "', targLevel,
          '".', sep=""))
        }
    yvar <- as.numeric(yvar2 == targLevel)
    sampResp <- sum(yvar)/length(yvar)
    print(sampResp)
    if(length(probVar)>0) {
        probVar <- unique(probVar)
        probModel <- modelList[!modAvail]
        warnString <- paste("The models", paste(probModel, collapse=", "),
          "are not in the lift chart because the variables",
          paste(probVar, collapse=", "), "are not available.")
        #warning(warnString, call.=FALSE, immediate.=TRUE)
        Message(message=gettextRcmdr(warnString), type="warning")
        modelList <- modelList[modAvail]
        if(length(modelList)==0) {
            Message(message=gettextRcmdr(paste(
              "All models are missing at least one of the variables: ",
              paste(probVar, collapse=", "), ".", sep="")), type="error")
            return()
            }
        }
    sampWt <- (sampResp*(1 - trueResp))/(trueResp*(1 - sampResp))
    nmodels <- length(modelList) # new
    colr <- rep(palette(), ceiling(nmodels/length(palette()))) # new
    if(type == "cumulative") {
        plot(seq(0.1, 1, 0.1), seq(0.1, 1, 0.1), main=
          "Weighted Cumulative Response Captured", sub=sub, xlab="Sample Proportion",
          ylab="Percent of Total Response Captured", type="l", lwd=2, xaxs="i",
          yaxs="i")
        for(i in 1:nmodels) {
            model1 <- eval(parse(text=modelList[i]))
            modtype <- class(model1)[1]
            model <- eval(model1$call)
            if(modtype == "glm" | modtype == "nnet.formula") {
                if(levels(yvar2)[1] == targLevel) {
                    var1 <- yvar[order(predict(model, newdata = data),
                      decreasing = FALSE)]
                    }
                else {
                    var1 <- yvar[order(predict(model, newdata = data),
                      decreasing = TRUE)]
                    }
               }
            else {
                var1 <- yvar[order(as.vector(predict(
                  model, newdata = data)[ , targLevel]), decreasing = TRUE)]
                }
            var.ind1 <- rep(1, length(var1))
            var.ind1[var1 == 0] <- sampWt
            var.ind <- cut(cumsum(var.ind1)/sum(var.ind1), seq(0, 1, 0.1),
              include.lowest = TRUE)
            var2 <- as.vector(by(var1, var.ind, sum))
            lines(seq(0.1, 1, 0.1), cumsum(var2)/sum(var2), col = colr[i], lwd = 2)
            points(seq(0.1, 1, 0.1), cumsum(var2)/sum(var2), col = colr[i], pch=i)
            }
        legend("bottomright", legend=modelList, col=colr[1:nmodels],
          pch=1:length(modelList), lty=1, lwd=2)
        }
    else {
        resp.matrix <- matrix(NA, nrow=10, ncol=length(modelList))
        for(i in 1:nmodels) {
            model1 <- eval(parse(text=modelList[i]))
            modtype <- class(model1)[1]
            model <- eval(model1$call)
            if(modtype == "glm" | modtype == "nnet.formula") {
                if(levels(yvar2)[1] == targLevel) {
                    var1 <- yvar[order(predict(model, newdata = data),
                      decreasing = FALSE)]
                    }
                else {
                    var1 <- yvar[order(predict(model, newdata = data),
                      decreasing = TRUE)]
                    }
               }
            else {
                var1 <- yvar[order(as.vector(predict(
                  model, newdata = data)[ , targLevel]), decreasing = TRUE)]
                }
            var.ind1 <- rep(1, length(var1))
            var.ind1[var1 == 0] <- sampWt
            var.ind <- cut(cumsum(var.ind1)/sum(var.ind1), seq(0, 1, 0.1),
              include.lowest = TRUE)
            var2 <- as.vector(by(var1, var.ind, sum))
            var3 <- as.vector(by(var.ind1, var.ind, sum))
            resp.matrix[, i] <- var2/var3
            }
        max.resp <- max(resp.matrix)
        plot(seq(0.1, 1, 0.1), seq(0, max.resp, length=10), type="n", main=
          "Weighted Incremental Response Rate", sub=sub, xlab="Sample Percentile",
          ylab="Resposne Rate", lwd=2, xaxs="i", yaxs="i")
        lines(seq(0.1, 1, 0.1), rep(trueResp, 10), lwd=2)
        for(j in 1:nmodels) {
            lines(seq(0.1, 1, 0.1), as.vector(resp.matrix[, j]), col=colr[j], lwd=2)
           points(seq(0.1, 1, 0.1), as.vector(resp.matrix[, j]), col=colr[j], pch=j)
            }
            print(length(modelList))
       legend("topright",legend=modelList,col=colr[1:nmodels],
          pch=1:nmodels, lty=1, lwd=2)
        }
    invisible()
    }

rankScore <- function(model, data, targLevel) {
    mod <- eval(parse(text=model))
    modtype <- class(mod)[1]
    if(modtype != "glm" & modtype != "rpart" & modtype != "nnet.formula") {
stop("Scoring can only be done for models estimated using glm, rpart, or nnet.")
        }
    yvar <- as.character(mod$call$formula)[2]
    origYs <- eval(parse(text=paste("unique(", ActiveDataSet(), "$", yvar, ")")))
    origYs <- as.character(origYs)
    origYs <- origYs[order(origYs)]
    xvars <- unlist(strsplit(as.character(mod$call$formula)[3]," + ",
      fixed=TRUE))
    if(!all(xvars %in% names(data))) {
        probVar <- c(xvars[!(xvars %in% names(data))])
        stop(paste("The model variables", paste(probVar, collapse=", "),
          "are not in the data set."))
        }
    modelReDo <- eval(mod$call)
    if(modtype == "glm" | modtype == "nnet.formula") {
        if(origYs[1] == targLevel) {
            scoreVar1 <- -1 * predict(modelReDo, newdata = data)
            }
        else {
            scoreVar1 <- predict(modelReDo, newdata = data)
            }
        }
    else {
        scoreVar1 <- predict(modelReDo, newdata = data)[ , targLevel]
        }
    score.df <- data.frame(scoreVar = 1:nrow(data), 
      oldOrd=order(scoreVar1, decreasing=TRUE))
    score.df <- score.df[order(score.df$oldOrd),]
    scoreVar <- score.df$scoreVar
    return(scoreVar)
    }

rawProbScore <- function(model, data, targLevel) {
    mod <- eval(parse(text=model))
    modtype <- class(mod)[1]
    if(modtype != "glm" & modtype != "rpart" & modtype != "nnet.formula") {
stop("Scoring can only be done for models estimated using glm, rpart, or nnet.")
        }
    yvar <- as.character(mod$call$formula)[2]
    origYs <- eval(parse(text=paste("unique(", ActiveDataSet(), "$", yvar, ")")))
    origYs <- as.character(origYs)
    origYs <- origYs[order(origYs)]
    xvars <- unlist(strsplit(as.character(mod$call$formula)[3]," + ",
      fixed=TRUE))
    if(!all(xvars %in% names(data))) {
        probVar <- c(xvars[!(xvars %in% names(data))])
        stop(paste("The model variables", paste(probVar, collapse=", "),
          "are not in the data set."))
        }
    modelReDo <- eval(mod$call)
    if(modtype == "glm" | modtype == "nnet.formula") {
        if(origYs[1] == targLevel) {
            scoreVar1 <- -1 * predict(modelReDo, newdata = data)
            }
        else {
            scoreVar1 <- predict(modelReDo, newdata = data)
            }
        }
    else {
        scoreVar1 <- predict(modelReDo, newdata = data)[ , targLevel]
        }
    if(modtype == "glm") {
        scoreVar <- exp(scoreVar1)/(exp(scoreVar1) + 1)
        }
    else {
        scoreVar <- scoreVar1
        }
    return(scoreVar)
    }

adjProbScore <- function(model, data, targLevel, trueResp) {
    mod <- eval(parse(text=model))
    modtype <- class(mod)[1]
    if(modtype != "glm" & modtype != "rpart" & modtype != "nnet.formula") {
stop("Scoring can only be done for models estimated using glm, rpart, or nnet.")
        }
    yvar <- as.character(mod$call$formula)[2]
    yvar1 <- eval(parse(text=paste("as.character(", ActiveDataSet(), "$", yvar, ")")))
    yvar2 <- as.numeric(yvar1 == targLevel)
    print(yvar1)
    sampResp <- sum(yvar2)/length(yvar2)
    adjWt <- trueResp/sampResp
    origYs <- eval(parse(text=paste("unique(", ActiveDataSet(), "$", yvar, ")")))
    origYs <- as.character(origYs)
    origYs <- origYs[order(origYs)]
    xvars <- unlist(strsplit(as.character(mod$call$formula)[3]," + ",
      fixed=TRUE))
    if(!all(xvars %in% names(data))) {
        probVar <- c(xvars[!(xvars %in% names(data))])
        stop(paste("The model variables", paste(probVar, collapse=", "),
          "are not in the data set."))
        }
    modelReDo <- eval(mod$call)
    if(modtype == "glm" | modtype == "nnet.formula") {
        if(origYs[1] == targLevel) {
            scoreVar1 <- -1 * predict(modelReDo, newdata = data)
            }
        else {
            scoreVar1 <- predict(modelReDo, newdata = data)
            }
        }
    else {
        scoreVar1 <- predict(modelReDo, newdata = data)[ , targLevel]
        }
    if(modtype == "glm") {
        scoreVar <- adjWt*(exp(scoreVar1)/(exp(scoreVar1) + 1))
        }
    else {
        scoreVar <- adjWt*scoreVar1
        }
    return(scoreVar)
    }

## The altered scatter plot functions
# fancy scatterplots  (J. Fox)

# 2010-09-05: J. Fox: changed color choice
# 2010-09-16: fixed point color when col is length 1
# 2010-12-19: J. Fox: added argument legend.coords to place legend.
# 2011-01-15: J. Fox: If x is a factor, calls Boxplot()
# 2011-02-24: Dan Putler: Altered color choice for the LS line

scatterplotBCA <- function(x, ...){
	UseMethod("scatterplotBCA", x)
}

scatterplotBCA.formula <- function (x, data, subset, xlab, ylab, legend.title, legend.coords, labels, ...) {
	na.save <- options(na.action=na.omit)
	on.exit(options(na.save))
	na.pass <- function(dframe) dframe
	m <- match.call(expand.dots=FALSE)
	if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
		m$data <- as.data.frame(data)
	m$na.action <- na.pass
	m$legend.coords <- m$legend.title <- m$labels <- m$xlab <- m$ylab <- m$... <- NULL
	m[[1]] <- as.name("model.frame")
	if (!inherits(x, "formula") | length(x) != 3) 
		stop("invalid formula")    
	x <- as.character(c(x))
	x <- as.formula(sub("\\|", "+", x))
	m$formula <- x
	if (missing(data)){ 
		X <- na.omit(eval(m, parent.frame()))
		if (missing(labels)) labels <- gsub("X", "", row.names(X)) 
	}
	else{
		if (!missing(labels)) row.names(data) <- labels
		X <- eval(m, parent.frame())
		labels <- row.names(X)
	}
	names <- names(X)
	if (missing(xlab)) xlab <- names[2]
	if (missing(ylab)) ylab <- names[1]
	if (ncol(X) == 2) scatterplotBCA(X[,2], X[,1], xlab=xlab, ylab=ylab, 
			labels=labels, ...)
	else {
		if (missing(legend.title)) legend.title <- names[3]
		scatterplotBCA(X[,2], X[,1], groups=X[,3], xlab=xlab, ylab=ylab,  
			legend.title=legend.title, legend.coords=legend.coords, labels=labels, ...)
	}
}

scatterplotBCA.default <- function(x, y, smooth=TRUE, spread=!by.groups, span=.5, loess.threshold=2, reg.line=lm, 
		boxplots=if (by.groups) "" else "xy",
		xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), las=par("las"),
		lwd=2, lwd.smooth=lwd, lwd.spread=lwd, lty=1, lty.smooth=lty, lty.spread=2,
		labels, id.method = "mahal", 
		id.n = if(id.method[1]=="identify") length(x) else 0, 
		id.cex = 1, id.col = palette()[1],
		log="", jitter=list(), xlim=NULL, ylim=NULL,
		cex=par("cex"), cex.axis=par("cex.axis"), cex.lab=par("cex.lab"), 
		cex.main=par("cex.main"), cex.sub=par("cex.sub"), 
		groups, by.groups=!missing(groups), legend.title=deparse(substitute(groups)), legend.coords,
		ellipse=FALSE, levels=c(.5, .95), robust=TRUE,
		col=if (n.groups == 1) palette()[c(2,1,3)] else rep(palette(), length=n.groups),
		pch=1:n.groups, 
		legend.plot=!missing(groups), reset.par=TRUE, grid=TRUE, ...){
		logged <- function(axis=c("x", "y")){
		axis <- match.arg(axis)
		0 != length(grep(axis, log))
	}
	err <- ""
	lowess.line <- function(x, y, col, span) {
		if (logged("x")) x <- log(x)
		if (logged("y")) y <- log(y)
		valid <- complete.cases(x, y)
		x <- x[valid]
		y <- y[valid]
		ord <- order(x)
		x <- x[ord]
		y <- y[ord]
#		if (length(unique(x)) < lowess.threshold || length(unique(y)) < lowess.threshold) return()
		if (length(unique(y)) < loess.threshold) return()
		warn <- options(warn=-1)
		if (!spread){
			fit <- try(loess.smooth(x, y, span=span), silent=TRUE)
			if (class(fit) == "try-error"){
				err <<- c(err, "smooth")
				options(warn)
				return()
			}
			x <-if (logged("x")) exp(fit$x) else fit$x  
			y <-if (logged("y")) exp(fit$y) else fit$y
			lines(x, y, lwd=lwd.smooth, col=col, lty=lty.smooth)
		}
		else{
			fit <- try(loess(y ~ x, degree=1, family="symmetric", span=span), silent=TRUE)
			if (class(fit) == "try-error"){
				err <<- c(err, "smooth")
				options(warn)
				return()
			}
			res <- residuals(fit)
			pos <- res > 0
			pos.fit <- try(loess(res^2 ~ x, span=span, degree=0, family="gaussian", subset=pos), silent=TRUE)
			neg.fit <- try(loess(res^2 ~ x, span=span, degree=0, family="gaussian", subset=!pos), silent=TRUE)
			if (class(pos.fit) == "try-error" || class(neg.fit) == "try.error"){
				err <<- c(err, "spread")
				options(warn)
				return()
			}
			if (logged("x")) x <- exp(x)
			y <- if (logged("y")) exp(fitted(fit)) else fitted(fit) 
			lines(x, y, lwd=lwd.smooth, col=col, lty=lty.smooth)
			y.pos <- if (logged("y")) exp(fitted(fit)[pos] + sqrt(fitted(pos.fit)))  
					else fitted(fit)[pos] + sqrt(fitted(pos.fit))
			lines(x[pos], y.pos, lwd=lwd.spread, lty=lty.spread, col=col)
			y.neg <- if (logged("y")) exp(fitted(fit)[!pos] - sqrt(fitted(neg.fit)))
					else fitted(fit)[!pos] - sqrt(fitted(neg.fit))
			lines(x[!pos], y.neg, lwd=lwd.spread, lty=lty.spread, col=col)
		}
		options(warn)
	}
	reg <- function(x, y, col){
		if (logged("x")) x <- log(x)
		if (logged("y")) y <- log(y)
		mod <- reg.line(y ~ x)
		y.hat <- fitted.values(mod)
		x <- model.matrix(mod)[, 2]
		min <- which.min(x)
		max <- which.max(x)
		if (!logged("x")){
			x1 <- x[min]
			x2 <- x[max]
		}
		else {
			x1 <- exp(x[min])
			x2 <- exp(x[max])
		}
		if (!logged("y")){
			y1 <- y.hat[min]
			y2 <- y.hat[max]
		}
		else {
			y1 <- exp(y.hat[min])
			y2 <- exp(y.hat[max])
		}
		lines(c(x1, x2), c(y1, y2), lwd=lwd, col=col, lty=lty)
	}
	hbox <- function(x){
		if (logged("x")){
			log.x <- "x"
			.x <- log(x)		
		}
		else {
			log.x <- ""
			.x <- x
		}
		plot(x, seq(0, 1, length=length(x)), type="n", axes=FALSE, xlab="", ylab="", log=log.x, xlim=xlim)
		res <- boxplot.stats(.x, coef = 1.5, do.conf=FALSE)
		if (logged("x")){
			res$stats <- exp(res$stats)
			if (!is.null(res$out)) res$out <- exp(res$out)
		}
		LW <- res$stats[1]
		Q1 <- res$stats[2]
		M <- res$stats[3]
		Q3 <- res$stats[4]
		UW <- res$stats[5]
		lines(c(Q1, Q1, Q3, Q3, Q1), c(0, 1, 1, 0, 0))
		lines(c(M, M), c(0, 1))
		lines(c(LW, Q1), c(.5, .5))
		lines(c(Q3, UW), c(.5, .5))
		if (!is.null(res$out)) points(res$out, rep(.5, length(res$out)), cex=cex)
	}
	vbox <- function(y){
		if (logged("y")){
			log.y <- "y"
			.y <- log(y)
		}
		else {
			log.y <- ""
			.y <- y
		}
		plot(seq(0, 1, length=length(y)), y, type="n", axes=FALSE, xlab="", ylab="", log=log.y, ylim=ylim)
		res <- boxplot.stats(.y, coef = 1.5, do.conf=FALSE)
		if (logged("y")){
			res$stats <- exp(res$stats)
			if (!is.null(res$out)) res$out <- exp(res$out)
		}
		LW <- res$stats[1]
		Q1 <- res$stats[2]
		M <- res$stats[3]
		Q3 <- res$stats[4]
		UW <- res$stats[5]
		lines(c(0, 1, 1, 0, 0), c(Q1, Q1, Q3, Q3, Q1))
		lines(c(0, 1), c(M, M))
		lines(c(.5, .5), c(LW, Q1))
		lines(c(.5, .5), c(Q3, UW))
		if (!is.null(res$out)) points(rep(.5, length(res$out)), res$out, cex=cex)
	}
	# force evaluation of some arguments
	by.groups
	legend.plot
	legend.title
	spread 
	if (missing(labels)){
		labels <- if (is.null(names(y)))
					seq(along=y)
				else names(y)
	}
	if (length(labels) != length(y)) stop("labels argument is the wrong length")
	if (is.factor(x)) {
		if (!(id.method %in% c("y", "identify", "none"))) id.method <- "y"
		return(Boxplot(y, x, id.method="y", labels=labels, xlab=xlab, ylab=ylab))
	}
	mar <- par("mar")
	mfcol <- par("mfcol")
	if (reset.par) on.exit(par(mar=mar, mfcol=mfcol))
	if( FALSE == boxplots) boxplots <- ""
	if (!missing(groups)){
		data <- na.omit(data.frame(groups, x, y, labels, stringsAsFactors=FALSE))
		groups <- data[,1]
		if (!is.factor(groups)) groups <- as.factor(groups)
		.x <- data[,2]
		.y <- data[,3]
		labels <- data[,4]
		top <- if (legend.plot && missing(legend.coords)) 
					# 4 + length(levels(as.factor(groups))) else mar[3]
					4 + nlevels(groups) else mar[3]
	}
	else {
		.x <- x
		.y <- y
		top <- mar[3]
		groups <- factor(rep(1, length(.x)))
	}
	xbox <- length(grep("x", boxplots)) > 0
	ybox <- length(grep("y", boxplots)) > 0
	# groups <- as.factor(if(missing(groups)) rep(1, length(.x)) else as.character(groups))
	if (xbox && ybox)
		layout(matrix(c(1, 0, 3, 2), 2, 2),
				widths = c(5, 95),
				heights= c(95, 5))
	else if (ybox)
		layout(matrix(c(1, 2),1, 2),
				widths = c(5, 95),
				heights= 100)
	else if (xbox)
		layout(matrix(c(2, 1), 2, 1),
				widths = 100,
				heights= c(95, 5))
	else layout (matrix(1, 1, 1),
				widths=100, heights=100)
	par(mar=c(mar[1], 0, top, 0))
	if (ybox > 0) vbox(.y) 
#	else plot(0, 0, xlab="", ylab="", axes=FALSE, type="n", xlim=xlim, ylim=ylim)
	par(mar=c(0, mar[2], 0, mar[4]))
	if (xbox > 0) hbox(.x) 
#	else plot(0, 0, xlab="", ylab="", axes=FALSE, type="n", xlim=xlim, ylim=ylim)
	par(mar=c(mar[1:2], top, mar[4]))
	plot(.x, .y, xlab=xlab, ylab=ylab, las=las, log=log, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab,
			cex.main=cex.main, cex.sub=cex.sub, type="n", xlim=xlim, ylim=ylim, ...)
	if(grid){
		grid(lty=1, equilogs=FALSE)
		box()}
	n.groups <- length(levels(groups))
	if (n.groups > length(col)) stop("number of groups exceeds number of available colors")
	if (length(col) == 1) col <- rep(col, 2)
	indices <- NULL
	range.x <- if (logged("x")) range(log(.x), na.rm=TRUE) else range(.x, na.rm=TRUE)
	for (i in 1:n.groups){
		subs <- groups == levels(groups)[i]
		points(if (is.null(jitter$x) || jitter$x == 0) .x[subs] else jitter(.x[subs], factor=jitter$x), 
				if (is.null(jitter$y) || jitter$y == 0) .y[subs] else jitter(.y[subs], factor=jitter$y), 
				pch=pch[i], col=col[if (n.groups == 1) 2 else i], cex=cex)
		if (by.groups){
			if (smooth) lowess.line(.x[subs], .y[subs], col=col[i], span=span)
			if (is.function(reg.line)) reg(.x[subs], .y[subs], col=col[i])
			if (ellipse) {
				X <- na.omit(data.frame(x=.x[subs], y=.y[subs]))
				if (logged("x")) X$x <- log(x)
				if (logged("y")) X$y <- log(y)
				with(X, dataEllipse(x, y, plot.points=FALSE, lwd=1, log=log,
								levels=levels, col=col[i], robust=robust))
			}
			if (id.method[1] != "identify") indices <- c(indices,
						showLabels(.x[subs], .y[subs], labels=labels[subs], id.method=id.method,
								id.n=id.n, id.cex=id.cex, id.col=col[i]))
		}}
	if (!by.groups){
		if (smooth) lowess.line(.x, .y, col=col[1], span=span)
		if (is.function(reg.line)) reg(.x, .y, col=col[3])
		if (ellipse) {
			X <- na.omit(data.frame(x=.x, y=.y))
			if (logged("x")) X$x <- log(X$x)
			if (logged("y")) X$y <- log(X$y)
			with(X, dataEllipse(x, y, plot.points=FALSE, lwd=1, log=log, levels=levels, col=col[1],
							robust=robust))
		}
		if (id.method[1] != "identify") indices <- showLabels(
					.x, .y, labels=labels, 
					id.method=id.method, id.n=id.n, id.cex=id.cex, id.col=id.col)
	}
	if (legend.plot) {
		xpd <- par(xpd=TRUE)
		on.exit(par(xpd=xpd), add=TRUE)
		usr <- par("usr")
		if (missing(legend.coords)){
			legend.x <- if (logged("x")) 10^(usr[1]) else usr[1]
			legend.y <- if (logged("y")) 10^(usr[4] + 1.2*top*strheight("x")) else usr[4] + 1.2*top*strheight("x")
			legend.coords <- list(x=legend.x, y=legend.y)
		}
		legend(legend.coords, legend=levels(groups), 
				pch=pch, col=col[1:n.groups], pt.cex=cex, cex=cex.lab, title=legend.title, bg="white")
	}
	if ("smooth" %in% err) warning("could not fit smooth")
	if ("spread" %in% err) warning("could not smooth spread")
	if (id.method[1] == "identify") indices <- showLabels(.x, .y, labels, 
				id.method=id.method, id.n=length(.x), id.cex=id.cex, id.col=id.col)
	if (is.null(indices)) invisible(indices) else if (is.numeric(indices)) sort(indices) else indices
} 

spBCA <- function(...) scatterplotBCA(...)

# fancy scatterplot matrices (J. Fox)

# 2010-09-04: J. Fox: changed color choice
# 2010-09-16: fixed point color when col is length 1
# 2011-02-25: Dan Putler: added a new color for the LS line

scatterplotMatrixBCA <- function(x, ...){
	UseMethod("scatterplotMatrixBCA")
}

scatterplotMatrixBCA.formula <- function (x, data=NULL, subset, labels, ...) {
	m <- match.call(expand.dots = FALSE)
	if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
		m$data <- as.data.frame(data)
	m$labels <- m$formula <- m$... <- NULL
	m[[1]] <- as.name("model.frame")
	if (!inherits(x, "formula") | length(x) != 2) 
		stop("invalid formula")
	rhs <- x[[2]]
	if ("|" != deparse(rhs[[1]])){
		groups <- FALSE
	}
	else{
		groups <- TRUE
		x<-as.character(c(x))
		x<-as.formula(sub("\\|", "+", x))   
	}
	m$formula <-x
	if (missing(data)){ 
		X <- na.omit(eval(m, parent.frame()))
		if (missing(labels)) labels <- gsub("X", "", row.names(X))
	}
	else{
		if (!missing(labels)) row.names(data) <- labels
		X <- eval(m, parent.frame())
		labels <- row.names(X)
	}
	if (!groups) scatterplotMatrixBCA(X, labels=labels, ...)
	else{
		ncol<-ncol(X)
		scatterplotMatrixBCA.default(X[, -ncol], groups=X[, ncol], labels=labels, ...)
	}
}

scatterplotMatrixBCA.default <- function(x, var.labels=colnames(x), 
	diagonal=c("density", "boxplot", "histogram", "oned", "qqplot", "none"), adjust=1, nclass,
#	plot.points=TRUE, smooth=TRUE, spread=smooth && !by.groups, span=.5, loess.threshold=5, reg.line=lm, 
	plot.points=TRUE, smooth=TRUE, spread=smooth && !by.groups, span=.5, loess.threshold=2, reg.line=lm, 
	transform=FALSE, family=c("bcPower", "yjPower"),
	ellipse=FALSE, levels=c(.5, .95), robust=TRUE,
	groups=NULL, by.groups=FALSE, 
	labels, id.method="mahal", id.n=0, id.cex=1, id.col=palette()[1],
#	col=if (n.groups == 1) palette()[2:1] else rep(palette(), length=n.groups),
	col=if (n.groups == 1) palette()[c(2,1,3)] else rep(palette(), length=n.groups),
	pch=1:n.groups, lwd=2, lwd.smooth=lwd, lwd.spread=lwd, lty=1, lty.smooth=lty, lty.spread=2,
	cex=par("cex"), cex.axis=par("cex.axis"), cex.labels=NULL, 
	cex.main=par("cex.main"), 
	legend.plot=length(levels(groups)) > 1, row1attop=TRUE, ...){
	spread # force evaluation
	if (id.method == "identify") stop("interactive point identification not permitted")
	family <- match.arg(family)
	err <- ""
	lowess.line <- function(x, y, col, span) {
		warn <- options(warn=-1)
		valid <- complete.cases(x, y)
		x <- x[valid]
		y <- y[valid]
		ord <- order(x)
		x <- x[ord]
		y <- y[ord]
#		if (length(unique(x)) < lowess.threshold || length(unique(y)) < lowess.threshold) return()
		if (length(unique(y)) < loess.threshold) return()
		if (!spread){
			fit <- try(loess.smooth(x, y, span=span), silent=TRUE)
			if (class(fit) == "try-error"){
				err <<- c(err, "smooth")
				options(warn)
				return()
			}
			lines(fit$x, fit$y, lty=lty.smooth, lwd=lwd.smooth, col=col)
		}
		else{
			fit <- try(loess(y ~ x, degree=1, family="symmetric", span=span), silent=TRUE)
			if (class(fit) == "try-error"){
				err <<- c(err, "smooth")
				options(warn)
				return()
			}
			res <- residuals(fit)
			pos <- res > 0
			pos.fit <- try(loess(res^2 ~ x, span=span, degree=0, family="gaussian", subset=pos), silent=TRUE)
			neg.fit <- try(loess(res^2 ~ x, span=span, degree=0, family="gaussian", subset=!pos), silent=TRUE)
			if (class(pos.fit) == "try-error" || class(neg.fit) == "try-error"){
				err <<- c(err, "spread")
				options(warn)
				return()
			}
			lines(x, fitted(fit), lty=lty.smooth, lwd=lwd.smooth, col=col)
			y.pos <- fitted(fit)[pos] + sqrt(fitted(pos.fit))
			lines(x[pos], y.pos, lty=lty.spread, lwd=lwd.spread, col=col)
			y.neg <- fitted(fit)[!pos] - sqrt(fitted(neg.fit))
			lines(x[!pos], y.neg, lty=lty.spread, lwd=lwd.spread, col=col)
		}
		options(warn)
	}
	if (missing(labels)){
		labels <- rownames(x)
		if (is.null(labels)) labels <- as.character(seq(length.out=nrow(x)))
	}
	if (!(missing(groups))){
		x <- na.omit(data.frame(groups, labels, x, stringsAsFactors=FALSE))
#		groups <- as.factor(as.character(x[, 1]))
		if (!is.factor(groups)) groups <- as.factor(as.character(x[,1]))
		labels <- x[, 2]
		x <- x[, -(1:2)]
	}
	else {
		x <- na.omit(data.frame(labels, x, stringsAsFactors=FALSE))
		labels <- x[, 1]
		x <- x[, -1]
	}
	if (missing(nclass)) nclass <- "FD"
	reg <- function(x, y, col){
		mod<-reg.line(y ~ x)
		y.hat <- fitted.values(mod)
		x <- model.matrix(mod)[,2]
		min <- which.min(x)
		max <- which.max(x)
		lines(c(x[min], x[max]), c(y.hat[min], y.hat[max]), lty=lty, lwd=lwd, col=col)
	}
	legendPlot <- function(){
		usr <- par("usr")
		legend("bottomleft", bg="white",
			legend=levels(groups), pch=pch, col=col[1:n.groups],
			cex=cex)
	}	
	do.legend <- legend.plot	
# The following panel function adapted from Richard Heiberger
	panel.density <- function(x, ...){
		dens.x <- density(x, adjust = adjust)
		lines(dens.x$x, min(x) + dens.x$y * diff(range(x))/diff(range(dens.x$y)))
		rug(x)
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	panel.histogram <- function(x, ...){
		par(new=TRUE)
		hist(x, main="", axes=FALSE, breaks=nclass, col=col[1])
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	panel.boxplot <- function(x, ...){
		par(new=TRUE)
		boxplot(x, axes=FALSE, main="", col=col[1])
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	# The following panel function adapted from Richard Heiberger
	panel.oned <- function(x, ...) {
		range <- range(x)
		delta <- diff(range)/50
		y <- mean(range)
		segments(x - delta, x, x + delta, x, col = col[1])
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	panel.qqplot <- function(x, ...){
		par(new=TRUE)
		qqnorm(x, axes=FALSE, xlab="", ylab="", main="", col=col[1])
		qqline(x)
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	panel.blank <- function(x, ...){
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	which.fn <- match(match.arg(diagonal),
		c("density", "boxplot", "histogram", "oned", "qqplot", "none"))
	diag <- list(panel.density, panel.boxplot, panel.histogram, panel.oned, panel.qqplot, panel.blank)[[which.fn]]
	groups <- as.factor(if(missing(groups)) rep(1, length(x[, 1])) else groups)
	n.groups <- length(levels(groups))
	if (n.groups > length(col)) stop("number of groups exceeds number of available colors")
	if (length(col) == 1) col <- rep(col, 2)
	if (transform != FALSE | length(transform) == ncol(x)){
		if (transform == TRUE & length(transform) == 1){
			transform <- if (by.groups) coef(powerTransform(as.matrix(x) ~ groups, family=family), round=TRUE)
				else coef(powerTransform(x, family=family), round=TRUE)
		}
		for (i in 1:ncol(x)){
			x[, i] <- if (family == "bcPower") 
					bcPower(x[, i], transform[i])
				else yjPower(x[, i], transform[i])
			var.labels[i] <- paste(var.labels[i], "^(", round(transform[i],2), ")", sep="")
		}
	}
	labs <- labels
	pairs(x, labels=var.labels, 
		cex.axis=cex.axis, cex.main=cex.main, cex.labels=cex.labels, cex=cex,
		diag.panel=diag, row1attop = row1attop,
		panel=function(x, y, ...){ 
			for (i in 1:n.groups){
				subs <- groups == levels(groups)[i]
				if (plot.points) points(x[subs], y[subs], pch=pch[i], col=col[if (n.groups == 1) 2 else i], cex=cex)
				if (by.groups){
					if (smooth) lowess.line(x[subs], y[subs], col=col[i], span)
					if (is.function(reg.line)) reg(x[subs], y[subs], col=col[i])
					if (ellipse) dataEllipse(x[subs], y[subs], plot.points=FALSE, 
							levels=levels, col=col[i], robust=robust, lwd=1)
					showLabels(x[subs], y[subs], labs[subs], id.method=id.method, 
					    id.n=id.n, id.col=col[i], id.cex=id.cex)
					#if (id.method != "none") 
					#	showLabelsScatter(x[subs], y[subs], labs[subs], id.var=id.var, id.method=id.method,
					#		id.n=id.n, id.col=col[i], id.cex=id.cex)
				}
			}
			if (!by.groups){
#				if (is.function(reg.line)) abline(reg.line(y ~ x), lty=lty, lwd=lwd, col=col[2])
				if (is.function(reg.line)) abline(reg.line(y ~ x), lty=lty, lwd=lwd, col=col[3])
				if (smooth) lowess.line(x, y, col=col[1], span)
				if (ellipse) dataEllipse(x, y, plot.points=FALSE, levels=levels, col=col[1],
						robust=robust, lwd=1)
				showLabels(x, y, labs, id.method=id.method, 
					    id.n=id.n, id.col=id.col, id.cex=id.cex)
				#if (id.method != "none") showLabelsScatter(x, y, labs, id.var=id.var, id.method=id.method, 

				#		id.n=id.n, id.col=col[1], id.cex=id.cex)
			}
		}, ...
	)
	if ("smooth" %in% err) warning("some smooths failed")
	if ("spread" %in% err) warning("some spreads failed")
}

spmBCA <- function(x, ...){
	scatterplotMatrixBCA(x, ...)
} 

