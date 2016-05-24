ftsm = function (y, order = 6, ngrid = max(500, ncol(y$y)), method = c("classical", 
    "M", "rapca"), mean = TRUE, level = FALSE, lambda = 3, weight = FALSE, 
    beta = 0.1, ...) 
{
	if(length(colnames(y$y)) > 0)
	{
        y$time = ts(as.numeric(colnames(y$y)), start = head(as.numeric(colnames(y$y)),1),
					end = tail(as.numeric(colnames(y$y)),1), frequency = 1/diff(as.numeric(colnames(y$y)))[1])
  	}
	else
	{
		y$time = ts(1:ncol(y$y), start=1, end = ncol(y$y))
	}
    method <- match.arg(method)
    if (!mean & !level & order < 1) 
        stop("No model to fit")
    if (order < 0) 
        stop("Order must not be negative")
    n <- ncol(y$y)
    m <- length(y$x)
    if (weight == FALSE) {
        y.pca <- fdpca(y$x, y$y, order = order, ngrid = ngrid, 
            method = method, mean = mean, level = level, lambda = lambda, 
            ...)
        mean.se = y.pca$mean.se
    }
    if (weight == TRUE) {
        q <- beta * (1 - beta)^(0:(n - 1))
        my <- apply(y$y, 1, weighted.mean, w = rev(q))
        newy <- sweep(t(y$y), 2, my)
        x = y$x
        y2 = y$y
        n = dim(y$y)[2]
        wq = diag(rev(q))
        newy2 = wq %*% newy
        dummy = svd(newy2)
        load = dummy$v[, 1:(order)]
        sco = newy %*% load
        yy <- matrix(NA, nrow = ngrid, ncol = n)
        xx = seq(min(x), max(x), l = ngrid)
        for (i in 1:n) {
            miss <- is.na(newy[i, ])
            yy[, i] <- spline(x[!miss], t(newy)[!miss, i], n = ngrid)$y
        }
        V <- rep(1, ngrid) * t(yy)
        dummy <- La.svd(V)
        varprop = (dummy$d[1:order])^2/(sum((dummy$d)^2))
        mean.se = approx(xx, sqrt(apply(yy, 1, var)/n), xout = x)$y
    }
    ytsp <- tsp(y$time)
    if (weight == TRUE) {
        colmeanrm = matrix(colMeans(sco), dim(sco)[2], 1)
        scomeanrm = sweep(sco, 2, colmeanrm)
        coeff <- ts(cbind(rep(1, dim(y$y)[2]), scomeanrm), start = ytsp[1], 
            frequency = ytsp[3])
    }
    else {
        coeff <- ts(y.pca$coeff, start = ytsp[1], frequency = ytsp[3])
    }
    if (weight == TRUE) {
        my = my + load %*% colmeanrm
        basis <- cbind(my, load)
        if(order == 0)
        {
            colnames(basis) = "mean"
            colnames(coeff) = "mean"
        }
        else
        {
	        colnames(basis) = c("mean", paste("phi", 1:order, sep = ""))
    	    colnames(coeff) = c("mean", paste("beta", 1:order, sep = ""))
    	}
        fits <- fts(y$x, sweep(load %*% t(scomeanrm), 1, my, 
            "+"), start = ytsp[1], frequency = ytsp[3], xname = y$xname, 
            yname = paste("Fitted", y$yname))
    }
    else {
        basis <- y.pca$basis
        fits <- fts(1:length(y$x), basis %*% t(coeff), start = ytsp[1], 
            frequency = ytsp[3], xname = y$xname, yname = paste("Fitted", 
                y$yname))
        fits$x = y$x
        if(order == 0)
        {
            colnames(basis) = "mean"
            colnames(coeff) = "mean"
        }
        else
        {
        	if(mean)
        	{
		        colnames(basis) = c("mean", paste("phi", 1:order, sep = ""))
    		    colnames(coeff) = c("mean", paste("beta", 1:order, sep = ""))
    		}
    		else
    		{
		   		colnames(basis) = c(paste("phi", 1:order, sep = ""))
    		    colnames(coeff) = c(paste("beta", 1:order, sep = ""))
    		}
    	}
    }
    # rownames(basis) <- paste(y$x)
    res <- fts(y$x, y$y - fits$y, start = ytsp[1], frequency = ytsp[3], 
        xname = y$xname, yname = paste("Residuals", y$yname))
    if (weight == FALSE) {
        out <- list(x1 = as.numeric(colnames(y$y)), y1 = as.numeric(rownames(y$y)), 
            y = fts(y$x, y$y, xname = y$xname, yname = y$yname), 
            basis = basis, coeff = coeff, fitted = fits, residuals = res, 
            varprop = y.pca$varprop, wt = ts(y.pca$weights, start = ytsp[1], 
                frequency = ytsp[3]), v = ts(y.pca$v, start = ytsp[1], 
                frequency = ytsp[3]), basis2 = y.pca$basis2, 
            coeff2 = y.pca$coeff2, mean.se = mean.se, call = match.call())
    }
    else {
        out <- list(x1 = as.numeric(colnames(y$y)), y1 = as.numeric(rownames(y$y)), 
            y = fts(y$x, y2, xname = y$xname, yname = y$yname), 
            basis = basis, coeff = coeff, fitted = fits, residuals = res, 
            varprop = varprop, wt = rev(q), mean.se = mean.se, 
            call = match.call())
    }
    return(structure(out, class = c("ftsm", "fm")))
}
