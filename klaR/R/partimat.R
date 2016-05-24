partimat <- function(x, ...) 
    UseMethod("partimat")


partimat.default <- function(x, grouping, method = "lda", prec = 100, 
    nplots.vert, nplots.hor, main = "Partition Plot", name, mar,
    plot.matrix = FALSE, plot.control = list(), ...){
  
    nvar <- ncol(x)
    if(nvar < 2) stop("at least 2 variables required")
    if(nlevels(grouping) < 2) stop("at least two classes required")
    nobs <- nrow(x)
    if(missing(name)) name <- colnames(x)
    # plot in scatterplot matrix
    if(plot.matrix){
        plot.new()
        if(missing(mar)) mar <- rep(0, 4)
        opar <- par(mfrow = c(nvar, nvar), mar = mar, oma = rep(3, 4), xpd = NA)
        on.exit(par(opar))
        for (i in 2:nvar)
            for (j in 1:(i-1))
            {
                par(mfg = c(i, j))
                drawparti(grouping, x[,j], x[,i], method = method, 
                    prec = prec, legend.err = plot.matrix, xlab="", ylab="", 
                    plot.control = c(xaxt="n", yaxt="n", plot.control), ...)
                if(j == 1) axis(2)
                if(i == nvar) axis(1)
                par(mfg = c(j, i))
                drawparti(grouping, x[,i], x[,j], method = method, 
                    prec = prec, legend.err = plot.matrix, xlab="", ylab="", 
                    plot.control = c(xaxt="n", yaxt="n", plot.control), ...)
                if(j == 1) axis(3) 
                if(i == nvar) axis(4)
            }
        for (i in 1:nvar)
            {
                par(mfg = c(i, i))
                plot(x[,i], x[,i], type = "n", xaxt="n", yaxt="n", xlab="", ylab="")
                if(i == 1){
                    axis(2); axis(3)
                }
                else if(i == nvar){
                    axis(1); axis(4)
                }
                mxi <- mean(range(x[,i]))
                do.call("text", c(list(mxi, mxi, name[i]), plot.control))
            }
    } # Plotting fuzzy in rows and columns to save space:
    else{
        ncomb <- round(0.5 * nvar * (nvar-1))
        if (missing(nplots.hor) && missing(nplots.vert)){
            nplots.hor<-ceiling(sqrt(ncomb))
            nplots.vert<-floor(sqrt(ncomb))
        }
        else if (missing(nplots.hor)) nplots.hor<-ceiling(ncomb/nplots.vert)
        else if (missing(nplots.vert)) nplots.vert<-ceiling(ncomb/nplots.hor)
        vars <- matrix(ncol=ncomb,nrow=2*nobs)
        varname <- matrix(ncol=ncomb,nrow=2)
        k <- 1
        for (i in 2:nvar)
            for (j in 1:(i-1))
            {
                vars[,k] <- c(x[,i], x[,j])
                varname[,k] <- c(name[i], name[j])
                k <- k + 1
            }

        if(missing(mar)) mar <- c(5.1, 4.1, 2.1, 1.1)
        opar <- par(mfrow = c(nplots.vert, nplots.hor), mar = mar, 
            oma = c(0, 0, !is.null(main), 0))
        on.exit(par(opar))
        
        sapply(1:ncomb, function(k) 
            drawparti(grouping = grouping, x = vars[(1:nobs), k], 
                y = vars[(nobs+1):(2*nobs), k], method = method, 
                xlab = varname[1,k], ylab = varname[2,k], prec = prec, 
                legend.err = plot.matrix, plot.control = plot.control, ...)
        )
        par(mfrow=c(1,1))
        title(main = main, outer = TRUE)
    }
    invisible()
}

partimat.formula <- function(formula, data = NULL, ..., subset, na.action = na.fail) 
{
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
        x <- x[, -xint, drop = FALSE]
    res <- partimat.default(x, grouping, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("partimat")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
    invisible()
}

partimat.matrix<-function (x, grouping, ..., subset, na.action = na.fail) 
{
    if (!missing(subset)) {
        x <- x[subset, , drop = FALSE]
        grouping <- grouping[subset]
    }
    if (!missing(na.action)) {
        dfr <- na.action(structure(list(g = grouping, x = x), 
            class = "data.frame"))
        grouping <- dfr$g
        x <- dfr$x
    }
    res <- partimat.default(x, grouping, ...)
    cl <- match.call()
    cl[[1]] <- as.name("partimat")
    res$call <- cl
    res
    invisible()
}

partimat.data.frame<-function (x, ...) 
{
   res <- partimat.matrix(structure(data.matrix(x), class = "matrix"), 
        ...)
    cl <- match.call()
    cl[[1]] <- as.name("partimat")
    res$call <- cl
    res
    invisible()    
}

drawparti <- function(grouping, x, y, method = "lda", prec = 100, 
    xlab=NULL, ylab=NULL, col.correct = "black", col.wrong = "red", 
    col.mean = "black", col.contour = "darkgrey", gs = as.character(grouping), 
    pch.mean = 19, cex.mean = 1.3, print.err = 0.7, legend.err = FALSE,
    legend.bg = "white", imageplot = TRUE, image.colors = cm.colors(nc), 
    plot.control = list(), ...){                       
    #grouping: class vec.
    #x: first data vec.
    #y: second data vec.
    #prec: nr. of hor/vert splits.
    
    z <- switch(method,
        lda = lda(grouping ~ x + y,...),
        qda = qda(grouping ~ x + y,...),
        svmlight = svmlight(grouping ~ x + y,...),
        rda = rda(grouping~ x + y, 
            data = cbind.data.frame("grouping" = grouping, "x" = x, "y" = y), ...),
        sknn = sknn(grouping ~ x + y,...),
        rpart = rpart::rpart(grouping~ x + y,...),
#### we need to move disco to CRAN!!!
#        disco = disco(grouping ~ x + y,
#            data = cbind.data.frame("grouping" = grouping, "x" = x, "y" = y), ...),
        naiveBayes = e1071::naiveBayes(grouping~ x + y, 
            data = cbind.data.frame("grouping" = grouping, "x" = x, "y" = y), ...),
        stop("method not yet supported"))

    # Build a grid on the 2 coordinates
    xg <- seq(min(x), max(x), length = prec)
    yg <- seq(min(y), max(y), length = prec)
    grd <- expand.grid(x = xg, y = yg)
    # Calcultate posterior Probabilities on grid points
    temp <- switch(method,
        lda = predict(z, grd,...)$post,
        qda = predict(z, grd,...)$post,
        svmlight = e.scal(predict(z, grd,...)$post)$sv,
        rda = predict(z, grd, posterior=TRUE, aslist=TRUE)$post,
        rpart = predict(z, grd, ...),
        sknn = predict(z, grd, ...)$post,
#        disco = predict(z, grd, ...)$post,
        naiveBayes = predict(z, grd , type="raw", ...),
        stop("method not yet supported"))
    khead <- switch(method,
        lda = predict(z, data.frame(cbind(x,y)),...)$class,
        qda = predict(z, data.frame(cbind(x,y)),...)$class,
        svmlight = predict(z, data.frame(cbind(x,y)),...)$class,
        rda = predict(z, data.frame(cbind(x,y)), posterior=TRUE, aslist=TRUE)$class,
        rpart = predict(z, data.frame(cbind(x,y)), type="class", ...),
        sknn = predict(z, data.frame(cbind(x,y)),...)$class,
#        disco = predict(z, data.frame(cbind(x,y)),...)$class,
        naiveBayes = predict(z, data.frame(cbind(x,y)), ...),
        stop("method not yet supported"))

    colorw <- grouping != khead
    err <- round(mean(colorw), 3)
    color <- ifelse(colorw, col.wrong, col.correct)
    if(is.character(gs) || is.factor(gs)) gs <- substr(gs, 1, 1)

    nc <- ncol(temp)
    if(imageplot){
        do.call("image", c(list(xg, yg, matrix(apply(temp, 1, which.max), ncol = prec), 
            main = NULL, col = image.colors, breaks = (0:nc) + .5, 
            xlab = xlab, ylab = ylab), plot.control))
        do.call("points", c(list(x, y, pch = gs, col = color), plot.control))
        box()
    }
    else 
        do.call("plot", c(list(x, y, pch = gs, col = color, main = NULL, xlab = xlab, ylab = ylab), plot.control))
    if((method=="lda") || (method=="qda")) 
        points(z$means, pch = pch.mean, cex = cex.mean, col = col.mean)
    
    # For each class calculate the difference between prob. and max(prob) for other class,
    # so, the obs is assigned to class iff diff>0
    if(!imageplot)
        for(i in 1:ncol(temp)){
            dummy <- temp[,i] - apply(temp[ , -i, drop = FALSE], 1, max)            
        # Draw contour line at hight=0, i.e. class border
            contour(xg, yg, matrix(dummy, ncol = prec), levels = 0, 
                add = TRUE, drawlabels = FALSE, col = col.contour)
        }

    if(print.err){
        if(legend.err)
            legend(par("usr")[1], par("usr")[4], 
                legend = paste("Error:", err), bg = legend.bg, cex = print.err)
        else
            mtext(paste("app. error rate:", err), 3, cex = print.err)
    }

}
