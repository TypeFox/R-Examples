plot.modwt <- function( x , levels = NULL, draw.boundary = FALSE,
                       type = "stack", col.plot = "black", col.boundary = "red",
                       X.xtick.at = NULL, X.ytick.at = NULL, Stack.xtick.at = NULL,
                       Stack.ytick.at = NULL, X.xlab = "t", y.rlabs = TRUE, plot.X = TRUE,
                       plot.W = TRUE, plot.V = TRUE, ...)
{
    stackplot.modwt <- function ( x , w.range, v.range, col.plot, col.boundary, draw.boundary, X.xtick.at, X.ytick.at, Stack.xtick.at, Stack.ytick.at, X.xlab = "t", plot.X = TRUE)
    {
        innerplot <- function(x, y, type = "l", xtick.at, ytick.at)
        {
            if(is.null(xtick.at) == FALSE || is.null(ytick.at) == FALSE) {
            plot(x, y, type = "l", axes = FALSE, frame.plot = TRUE)
                if(is.null(xtick.at) == FALSE) {
                    axis(1, at = axTicks(1, xtick.at))
                    xtickrate <- xtick.at
                }
                else {
                    axis(1)
                    xtickrate <- par("xaxp")
                }
                if(is.null(ytick.at) == FALSE) {
                    axis(2, at = axTicks(2, ytick.at))                
                    ytickrate <- ytick.at
                }
                else {
                    axis(2)
                    ytickrate <- par("yaxp")
                }
            }
            else {
                plot(x, y, type = "l")
                xtickrate <- par("xaxp")
                ytickrate <- par("yaxp")
            }
            tickrate <- list(xtick = xtickrate, ytick = ytickrate)
            tickrate
        }

        if(plot.X) {    
            nf <- layout(matrix(c(2,2,1,1), 2, 2, byrow=TRUE), c(1,2), c(2,1), TRUE)
            par(mai = c(.6, .4, .1, .6))
            
            if( x @class.X == "ts" ||  x @class.X == "mts") {
                x.range <-  x @attr.X$tsp[1]: x @attr.X$tsp[2]
            }
            else{
                x.range <- 1:dim( x @series)[1]
            }
            tickrate <- innerplot(x.range,  x @series[,1], type = "l", X.xtick.at, X.ytick.at)
            right.usrplotrange <- par()$usr[2] - par()$usr[1]
            NDCplotrange <- par()$plt[2] - par()$plt[1] 
            marginpos <- (1-par()$plt[2])/2
            right.usrlabelpos <- ((marginpos*right.usrplotrange)/NDCplotrange) + par()$usr[2]
            text(right.usrlabelpos, 0, "X", xpd = TRUE)
            mtext(X.xlab, side = 1, line = 2)

            par(mai = c(0, .4, .1, .6))
        }
        if(plot.X == FALSE) {
            par(mai = c(.4, .4, .1, .6))
            if(is.null(Stack.xtick.at) == FALSE) {
                xtickrate <- Stack.xtick.at
            } 
            else {
                xtickrate <- NULL
            }
            if(is.null(Stack.ytick.at) == FALSE) {
                ytickrate <- Stack.ytick.at
            }
            else {
                 ytickrate <- NULL
            }
            tickrate <- list(xtick = xtickrate, ytick = ytickrate)
        }
        if(is.null(w.range) == FALSE) {    
            gammawave = wt.filter.shift( x @filter, w.range, wavelet = TRUE, modwt = TRUE)
        }

        if(is.null(v.range) == FALSE) {    
            gammascale = wt.filter.shift( x @filter, v.range, wavelet = FALSE, modwt = TRUE)
        }

        if(y.rlabs) {
            rightlabels <- labels.modwt(w.range = w.range, v.range = v.range, gammah = gammawave, gammag = gammascale)
        }
        else {
            rightlabels <- NULL
        }

        if (draw.boundary) {
            matrixlist <- list(modwt = as.matrix.modwt( x , w.range, v.range), posbound = boundary.as.matrix.modwt( x , w.range, v.range, positive = TRUE), negbound = boundary.as.matrix.modwt( x ,  w.range, v.range, positive = FALSE))
            col <- c(col.plot, col.boundary, col.boundary)
            stackplot(matrixlist, y = NULL, y.rlabs = rightlabels, type = c("l", "h", "h"), col = col, xtick.at = tickrate$xtick, ytick.at = tickrate$ytick) 
        }
        else {
            matrixlist <- list(modwt = as.matrix.modwt( x , w.range, v.range))
            col <- col.plot
            stackplot(matrixlist, y = NULL, y.rlabs = rightlabels, type = "l", col = col, xtick.at = tickrate$xtick, ytick.at = tickrate$ytick) 
        }
    }

    boundary.as.matrix.modwt <- function( x , w.range, v.range, positive = TRUE) 
    {
        if(is.null(w.range) == FALSE) {    
            wavecoefmatrix <- array(NA, c(2*dim( x @series)[1], length(w.range)))
            Wjplot <- rep(NA, 2*dim( x @series)[1])
            wavecoefmatrix.index <- 0
            W.Ljs <- ((2^w.range) - 1)*( x @filter@L - 1) + 1
    
            for (j in w.range) 
            {
                wavecoefmatrix.index <- wavecoefmatrix.index + 1

                if(positive) {
                    boundaryheight <- max( x @W[[j]])
                }
                else {
                    boundaryheight <- min( x @W[[j]])
                }
    
                leftspace <- rep(NA, 2*(W.Ljs[wavecoefmatrix.index] - 2 - vjH.modwt( x @filter@L, j, dim( x @series)[1])) - 1)
                rightspace <- rep(NA, 2*(vjH.modwt( x @filter@L, j, dim( x @series)[1])))
                middlespace <- rep(NA, 2*dim( x @series)[1] - 2 - length(leftspace) - length(rightspace))
                Wjplot <- c(leftspace, boundaryheight, middlespace, boundaryheight, rightspace)
                wavecoefmatrix[,wavecoefmatrix.index] <- Wjplot
            }

            rownames(wavecoefmatrix) <- seq(.5, dim( x @series)[1], by = .5)
        }


        if(is.null(v.range) == FALSE) {
            scalecoefmatrix <- array(NA, c(2*dim( x @series)[1], length(v.range)))
            Vjplot <- rep(NA, 2*dim( x @series)[1])
            scalecoefmatrix.index <- 0
            V.Ljs <- ((2^v.range) - 1)*( x @filter@L - 1) + 1

            for(j in v.range) 
            {
                scalecoefmatrix.index <- scalecoefmatrix.index + 1

                Vj <-  x @V[[j]][,1] - mean( x @V[[j]][,1])

                if(positive) {
                    boundaryheight <- max(Vj)
                }
                else {
                    boundaryheight <- min(Vj)
                }
     
                leftspace <- rep(NA, 2*(V.Ljs[scalecoefmatrix.index] - 2 - vjG.modwt( x @filter@L, j, dim( x @series)[1])) - 1)
                rightspace <- rep(NA, 2*(vjG.modwt( x @filter@L, j, dim( x @series)[1])))
                middlespace <- rep(NA, 2*dim( x @series)[1] - 2 - length(leftspace) - length(rightspace))
                Vjplot <- c(leftspace, boundaryheight, middlespace, boundaryheight, rightspace)
                scalecoefmatrix[,scalecoefmatrix.index] <- Vjplot
            }
            rownames(scalecoefmatrix) <- seq(.5, dim( x @series)[1], by = .5)
        }

        if(is.null(w.range) == FALSE && is.null(v.range) == FALSE) {  
            results <- cbind(wavecoefmatrix, scalecoefmatrix)    
        }
        if(is.null(w.range) == FALSE && is.null(v.range)) {  
            results <- wavecoefmatrix
        }
        if(is.null(w.range) && is.null(v.range) == FALSE) {  
            results <- scalecoefmatrix
        }

        results
    }

    as.matrix.modwt <- function ( x , w.range, v.range)
    {
        if(is.null(w.range) == FALSE) {
            wavecoefmatrix <- array(NA, c(dim( x @series)[1], length(w.range)))
            wavecoefmatrix.index <- 0
                for (j in w.range) {
                    wavecoefmatrix.index <- wavecoefmatrix.index + 1
                    Wjplot <-  x @W[[j]][,1]
                    Wjplot <- levelshift.modwt(Wjplot, wt.filter.shift( x @filter, j, wavelet = TRUE, modwt=TRUE))
                    wavecoefmatrix[,wavecoefmatrix.index] <- Wjplot
                }
            rownames(wavecoefmatrix) <- 1:dim( x @series)[1]
        }

        if(is.null(v.range) == FALSE) {
            scalecoefmatrix <- array(NA, c(dim( x @series)[1], length(v.range)))
            scalecoefmatrix.index <- 0

            for(k in v.range) {
                scalecoefmatrix.index <- scalecoefmatrix.index + 1
                Vjplot <-  x @V[[k]][,1] - mean( x @V[[k]][,1])
                Vjplot <- levelshift.modwt(Vjplot, wt.filter.shift( x @filter, k, wavelet = FALSE, modwt=TRUE))
                scalecoefmatrix[,scalecoefmatrix.index] <- Vjplot                   
            }
            rownames(scalecoefmatrix) <- 1:dim( x @series)[1]    
        }

        if(is.null(w.range) == FALSE && is.null(v.range) == FALSE) {  
            results <- cbind(wavecoefmatrix, scalecoefmatrix)    
        }
        if(is.null(w.range) == FALSE && is.null(v.range)) {  
            results <- wavecoefmatrix
        }
        if(is.null(w.range) && is.null(v.range) == FALSE) {  
            results <- scalecoefmatrix
        }
        
        results
    }        

    labels.modwt <- function (w.range = NULL, v.range = NULL, gammah = NULL, gammag = NULL)
    {
        verticallabel <- list()
 
        if(is.null(w.range) == FALSE && is.null(gammah) == FALSE) {
            for (j in 1:length(w.range)) {
                label <- substitute(paste(T^-gamma,W[level]), list(gamma = gammah[j], level = w.range[j]))
                verticallabel <- c(verticallabel, label)
            }
        }
       
        if(is.null(v.range) == FALSE && is.null(gammag) == FALSE) {
            for (i in 1:length(v.range)) {
                label <- substitute(paste(T^-gamma,V[level]), list(gamma = gammag[i], level = v.range[i]))
                verticallabel <- c(verticallabel, label)
            }
        }
        results <- verticallabel
    
        results
    }

    levelshift.modwt <- function(level, shift)
    {
        if(shift != 0) {
            level <- c(level[(round(shift)+1):length(level)], level[1:round(shift)])
        }

        level
    }

    shift.modwt <- function(L, j, N)
    {
        Lj <- ((2^j)-1)*(L-1) + 1
        shift <- min(Lj - 2, N-1)

        shift
    }

    vjH.modwt <- function (L, j, N)
    {
        Lj <- ((2^j)-1)*(L-1) + 1
        if (L == 10 || L == 18) {
            vjH <- (-Lj/2) + 1
        }
        else if (L == 14) {
            vjH <- (-Lj/2) - 1
        }
        else {
            vjH <- -Lj/2
        }
        vjH  <- abs(vjH)
    }

    vjG.modwt <- function (L, j, N)
    {
        Lj <- ((2^j)-1)*(L-1) + 1
        if (L == 10 || L == 18) {
            vjG <- -((Lj-1)*L)/(2*(L-1))
        }
        else if (L == 14) {
            vjG <- -((Lj-1)*(L-4))/(2*(L-1))
        }
        else {
            vjG <- -((Lj-1)*(L-2))/(2*(L-1))
        }
        vjG <- abs(vjG)

        vjG
    }

    if (type == "stack") {  
        if(class( x ) != "modwt") {
            stop("Invalid argument: 'modwt' object must be of class modwt.")
        }
        if(is.null(levels)) {
            w.range <- 1: x @level
            v.range <- max(w.range)  
        }
        if(class(levels) == "numeric") {
            if(length(levels) == 1) {
                w.range <- 1:levels
                v.range <- max(w.range)  
            }
            else {
                w.range <- levels
                v.range <- max(w.range)
            }
        }
        if(class(levels) == "list") {
            if(length(levels) < 1) {
                w.range <- 1:x@level
                v.range <- max(w.range)
            }
            if(length(levels) == 1) {
                w.range <- levels[[1]]
                v.range <- max(w.range)
            }
            else {
            w.range <- levels[[1]]
            v.range <- levels[[2]]
            }
        }        
        if(class(levels) != "list" && class(levels) != "vector" && class(levels) != "numeric" && is.null(levels) == FALSE) {
            stop("Invalid argument: 'levels' must be numeric, vector, or list.")
        }
        if(plot.W == FALSE) {
            w.range <- NULL
        }
        if(plot.V == FALSE) {
            v.range <- NULL
        } 
        if(plot.W == FALSE && plot.V == FALSE) {
            stop("At least one of plot.W or plot.V must be TRUE")
        }
        if(is.null(w.range) == FALSE) {
            if(min(w.range) < 1 ||  x @level < max(w.range)) {
                stop("Invalid argument: elements of 'levels' must be compatible with the level of decomposition of the 'modwt' object.")
            } 
        }
        if(is.null(v.range) == FALSE) {
            if(min(v.range) < 1 ||  x @level < max(v.range)) {
                stop("Invalid argument: elements of 'levels' must be compatible with the level of decomposition of the 'modwt' object.")
            }  
        }
        stackplot.modwt( x , w.range, v.range, col.plot, col.boundary, draw.boundary = draw.boundary, X.xtick.at = X.xtick.at, X.ytick.at = X.ytick.at, Stack.xtick.at = Stack.xtick.at, Stack.ytick.at = Stack.ytick.at, X.xlab = X.xlab, plot.X = plot.X)
    }
    else {
        stop("Only the stackplot is currently implemented.")
    }	
}
