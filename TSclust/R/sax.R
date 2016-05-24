#################################
#############   SAX   ###########
#################################

#create a table with the breakpoints of dividing the gaussian(0,1)
#in n equiprobable parts
#used for generating equiprobable symbols in the SAX representation
SAX.breakpoints.table <- function (n) {
    qnorm(0:n/n)
}

#Dimensionality reduction
# x the series
# w the number of parts desired
PAA <- function( x, w) {
    if ((w - floor(w)) > 0) {
        stop("w (number of frames) must be an integer")
    }
    n <- length(x)
    if (w > n) {
        stop("cannot have more parts than the length of the series")
    }

    PAA <- rep(0,w)
    d = n/w
    
    breakpoints <- seq(0,n, d)
    
    for ( i in 1:w) {
        init <- breakpoints[i] +1
        end <- breakpoints[i+1]
        frac_first <-  ceiling(init) - init
        frac_end <- end - floor(end)    
        
        interv = floor(init):ceiling(end)
        sec <- x[interv]
        if( frac_first > 0) {
            sec[1] = sec[1]*frac_first
        }
        if (frac_end > 0) {
            sec[length(sec)] = sec[length(sec)]*frac_end
        }
        PAA[i] = sum(sec)/d
    }
    PAA
}

#creates the SAX symbolic representation
#x is the series, usually in PAA reduced form
#alpha the amount of symbols desired
convert.to.SAX.symbol <- function( x, alpha) {
    symb <- SAX.breakpoints.table(alpha) #divide the gaussian in n parts
    saxstring <- NULL
    for (s in x) {
        saxstring <- c(saxstring, sum(symb < s)) #basically in which part of the gaussian
    }            
    saxstring 
}

#distance between SAX representations
#this distance again uses breakpoints of dividing the gaussian in n parts
#whith neighbouring parts having a distance of 0
MINDIST.SAX <- function(x, y, alpha, n) {
    w <- length(x)
    symb <- SAX.breakpoints.table(alpha)
    d <- 0
    for (i in 1:w) {
        xi <- x[i]
        yi <- y[i]
        if (abs(xi-yi)> 1) {
          d <- d + (symb[max(xi,yi)] - symb[ min(xi,yi) +1])^2
        }
    }
    sqrt( (n/w)*d )
}

to.char.representation <- function(x) {
    char_table <- "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    char_table <- unlist( strsplit(char_table, "") )
    if (min(x) < 1 || max(x) > length(char_table)) {
        stop("Error: cannot translate symbols to characters, alphabet size not supported for plotting")
    }
    char_table[x]
}

diss.MINDIST.SAX <- function(x, y, w, alpha=4, plot=FALSE) {
    .ts.sanity.check(x, y)
    n <- length(x)
    x <- (x - mean(x)) / sd(x) #series z-normalization
    y <- (y - mean(y)) / sd(y)
    PAAx <- PAA(x, w) #PAA dimensionality reduction
    PAAy <- PAA(y, w)
    SAXx <- convert.to.SAX.symbol(PAAx, alpha) #discretization
    SAXy <- convert.to.SAX.symbol(PAAy, alpha)
    dSAX <- MINDIST.SAX(SAXx, SAXy, alpha, n)
    
    if (plot) {
        def.par <- par(no.readonly = TRUE) # save default, for resetting...
        layout( matrix(c(1,2), nrow=1 ) ,c(1,5) )
        
        mai <- par("mai") #change margins to get plots together
        mai1 <- mai #margins for plot1 (gaussian)
        mai1[4] <- 0 #remove right margin
        mai2 <- mai #margins for plot2 (series)
        mai2[2] <- 0 #remove left margin
        
        r <- range(c(PAAx, PAAy))  
        plottitle <- paste("MINDIST.SAX DISSIMILARITY IS:", dSAX)
        
        par(mai=mai1)
        b <- seq(min(r), max(r), 0.01)
        plot(dnorm(b), b, type="l", xlab="", ylab="", ann=FALSE, xaxt="n", bty="n", ylim=r)
        abline(h=SAX.breakpoints.table(alpha), col=rainbow(alpha), lty=2) #horizontal lines are gaussian breakpoints
        
        n <- length(x)
        k <- 1
        linetype <- "b"
        if ( n %% w == 0) {
            k <- n/w
            PAAx <- rep(PAAx,each=k)
            PAAy <- rep(PAAy,each=k)
            linetype <- "l"
        } else {
            warning("Series length is not multiple of w (number of frames), plotting the dimensionality reduced series")
        }
        
        par(mai=mai2)
             
        plot(PAAx, type=linetype, pch=to.char.representation(SAXx), ylim=r, col="red",
             main=plottitle, xlab="", yaxt="n", bty="n", ylab="")
        if (linetype=="l") {
            points(seq(1,n,k) +k/2, PAAx[seq(1,n,k)], pch=to.char.representation(SAXx), col="red")
        }
        
        lines(PAAy, type=linetype, pch=to.char.representation(SAXy), col="blue")
        abline(h=SAX.breakpoints.table(alpha), col=rainbow(alpha), lty=2) #horizontal lines are gaussian breakpoints
        if (linetype=="l") {
            points(seq(1,n,k) +k/2, PAAy[seq(1,n,k)], pch=to.char.representation(SAXy), col="blue")
        }
        legend("topright", pch=16, legend=c("x","y"), col=c("red", "blue"))
        par(def.par)  #- reset to default
    }
    
    dSAX
}


SAX.plot <- function(series, w, alpha, col.ser = rainbow(ncol(series))) {
    x <- series
    stopifnot(is.ts(x) || is.mts(x))
    n = length(x)
    x <- as.matrix(x)
    #calc the SAX representation of the series
    PAAser <- apply( x, 2, function(ser) {
        ser <- (ser - mean(ser)) / sd(ser) #series z-normalization
        PAAser <- PAA(ser, w) #PAA dimensionality reduction
        })
    
    SAXser <- apply( PAAser, 2, function(ser) {
        convert.to.SAX.symbol(ser, alpha) #discretization
        })

           
    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    layout( matrix(c(1,2), nrow=1 ) ,c(1,5) )
        
    mai <- par("mai") #change margins to get plots together
    mai1 <- mai #margins for plot1 (gaussian)
    mai1[4] <- 0 #remove right margin
    mai2 <- mai #margins for plot2 (series)
    mai2[2] <- 0 #remove left margin
        
    r <- range(SAX.breakpoints.table(alpha)[c(-1, -(alpha+1))], PAAser)
    par(mai=mai1)
    b <- seq(min(r), max(r), 0.01)
    plot(dnorm(b), b, type="l", xlab="", ylab="", ann=FALSE, xaxt="n", bty="n", ylim=r)
    abline(h=SAX.breakpoints.table(alpha), col=rainbow(alpha), lty=2) #horizontal lines are gaussian breakpoints
        
    n <- length(x)
    k <- 1
    linetype <- "b"
    if ( n %% w == 0) {
        k <- n/w
        PAAser <- matrix( rep(PAAser,each=k), byrow=F, ncol = ncol(x) )
        linetype <- "l"
    } else {
        warning("Series length is not multiple of w (number of frames), plotting the dimensionality reduced series")
    }
        
    par(mai=mai2)
        
    plot(PAAser[,1], type=linetype, pch=to.char.representation(SAXser[,1]), ylim=r, col=col.ser[1], xlab="", yaxt="n", bty="n")
    abline(h=SAX.breakpoints.table(alpha), col=rainbow(alpha), lty=2) #horizontal lines are gaussian breakpoints

    if (linetype=="l") {
        points(seq(1,n,k) +k/2, PAAser[seq(1,n,k),1], pch=to.char.representation(SAXser[,1]), col=col.ser[1])
    }
    
    if (ncol(x) > 1) {
        for (i in 2:ncol(x)) {
            lines(PAAser[,i], type=linetype, pch=to.char.representation(SAXser[,i]), ylim=r, col=col.ser[i])
            points(seq(1,n,k) +k/2, PAAser[seq(1,n,k), i], pch=to.char.representation(SAXser[,i]), col=col.ser[i])
        }
    }
    legnames <- colnames(x)
    if (is.null(legnames)) {
        legnames = 1:ncol(x)
    }
    legend("topright", pch=16, legend=legnames, col=col.ser)
    par(def.par)  #- reset to default
}
