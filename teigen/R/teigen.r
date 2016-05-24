## This file is part of teigen.
## Copyright (C) 2012-2015  Jeffrey L. Andrews
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA  02110-1301, USA.



##x is a matrix or data frame, wt is a vector, center is a vector
covwtC <- function(x, wt, center){
    if (is.data.frame(x)) {
        x <- as.matrix(x)
    }
    else if (!is.matrix(x)) {
        stop("'x' must be a matrix or a data frame")
    }
    if (!all(is.finite(x))) {
        stop("'x' must contain finite values only")
    }
    if(length(wt) != nrow(x)){
        stop("length of 'wt' must equal the number of rows in 'x'")
    }
    if (any(wt < 0) || sum(wt) == 0) {
        stop("weights must be non-negative and not all zero")
    }
    if (length(center) != ncol(x)) {
        stop("length of 'center' must equal the number of columns in 'x'")
    }
    rv <- .C("mycovwt", as.double(as.vector(x)), as.integer(dim(x)[1]), as.integer(dim(x)[2]), 
             as.double(wt),
             as.double(center), as.integer(length(center)), 
             double((dim(x)[2]^2)))[[7]]
    dim(rv) <- c(dim(x)[2], dim(x)[2])
    rv
}

maha <- function(x, cm, Sx){
    x <- if (is.vector(x))
        matrix(x, ncol = length(x))
    else as.matrix(x)
    
    .C("mymaha",
       as.double(as.vector(x)), as.integer(dim(x)[2]), as.integer(dim(x)[1]),
       as.double(cm),
       as.double(as.vector(Sx)), as.integer(dim(Sx)[1]),
       double(dim(x)[1]))[[7]]    
}

deltaup <-
    function(x,mug,sigma,sigmainv,G,n,univar){ 
        delta <- matrix(0,n,G)
        if(univar){
            for(g in 1:G){
                delta[,g] <- maha(x, mug[g,],solve(1/sigma[,,g])) 
            }
        }
        else{
            for(g in 1:G){
                delta[,g] <- maha(x, mug[g,], solve(sigmainv[,,g]))
            }
        }
        delta
    }


## estimateTime is run if teigen is run in verbose mode, not in parallel, and will 
## output the progress of the program. Updated progress is printed to the screen 
## after a new model has been calculated.
estimateTime <- function(stage, start.time, totmod, backspace){
    curwidth <- getOption("width")
    ##Output enough backspace characters to erase the previous output progress information
    ##\b is used instead of \r due to RStudio's difficulty with handling carriage returns
    cat(paste(rep("\b", backspace), collapse = ""))
    
    ##the string that will eventually be output to the screen
    output.string <- "" 
    
    ##[short,med,long]string all are strings that will be output depending on size of the
    ##console. These strings will be one of two things, depending on if stage is "init" or not
    if(stage=="init"){
        medstring <- longstring <- "???"
        shortstring <- 0
        modelcount <- NA 
    }
    else{
        modelcount <- stage+1 ##number of models calculated
        modsleft <- totmod-modelcount
        timerun <- difftime(Sys.time(),start.time, units="secs")
        timeremain <- (timerun/modelcount)*modsleft
        
        if(timeremain>60){
            units(timeremain) <- "mins"
        } else if(timeremain>3600){
            units(timeremain) <- "hours"
        } else if(timeremain>86400){
            units(timeremain) <- "days"
        } else if(timeremain>604800){
            units(timeremain) <- "weeks"
        }
        
        ## % complete
        shortstring <- format(round((1-modsleft/totmod)*100), width=4, justify="right")
        ##approx. time remaining
        medstring <- format(round(timeremain,1), width=6, justify="right")
        ##Time taken
        longstring <- format(round(timerun,1), width=6, justify="right")
    }
    
    ##start to build up the output string depending on console size
    if(curwidth>=16){
        output.string <- paste(shortstring, "% complete", sep="")
        
        ##if larger console, output more information
        if(curwidth>=50){ 
            output.string <- paste("Approx. remaining:", medstring, "  |  ", output.string, sep = "")
            
            ##if even larger console, output even more information
            if(curwidth>=80){
                output.string <- paste("Time taken:", longstring , "  |  ", output.string, sep = "")
            }
        }
    }
    
    cat(output.string)
    flush.console()
    ## return model count and output string size
    c(modelcount, nchar(output.string)) 
    
    
}

modelgen <- function(){
    teigenModels <- list()
    teigenModels[["altnames"]] <- c("VVVV", "VVVE","EEVV", "EEVE","EVVV", "EVVE","EEEV","EEEE",
                                    "EVIV", "EVIE","EEIV", "EEIE","VIIV","VIIE","EIIV","EIIE",
                                    "VVIV","VVIE","VEEV", "VEEE","VEVV", "VEVE","VEIV", "VEIE",
                                    "VVEV","VVEE","EVEV","EVEE")
    teigenModels[["altunivariate"]] <- c("univVV", "univVE", "univEV", "univEE")
    teigenModels[["multivariate"]] <- c("UUUU", "UUUC","CUCU", "CUCC","CUUU", "CUUC","CCCU","CCCC",
                                        "CIUU", "CIUC","CICU", "CICC","UIIU","UIIC","CIIU","CIIC",
                                        "UIUU","UIUC","UCCU", "UCCC","UUCU", "UUCC","UICU", "UICC",
                                        "UCUU","UCUC","CCUU","CCUC")
    teigenModels[["univariate"]] <- c("univUU", "univUC", "univCU", "univCC")
    teigenModels[["dfconstrained"]] <- c("UUUC","CUCC","CUUC","CCCC",
                                         "CIUC", "CICC","UIIC","CIIC",
                                         "UIUC","UCCC","UUCC", "UICC",
                                         "UCUC","CCUC")
    teigenModels[["altdfconstrained"]] <- c("VVVE", "EEVE", "EVVE","EEEE",
                                            "EVIE","EEIE","VIIE","EIIE",
                                            "VVIE", "VEEE", "VEVE", "VEIE","VVEE","EVEE")
    teigenModels[["dfunconstrained"]] <- c("UUUU","CUCU","CUUU","CCCU",
                                           "CIUU", "CICU","UIIU","CIIU",
                                           "UIUU","UCCU","UUCU", "UICU",
                                           "UCUU","CCUU")
    teigenModels[["altdfunconstrained"]] <- c("VVVV", "EEVV", "EVVV","EEEV",
                                              "EVIV","EEIV","VIIV","EIIV",
                                              "VVIV", "VEEV", "VEVV", "VEIV","VVEV","EVEV")
    teigenModels
}

## Display uncertainty plot or marginal contour plot, specified by the "what" argument.
## If what is both plots (default) then user has option to switch between displaying one or the other.
plot.teigen <- function(x, xmarg = 1, ymarg = 2, res = 200, what = c("contour", "uncertainty"), alpha = 0.4,
                        col = rainbow(x$G), pch = 21, cex = NULL, bg = NULL, lty = 1, uncmult = 0,
                        levels=c(seq(0.01, 1, by = 0.025), 0.001), main=NULL, xlab=NULL, draw.legend = TRUE, ...){
    
    teigen <- x
    G <- teigen$G ## number of groups
    margs <- par()$mar
    
    if(length(col) < G)
        stop("Not enough colors provided for all the groups (", G, " groups). Please provide more colors.")
    
    basecols <- col ##to be potentially used for legend later
    ## bg will only be used if pch is one of 21:25
    basebg <- bg ##to be potentially used for legend later
    
    bg <- bg[teigen$clas] ## will still be NULL if was originally NULL
    groups <- paste("Group", 1:G)
        
    if((ncol(teigen$x) > 1) && (!is.null(ymarg))){
        cols <- col[teigen$clas] ## colors corresponding to data
        
        if(pch %in% 21:25){
            if(is.null(bg)){ ##if user did not provide a bg argument
                bg <- adjustcolor(cols, alpha.f = alpha)
            }
            else{
                bg <- adjustcolor(bg, alpha.f = alpha)
            }
            col <- adjustcolor("black", alpha.f = alpha)
        }
        else{
            col <- adjustcolor(cols, alpha.f = alpha)
        }
        
        ds <- ceiling(par()$pin[2]/(par()$cin[2] + 0.1))
        columns <- ceiling(G/ds)
        rmargin <- 8 + (4.5*(columns-1))
        
        mx <- max(teigen$x[,c(1)]) ## used in placement of legend
        my <- max(teigen$x[,c(2)])
        
        myuncertainty <- function(teigen, i = 10,j = 1.1,  pch, cex, bg, col, alpha = 0.4, uncmult,
                                  xmarg = 1, ymarg = 2, res=200, ...) {
            
            mycex <- (i-3)*(j-apply(teigen$fuzzy, 1, max))
            if(!is.null(cex)){ ##if user provided a cex argument (minimum point size)
                dif <- cex - min(mycex) ##amount to change mycex by
                mycex <- mycex + dif
            }
            mycex <- (1 + uncmult/20)*mycex^(1 + uncmult/20)
            
            plot(teigen$x[,1:2], col= col, bg = bg, pch = pch, cex = mycex, ...)
            
            if(is.null(main)){title("Uncertainty Plot", ...)}else{title(main, ...)}
        }
        
        mycontour <- function(teigen, i = 10, j = 1.1, pch, bg, col, alpha = 0.4, 
                              xmarg = 1, ymarg = 2, res=200, levels=c(seq(0.01, 1, by = 0.025), 0.001), ...) { 
                        
            plot(teigen$x[,c(xmarg,ymarg)], col=col, pch=pch, bg=bg, main="", xlab="",...)
            if(is.null(main)){title("Marginal Contour Plot", ...)}else{title(main, ...)}
            if(is.null(xlab)){title(xlab=colnames(teigen$x)[xmarg], ...)}else{title(xlab=xlab, ...)}
          
            lims <- par()$usr
            xseq <- seq(lims[1], lims[2], length.out=res)
            yseq <- seq(lims[3], lims[4], length.out=res)
            seqmat <- matrix(NA, res^2, 2)
            seqmat[,1] <- rep(xseq, each=res)
            seqmat[,2] <- rep(yseq, res)
            val <- matrix(NA,res,res)
            if(!teigen$info$univar){
                sigmainv <- array(NA, dim=c(2,2,teigen$G))
                for(g in 1:teigen$G){
                    sigmainv[,,g] <- solve(teigen$par$sigma[c(xmarg,ymarg),c(xmarg,ymarg),g])
                }
            }
            delt <- deltaup(seqmat, matrix(teigen$par$mean[,c(xmarg,ymarg)],nrow=teigen$G,ncol=2), teigen$par$sigma[c(xmarg,ymarg),c(xmarg,ymarg),],sigmainv, teigen$G, res^2, teigen$info$univar)
            dens <- rowSums(exp(tft(seqmat,teigen$G,colSums(teigen$fuzzy)/nrow(teigen$x),teigen$par$df,2,matrix(teigen$par$mean[,c(xmarg,ymarg)],nrow=teigen$G,ncol=2),sigmainv,res^2,array(teigen$par$sigma[c(xmarg,ymarg),c(xmarg,ymarg),], dim=c(2,2,teigen$G)),teigen$info$univar,delt,teigen$info$gauss)))
            val <- matrix(dens, res, res, byrow=TRUE)
            contour(x=xseq, y=yseq, z=val, add=TRUE, levels=levels, col=rgb(0.5,0.5,0.5,alpha=0.7))
            
        }
        
        if (interactive() || length(what) == 1) {
            title <- "Choose which graph you would like to view. Enter 0 to exit."
            choice <- if(length(what) != 1)
                choice <- menu(what, graphics = FALSE, title = title)
            else
                choice <- 1
            
            while (choice != 0) {
                
                margs <- par()$mar ## make space on right margin for legend
                if(draw.legend){
                    par(mar = c(margs[1:3], rmargin))
                }
                
                if (what[choice] == "contour") {
                    mycontour(teigen, xmarg = xmarg, ymarg = ymarg, res = res, levels = levels, 
                              col = col, pch = pch, cex = cex, bg = bg, alpha = alpha, ...) 
                }
                else if (what[choice] == "uncertainty") {
                    myuncertainty(teigen, xmarg = xmarg, ymarg = ymarg, res = res, 
                                  uncmult = uncmult, col = col, bg = bg, cex = cex, pch = pch, alpha = alpha, ...)
                }
                else
                    stop("Unknown type of graph. Please specify one (or both) of 'uncertainty' or 'contour'. ")
                
#                 if(!identical(main,"")){
#                     text <- paste("Graphing teigen object: ' ", deparse(substitute(x)), " '")
#                     mtext(text, 3, font = 3)
#                 }
                
                if(draw.legend){
                                        
                    ptcols <- if((pch %in% 21:25) && !(is.null(basebg)))
                        basebg
                    else
                        basecols
                    
                    legend(groups, pch = 21, pt.bg=adjustcolor(ptcols, alpha.f = alpha), x = mx + 0.25, y = my + 0.25, 
                                      bty = "n", x.intersp = 1, y.intersp = 1.25, ncol = columns, xpd = NA)
                }
                
                choice <- if(length(what) != 1)
                    choice <- menu(what, graphics = FALSE, title = title)
                else
                    choice <- 0
                
                if (choice != 0){
                    dev.off() ## without this, graphics display incorrectly on Windows
                }
            }
        }
    }
    else{
        
        ds <- ceiling(par()$pin[2]/(par()$cin[2] + 0.1))
        columns <- ceiling(G/ds)
        rmargin <- 8 + (6.2*(columns-1))
        
        objname <- deparse(substitute(x))
        if((ncol(teigen$x) > 1) && (is.null(ymarg))){ ##if this is not conventional univariate data, modify the teigen object
            ##to interpret it as if it were univariate
            if(xmarg > ncol(teigen$x))
                stop("Invalid xmarg value for univariate plot: xmarg too large.")
            
            colname <- colnames(teigen$x)[xmarg] ##get name of the column
            teigen$x <- as.matrix(teigen$x[,xmarg], ncol = 1) ##adjust data set
            colnames(teigen$x) <- colname
            
            teigen$par$mean <- matrix(teigen$par$mean[,xmarg], ncol = 1) ##adjust mean matrix
            a <- array(dim = c(1,1,G)) 
            for(i in 1:G){ ##adjust sigma matrix 
                a[1,1,i] <- teigen$par$sigma[xmarg,xmarg,i]
            }
            teigen$par$sigma <- a
        }
                            
        dunivt <- function(xdum,df,sig,mean,pig,gauss){ ##group curves
            if(!gauss){
                exp(log(pig)+lgamma((df+1)/2)-(1/2)*log(sig)-((1/2)*(log(pi)+log(df))+lgamma(df/2)+((df+1)/2)*(log(1+ mahalanobis(matrix(xdum,ncol=1), mean, 1/sig, inverted=TRUE)/df))))
            }
            else{
                dnorm(xdum, mean=mean, sd=sqrt(sig))
            }
        }
        mixuniv <- function(x, teigen){ ## mixture curve
            summat <- rep(0, length(x))
            for(g in 1:G){
                summat <- summat + dunivt(x, teigen$par$df[g], teigen$par$sigma[,,g], teigen$par$mean[g,], teigen$par$pig[g], teigen$info$gauss)
            }
            summat
        }
        bigdens <- NA ##vector to store the estimated max values of the curves to be plotted (used for plot dimensions)
        
        for(g in 1:G){ ##find max for each group's curve
            bigdens[g] <- dunivt(teigen$par$mean[g,],teigen$par$df[g],teigen$par$sigma[,,g],teigen$par$mean[g,],teigen$par$pig[g], teigen$info$gauss)
        }
        
        teig.dens <- density(teigen$x)
        ##mixcurve curve has interval c(from, to) multiplied by 1.1 to account for slight increase
        ##in the x dimensions of the plot that par()$usr adds. Do not plot it yet
        mixcurve <- curve(mixuniv(x, teigen), from = 1.1*min(teig.dens$x), to = 1.1*max(teig.dens$x),
                          n=res, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
        
        par(mfg = par()$mfg[1:2])
        margs <- par()$mar ## make space on right margin for legend
        if(draw.legend){
            par(mar = c(margs[1:3], rmargin))
        }
        
        xlim <- c(min(teig.dens$x), max(teig.dens$x)) ## use these numbers for base plot dimensions
        
        bigdens[G + 1] <- max(mixcurve$y)
        bigdens[G + 2] <- max(teig.dens$y)
        
        ylim <- c(0,max(bigdens)+0.1*max(bigdens)) ## find max y limit
        
        plot(teig.dens, ylim=ylim, lty=4, main="", xlab="", ...) ##plot first curve, density curve
        if(is.null(main)){title("Univariate Density Plot", ...)}else{title(main, ...)}
        if(is.null(xlab)){title(xlab=colnames(teigen$x)[xmarg], ...)}else{title(xlab=xlab, ...)}
        
        ##now plot the remaining curves
        lines(mixcurve$x, mixcurve$y, col="black", ...)
        for(g in 1:G){
            curve(dunivt(x,teigen$par$df[g],teigen$par$sigma[,,g],teigen$par$mean[g,],teigen$par$pig[g], teigen$info$gauss),add=TRUE,
                  col=col[g], n=res, lty = lty, ...)
        }
        
#         if(!identical(main,"")){
#             text <- paste("Graphing teigen object: ' ", objname, " '")
#             mtext(text, 3, font = 3)
#         }
        
        if(draw.legend){
            mx <- par()$usr[2]
            my <- par()$usr[4]

            legendvec <- c("density()", "mixture", paste("Group",1:G))
                        
            legend(legendvec, x = mx, y = my, col=c("black", "black", col),lty=c(4,1,rep(lty,G)),
                   lwd=c(1,2,rep(1,G)), bty = "n", x.intersp = 1, y.intersp = 1.25, 
                   ncol = columns, xpd = NA)
        }
    }
    par(mar = margs)
}

print.teigen <- function(x, ...){
    if(x$G==x$iclresults$G & x$modelname==x$iclresults$modelname){
        s <- strsplit(x$bestmodel, split = "\\)")
        cat("BIC and ICL select the same model and groups.\n",
            paste(s[[1]][1], ", ICL of ", round(x$iclresults$icl, 4), ")", s[[1]][2], sep = ""), "\n", sep = "")
    }
    else
        cat(cat(x$bestmodel), "\n", x$iclresults$bestmodel, "\n", sep = "")
}

summary.teigen <- function(object, ...){
    x <- object
    mysummary <- list()
    mysummary$bicloglik <- x$logl
    mysummary$bic       <- x$bic
    mysummary$icl       <- x$iclresults$icl
    mysummary$bicmodel  <- x$modelname
    mysummary$bicgroups <- x$G
    mysummary$bicconv <- x$parameters$conv
    mysummary$iclconv <- x$iclresults$parameters$conv
    mysummary$class <- x$classification
    mysummary$disagree <- FALSE
    
    if(!(x$G==x$iclresults$G & x$modelname==x$iclresults$modelname)){
        mysummary$disagree <- TRUE
        mysummary$iclloglik <- x$iclresults$logl
        mysummary$iclmodel  <- x$iclresults$modelname
        mysummary$iclgroups <- x$iclresults$G
        mysummary$iclclass <- x$iclresults$classification
    }
    
    class(mysummary) <- "summary.teigen"
    return(mysummary)
}

print.summary.teigen <- function(x, ...){
    teigsummary <- x
    convb <- teigsummary$bicconv
    convi <- teigsummary$iclconv
    
    cat("------------- Summary for teigen -------------", "\n")
    
    if(!(is.null(convb) || is.null(convi))){
        if(!convb && !convi){
            cat("BIC and ICL results both obtained from a model\nthat did not converge in the imposed maximum\n         iterations of the algorithm.\n\n")
        }
        else if(!convb){
            cat("  BIC results obtained from a model that did\nnot converge in the imposed maximum iterations\n               of the algorithm.\n\n")
        }
        else if(!convi){
            cat("  ICL results obtained from a model that did\nnot converge in the imposed maximum iterations\n               of the algorithm.\n\n")
        }
    }
    
    if(!teigsummary$disagree){
        cat("            ------   RESULTS   ------    ", "\n")
        cat("            Loglik:        ", round(teigsummary$biclogl, 4), "\n")
        cat("            BIC:           ", round(teigsummary$bic, 4), "\n")
        cat("            ICL:           ", round(teigsummary$icl, 4), "\n")
        cat("            Model:         ", teigsummary$bicmodel, "\n")
        cat("            # Groups:      ", teigsummary$bicgroups, "\n\n\n")
        cat("Clustering Table:  ", "\n")
        print(table(teigsummary$class))
        
    }
    else{
        cat("         ---- BIC RESULTS ----    ", "\n")
        cat("         Loglik:     ", round(teigsummary$biclogl, 4), "\n")
        cat("         BIC:        ", round(teigsummary$bic, 4), "\n")
        cat("         Model:      ", teigsummary$bicmodel, "\n")
        cat("         # Groups:   ", teigsummary$bicgroups, "\n\n\n")
        cat("Clustering Table:  ", "\n")
        print(table(teigsummary$class))
        cat("\n\n\n\n")
        
        cat("         ---- ICL RESULTS ----    ", "\n")
        cat("         Loglik:     ", round(teigsummary$icllogl, 4), "\n")
        cat("         ICL:        ", round(teigsummary$icl, 4), "\n")
        cat("         Model:      ", teigsummary$iclmodel, "\n")
        cat("         # Groups:   ", teigsummary$iclgroups, "\n\n\n")
        cat("Clustering Table:  ", "\n")
        print(table(teigsummary$iclclass))
    }
    
    invisible(teigsummary)
}

tBICcalc <-
    function(conv,G,p,mod,logl,n,gauss,univar,submod13){
        if(conv==0){
            BIC <- -Inf
        }
        if(conv %in% 1:2){
            freepar <- G-1 + G*p
            if(univar){
                if(mod=="univUU"){
                    freepar <- freepar + G + G
                }
                if(mod=="univUC"){
                    freepar <- freepar + G + 1
                }
                if(mod=="univCU"){
                    freepar <- freepar + 1 + G
                }
                if(mod=="univCC"){
                    freepar <- freepar + 1 + 1
                }
            }
            if(submod13=="UUU"){
                freepar <- freepar + G*p*(p+1)/2
            }
            if(submod13=="CCC"){
                freepar <- freepar + p*(p+1)/2
            }
            if(submod13=="CUC"){
                freepar <- freepar + G*p*(p+1)/2 - (G-1)*p
            }
            if(submod13=="CUU"){
                freepar <- freepar + G*p*(p+1)/2 - (G-1)
            }
            if(submod13=="CIC"){
                freepar <- freepar + p
            }
            if(submod13=="CIU"){
                freepar <- freepar + G*p - (G-1)
            }
            if(submod13=="UIU"){
                freepar <- freepar + G*p
            }
            if(submod13=="UII"){
                freepar <- freepar + G
            }
            if(submod13=="CII"){
                freepar <- freepar + 1
            }
            ##ITERATIVE PROCEDURE MODELS
            if(submod13=="UCC"){
                freepar <- freepar + p*(p+1)/2 + G-1
            }
            if(submod13=="UUC"){
                freepar <- freepar + G*p*(p+1)/2 - (G - 1)*(p-1)
            }
            if(submod13=="CCU"){
                freepar <- freepar + p*(p+1)/2 + (G - 1)*(p-1)
            }
            if(submod13=="UCU"){
                freepar <- freepar + p*(p+1)/2 + (G - 1)*p
            }
            if(submod13=="UIC"){
                freepar <- freepar + p + G-1
            }
            if(substring(mod,4,4)=="U"){
                freepar <- freepar + G
            }
            if(substring(mod,4,4)=="C"){
                freepar <- freepar + 1
            }
            if(gauss){freepar <- freepar-1}
            BIC <- 2 * max(logl) - freepar*log(n)
        }
        BIC
    }

tICLcalc <-
    function(conv,n,zmat,bic,modnum,G){
        if(conv==0){
            ICL <- -Inf
        }
        if(conv %in% 1:2){
            ent <- apply(zmat,1,max)
            ICL <- bic[modnum,G] + 2*sum(log(ent))
        }
        ICL
    }

taginit <-
    function(p,G,sg,mod,sgc){
        ag <- array(0, dim=c(p, p, G))
        if(substring(mod,3,3)=="U"){
            for(g in 1:G){
                diag(ag[,,g]) <- (eigen(sg[,,g])$values)/(det(sg[,,g])^(1/p))
            }
        }
        if(substring(mod,3,3)=="C"){
            for(g in 1:G){
                diag(ag[,,g]) <- (eigen(sgc)$values)/(det(sgc)^(1/p))
            }
        }
        ag
    }

tagupdate <-
    function(p,G,mod,sg,lambdag,ng,n,dg,submod13){
        ag <- array(0, dim=c(p, p, G))
        negdet <- FALSE
        if(submod13=="CUC"){
            omegag <- matrix(0,p,p)
            for(g in 1:G){
                diag(omegag) <- diag(omegag) + eigen(ng[g]*sg[,,g])$values
            }
            ag[,,] <- omegag/(det(omegag)^(1/p))
        }
        if(submod13=="CUU"){
            for(g in 1:G){
                diag(ag[,,g]) <- eigen((ng[g]*sg[,,g])/(det(ng[g]*sg[,,g])^(1/p)))$values
            }
        }
        if(submod13=="CIC"){
            for(g in 1:G){
                diag(ag[,,g]) <- diag(n*sg[,,g])/(prod(diag(n*sg[,,g]))^(1/p))
            }
        }
        if(any(submod13==c("CIU","UIU"))){
            for(g in 1:G){
                diag(ag[,,g]) <- diag(ng[g]*sg[,,g])/(prod(diag(ng[g]*sg[,,g]))^(1/p))
            }
        }
        if(substring(mod,3,3)=="I"){
            ag[,,] <- diag(p)
        }
        if(submod13=="CCU"){
            for(g in 1:G){
                dwd <- diag(t(dg[,,g]) %*% (ng[g]*sg[,,g]) %*% dg[,,g])
                diag(ag[,,g]) <- dwd/(prod(dwd)^(1/p))
            }
        }
        if(submod13=="UCU"){
            for(g in 1:G){
                dwd <- diag(t(dg[,,g]) %*% (ng[g]*sg[,,g]) %*% dg[,,g])
                diag(ag[,,g]) <- dwd/(prod(dwd)^(1/p))
            }
        }
        ##IP
        if(any(submod13==c("UCC","UUC"))){
            dum <- matrix(0,p,p)
            for(g in 1:G){
                dum <- dum + ng[g]*sg[,,g]/lambdag[g]
            }
            detdum <- det(dum)
            if(is.finite(detdum)){
                if(detdum<=0){
                    negdet <- TRUE
                }
                else{
                    dum2 <- dum/(detdum^(1/p))
                    dum3 <- eigen(dum2)$values
                    for(g in 1:G){
                        diag(ag[,,g]) <- dum3
                    }
                }
            }
            else{
                negdet <- TRUE
            }
        }
        if(submod13=="UIC"){
            dum <- matrix(0,p,p)
            for(g in 1:G){
                dum <- dum + ng[g]*sg[,,g]/lambdag[g]
            }
            dum2 <- diag(dum)/(prod(diag(dum))^(1/p))
            for(g in 1:G){
                diag(ag[,,g]) <- dum2
            }
        }
        for(g in 1:G){
            if(any(diag(ag[,,g]) <= sqrt(.Machine$double.eps))){
                negdet <- TRUE
            }
        }
        if(negdet){
            ag <- FALSE
        }
        ag
    }

tsoftrandz <-
    function(n,G,clas,kno,known){
        zmat <- matrix(runif(n*G,0,1), n, G)
        zmat <- zmat/rowSums(zmat)
        
        if(clas>0){
            known.used <- as.numeric(known)[kno==1] #the elements of the known vector that we are using for classification
            zreplace <- matrix(0, nrow = length(known.used), G)
            for(i in 1:G){
                zreplace[known.used==i, i] <- 1
            }
            zmat[kno==1, ] <- zreplace
        }
        zmat
    }

tdginit <-
    function(p,G,sg,mod,sgc){
        dg <- array(0, dim=c(p, p, G))
        if(substring(mod,3,3)=="U"){
            for(g in 1:G){
                dg[,,g] <- eigen(sg[,,g])$vectors
            }
        }
        if(substring(mod,3,3)=="C"){
            dg[,,] <- eigen(sgc)$vectors
        }
        dg
    }

tdgupdate <-
    function(p,G,mod,sg,lambdag,ng,ag,submod13, dg, flipswitch=NULL){
        dg <- array(0, dim=c(p, p, G))
        negdet <- FALSE
        if(any(submod13==c("CUC","UUC"))){
            for(g in 1:G){
                dg[,,g] <- eigen(ng[g]*sg[,,g])$vectors
            }
        }
        else{
            if(submod13=="CUU"){
                for(g in 1:G){
                    detngsg <- det(ng[g]*sg[,,g])
                    if(detngsg<=0){
                        negdet <- TRUE
                        break
                    }
                    dg[,,g] <- eigen((ng[g]*sg[,,g])/(detngsg^(1/p)))$vectors
                }
            }
            else{
                if(submod13=="UCC"){
                    dum <- matrix(0,p,p)
                    for(g in 1:G){
                        dum <- dum + ng[g]*sg[,,g]/lambdag[g]
                    }
                    detdum <- det(dum)
                    if(is.finite(detdum)){
                        if(detdum<=0){
                            negdet <- TRUE
                        }
                        else{
                            dum2 <- dum/(detdum^(1/p))
                            dum3 <- eigen(dum2)$vectors
                            dg[,,] <- dum3
                        }
                    }
                    else{
                        negdet <- TRUE
                    }
                }
                else{
                    if(any(submod13==c("UCU","CCU"))){
                        if(flipswitch){
                            ##MM1
                            dum <- matrix(0,p,p)
                            ainv <- dum
                            for(g in 1:G){
                                wg <- ng[g]*sg[,,g]
                                wk <- eigen(wg, symmetric=TRUE, only.values=TRUE)$values[1]
                                diag(ainv) <- diag(1/(lambdag[g]*ag[,,g]))
                                dum <- dum + ainv %*% t(dg[,,g]) %*% wg - wk*ainv %*% t(dg[,,g])
                            }
                            ##print(dum)
                            if(all(is.finite(dum))){
                                svdum <- svd(dum)
                                dg[,,] <- svdum$v %*% t(svdum$u)
                            }
                            else{
                                negdet <- TRUE
                            }
                        }
                        else{
                            ##MM2
                            dum <- matrix(0,p,p)
                            ainv <- dum
                            for(g in 1:G){
                                wg <- ng[g]*sg[,,g]
                                diag(ainv) <- diag(1/(lambdag[g]*ag[,,g]))
                                ##dum <- dum + wg %*% dg[,,g] %*% ainv  - (lambdag[g]*ag[p,p,g])*wg %*% dg[,,g]
                                dum <- dum + wg %*% dg[,,g] %*% ainv  - (max(ainv))*wg %*% dg[,,g]
                            }
                            svdum <- svd(dum)
                            dg[,,] <- svdum$v %*% t(svdum$u)
                            
                        }
                    }
                    else{
                        dg[,,] <- diag(p)
                    }
                }
            }
        }
        if(negdet){
            dg <- FALSE
        }
        dg
    }

thardrandz <-
    function(n,G, clas,kno,known){
        known.used <- as.numeric(known)[kno==1] #the elements of the known vector that we are using for classification
        zmat <- matrix(0, n, G)
        zreplace <- matrix(0, nrow = length(known.used), G)
        
        dum <- sample(1:G,n,replace=TRUE)
        for(i in 1:G){
            zmat[dum==i,i] <- 1
            zreplace[known.used==i, i] <- 1
        }
        if(clas>0){
            zmat[kno==1, ] <- zreplace
        }
        zmat
    }

teigen <- 
    function(x, Gs=1:9, models="all", init="kmeans", scale=TRUE, dfstart=50, known=NULL, training=NULL, 
             gauss=FALSE, dfupdate="approx", eps=c(0.001,0.1), verbose=TRUE, maxit=c(Inf,Inf), 
             convstyle = "aitkens", parallel.cores=FALSE, ememargs=NULL){
        fcall <- match.call() ##get the call of this function
        
        ##error checking
        if(class(init)!="list" && !(init %in% c("kmeans", "hard", "soft", "uniform", "emem"))){
            stop("'init' must be one of 'kmeans', 'hard', 'soft', 'uniform', 'emem' or a list. See ?teigen.")
        }
        if(!(dfupdate %in% c("approx", "numeric", FALSE, TRUE))){
            stop("'dfupdate' must be either 'approx', 'numeric', or FALSE. See ?teigen.")
        }
        
        origx <- x
        if(is.vector(x)){
          x <- matrix(origx,length(origx),1)
        }
        if(scale){
          x <- scale(x, center = TRUE, scale = TRUE)
        }        
        n <- nrow(x)
        if(!is.null(known)){
          origknown <- known
          known <- factor(known)
          fcall$"origknown" <- origknown
        }
        
        if(is.null(known)){
          clas <- 0
          if(!is.null(training)){
            stop("Known classifications vector not given, or not the same length as the number of samples (see help file)")
          }
          kno <- NULL
          testindex <- NULL
        }
        else{
          if(length(known)!=n){
            stop("Known classifications vector not given, or not the same length as the number of samples (see help file)")
            return(NULL)
          }
          if(is.null(training)){
            if(!anyNA(known)){
               #Gs <- length(unique(known))
               warning("No NAs in 'known' nor training vector supplied, all values have known classification (parameter estimation only)")
               testindex <- 1:n
               kno <- rep(1, n)
               unkno <- (kno-1)*(-1)
               Gs <- length(unique(known))
               clas <- length(training)/nrow(x)
               init <- list()
               init[[Gs]] <- as.numeric(known)
            }
            else{
              training <- which(!is.na(known))
              testindex <- training
              kno <- vector(mode="numeric", length=n)
              kno[testindex] <- 1
              unkno <- (kno-1)*(-1)
              Gs <- length(unique(known[training]))
              clas <- length(training)/nrow(x)
            }
          }
          else{
            if(anyNA(known)){
              warning("training vector ignored, NAs in known vector assumed unknown")
              training <- NULL
            }
            else{
              testindex <- training
              newknown <- rep(NA, n)
              origG <- length(unique(known))
              newknown[testindex] <- known[testindex]
              known <- newknown
              kno <- vector(mode="numeric", length=n)
              kno[testindex] <- 1
              unkno <- (kno-1)*(-1)
              Gs <- length(unique(known[training]))
              if(!(origG==Gs)){warning("Not all groups represented in training set, fitting for fewer groups...")}
            }
          }
          #
#           if(anyNA(known)){
#             training <- which(!is.na(known))
#             testindex <- training
#             kno <- vector(mode="numeric", length=n)
#             kno[testindex] <- 1
#             unkno <- (kno-1)*(-1)
#             Gs <- length(unique(known[training]))
#             clas <- length(training)/nrow(x)
#           }
#           else{
#               
#           }
          clas <- 1
        }
        
        
 
  ####      
#        if(clas>0){
            
#             testindex <- sample(1:n, ceiling(n*(clas/100)))
#             kno <- vector(mode="numeric", length=n)
#             kno[testindex] <- 1
#             unkno <- (kno-1)*(-1)
#             Gs <- length(unique(known))
#        }
#         else{
#           kno <- NULL
#         }
#         if(length(training)>0){
#             if(length(known)!=n){
#                 stop("Known classifications vector not given, or not the same length as the number of samples (see help file)")
#                 return(NULL)
#             }
#             
#         }
   ####
        
        zmatin <- list()
        ##build up the zmatin list of matrices based on the provided init variable
       
          for(G in sort(Gs, decreasing=TRUE)){
              if(G==1){
                  zmatin[[G]] <- matrix(1,n,1)
              }
              else{
                  if(is.character(init)){
                      if(init == "hard"){
                          zmatin[[G]] <- thardrandz(n,G, clas,kno,known)
                      }
                      if(init == "soft"){
                          zmatin[[G]] <- tsoftrandz(n,G,clas,kno,known)
                      }
                      if(init == "uniform"){
                          if(clas>0){
                              zmatin[[G]] <- tuniformz(n,G,clas,kno,known)
                          }
                          else{
                              stop("Uniform initialization not available for clustering.")
                              return(NULL)
                          }
                      }
                      if(init == "kmeans"){
                          zmatin[[G]] <- tkmeansz(x,n,G,known,kno,testindex,clas)
                      }
                      if(init == "emem"){
                        zmatin[[G]] <- ememinit(ememargs, x, G, scale, dfstart, clas, known, training, gauss, dfupdate, eps, n)
                      }
                  }
                  else{
                      zmatin[[G]] <- tgivenz(n,G,known,init[[G]],testindex,clas)
                  }
              }
          }
        
        fcall$"training" <- testindex
        fcall$"Gs" <- Gs
        fcall$"kno" <- kno
        fcall$"zmatin" <- zmatin
        fcall$"init" <- NULL
        fcall$"clas" <- clas
        
        ##if parallel.cores was provided as an argument
        if("parallel.cores" %in% names(fcall)){ 
            
            ##One of two conditions for teigen.parallel to run are:
            ##if parallel.cores is logical (TRUE or T)
            if(TRUE %in% (fcall$"parallel.cores" == c(TRUE, "T"))){
                fcall$"parallel.cores" <- NULL  ##parallel will use default cores
                ##set teigen.parallel to be the function that is run instead
                fcall[[1]] <- teigen.parallel
            }
            else if (TRUE %in% (fcall$"parallel.cores" == c(FALSE, "F"))){
                fcall$"parallel.cores" <- NULL  ##we will run EM
                fcall[[1]] <- teigen.EM    
            }
            ##else if parallel.cores is integer 
            ##here's a safer way to check if integer, since is.integer(5) is FALSE
            else if ((is.numeric(fcall$"parallel.cores")) && (fcall$"parallel.cores"%%1==0)) {
                ##set teigen.parallel to be the function that is run instead
                fcall[[1]] <- teigen.parallel
                ##change argument name to match the argument in teigen.parallel and use the value provided
                names(fcall)[match("parallel.cores", names(fcall))] <- "numcores"
            } 
            
            else{ 
                stop("Wrong format for parallel.cores. Please provide a logical or integer indicating the number of cores you want to use.")
            }
            
        }
        else{ ##else we want to run the non parallel version
            fcall[[1]] <- teigen.EM
        }
        ##run teigen
        
        eval(fcall, parent.frame())
    }

teigen.EM <-
    function(x, zmatin, Gs=1:9, models="all", scale=TRUE, dfstart=50, clas=0, known=NULL, training=NULL, gauss=FALSE, 
             dfupdate="approx", eps=c(0.001,0.1), verbose=TRUE, maxit=c(20,1000), convstyle = "aitkens", 
             kno = NULL, ememargs = NULL, origknown = NULL){
        if(verbose){
            backspace.amount <- estimateTime("init")[2]
        }
        teigenModels <- modelgen()
        alternatenames <- FALSE
        origx <- x
        if(is.vector(x)){
            x <- matrix(origx,length(origx),1)
            univar <- TRUE
        }
        else{
            if(ncol(x)<2){
                univar <- TRUE
            }
            else{
                if(length(ncol(x))<0){
                    univar <- TRUE
                }
                else{
                    univar <-  FALSE
                }
            }
        }
        if(univar){
            if(length(models)>1){
                if(gauss){
                    dfstart <- Inf
                    gauss <- TRUE
                    dfupdate <- FALSE
                    models <- c("univUC","univCC")
                }
                else{
                    if(any(!models%in%teigenModels$univariate)){
                        models <- teigenModels$univariate
                    }
                }
            }
            else{
                if(models=="gaussian"){
                    dfstart <- Inf
                    gauss <- TRUE
                    dfupdate <- FALSE
                    models <- c("univUC","univCC")
                }
                else{
                    if((!models%in%teigenModels$univariate) & gauss){
                        models <- teigenModels$univariate
                    }
                }
            }
        }
        
        known <- factor(known)
        testindex <- training
        if(nrow(x)<ncol(x)){
            warning("Dimensionality of data exceeds the number of samples, may result in model-fitting failure")
        }
        if(!all(models %in% c(teigenModels$multivariate, teigenModels$univariate, teigenModels$altnames,"all", "gaussian", "dfconstrained", "dfunconstrained", "univariate", "altall"))){
            stop("You have specified at least one unknown model abbreviation...please select a different set of models")
            return(NULL)
        }
        if(gauss){
            dfstart <- Inf
            dfupdate <- FALSE
        }
        if(scale){
            x <- scale(x, center=TRUE, scale=TRUE)
        }
        p <- ncol(x)
        n <- nrow(x)
        if(length(models)==1){
            if(models=="dfunconstrained"|models=="altdfunconstrained"){
                if(models=="altdfunconstrained"){
                    origmodels <- teigenModels$altdfunconstrained
                    alternatenames <- TRUE
                }
                models <- teigenModels$dfunconstrained
            }
            else{
                if(models=="all"){
                    if(univar){
                        models <- teigenModels$univariate
                    }
                    else{
                        models <- teigenModels$multivariate
                    }
                }
                else{
                    if(models=="gaussian"){
                        models <- teigenModels$dfconstrained
                        dfstart <- Inf
                        gauss <- TRUE
                        dfupdate <- FALSE
                    }
                    else{
                            if(models=="dfconstrained"|models=="altdfconstrained"){
                                if(models=="altdfconstrained"){
                                    origmodels <- teigenModels$altdfconstrained
                                    alternatenames <- TRUE
                                }
                                models <- teigenModels$dfconstrained
                            }
                            else{
                                if(models=="univariate"){
                                    models <- teigenModels$univariate
                                }
                                else{
                                    if(models=="altall"){
                                        if(univar){
                                            origmodels <- teigenModels$altunivariate
                                            models <- teigenModels$univariate
                                        }
                                        else{
                                            origmodels <- teigenModels$altnames
                                            models <- teigenModels$multivariate
                                        }
                                        alternatenames <- TRUE
                                    }
                                    else{
                                        if(models %in% teigenModels$altnames){
                                            origmodels <- models
                                            models <- teigenModels$multivariate[teigenModels$altnames %in% models]
                                            alternatenames <- TRUE
                                        }
                                    }
                                }
                            }
                    }
                }
            }
        }
        else{
            if(univar){
                if(any(models %in% teigenModels$altunivariate)){
                    origmodels <- models
                    models[models %in% teigenModels$altunivariate] <- teigenModels$univariate[teigenModels$altunivariate %in% models]
                    alternatenames <- TRUE
                }
            }
            else{
                if(any(models %in% teigenModels$altnames)){
                    origmodels <- models
                    models[models %in% teigenModels$altnames] <- na.omit(teigenModels$multivariate[match(models, teigenModels$altnames)])
                    alternatenames <- TRUE
                }
            }
        }
        hh8 <- dfupdate
        if(clas>0){
#             testindex <- sample(1:n, ceiling(n*(clas/100)))
#             kno <- vector(mode="numeric", length=n)
#             kno[testindex] <- 1
            unkno <- (kno-1)*(-1)
            # Gs <- length(unique(known))
        }
#         if(length(training)>0){
# #             testindex <- training
# #             kno <- vector(mode="numeric", length=n)
# #             kno[testindex] <- 1
#             unkno <- (kno-1)*(-1)
#             #Gs <- length(unique(known[training]))
#             #clas <- length(training)/nrow(x)
#         }
        gvec <- 1:max(Gs)
        gstuff <- paste("G=",gvec,sep="")
        
        ##initialize main bic, icl, logls matrices
        bic <- icl <- logls <- matrix(-Inf, length(models), max(Gs))
        
        unc <- matrix(Inf,length(models),max(Gs))

         if(verbose){
             start.time <- Sys.time()
             totmod <- length(Gs)*length(models)
             modelcount <- 0
         }
        
        store <- list() ##what will be returned at the end of teigen.EM
        icllist <- list() 
        ##these lists will hold the best bic and icl parameters
        parameters <- list()
        parametersicl <- list()
        
        ##MAIN LOOP FOR EACH MODEL
        for(modnum in 1:length(models)){
            mod <- models[modnum]
            submod13 <- substring(mod,1,3)
            
            for(G in Gs){ ##LOOP THROUGH EACH GROUP
                singular <- 0
                killit <- FALSE
                if(G==1){
                    if(length(models)>1){
                        CCCCgroup <- c("CCCC", "CCCU", "CUCC", "CUCU", "CUUC","CUUU","UCCC","UCCU","UUCU","UUCC", "UUUC","UUUU","UCUU","UCUC","CCUU","CCUC")
                        if(any(mod==CCCCgroup)){
                            cccdum <- models[models %in% CCCCgroup]
                            if(length(cccdum)>0){
                                if(mod!=cccdum[1]){
                                    killit <- TRUE
                                }
                            }
                        }
                        CICCgroup <- c("CICC","CICU","UICC","UICU","CIUC","CIUU","UIUC","UIUU")
                        if(any(mod==CICCgroup)){
                            cicdum <- models[models %in% CICCgroup]
                            if(length(cicdum)>0){
                                if(mod!=cicdum[1]){
                                    killit <- TRUE
                                }
                            }
                        }
                        CIICgroup <- c("CIIC", "CIIU", "UIIC","UIIU")
                        if(any(mod==CIICgroup)){
                            ciidum <- models[models %in% CIICgroup]
                            if(length(ciidum)>0){
                                if(mod!=ciidum[1]){
                                    killit <- TRUE
                                }
                            }
                        }
                        univgroup <- c("univCC","univCU","univUC","univUU")
                        if(any(mod==univgroup)){
                            unidum <- models[models %in% univgroup]
                            if(length(unidum)>0){
                                if(mod!=unidum[1]){
                                    killit <- TRUE
                                }
                            }
                        }
                    }
                }
                zmat <- zmatin[[G]]
                vg <- tvginit(dfstart,G)
                ng <- tngupdate(zmat)
                if(any(ng<1.5)){killit <- TRUE}
                pig <- tpigupdate(ng,n)
                if(!killit){
                    mug <- matrix(0,G,p)
                    sg <- array(0,dim=c(p,p,G))
                    for(g in 1:G){
                        wtcov <- cov.wt(x,wt=zmat[,g],method="ML")
                        mug[g,] <- wtcov$center
                        sg[,,g] <- wtcov$cov
                    }
                    sgc <- tsginitc(G,sg,pig,p,n,x)
                    if(any(submod13==c("CCC","CIC")) | mod=="univCC" | mod=="univCU"){
                        sg[,,] <- sgc
                    }
                    if(!univar){
                        for(g in 1:G){
                            test <- rcond(sg[,,g])
                            if(test <= sqrt(.Machine$double.eps)){
                                singular <- 1
                            }
                        }
                        if(singular==1){
                            killit <- TRUE
                        }
                        if(!killit){
                            if(all(submod13!=c("CCC","UUU"))){
                                if(any(submod13==c("UCC","CCU","UCU","UUC","UIC","UCU"))){
                                    lambdag <- tlambdaginit(p,G,sg,mod)
                                    dg <- tdginit(p,G,sg,mod,sgc)
                                    ag <- taginit(p,G,sg,mod,sgc)
                                }
                                else{
                                    lambdag <- tlambdagupdate(G,mod,sg,dg,ag,p,n,ng,submod13)
                                    dg <- tdgupdate(p,G,mod,sg,lambdag,ng,ag,submod13)
                                    if(any(is.logical(dg))){
                                        killit <- TRUE
                                    }
                                    else{
                                        ag <- tagupdate(p,G,mod,sg,lambdag,ng,n,dg,submod13)
                                        if(any(is.logical(ag))){
                                            killit <- TRUE
                                        }
                                    }
                                }
                            }
                            else{
                                dg <- Inf
                                ag <- Inf
                                lambdag <- Inf
                            }
                        }
                    }
                    if(!killit){
                        sigma <- tsigmaup(p,G,sg,lambdag,dg,ag,mod,univar,submod13)
                        if(!univar){
                            for(g in 1:G){
                                test <- rcond(sigma[,,g])
                                if(test <= sqrt(.Machine$double.eps)){
                                    singular <- 1
                                }
                            }
                            if(singular==1){killit <- TRUE}
                            if(!killit){
                                sigmainv <- tsigmainvup(p,G,sigma)
                            }
                        }
                        if(!killit){
                            if(!gauss){
                                w <- twinit(x,n,G,mug,sigmainv,vg,p,sg,zmat,univar,sigma)
                            }
                            else{ w <- matrix(1,n,G) }
                        }
                    }
                }
                cycle <- 0
                dhfgs78 <- vg
                conv <- 0
                num <- ft <- matrix(0,n,G)
                logl <- NaN
                if(!killit){
                     delta <- deltaup(x,mug,sigma,sigmainv,G,n,univar)
                 }
                while(conv != 1){
                    if(killit){break}
                    ng <- tngupdate(zmat)
                    if(any(ng<1.5)){break}
                    pig <- tpigupdate(ng,n)
                    mug <- tmugupdate(G,zmat,w,x,p,univar)
                    if(hh8=="approx"){
                        testing <- try(jk861 <- yxf8(mod,dhfgs78,ng,zmat,w,G,p,n,x,mug,sigmainv),silent=TRUE)
                        if(all(is.finite(testing))){
                            dhfgs78 <- jk861
                        }
                        else{break}
                    }
                    else{
                        if(hh8=="numeric"){
                            testing <- try(jk861 <- yxf7(mod,dhfgs78,ng,zmat,w,G,p,n,x,mug,sigmainv),silent=TRUE)
                            if(all(is.finite(testing))){
                                dhfgs78 <- jk861
                            }
                            else{break}
                        }
                    }
                    if(any(is.nan(zmat))){
                        break
                    }
                    if(!gauss){
                        w <- twupdate(x,n,G,mug,sigmainv,dhfgs78,p,univar,sigma,delta)
                    }
                    ng <- tngupdate(zmat)
                    if(any(ng<1.5)){break}
                    sg <- tsgupdate(p,G,n,x,mug,zmat,w,ng,mod,pig,submod13)
                    if(!univar){
                        if(all(submod13!=c("CCC","UUU"))){
                            if(any(submod13==c("UCC","UUC","UIC","UCU","CCU"))){
                                mcyc <- 0
                                fmin <- 0
                                mtest <- Inf
                                while(mtest > eps[1] & mcyc<maxit[1]){
                                    mcyc <- mcyc + 1
                                    flipswitch <- mcyc%%2 == 1
                                    lambdag <- tlambdagupdate(G,mod,sg,dg,ag,p,n,ng,submod13)
                                    dg <- tdgupdate(p,G,mod,sg,lambdag,ng,ag,submod13,dg,flipswitch)
                                    if(any(is.logical(dg))){
                                        killit <- TRUE
                                        break
                                    }
                                    ag <- tagupdate(p,G,mod,sg,lambdag,ng,n,dg,submod13)
                                    if(any(is.logical(ag))){
                                        killit <- TRUE
                                        break
                                    }
                                    fmin[mcyc] <- tfminup(mod,G,sg,dg,ag,p,ng,lambdag,submod13)
                                    if(mcyc > 1){
                                        mtest <- fmin[mcyc-1] - fmin[mcyc]
                                    }
                                }
                                if(mcyc>=maxit[1]){
                                    warning(paste("Max M-step iteration of ", maxit[1], " met for model ", mod, " and G=", G, ", increase 'eps[1]' or 'maxit[1]'", sep=""))
                                    conv <- 2
                                }
                            }
                            else{
                                lambdag <- tlambdagupdate(G,mod,sg,dg,ag,p,n,ng,submod13)
                                dg <- tdgupdate(p,G,mod,sg,lambdag,ng,ag,submod13)
                                if(any(is.logical(dg))){
                                    break
                                }
                                ag <- tagupdate(p,G,mod,sg,lambdag,ng,n,dg,submod13)
                                if(any(is.logical(ag))){
                                    break
                                }
                            }
                            if(killit){break}
                        }
                    }
                    sigma <- tsigmaup(p,G,sg,lambdag,dg,ag,mod,univar,submod13)
                    if(!univar){
                        for(g in 1:G){
                            test <- rcond(sigma[,,g])
                            if(test <= sqrt(.Machine$double.eps)){
                                singular <- 1
                            }
                        }
                        if(singular==1){break}
                        sigmainv <- tsigmainvup(p,G,sigma)
                    }
                    delta <- deltaup(x,mug,sigma,sigmainv,G,n,univar)
                    ft <- exp(tft(x,G,pig,dhfgs78,p,mug,sigmainv,n,sigma,univar,delta,gauss))
                    ft_den <- rowSums(ft)
                    zmat <- ft/ft_den
                    if(clas>0){
                        zmat <- unkno*zmat
                        for(i in 1:n){
                            if(kno[i]==1){
                                zmat[i, known[i]] <- 1
                            }
                        }
                    }
                    if(any(is.nan(zmat))){
                        break
                    }
                    if(!gauss){
                        w <- twupdate(x,n,G,mug,sigmainv,dhfgs78,p,univar,sigma,delta)
                    }
                    ng <- tngupdate(zmat)
                    if(any(ng<1.5)){break}
                    cycle <- cycle + 1
                    logl[cycle]<- sum(log(ft_den))
                    if(is.na(logl[cycle])){break}
                    if(cycle>3){
                        #############Test for convergence##############
                        if(convstyle == "aitkens"){
                            denomi <- (logl[cycle-1]-logl[cycle-2])
                            if(denomi==0){
                                conv <- 1
                            }
                            else{
                                    ak <- (logl[cycle]-logl[cycle-1])/denomi
                                    linf <- logl[cycle-1] + (logl[cycle]-logl[cycle-1])/(1-ak)
                                    if(abs(linf-logl[cycle-1]) < eps[2]){
                                        conv<-1
                                    }
                                    if((logl[cycle]-logl[cycle-1])<0){
                                        break
                                    }
                            }
                        }
                        else if (convstyle == "lop"){
                            if((logl[cycle] - logl[cycle-1]) < eps[2]){
                                conv <- 1
                            }
                            if((logl[cycle]-logl[cycle-1])<0){
                                break
                            }
                        }
                        ################################################
                    }
                    if(cycle>(maxit[2])){
                        warning(paste("Max EM iteration of ", maxit[2], " met for model ", mod, " and G=", G, ", increase 'eps[2]' or 'maxit[2]'", sep=""))
                        conv <- 2
                        break
                    }
                }
                if(verbose){
                    rv <- estimateTime(modelcount, start.time, totmod, backspace.amount)
                    modelcount <- rv[1]
                    backspace.amount <- rv[2]
                }
                if(conv != 0){ ##if model converged (done running) then take all this info, cause we might need it later
                    maxbic <- max(bic)
                    maxicl <- max(icl)
                    
                    logls[modnum,G] <- max(logl) ##log likelihoods
                    bic[modnum,G] <- bicnum <- tBICcalc(conv,G,p,mod,logl,n,gauss,univar,submod13)
                    icl[modnum,G] <- iclnum <- tICLcalc(conv,n,zmat,bic,modnum,G)
                    
                    ##check if its the best model, then overwrite the best model object to get the new best
                    if(bicnum > maxbic){
                        store$iter <- cycle
                        parameters$df <- dhfgs78
                        parameters$mean <- mug
                        if(!univar){
                            parameters$lambda <- lambdag
                            parameters$d <- dg
                            parameters$a <- ag
                        }
                        parameters$weights <- w
                        parameters$sigma <- sigma
                        store$fuzzy <- zmat
                        parameters$pig <- pig
                        parameters$conv <- if(conv == 1)
                            TRUE
                        else
                            FALSE
                        
                        
                        
                    } 
                    if(iclnum > maxicl){
                        icllist$iter <- cycle
                        parametersicl$df <- dhfgs78
                        parametersicl$mean <- mug
                        if(!univar){
                            parametersicl$lambda <- lambdag
                            parametersicl$d <- dg
                            parametersicl$a <- ag
                        }
                        parametersicl$weights <- w
                        parametersicl$sigma <- sigma
                        icllist$fuzzy <- zmat
                        parametersicl$pig <- pig
                        parametersicl$conv <- if(conv == 1)
                            TRUE
                        else
                            FALSE
                    } 
                }
            } ##end of for(G in GS) loop
        }
        ##end of algorithm
        ##everything after here only runs ONCE
        
        ##if at least one model converged
        if(!all(is.infinite(bic))){
            if(alternatenames){
                models <- origmodels
            }
            ##gstuff is character vector with "G=1" up to "G=max(Gs)"
            dimnames(bic) <- list(models,gstuff[1:max(Gs)])
            dimnames(logls) <- list(models,gstuff[1:max(Gs)])
            dimnames(icl) <- list(models,gstuff[1:max(Gs)])
            dimnames(unc) <- list(models,gstuff[1:max(Gs)])
            unc[,1] <- Inf
            maxes <- which(bic==max(bic), arr.ind=TRUE)
            maxicl <- which(icl==max(icl), arr.ind=TRUE)
            minunc <- which(unc==min(unc), arr.ind=TRUE)
            known <- as.character(known)
            known[is.na(known)] <- "unknown"
            known <- factor(as.character(known))
            if(nrow(maxes)>1){
                warning("Maximum BIC tie between two or more models")
                bestmodnum <- maxes[1:nrow(maxes),1]
                bestmod <- models[bestmodnum]
                bestg <- maxes[1:nrow(maxes),2]
                bestzmap <- adjrand <- tab <- "MULTIPLE"
                
                parnames <- names(parameters)
                parameters <- as.list(rep("MULTIPLE", length(parameters)))
                names(parameters) <- parnames
                parameters$conv <- FALSE
                store$fuzzy <- store$iter <- "MULTIPLE"
                
            }
            if(nrow(maxes)==1){ ##if bic matrix has a highest converged number
                bestmodnum <- maxes[1]
                bestmod <- models[bestmodnum]
                bestg <- maxes[2]
                bestzmap <- apply(store$fuzzy,1,which.max)
                if(clas>0){
                    newmap <- bestzmap
                    newmap <- as.factor(newmap)
                    levels(newmap) = levels(known)[levels(known) != "unknown"]
                    newmap[testindex] <- NA
                    newknown <- known
                    newknown[testindex] <- NA
                    tab <- table(factor(known),newmap)
                }
                else{
                    if(length(known)>0){
                        bestzmap <- as.factor(bestzmap)
                        levels(bestzmap) <- levels(known)[levels(known) != "unknown"]
                        tab <- table(known,bestzmap)
                    }
                }
                if(!univar){
                    if(models[bestmodnum]%in%c("UUUU","UUUC","CCCC","CCCU")){
                        decom <- list()
                        parameters$d <- array(0, dim=c(p, p, bestg))
                        parameters$lambda <- NA
                        parameters$a  <- parameters$d
                        for(g in 1:bestg){
                            decom[[g]] <- eigen(parameters$sigma[,,g])
                            parameters$d[,,g] <- decom[[g]]$vectors
                            eigvals <- decom[[g]]$values
                            parameters$lambda[g] <- prod(eigvals)^(1/p)
                            parameters$a[,,g] <-  diag(eigvals)/parameters$lambda[g]
                        }
                    }
                }
            }
            if(nrow(maxicl)>1){
                warning("Maximum ICL tie between two or more models")
                bestmodnumicl <- maxicl[1:nrow(maxicl),1]
                bestmodicl <- models[bestmodnumicl]
                bestgicl <- maxicl[1:nrow(maxicl),2]
                bestzmapicl <- adjrandicl <- tabicl <- "MULTIPLE"
                
                parnames <- names(parametersicl)
                parametersicl <- as.list(rep("MULTIPLE", length(parametersicl)))
                names(parametersicl) <- parnames
                parametersicl$conv <- FALSE
                icllist$iter <- icllist$fuzzy <- "MULTIPLE"
            }
            if(nrow(maxicl)==1){
                bestmodnumicl <- maxicl[1]
                bestmodicl <- models[bestmodnumicl]
                bestgicl <- maxicl[2]
                bestzmapicl <- apply(icllist$fuzzy,1,which.max) 
                if(clas>0){
                    newmapicl <- bestzmapicl
                    newmapicl[testindex] <- NA
                    newmapicl <- as.factor(newmapicl)
                    levels(newmapicl) = levels(known)[levels(known) != "unknown"]
                    newknown <- known
                    newknown[testindex] <- NA
                    tabicl <- table(factor(known),newmapicl)
                }
                else{
                    if(length(known)>0){
                        bestzmapicl <- as.factor(bestzmapicl)
                        levels(bestzmapicl) <- levels(known)[levels(known) != "unknown"]
                        tabicl <- table(known,bestzmapicl)
                    }
                }
            }
            
            store <- c(store, list(parameters = parameters, allbic = bic[,Gs], bic = max(bic), 
                                   bestmodel = paste("The best model (BIC of ",round(max(bic),2),") is ",bestmod," with G=",bestg,sep=""),
                                   modelname = bestmod, classification = bestzmap, G = bestg, x = x,
                                   logl = logls[which(bic==max(bic), arr.ind=TRUE)[1],which(bic==max(bic), arr.ind=TRUE)[2]]))
            icllist <- c(icllist, list(allicl = icl[,Gs], parameters = parametersicl, icl = max(icl), 
                                       bestmodel = paste("The best model (ICL of ",round(max(icl),2),") is ",bestmodicl," with G=",bestgicl,sep=""),
                                       classification = bestzmapicl, modelname = bestmodicl, G = bestgicl,
                                       logl = logls[which(icl==max(icl), arr.ind=TRUE)[1],which(icl==max(icl), arr.ind=TRUE)[2]]))
            
            if(length(known)>0){
                store[["tab"]] <- tab
                icllist[["tab"]] <- tabicl
            }
            if(clas>0){
                store[["index"]] <- testindex
            }
            store[["iclresults"]] <- icllist
        }
        else{ ##no models convgerged, BIC matrix only contains -Inf values
            warning("No models converged. Try different models, number of components, or different initialization.")
            
            store <- c(store, list(bic = -Inf, bestmodel = "No models converged...", modelname = models, G = Gs))
            iclresults <- list(icl = -Inf)
            store[["iclresults"]] <- iclresults
        }
        info <- list(univar = univar, gauss = gauss)
        store[["info"]] <- info
        class(store) <- "teigen"
        if(verbose){cat("\n")}
        store
    }

teigen.parallel <- function(x, zmatin, Gs=1:9, numcores=NULL, models="all", scale=TRUE, dfstart=50, clas=0, known=NULL, training=NULL, 
                            gauss=FALSE, dfupdate="approx", eps=c(0.001,0.1), maxit=c(20,1000), convstyle = "aitkens", kno = NULL, ememargs = NULL, origknown = NULL){
    if(!requireNamespace("parallel", quietly = TRUE)){
        stop("Running in parallel uses the parallel package - if unavailable, use teigen(parallel.cores = FALSE) instead")
    }
    
    if(is.null(numcores)){
        numcores <- detectCores()
    }
    teigenModels <- modelgen()
    if(length(models)==1){
        if(models=="dfunconstrained"){
            models <- teigenModels$dfunconstrained
        }
        else{
            if(models=="all"){
                if(ncol(x)==1){
                    models <- teigenModels$univariate
                }
                else{
                    models <- teigenModels$multivariate
                }
            }
            else{
                if(models=="gaussian"){
                    models <- teigenModels$dfconstrained
                }
                else{
                        if(models=="dfconstrained"){
                            models <- teigenModels$dfconstrained
                        }
                        else{
                            if(models=="univariate"){
                                models <- teigenModels$univariate
                            }
                            else{
                                if(models=="altall"){
                                    if(ncol(x)==1){
                                        models <- teigenModels$altunivariate
                                    }
                                    else{
                                        models <- teigenModels$altnames
                                    }
                                }
                            }
                        }
                }
            }
        }
    }
    modrep <- rep(models,length(Gs))
    backwards <- sort(Gs, decreasing=TRUE)
    grep <- rep(backwards, each=length(models))
    if(length(grep[grep==1])>0){
        mod1 <- NA
        cuont <- 1
        CCCCgroup <- c("CCCC", "CCCU", "CUCC", "CUCU", "CUUC","CUUU","UCCC","UCCU","UUCU","UUCC", "UUUC","UUUU")
        CCCCgroup <- c(CCCCgroup,teigenModels$altnames[teigenModels$multivariate%in%CCCCgroup])
        cccdum <- models[models %in% CCCCgroup]
        if(length(cccdum)>0){
            mod1[cuont] <- cccdum[1]
            cuont <- cuont+1
        }
        CICCgroup <- c("CICC","CICU","UICC","UICU","CIUC","CIUU","UIUC","UIUU")
        CICCgroup <- c(CICCgroup,teigenModels$altnames[teigenModels$multivariate%in%CICCgroup])
        cicdum <- models[models %in% CICCgroup]
        if(length(cicdum)>0){
            mod1[cuont] <- cicdum[1]
            cuont <- cuont+1
        }
        CIICgroup <- c("CIIC", "CIIU", "UIIC","UIIU")
        CIICgroup <- c(CIICgroup,teigenModels$altnames[teigenModels$multivariate%in%CIICgroup])
        ciidum <- models[models %in% CIICgroup]
        if(length(ciidum)>0){
            mod1[cuont] <- ciidum[1]
            cuont <- cuont+1
        }
        univgroup <- c("univCC","univCU","univUC","univUU","univVV", "univVE", "univEV", "univEE")
        unidum <- models[models %in% univgroup]
        if(length(unidum)>0){
            mod1[cuont] <- unidum[1]
            cuont <- cuont+1
        }
        modrep <- c(mod1, modrep[!grep==1])
        grep <- c(rep(1,length(mod1)), grep[!grep==1])
    }
    runvec <- 1:length(modrep)
    clus <- makeCluster(numcores)
    clusterEvalQ(clus, library(teigen))
    clusterExport(clus, ls(environment()), envir=environment())
    testparallel <- try(runlist <- parallel::clusterApplyLB(clus, runvec, function(g) teigen.EM(x, zmatin, Gs = grep[g],
                        models=modrep[g], verbose=FALSE, scale=scale, dfstart=dfstart, clas=clas, known=known,
                        training = training, gauss=gauss, dfupdate=dfupdate, eps=eps, maxit=maxit,
                        convstyle = convstyle, kno=kno, ememargs=ememargs, origknown=origknown)), silent=TRUE)
    
    stopCluster(clus)
    if(class(testparallel)=="try-error"){
        stop(testparallel)
    }
    ##runlist contains return values from each parallel node
    gmap <-  unlist(lapply(runlist, function(x) x$G))
    bicmap <-  unlist(lapply(runlist, function(x) x$bic))
    iclmap <-  unlist(lapply(runlist, function(x) x$iclresults$icl))
    modmap <-  unlist(lapply(runlist, function(x) x$modelname))
    ##	modmapicl <- unlist(lapply(runlist, function(x) x$iclresults$modelname))
    ##dim(nodemap) <- c(3, length(models)*length(Gs))
    ##	nodemap <- t(nodemap)
    bestbic <- runlist[[which.max(bicmap)]]
    besticl <- runlist[[which.max(iclmap)]]
    bigbic <- matrix(-Inf, nrow=length(models), ncol=length(Gs))
    colnames(bigbic) <- paste("G=", Gs, sep="")
    rownames(bigbic) <- models
    bigicl <- bigbic
    for(i in 1:length(gmap)){
        bigbic[models==modmap[i],Gs==gmap[i]] <- bicmap[i]
        bigicl[models==modmap[i],Gs==gmap[i]] <- iclmap[i]
    }
    store <- list()
    store <- bestbic
    ##	store[["runlist"]] <- runlist
    store[["iclresults"]] <- besticl$iclresults
    store[["allbic"]] <- bigbic
    store$iclresults[["allicl"]] <- bigicl
    store
}

tfminup <-
    function(mod,G,sg,dg,ag,p,ng,lambdag,submod13){
        fmin <- 0
        if(submod13=="UCC"){
            for(g in 1:G){
                fmin <- fmin + sum(diag((ng[g]*sg[,,g])%*%solve(dg[,,g]%*%ag[,,g]%*%t(dg[,,g]))))/lambdag[g] + p*ng[g]*log(lambdag[g])
            }
        }
        if(submod13=="UUC"){
            for(g in 1:G){
                fmin <- fmin + sum(diag((ng[g]*sg[,,g])%*%solve(dg[,,g]%*%ag[,,g]%*%t(dg[,,g]))))/lambdag[g] + p*ng[g]*log(lambdag[g])
            }
        }
        if(submod13=="CCU"){
            for(g in 1:G){
                fmin <- fmin + sum(diag(dg[,,g] %*% solve(ag[,,g]) %*% t(dg[,,g]) %*% (ng[g]*sg[,,g])))/lambdag[g] + p*ng[g]*log(lambdag[g])
            }
        }
        fmin
    }

tft <-
    function(x,G,pig,dhfgs78,p,mug,sigmainv,n,sigma,univar,delta,gauss){
        log_num <- matrix(0,n,G)
        if(gauss){
            if(univar){
                for(g in 1:G){
                    log_num[,g] <- log(pig[g])-(p/2)*log(2*pi)-(1/2)*log(sigma[,,g])-
                        (1/2)*delta[,g]
                }
            }
            else{
                for(g in 1:G){
                    log_num[,g] <- log(pig[g])-(p/2)*log(2*pi)-(1/2)*log(det(sigma[,,g]))-
                        (1/2)*delta[,g]
                }
            }
        }
        else{
            if(univar){
                for(g in 1:G){
                    log_num[,g]<-log(pig[g])+lgamma((dhfgs78[g]+1)/2)-(1/2)*log(sigma[,,g])-
                        ((p/2)*(log(pi)+log(dhfgs78[g]))+lgamma(dhfgs78[g]/2)+((dhfgs78[g]+p)/2)*
                             (log(1+ delta[,g]/dhfgs78[g])))
                }
            }
            else{
                for(g in 1:G){
                    log_num[,g]<-log(pig[g])+lgamma((dhfgs78[g]+p)/2)-(1/2)*log(det(sigma[,,g]))-
                        ((p/2)*(log(pi)+log(dhfgs78[g]))+lgamma(dhfgs78[g]/2)+((dhfgs78[g]+p)/2)*(log(1+
                                                                                                          delta[,g]/dhfgs78[g])))
                }
            }
        }
        log_num
    }

tgivenz <-
    function(n,G,known,initg,testindex,clas){
        zmat <- matrix(0,n,G)
        if(clas>0){
            matchtab <- table(known[testindex],initg[testindex])
            matchit <- NA
            for(i in 1:G){
                matchit[i] <- which.max(matchtab[,i])
                
            }
            if(length(unique(matchit))==length(matchit)){
                initnew <- initg
                for(i in 1:G){
                    initnew[initg==i] <- matchit[i]
                }
                initg <- initnew
            }
        }
        for(i in 1:G){
            zmat[initg==i, i]<-1
        }
        zmat
    }

tkmeansz <-
    function(x,n,G,known,kno,testindex,clas){
        kclus <- kmeans(x,G, nstart = 50)$cluster
        zmat <- matrix(0,n,G)
        if(clas>0){
            matchtab <- table(known[testindex],kclus[testindex])
            matchit <- NA
            for(i in 1:G){
                matchit[i] <- which.max(matchtab[,i])
                
            }
            if(length(unique(matchit))==length(matchit)){
                knew <- kclus
                for(i in 1:G){
                    knew[kclus==i] <- matchit[i]
                }
                kclus <- knew
            }
        }
        for(i in 1:G){
            zmat[kclus==i, i]<-1
        }
        zmat
    }

tlambdaginit <-
    function(p,G,sg,mod){
        lambdag <- rep.int(0,G)
        if(substring(mod,1,1)=="U"){
            for(g in 1:G){
                lambdag[g] <- det(sg[,,g])^(1/p)
            }
        }
        if(substring(mod,1,1)=="C"){
            thing <- 0
            for(g in 1:G){
                thing <- thing + det(sg[,,g])^(1/p)
            }
            lambdag <- rep(thing,G)
        }
        lambdag
    }

tlambdagupdate <-
    function(G,mod,sg,dg,ag,p,n,ng,submod13){
        lambdag <- rep.int(0,G)
        if(submod13=="CUC"){
            omegag <- matrix(0,p,p)
            for(g in 1:G){
                diag(omegag) <- diag(omegag) + eigen(ng[g]*sg[,,g])$values
            }
            lambdag <- rep((det(omegag)^(1/p))/n,G)
        }
        if(submod13=="CUU"){
            thing <- 0
            for(g in 1:G){
                thing <- thing + (det(ng[g]*sg[,,g])^(1/p))/n
            }
            lambdag <- rep(thing,G)
        }
        if(submod13=="CIC"){
            for(g in 1:G){
                lambdag[g] <- (prod(diag(n*sg[,,g]))^(1/p))/n
            }
        }
        if(submod13=="CIU"){
            thing <- 0
            for(g in 1:G){
                thing <- thing + (prod(diag(ng[g]*sg[,,g]))^(1/p))/n
            }
            lambdag <- rep(thing,G)
        }
        if(submod13=="UIU"){
            for(g in 1:G){
                lambdag[g] <- (prod(diag(ng[g]*sg[,,g]))^(1/p))/ng[g]
            }
        }
        if(any(submod13==c("CII","UII"))){
            for(g in 1:G){
                lambdag[g] <- sum(diag(sg[,,g]))/p
            }
        }
        ##IP
        if(submod13=="UCC"){
            for(g in 1:G){
                lambdag[g] <- sum(diag(sg[,,g] %*% solve(dg[,,g] %*% ag[,,g] %*% t(dg[,,g]))))/p
            }
        }
        if(submod13=="UUC"){
            for(g in 1:G){
                lambdag[g] <- sum(diag(sg[,,g] %*% dg[,,g] %*% diag(1/diag(ag[,,g])) %*% t(dg[,,g])))/p
            }
        }
        if(submod13=="UIC"){
            for(g in 1:G){
                lambdag[g] <- sum(diag(sg[,,g] %*% diag(1/diag(ag[,,g]))))/p
            }
        }
        if(submod13=="CCU"){
            for(g in 1:G){
                lambdag <- lambdag + sum(diag(ng[g]*sg[,,g] %*% dg[,,g] %*% diag(1/diag(ag[,,g])) %*% t(dg[,,g])))/(n*p)
            }
        }
        if(submod13=="UCU"){
            for(g in 1:G){
                lambdag[g] <- sum(diag(sg[,,g] %*% dg[,,g] %*% diag(1/diag(ag[,,g])) %*% t(dg[,,g])))/p
            }
        }
        lambdag
    }

tmuginit <-
    function(G,p,x,zmat,ng,univar){
        mug <- matrix(0,G,p)
        for(g in 1:G){
            mug[g,] <- colSums(zmat[,g]*x)/ng[g]
        }
        mug
    }

tmugupdate <-
    function(G,zmat,w,x,p,univar){
        mug <- matrix(0,G,p)
        if(univar){
            for(g in 1:G){
                mug[g,1] <- colSums(zmat[,g]*w[,g]*x)/sum(zmat[,g]*w[,g])
            }
        }
        else{
            for(g in 1:G){
                mug[g,] <- colSums(zmat[,g]*w[,g]*x)/sum(zmat[,g]*w[,g])
            }
        }
        mug
    }

tngupdate <-
    function(zmat){
        ng <- colSums(zmat)
        ng
    }

tpigupdate <-
    function(ng,n){
        pig <- ng/n
        pig
    }

tsginit <-
    function(p,G,x,mug,zmat,n,ng){
        sg <- array(0, dim=c(p, p, G))
        for(g in 1:G){
            sg[,,g] <- covwtC(x, zmat[,g], mug[g,])#cov.wt(x,wt=zmat[,g],center=mug[g,],method="ML")$cov
        }
        sg
    }

tsginitc <-
    function(G,sg,pig,p,n,x){
        sgc <- matrix(0,p,p)
        for(g in 1:G){
            sgc <- sgc + pig[g]*sg[,,g]
        }
        sgc
    }

tsgupdate <-
    function(p,G,n,x,mug,zmat,w,ng,mod,pig,submod13){
        sg <- array(0, dim=c(p, p, G))
        for(g in 1:G){
            sg[,,g] <- covwtC(x, zmat[,g]*w[,g], mug[g,])*(sum(zmat[,g]*w[,g])/ng[g])
        }
        if(any(submod13==c("CCC","CIC","CII")) | mod=="univCC"|mod=="univCU"){
            sgc <- matrix(0,p,p)
            for(g in 1:G){
                sgc <- sgc + pig[g]*sg[,,g]
            }
            sg[,,] <- sgc
        }
        sg
    }

tsigmainvup <-
    function(p,G,sigma){
        sigmainv <- array(0, dim=c(p, p, G))
        for(g in 1:G){
            sigmainv[,,g] <- solve(sigma[,,g])
        }
        sigmainv
    }

tsigmaup <-
    function(p,G,sg,lambdag,dg,ag,mod,univar,submod13){
        sigma <- array(0, dim=c(p, p, G))
        if(any(submod13==c("CCC","UUU")) | univar){
            sigma <- sg
        }
        else{
            if(substring(mod,2,2)=="I"){
                for(g in 1:G){
                    sigma[,,g] <- lambdag[g] * ag[,,g]
                }
            }
            else{
                for(g in 1:G){
                    sigma[,,g] <- lambdag[g] * (dg[,,g] %*% ag[,,g] %*% t(dg[,,g]))
                }
            }
        }
        sigma
    }

tuniformz <-
    function(n,G,clas,kno,known){
        zmat <- matrix(0, n, G)
        if(clas>0){
            for(i in 1:G){
                zmat[as.numeric(known)==i,i] <- 1
            }
            zmat[kno==0,] <- 1/G
        }
        zmat
    }

tvginit <-
    function(dfstart,G){
        vg <- c(rep.int(dfstart, G))
        vg
    }

twinit <-
    function(x,n,G,mug,sigmainv,vg,p,sg,zmat,univar,sigma){
        w <- matrix(0,n,G)
        delta <- matrix(0,n,G)
        if(!univar){
            for(g in 1:G){
                delta[,g] <- mahalanobis(x, mug[g,], sigmainv[,,g], inverted=TRUE)
                w[,g] <- (vg[g]+p)/(vg[g]+delta[,g])
            }
        }
        else{
            for(g in 1:G){
                delta[,g] <- mahalanobis(x, mug[g,], 1/sigma[,,g],inverted=TRUE)
                w[,g] <- (vg[g]+p)/(vg[g]+delta[,g])
            }
        }
        w
    }

twupdate <-
    function(x,n,G,mug,sigmainv,dhfgs78,p,univar,sigma,delta){
        w <- matrix(0,n,G)
        if(!univar){
            for(g in 1:G){
                w[,g] <- (dhfgs78[g]+p)/(dhfgs78[g]+delta[,g])
            }
        }
        else{
            for(g in 1:G){
                w[,g] <- (dhfgs78[g]+p)/(dhfgs78[g]+delta[,g])
            }
        }
        w
    }

tzupdate <-
    function(x,G,pig,dhfgs78,p,mug,sigmainv,n,sigma,clas,kno,known,unkno,univar,delta,gauss,cycle){
        log_num <- matrix(0,n,G)
        if(gauss){
            if(univar){
                for(g in 1:G){
                    log_num[,g] <- log(pig[g])-(p/2)*log(2*pi)-(1/2)*log(sigma[,,g])-
                        (1/2)*delta[,g]
                }
            }
            else{
                for(g in 1:G){
                    log_num[,g] <- log(pig[g])-(p/2)*log(2*pi)-(1/2)*log(det(sigma[,,g]))-
                        (1/2)*delta[,g]
                }
            }
        }
        else{
            if(univar){
                for(g in 1:G){
                    log_num[,g]<-log(pig[g])+lgamma((dhfgs78[g]+1)/2)-(1/2)*log(sigma[,,g])-
                        ((p/2)*(log(pi)+log(dhfgs78[g]))+lgamma(dhfgs78[g]/2)+((dhfgs78[g]+p)/2)*
                             (log(1+ delta[,g]/dhfgs78[g])))
                }
            }
            else{
                for(g in 1:G){
                    log_num[,g]<-log(pig[g])+lgamma((dhfgs78[g]+p)/2)-(1/2)*log(det(sigma[,,g]))-((p/2)*(log(pi)+log(dhfgs78[g]))
                                                                                                  +lgamma(dhfgs78[g]/2)+((dhfgs78[g]+p)/2)*(log(1+ delta[,g]/dhfgs78[g])))
                }
            }
        }
        num <- exp(log_num)
        zmat <- num/rowSums(num)
        if(clas>0){
            zmat <- unkno*zmat
            for(i in 1:n){
                if(kno[i]==1){
                    zmat[i, known[i]] <- 1
                }
            }
        }
        zmat
    }

yxf7 <-
    function(mod,dhfgs78,ng,zmat,w,G,p,n,x,mug,sigmainv){
        if(substring(mod,4,4)=="U"|mod=="univUU"|mod=="univCU"){
            dfoldg <- dhfgs78
            for(g in 1:G){
                constn <- 1 + (1/ng[g]) * sum(zmat[,g]*(log(w[,g])-w[,g])) +
                    digamma((dfoldg[g]+p)/2) - log((dfoldg[g]+p)/2)
                dhfgs78[g] <- uniroot( function(v) log(v/2) - digamma(v/2) + constn ,
                                       lower=0.0001, upper=1000, tol=0.00001)$root
                ##constn <- -constn
                ##print(constn)
                ##print(paste("fix", dhfgs78[g] - ((-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg[g]/2))-(dfoldg[g]/2-1/2)))/(1-exp(constn)))))
                ##print(paste("old", dhfgs78[g] - (-exp(constn)/(1-exp(constn)))))
                ##print(c((2-sqrt(4+4*4*exp(constn)))/2, (2+sqrt(4+4*4*exp(constn)))/2))
                ##print((-exp(constn))/(1-exp(constn)))
                if(dhfgs78[g]>200){
                    dhfgs78[g]<-200
                }
                if(dhfgs78[g]<2){
                    dhfgs78[g]<-2
                }
            }
        }
        if(substring(mod,4,4)=="C"|mod=="univUC"|mod=="univCC"){
            dfoldg <- dhfgs78[1]
            constn <- 1 + (1/n) * sum(zmat * (log(w)-w)) +
                digamma((dfoldg+p)/2) - log((dfoldg+p)/2)
            dfsamenewg <- uniroot( function(v) log(v/2) - digamma(v/2) + constn,
                                   lower=0.0001, upper=1000, tol=0.01)$root
            if(dfsamenewg>200){
                dfsamenewg<-200
            }
            if(dfsamenewg<2){
                dfsamenewg<-2
            }
            dhfgs78 <- c(rep(dfsamenewg,G))
        }
        dhfgs78
    }

yxf8 <- function(mod,dhfgs78,ng,zmat,w,G,p,n,x,mug,sigmainv){
    if(substring(mod,4,4)=="U"|mod=="univUU"|mod=="univCU"){
        dfoldg <- dhfgs78
        for(g in 1:G){
            constn <- 1 + (1/ng[g]) * sum(zmat[,g]*(log(w[,g])-w[,g])) + digamma((dfoldg[g]+p)/2) - log((dfoldg[g]+p)/2)
            ## 			dhfgs78[g] <- uniroot( function(v) log(v/2) - digamma(v/2) + constn ,
            ## 					lower=0.0001, upper=1000, tol=0.00001)$root
            constn <- -constn
            ##print(constn)
            ##print(paste("fix", dhfgs78[g] - ((-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg[g]/2))-(dfoldg[g]/2-1/2)))/(1-exp(constn)))))
            ##print(paste("old", dhfgs78[g] - (-exp(constn)/(1-exp(constn)))))
            ##print(c((2-sqrt(4+4*4*exp(constn)))/2, (2+sqrt(4+4*4*exp(constn)))/2))
            ##print((-exp(constn))/(1-exp(constn)))
            dhfgs78[g] <- (-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg[g]/2))-(dfoldg[g]/2-1/2)))/(1-exp(constn))
            ##     dhfgs78[g] <- (-exp(constn))/(1-exp(constn))
            if(dhfgs78[g]>200){
                dhfgs78[g]<-200
            }
            if(dhfgs78[g]<2){
                dhfgs78[g]<-2
            }
        }
    }
    if(substring(mod,4,4)=="C"|mod=="univUC"|mod=="univCC"){
        dfoldg <- dhfgs78[1]
        constn <- 1 + (1/n) * sum(zmat * (log(w)-w)) + digamma((dfoldg+p)/2) - log((dfoldg+p)/2)
        constn <- -constn
        ##		dfsamenewg <- uniroot( function(v) log(v/2) - digamma(v/2) + constn,
        ##				lower=0.0001, upper=1000, tol=0.01)$root
        dfsamenewg <- (-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg/2))-(dfoldg/2-1/2)))/(1-exp(constn))
        if(dfsamenewg>200){
            dfsamenewg<-200
        }
        if(dfsamenewg<2){
            dfsamenewg<-2
        }
        dhfgs78 <- c(rep(dfsamenewg,G))
    }
    dhfgs78
}

##Initialization using emEM approach from Biernacki et al (2003)
ememinit <- function(ememargs, x, g, scale, dfstart, clas, known, training, gauss, dfupdate, eps, n){
  ##ememargs list with
  ##numstart - number of random starts
  ##iter -  number of EM iterations
  ##model -  model name
  ##init -  initialization for emEM - hard, soft, kmeans
  emruns <- list()
  for(i in 1:ememargs$numstart){
    suppressWarnings(emruns[[i]] <- teigen(x,Gs=g, models=ememargs$model,init=ememargs$init,scale,dfstart,known,training,gauss,dfupdate,eps,verbose=FALSE,maxit=c(ememargs$iter,5))[c("logl","classification")])
  }
  logls <- unlist(lapply(emruns, function(v) v$logl))
  clindex <- which.max(logls)
  cl <- emruns[[clindex]]$classification
  zmat <- matrix(0,n,g)
  for(i in 1:g){
    zmat[cl==i, i]<-1
  }
  zmat
}