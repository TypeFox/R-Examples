### R code from vignette source 'SII.Rnw'

###################################################
### code chunk number 1: load.package
###################################################
library(SII)


###################################################
### code chunk number 2: load.data
###################################################
## Table 1: Critical band SII procedure constants
data("critical")
head(critical)

## Table 2:Equally contributing (17 band) critical band SII 
## procedure constants
data("equal")
  
## Table 3: One-third octave band SII procedure constants
data("onethird")

## Table 4: Octave band SII procedure constants
data("octave")
      
## Overall SPL constants
data("overall.spl")
overall.spl


###################################################
### code chunk number 3: load.alternative.data
###################################################
## Table B.1: Critical band importance functions for various speech tests.
data(sic.critical)
head(sic.critical)

## Table B.2: One-third octave band importance functions for various speech tests.
data(sic.onethird)

## Table B.3: Octave band importance functions for various speech tests.
data(sic.octave)


###################################################
### code chunk number 4: sicplot
###################################################
data(sic.critical)
ngroup <- ncol(sic.critical)
matplot(x=sic.critical[,1], y=sic.critical[,-1],
        type="o",
        xlab="Frequency, Hz",
        ylab="Weight",
        log="x",
        lty=1:ngroup,
        col=rainbow(ngroup)
)
legend(
       "topright",
       legend=names(sic.critical)[-1],
       pch=as.character(1:ngroup),
       lty=1:ngroup,
       col=rainbow(ngroup)
       )



###################################################
### code chunk number 5: header
###################################################
args(sii)


###################################################
### code chunk number 6: ExampleC.1
###################################################
sii.C1 <- sii(
              speech   = c(50.0, 40.0, 40.0, 30.0, 20.0,  0.0),
              noise    = c(70.0, 65.0, 45.0, 25.0,  1.0,-15.0),
              threshold= c( 0.0,  0.0,  0.0,  0.0,  0.0,  0.0),
              method="octave"
              #, importance="SII"
              #, importance=octave$Ii
              , importance="CST"
              )
round(sii.C1$table[,-c(5:7,13)],2)
sii.C1$sii


###################################################
### code chunk number 7: ExampleC.1
###################################################
sii.C2 <- sii(
              speech   = rep(54.0, 18),
              noise    = c(40.0, 30.0, 20.0, rep(0, 18-3) ),
              threshold= rep(0.0,  18),
              method="one-third"
              )
sii.C2$table[1:3,1:8]
sii.C2$sii


###################################################
### code chunk number 8: ExampleC.1
###################################################
sii.C1 <- sii(
              speech   = c(50.0, 40.0, 40.0, 30.0, 20.0,  0.0),
              noise    = c(70.0, 65.0, 45.0, 25.0,  1.0,-15.0),
              threshold= c( 0.0,  0.0,  0.0,  0.0,  0.0,  0.0),
              method="octave"
              #, importance="SII"
              #, importance=octave$Ii
              , importance="CST"
              )
round(sii.C1$table[,-c(5:7,13)],2)
sii.C1$sii


###################################################
### code chunk number 9: LoadData (eval = FALSE)
###################################################
## library(gdata)
## patInfo <- read.xls("../AI subject list.xls")


###################################################
### code chunk number 10: ExamineData (eval = FALSE)
###################################################
## ## Information about variables
## str(patInfo)


###################################################
### code chunk number 11: ExamineData (eval = FALSE)
###################################################
## ## First 6 rows
## head(patInfo)


###################################################
### code chunk number 12: vars (eval = FALSE)
###################################################
## ## measured frequencies
## freq <- c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000)
## 
## ## columns containing frequencies for the right/left ear
## rt.cols <- paste("PTR",freq, sep="")
## lt.cols <- paste("PTL",freq, sep="")
## 
## rt.vals <- patInfo[,rt.cols]
## lt.vals <- patInfo[,lt.cols]


###################################################
### code chunk number 13: missing (eval = FALSE)
###################################################
## rt.vals[rt.vals==-888] <- NA
## lt.vals[rt.vals==-888] <- NA


###################################################
### code chunk number 14: fun (eval = FALSE)
###################################################
## fun <- function(X)
##   {
##     ret <- try( 
##                sii(X, 
##                    speech="raised",
##                    threshold=c(15,15,20,25,35,35,45,50), 
##                    freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
##                    importance="SII",
##                    interpolate=TRUE
##                    )$sii
##                )
##   if("try-error" %in% class(ret))
##     return(NA)
##   else
##     return(ret)
##   }
## 
## ## Test it
## fun( rt.vals[1,] )


###################################################
### code chunk number 15: CalculateSII (eval = FALSE)
###################################################
## sii.right <- apply(rt.vals, 1, fun )
## 
## sii.left  <- apply(lt.vals, 1, fun )


###################################################
### code chunk number 16: AddSIIToTable (eval = FALSE)
###################################################
## patInfo$"SII.right" <- sii.right
## patInfo$"SII.left" <- sii.left
## 
## tail(patInfo)


###################################################
### code chunk number 17: SavePatInfo (eval = FALSE)
###################################################
## write.table(patInfo, 
##             file="../AI subject list-SII.xls",
##             row.names=FALSE, 
##             sep="\t"
##             )


###################################################
### code chunk number 18: SIIFileFun (eval = FALSE)
###################################################
## sii.dina <- function(infile, outfile, verbose=TRUE) 
##   {
##     if(verbose)
##       cat("\nLoading data file '", infile, "'...\n", sep="")
##     ## Load the data
##     library(gdata)
##     patInfo <- read.xls(infile)
##     
##     ## measured frequencies
##     freq <- c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000)
## 
##     if(verbose)
##       cat("\nExtracting hearing thresholds...\n")
## 
##     ## columns containing frequencies for the right/left ear
##     rt.cols <- paste("PTR",freq, sep="")
##     lt.cols <- paste("PTL",freq, sep="")
##     
##     rt.vals <- patInfo[,rt.cols]
##     lt.vals <- patInfo[,lt.cols]
##     
##     ## Handle missing code '-888'
##     rt.vals[rt.vals==-888] <- NA
##     lt.vals[rt.vals==-888] <- NA
##     
##     ## define function to compute SII with our defaults
##     fun <- function(X)
##       {
##         ret <- try( 
##                    sii(X, 
##                        speech="raised",
##                        threshold=c(15,15,20,25,35,35,45,50), 
##                        freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
##                        importance="SII",
##                        interpolate=TRUE
##                        )$sii #$
##                    )
##       if("try-error" %in% class(ret))
##         return(NA)
##       else
##         return(ret)
##       }
##     
##     ## Calculate SII
##     if(verbose)
##       cat("\nCalculating right ear SII...\n")
##     sii.right <- apply(rt.vals, 1, fun )
##     if(verbose)
##       cat("\nCalculating left ear SII...\n")
##     sii.left  <- apply(lt.vals, 1, fun )
##     
##     ## Add back onto the table
##     patInfo$"SII.right" <- sii.right
##     patInfo$"SII.left" <- sii.left
## 
##     if(verbose)
##       cat("\nWriting new data table as '", outfile, "'...\n", sep="")
## 
##     ## Save file
##     write.table(patInfo, 
##                 file=outfile,
##                 row.names=FALSE, 
##                 sep="\t"
##                 )
## 
##     if(verbose)
##       cat("\nDone.\n\n")
##   }


###################################################
### code chunk number 19: RunFun (eval = FALSE)
###################################################
## sii.dina(infile="../AI subject list.xls",
##          outfile="../AI subject list-SII.xls") 


###################################################
### code chunk number 20: CompareAItoSII (eval = FALSE)
###################################################
## library(xtable)
## xt <- xtable(patInfo[,c("AI","SII.right","SII.left")], 
##              caption="Comparison of original AI and new right and left SII values \\label{table}", 
##              digits=2)
## print(xt)


###################################################
### code chunk number 21: PlotAIandSII (eval = FALSE)
###################################################
## ## put histograms on the diagonal
## panel.hist <- function(x, ...)
##   {
##     usr <- par("usr"); on.exit(par(usr))
##     par(usr = c(usr[1:2], 0, 1.5) )
##     h <- hist(x, plot = FALSE)
##     breaks <- h$breaks; nB <- length(breaks)
##     y <- h$counts; y <- y/max(y)
##     rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
##   }
## pairs.2 <- function(x)
##   pairs(x, 
##       panel=panel.smooth,
##       cex = 1.5, 
##       pch = 24, 
##       bg="light blue",
##       diag.panel=panel.hist, 
##       cex.labels = 2, 
##       font.labels=2
##       )
## 
## 
## pairs.2(patInfo[,c("AI","SII.right","SII.left")])


###################################################
### code chunk number 22: SII.Rnw:438-479
###################################################
latexhelp <- function(topics, package=NULL)
  {
    if(missing(topics) && missing(package))
      stop("Help on what?")
    else if (missing(topics))
      {
        callObj <- call("library", help=package)
        topicStrings <- eval(callObj)$info[[2]]
        topicStrings <- lapply( topicStrings, strsplit, split=" ")
        topicStrings <- sapply( topicStrings, function(x) x[[1]][1])
        topicStrings <- topicStrings[topicStrings > " "]
        topics <- topicStrings
        }
    
    for(topic in topics)
      {
        if(!is.character(topic))
          topic <- deparse(substitute(topic))
       
        callObj <- call("help", topic=topic, package=package)
        
        ## Get the file path
        help.file<- as.character(eval(callObj))

        ## Get R to translate it to latex and display it
        tools::Rd2latex(utils:::.getHelpFile(help.file))
        
        ## new page
        cat("\n")
        cat("\\clearpage")
        cat("\n")
      }
  }


tmp <- function (x, ...) 
{
}

latexhelp(package="SII")



###################################################
### code chunk number 23: splineFit
###################################################

THDI=c(25,25,30,35,45,45,55,60)
freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000)

sii.freqs <- SII:::sii.constants[,1]

xlist <- sort(c(SII:::sii.constants[, 1], freq))
ylist <- rep(NA, length=length(xlist))
names(ylist) <- xlist
ylist[ as.character(freq) ] <- THDI


###################################################
### code chunk number 24: interSpline
###################################################
library(splines)
ispl <- interpSpline( THDI ~ freq )
ispl <- predict(ispl, sii.freqs)$y
      
ispl.l <- interpSpline( THDI ~ log(freq) )
ispl.l <- predict(ispl.l, log(sii.freqs) )$y

approx.l <- function(x,y,xout,...)
  { 
    retval <- approx(log(x), y, log(xout), ...)
    retval$x <- xout
    retval
  }



###################################################
### code chunk number 25: Spline_comparison_function
###################################################
doplot <- function(logx=FALSE, separate=FALSE)
  {
    if(separate)
      layout( cbind( c(1,2,3,7), c(4,5,6,7) ))
    tmp <- function() {
      plot(x=freq, y=THDI, col="black", cex=2, log=if(logx) "x" else "", 
           xlab="Frequency (Herz)", ylab="Threshhold of Detection (dB)",
           xlim=range(sii.freqs), ylim=c(20,60) )
      lines(SII:::sii.constants[, 1], SII:::sii.constants[, "Ti'.THDN"], col="red", lwd=3)
      abline(v=sii.freqs,lty=3)
    }
    
    tmp()
    lines( x=SII:::sii.constants[, 1], y=ispl, col="blue", lwd=2)
    
    if(separate) tmp()
    lines( x=SII:::sii.constants[, 1], y=ispl.l, col="cyan", lwd=2)
    
    if(separate) tmp()
    lines( spline  (x=freq, y=THDI), col="green", lwd=2)
    
    if(separate) tmp()
    lines( spline  (x=freq, y=THDI, method="natural", xout=sii.freqs), col="orange", lwd=2)
    
    if(separate) tmp()
    lines( approx  (x=freq, y=THDI, method="linear",  xout=sii.freqs, rule=2), col="magenta", lwd=2)
    
    if(separate) tmp()
    lines( approx.l(x=freq, y=THDI, method="linear",  xout=sii.freqs, rule=2), col="brown",   lwd=3)

    if(separate) plot.new()
    legend(if(separate) "center" else { if(logx)  "topleft" else "bottomright"  },
           legend=c(
             "Measured data",
             "Matlab INTERP1(X, Y, XI, 'spline', 'extrap')",
             "R's predict( interSpline(X,Y), XI )",
             "R's predict( interSpline(log(X),Y), log(XI) )",
             "R's spline(x, y)",
             "R's spline(x, y, method='natural')",
             "R's approx(x, y, method='linear')",
             "R's approx(x, y, method='linear') (on log scale)",
             "Frequencies used for SII calculation"
             ),
           col=c("black", "red", "blue","cyan","green","orange","magenta","brown", "black"),
           pch=c(      1,    NA,    NA,    NA,     NA,      NA,       NA,      NA,     NA),
           lty=c(     NA,     1,     1,     1,      1,       1,        1,       1,      3),
           lwd=c(     NA,     3,     2,     2,      2,       2,        2,       2,      1),
           bg="white",
           cex=if(separate) 1.25 else 0.75,
           ncol=if(separate)   2 else 1
           )
    if(separate)
      layout( 1 )
    title(paste("Spline method comparison\n", if(logx)"(log scale)"else"(natural scale)", "\n"))
    
  }


###################################################
### code chunk number 26: natural_scale_one_plot
###################################################
doplot(FALSE, FALSE)


###################################################
### code chunk number 27: natural_scale_separate_plots
###################################################
doplot(FALSE, TRUE)


###################################################
### code chunk number 28: log_scale_one_plot
###################################################
doplot(TRUE, FALSE)


###################################################
### code chunk number 29: log_scale_separate_plots
###################################################
doplot(TRUE, TRUE)


###################################################
### code chunk number 30: SII:::sii.excel
###################################################
args(SII:::sii.excel)


###################################################
### code chunk number 31: sii
###################################################
args(sii)


###################################################
### code chunk number 32: sii.print
###################################################



###################################################
### code chunk number 33: show.reload.constants
###################################################
SII:::reload.constants 


###################################################
### code chunk number 34: run.reload.constants (eval = FALSE)
###################################################
## SII:::reload.constants(xls.path="./SII/extdata")


###################################################
### code chunk number 35: ComparisonTable1
###################################################
sii.left <- sii(
                speech="raised",
                threshold=c(25,25,30,35,45,45,55,60),
                freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                importance="NU6",
                interpolate=TRUE
                )


## comparison of our interpolation and matlab's
tab <- cbind(
             matlab=SII:::sii.constants[,"Ti'.THDN"],
             R     =sii.left$table[,"T'i"],
             delta =NA
             )
tab[,3] <- tab[,1] - tab[,2]
rownames(tab) <- SII:::sii.constants[,".NFreqLin"]
t(round(tab,2))


###################################################
### code chunk number 36: <plot.sii
###################################################
compare.plot <- function(x, matlab, title)
  {
    plot(x)
    lines(SII:::sii.constants[,".NFreqLin"],
          matlab,
          type="l", col="red", lwd=3)

    legend("topleft",
           legend=c(
             "Measured data",
             "Matlab INTERP1(X, Y, XI, 'spline', 'extrap')",
             "R's approx(X,Y, XI, xout=XI,\n method='linear', rule=2)"  
             ),
           col=c("black", "red","blue","green","orange","magenta"),
           pch=c(      1,    NA,    2,     NA,      NA,       NA),
           lty=c(     NA,     1,    1,      1,       1,        1),
           lwd=c(     NA,     3,    2,      2,       2,        2),
           bg="white"
           )
    title(title)
  }


###################################################
### code chunk number 37: ComparisonFigure1
###################################################

compare.plot(sii.left, matlab=SII:::sii.constants[, "Ti'.THDN"], title="Spline method comparison, Left Ear")
      


###################################################
### code chunk number 38: ComparisonTable2
###################################################
# comparison of our interpolation and matlab's

matlab <- c(15.00, 15.00, 13.34, 14.27, 16.06, 17.82, 19.18, 20.00,
            20.26, 20.54, 21.44, 23.36, 26.93, 31.38, 34.63, 35.13, 
            35.00, 38.29, 43.98, 48.80, 49.58)

sii.right <- sii( 
                 speech="raised",
                 threshold=c(15,15,20,25,35,35,45,50), 
                 freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                 importance="NU6",
                 interpolate=TRUE
                 )

tab <- cbind(
             matlab=matlab,
             R     =sii.right$table    [,"T'i"],
             delta =NA
             )
tab[,3] <- tab[,1] - tab[,2]
rownames(tab) <- SII:::sii.constants[,".NFreqLin"]
t(round(tab,2))


###################################################
### code chunk number 39: ComparisonFigure2
###################################################
compare.plot(sii.right, matlab=matlab, title="Spline method comparison, Right Ear")


###################################################
### code chunk number 40: <sii.left.old
###################################################
SII:::sii.excel( 
          c(25,25,30,35,45,45,55,60),
          c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000) 
          )


###################################################
### code chunk number 41: sii.left
###################################################
sii.left <- sii(
                speech="raised",
                threshold=c(25,25,30,35,45,45,55,60),
                freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                importance="NU6",
                interpolate=TRUE
                )
sii.left


###################################################
### code chunk number 42: <sii.right.old
###################################################
SII:::sii.excel( 
          c(15,15,20,25,35,35,45,50), 
          c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000) 
          )


###################################################
### code chunk number 43: sii.right
###################################################
sii.right <- sii( 
                 speech="raised",
                 threshold=c(15,15,20,25,35,35,45,50), 
                 freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                 importance="NU6",
                 interpolate=TRUE
                 )
sii.right


###################################################
### code chunk number 44: <sii.best.old
###################################################
SII:::sii.excel( 
          rep(0,8),
          freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000) 
          )


###################################################
### code chunk number 45: sii.best
###################################################
sii.best <- sii( 
                threshold=rep(0,8),
                freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                interpolate=TRUE 
                )
sii.best


###################################################
### code chunk number 46: <sii.worst.old
###################################################
SII:::sii.excel( 
          rep(100,8),
          c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000) 
          )


###################################################
### code chunk number 47: sii.worst
###################################################
sii.worst <- sii( 
                 threshold=rep(100,8),
                 freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                 interpolate=TRUE
                )
sii.worst


###################################################
### code chunk number 48: sii.missing
###################################################
sii.worst <- sii( 
                 threshold=c(NA, rep(100,7)),
                 freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                 interpolate=TRUE
                )
sii.worst


###################################################
### code chunk number 49: sii.missing
###################################################
sii.right <- sii( 
                 speech="raised",
                 threshold=c(0,15,15,20,25,35,35,45,50), 
                 freq=c(NA, 250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                 importance="NU6",
                 interpolate=TRUE
                 )
sii.right


###################################################
### code chunk number 50: sii.all.missing
###################################################
## This should fail, because there is no data!
sii.NONE <- try(
             sii( 
                 threshold=rep(NA,8),
                 freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                 interpolate=TRUE
                )
            )
sii.NONE


###################################################
### code chunk number 51: sii.missing
###################################################
sii.right <- sii( 
                 speech="raised",
                 threshold=c(15,15,20,NA,35,35,45,50), 
                 freq=c(250, 500, 1000, 2000, 3000, 4000, 6000, 8000),
                 importance="NU6",
                 interpolate=TRUE
                 )
sii.right


###################################################
### code chunk number 52: missing.C.1
###################################################
sii.C1.NA <- sii(
              speech   = c(50.0, 40.0, 40.0, NA,   20.0,  0.0),
              noise    = c(70.0, 65.0, 45.0, 25.0,  1.0,-15.0),
              threshold= c( 0.0,  0.0,  0.0,  0.0,   NA,  0.0),
              freq     = c( 250,  500, 1000, 2000, 4000, 8000),
              method="octave",
              importance="CST",
              interpolate=TRUE
              )
sii.C1
sii.C1.NA


###################################################
### code chunk number 53: save data
###################################################
save.image("SII-Code.Rda")


