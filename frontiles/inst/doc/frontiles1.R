### R code from vignette source 'frontiles1.rnw'

###################################################
### code chunk number 1: frontiles1.rnw:28-33
###################################################
owidth <- getOption("width")
options("width"=70)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0


###################################################
### code chunk number 2: bfig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 7, height = 7, pointsize = 18, bg = "white")


###################################################
### code chunk number 3: bfig2 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 7, height = 7, pointsize = 14, bg = "white")


###################################################
### code chunk number 4: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.65\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 5: zfig2 (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.65\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 6: frontiles1.rnw:64-68
###################################################
library(frontiles)
data(burposte)
smallest<-sort(burposte$xinput,index.return=TRUE)$ix[1:4000]
sample.burposte<-burposte[1:9521%in%smallest,]


###################################################
### code chunk number 7: plot1 (eval = FALSE)
###################################################
## xtab<-sample.burposte$xinput
## ytab<-sample.burposte$yprod
## plot(ytab~xtab, pch=16,col='blue2')


###################################################
### code chunk number 8: frontiles1.rnw:81-84
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 7, height = 7, pointsize = 18, bg = "white")
xtab<-sample.burposte$xinput
ytab<-sample.burposte$yprod
plot(ytab~xtab, pch=16,col='blue2')
dev.null <- dev.off()
cat("\\includegraphics[width=0.65\\textwidth]{", file, "}\n\n", sep="")


