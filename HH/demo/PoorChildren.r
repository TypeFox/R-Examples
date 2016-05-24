## Poor Children, Working Parents

## source("PoorChildren.r", echo=TRUE)
## R code by Richard Heiberger, December 20, 2011
## Revised June 18, 2013
##
## This file requires the HH package, version 2.3-38 or newer.

require(HH)
if (packageVersion("HH") < "2.3-38") stop("Please update HH to version 2.3-38 or newer.")

try(windows.options(record=TRUE)) ## windows.
## The try() will allow this to go through if a different device is in use.

## Colors taken from NY Times figure using Snagit editor
PCWPpalette <- c("#A8C2C0", "#C9E9E6", "#D0E39A", "#A3B37B")

## for the Forbes online article we needed 640x480 pixels
## http://www.forbes.com/sites/naomirobbins/2011/12/20/alternative-to-charles-blows-figure-in-newts-war-on-poor-children-2/
## try(windows(width=6.6666666666666661, height=5.000)) ## this gives PNG with 640 x 480 pixels
## The try() will allow this to go through if a different device is in use.
## The default is
## windows(width=7, height=7) ## this gives PNG with 672 x 772 pixels

data(PoorChildren)

## > PoorChildren
##            Extremely.poor.NWP Moderately.poor.NWP Moderately.poor.1+WP Extremely.poor.1+WP
## 10 or less             606352              354007               930909              337828
## 10 to 15              1001894              641736              1548697              566705
## 15 to 20               998587              629566              1417142              535242
## 20 to 25               533108              308569               637017              251701
## 25 to 30               406372              242874               457117              197478
## 30 or more             391888              197773               348115              183158
## >


## Figure 2a
PCcount <- likert(PoorChildren, col=PCWPpalette,
                  ylab="Percent of poor households in area",
                  xlab="Number of Children",
                  xlab.top=c("No Working Parents", "One or more Working Parents"),
                  ylab.right="Row Count Totals",
                  main="Poor Children, Working Parents",
                  scales=list(
                    x=list(
                      at=seq(-2,2,1)*1000000,
                      labels=c("2,000,000", "1,000,000","0","1,000,000","2,000,000"))),
                  rightAxisLabels=format(rowSums(PoorChildren[,1:4]), big.mark=","),
                  par.settings=list(
                    axis.line=list(col="transparent"),
                    layout.widths=list(ylab.right=1, axis.key.padding=0)
                    )
                  )
PCcount
## update(PCcount, sub="Figure 2a")

## Population pyramid style, with y-axis in the center
PCcountPP <- likert(PoorChildren, col=PCWPpalette,
                  ylab="Percent of\npoor households\nin area",
                  xlab="Number of Children",
                  xlab.top=c("\n\nNo Working Parents", "\n\nOne or more Working Parents"),
                  ylab.right=NULL,
                  scales=list(x=list(
                      at=seq(-2,2,1)*1000000,
                      labels=c("2,000,000", "1,000,000","0","1,000,000","2,000,000"))),
                  rightAxis=FALSE)
as.pyramidLikert(PCcountPP, panel.width=.44)


## Figure 2b
PCpercent <- likert(PoorChildren, col=PCWPpalette, as.percent=TRUE,
                    ylab="Percent of poor households in area",
                    xlab="Percent of Children",
                    xlab.top=c("No Working Parents", "One or more Working Parents"),
                    ylab.right="Row Count Totals",
                    main="Poor Children, Working Parents",
                    scales=list(
                      x=list(
                        at=seq(-40,60,20),
                        labels=c(40,20,0,20,40,60))),
                    rightAxisLabels=format(rowSums(PoorChildren[,1:4]), big.mark=","))
## PCpercent
PCpercent2 <- update(PCpercent, par.settings=list(axis.line=list(col="transparent")))
PCpercent2
## update(PCpercent2, sub="Figure 2b")


## Figure 2
PL3 <-
likert(as.listOfNamedMatrices(PoorChildren),
       col=PCWPpalette, as.percent=TRUE,
       ylab="Percent of poor households in area",
       xlab="Percent of Children",
       xlab.top=c("No Working Parents", "One or more Working Parents"),
       ylab.right="Row Count Totals",
       main="Poor Children, Working Parents",
       strip=FALSE,
       strip.left=FALSE,
       rightAxisLabels=format(rowSums(PoorChildren[,1:4]), big.mark=","),
       resize.height="rowSums", par.settings=list(axis.line=list(col="transparent")))
PL3
## update(PL3, sub="Figure 2")


## Figure 3
PCpercentAvg <- likert(colSums(PoorChildren[,1:4]), col=PCWPpalette, as.percent=TRUE,
                       xlab="Percent of Children",
                       ylab=list("All Areas", rot=0),
                       xlab.top=c("No Working Parents", "One or more Working Parents"),
                       ylab.right="Total All Children",
                       main="Poor Children, Working Parents",
                       box.ratio=1000,
                       rightAxisLabels=format(sum(PoorChildren[,1:4]), big.mark=","),
                       par.settings=list(axis.line=list(col="transparent"),
                          layout.widths=list(axis.right=1, axis.key.padding=0,
                            ylab.right=1, right.padding=1)))
PCpercentAvg  ## this plot needs at least 8 inch wide plotting device to get the full 13,723,835 visible
## update(PCpercentAvg, sub="Figure 3")


## Figure 4
tmp4 <- colSums(PoorChildren[,c(2,1,3,4)])
ByWP <- rbind(matrix(0,4,2,dimnames=list(letters[1:4],NULL)),
              cbind(NWP=tmp4, '1+WP'=tmp4))
ByWP[cbind(5:8, c(2,2,1,1))] <- 0
ByWP
##
PL4vert <-
likert(as.listOfNamedMatrices(t(ByWP)), as.percent=TRUE, horizontal=FALSE,
       col=c(rep("transparent",4), PCWPpalette[c(2,1,3,4)]),
       strip.left=FALSE, strip=FALSE,
       ylab.right=list(c("Moderately\nPoor", "Extremely\nPoor"), rot=0),
       rightAxisLabels=format(colSums(ByWP), big.mark=","),
       resize.width=1,
       resize.height=colSums(ByWP), ## horizontal=FALSE applies this to width
       layout=c(2,1),
       box.ratio=1000,
       xlab="Percent of Children, All Areas",
       xlab.top=c("No Working Parents", "One or more Working Parents"),
       ylab="Number of Children",
       main="Poor Children, Working Parents",
       par.settings=list(axis.line=list(col="transparent"),
         layout.widths=list(axis.right=0, ylab.right=1)))
## PL4vert
PL4vert$legend$bottom$args$text <- PL4vert$legend$bottom$args$text[c(6,5,7,8)]
PL4vert$legend$bottom$args$rect$col <- PL4vert$legend$bottom$args$rect$col[c(6,5,7,8)]
PL4vert$legend$bottom$args$columns <- 4
## for (i in seq(along=PL4vert$x.limits))
##   PL4vert$x.limits[[i]] <- names(PL4vert$x.limits[[i]])
PL4vert$x.limits <- as.list(format(colSums(ByWP), big.mark=","))
PL4vert
## update(PL4vert, sub="Figure 4")



## Figure 5
tmp5 <- colSums(PoorChildren[,c(2,3,1,4)])
ByPr <- rbind(matrix(0,4,2,dimnames=list(letters[1:4],NULL)),
              cbind(MP=tmp5, EP=tmp5))
ByPr[cbind(5:8, c(2,2,1,1))] <- 0
ByPr
##
PL5vert <-
likert(as.listOfNamedMatrices(t(ByPr)), as.percent=TRUE, horizontal=FALSE,
       col=c(rep("transparent",4), PCWPpalette[c(2,3,1,4)]),
       strip.left=FALSE, strip=FALSE,
       ylab.right=list(c("No Working\nParents", "One or more\nWorking\nParents"), rot=0),
       rightAxisLabels=format(colSums(ByPr), big.mark=","),
       resize.width=1,
       resize.height=colSums(ByPr), ## horizontal=FALSE applies this to width
       layout=c(2,1),
       box.ratio=1000,
       xlab="Percent of Children, All Areas",
       xlab.top=list(c("Moderately Poor", "Extremely Poor"), rot=0),
       ylab="Number of Children",
       main="Poor Children, Working Parents",
       par.settings=list(axis.line=list(col="transparent"),
         layout.widths=list(axis.right=0, ylab.right=1)))
## PL5vert
PL5vert$legend$bottom$args$text <- PL5vert$legend$bottom$args$text[c(7,5,6,8)]
PL5vert$legend$bottom$args$rect$col <- PL5vert$legend$bottom$args$rect$col[c(7,5,6,8)]
PL5vert$legend$bottom$args$columns <- 4
## for (i in seq(along=PL5vert$x.limits))
##   PL5vert$x.limits[[i]] <- names(PL5vert$x.limits[[i]])
PL5vert$x.limits <- as.list(format(colSums(ByPr), big.mark=","))
PL5vert
## update(PL5vert, sub="Figure 5")
