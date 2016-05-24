## From: VINOD@FORDHAM.EDU
## To: maechler@stat.math.ethz.ch
## X-Spam-Level: *
## Subject: fracdiff in R  does not work for gnp series  "insufficient workspace"
## Date: Sun, 15 May 2005 13:24:46 -0400

## Dear Martin Maechler

## I teach econometrics at Fordham.  For some reason the fracdiff
## does not work for the basic gnp series.

library(fracdiff)

if(FALSE) {
    ##MM library(urca)
    ##MM data(npext)
    data(npext, package = "urca")       # Nelson Plosser data
    ## "bad practice": attach(npext)
    realgnp2 <- npext[50:129, "realgnp"] # to exclude missing data
} else { ## keep test independent:
    realgnp2 <-
        c(4.7604631, 4.7883247, 4.8138091, 4.8690717, 4.8782461, 4.8331023,
          4.8243057, 4.9000761, 4.9067552, 5.0225639, 4.9863426, 4.9416424,
          4.8504665, 4.9972123, 5.1113852, 5.1089712, 5.1896179, 5.2470241,
          5.2459709, 5.2517497, 5.3161573, 5.2122147, 5.1316723, 4.9712012,
          4.9522997, 5.0388988, 5.1328529, 5.2626902, 5.3141907, 5.2621719,
          5.3442463, 5.4258307, 5.5748121, 5.6964221, 5.8203796, 5.8897086,
          5.872681, 5.7449244, 5.7362497, 5.7798172, 5.7810521, 5.8729625,
          5.9490788, 5.9791389, 6.0229632, 6.0088132, 6.0822189, 6.1005431,
          6.1147878, 6.1032295, 6.1652077, 6.1897005, 6.2089924, 6.2724996,
          6.3117348, 6.3649229, 6.4261648, 6.4893569, 6.5150089, 6.5604647,
          6.5857578, 6.5792512, 6.6067, 6.6552832, 6.705961, 6.700553,
          6.687906, 6.7356178, 6.781224, 6.8328012, 6.8572808, 6.8556192,
          6.8747935, 6.8489768, 6.8840768, 6.9496707, 6.9826227, 7.0096668,
          7.0455416, 7.0888837)
}

fr1 <-  fracdiff(realgnp2, nar = 0, nma = 0,  M = 100)

## COMPUTER SAYS
## Error in switch(result$info, stop("insufficient workspace"), stop("error in
## gamma function"),  :
##         insufficient workspace

fr1

## ...

## Hrishikesh D. Vinod
## Professor of Economics, Fordham University
## E-Mail: Vinod@fordham.edu
## Web page:  http://www.fordham.edu/economics/vinod

summary(fr1)
