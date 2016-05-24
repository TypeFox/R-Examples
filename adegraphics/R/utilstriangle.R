## projection dans triangle
## de la base e1=c(1,0,0), e2=c(0,1,0), e3=c(0,0,1)
## a c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)), c(-1/sqrt(2),1/sqrt(2),0), c(-1/sqrt(6),-1/sqrt(6),2/sqrt(6))
.coordtotriangleUnity <- function(mdata3) {
  x <- mdata3[, 1]
  y <- mdata3[, 2]
  z <- mdata3[, 3]
  return(cbind(0, (y - x) / sqrt(2), (2 * z - x - y) / sqrt(6)))
}


## projection depend also on scale defined by min and max
## need to rescale coordinates to maintain distances
## in the new space
.coordtotriangleM <- function(ta, mini3, maxi3) {
  data3d <- t(apply(ta, 1, FUN = function(x) {
    x <- (x - mini3) / (maxi3 - mini3)
    return(x / sum(x))}))
  return(.coordtotriangleUnity(data3d))
}


## TODO: redo this, from ade4
.trranges <- function(df, adjust = TRUE, min3 = NULL, max3 = NULL) {
  ta <- sweep(df, 1, rowSums(df), "/")

  if(ncol(ta) != 3)
    stop("Non convenient data")
  if(min(ta) < 0) 
    stop("Non convenient data")
  if((!is.null(min3)) & (!is.null(max3))) 
    adjust <- TRUE
  cal <- matrix(0, 9, 3)
  tb <- t(apply(ta, 1, FUN = function(x) {x / sum(x)}))
  mini <- apply(tb, 2, min)
  maxi <- apply(tb, 2, max)
  mini <- (floor(mini / 0.1)) / 10
  maxi <- (floor(maxi / 0.1) + 1) / 10
  mini[mini < 0] <- 0
  maxi[maxi > 1] <- 1
  if(!is.null(min3))
    mini <- min3
  if(!is.null(max3))
    maxi <- min3
  ampli <- maxi - mini
  amplim <- max(ampli)
  if(!all(ampli == amplim)) {
    for (j in 1:3) {
      k <- amplim - ampli[j]
      while (k > 0) {
        if((k > 0) & (maxi[j] < 1)) {
          maxi[j] <- maxi[j] + 0.1
          k <- k - 1
        }
        if((k > 0) & (mini[j] > 0)) {
          mini[j] <- mini[j] - 0.1
          k <- k - 1
        }
      }
    }
  }
  cal[1, 1] <- mini[1]
  cal[1, 2] <- mini[2]
  cal[1, 3] <- 1 - cal[1, 1] - cal[1, 2]
  cal[2, 1] <- mini[1]
  cal[2, 2] <- maxi[2]
  cal[2, 3] <- 1 - cal[2, 1] - cal[2, 2]
  cal[3, 1] <- maxi[1]
  cal[3, 2] <- mini[2]
  cal[3, 3] <- 1 - cal[3, 1] - cal[3, 2]
  cal[4, 1] <- mini[1]
  cal[4, 3] <- mini[3]
  cal[4, 2] <- 1 - cal[4, 1] - cal[4, 3]
  cal[5, 1] <- mini[1]
  cal[5, 3] <- maxi[3]
  cal[5, 2] <- 1 - cal[5, 1] - cal[5, 3]
  cal[6, 1] <- maxi[1]
  cal[6, 3] <- mini[3]
  cal[6, 2] <- 1 - cal[6, 1] - cal[6, 3]
  cal[7, 2] <- mini[2]
  cal[7, 3] <- mini[3]
  cal[7, 1] <- 1 - cal[7, 2] - cal[7, 3]
  cal[8, 2] <- mini[2]
  cal[8, 3] <- maxi[3]
  cal[8, 1] <- 1 - cal[8, 2] - cal[8, 3]
  cal[9, 2] <- maxi[2]
  cal[9, 3] <- mini[3]
  cal[9, 1] <- 1 - cal[9, 2] - cal[9, 3]
  mini <- apply(cal, 2, min)
  mini <- round(mini, digits = 4)
  maxi <- apply(cal, 2, max)
  maxi <- round(maxi, digits = 4)
  ampli <- maxi - mini
  if(!adjust) {
    mini <- c(0, 0, 0)
    maxi <- c(1, 1, 1)
  }
  return(list(mini = mini, maxi = maxi))
}


## ## calcul maximum et minimum pour triangle
## ## data as list
## .trranges <- function(data, mini, maxi, adjust){
##   if(is.null(mini))mini <-c(0,0,0)
##   if(is.null(maxi))maxi <-c(1,1,1)
##   if(adjust){
##     if(!is.null(data$frame))
##       ta <- t(apply(eval(data$ta, envir = sys.frame(data$frame)), 1, function(x) x/sum(x)))
##     else
##      ta <- t(apply(data$ta, 1, function(x)x/sum(x)))
##     tb <- t(apply(ta, 1, function(x) x/sum(x)))
##     mini <- apply(tb, 2, min)
##     maxi <- apply(tb, 2, max)
##     mini <- (floor(mini/0.1))/10
##     maxi <- (floor(maxi/0.1) + 1)/10
##     mini[mini < 0] <- 0
##     maxi[maxi > 1] <- 1
##   }
##   ampli <- maxi-mini
##   amplim <- max(ampli)
##   if(! all(ampli == amplim)){#on doit avoir la meme chose.
##     for(i in 1:3){
##       diffv <- amplim -ampli[i]/2
##       mini[i] <- mini[i]-diffv
##       maxi[i] <- maxi[i]+diffv
##       if(mini[i]<0){
##         maxi[i] <- maxi[i]-mini[i]
##         mini[i] <- 0
##        }
##       if(maxi[i]>1){
##         mini[i] <- mini[i]-(maxi[i]-1)
##         maxi[i] <- 1
##       }
##     }
##   }
##   if(any(mini<0) | any(maxi>1))
##     stop("wrong calculus for limits", call.  = FALSE)
##   ##"ici partie cal non reprise. a voir ensuite
##   return(list(mini=mini, maxi=maxi))
## }


.showpos <- function(object) {
  ## from ade4
  mini <- object@g.args$min3d
  maxi <- object@g.args$max3d
  w <- matrix(0, 3, 3)
  w[1, 1] <- mini[1]
  w[1, 2] <- mini[2]
  w[1, 3] <- maxi[3]
  w[2, 1] <- maxi[1]
  w[2, 2] <- mini[2]
  w[2, 3] <- mini[3]
  w[3, 1] <- mini[1]
  w[3, 2] <- maxi[2]
  w[3, 3] <- mini[3]
  smallT <- .coordtotriangleM(matrix(c(0, 0, 1, 1, 0, 0, 0, 1, 0), byrow = TRUE, ncol = 3), mini3 = rep(0, 3), maxi3 = rep(1, 3))[, -1]
  A <- smallT[1, ]
  B <- smallT[2, ]
  C <- smallT[3, ]
  shadowS <- .coordtotriangleM(w, c(0, 0, 0), c(1, 1, 1))[, -1]
  a <- shadowS[1, ]
  b <- shadowS[2, ]
  c <- shadowS[3, ]
  aa <- xyplot(0 ~ 0, xlim = c(-0.7, 0.7), ylim = c(-0.55, 0.9), aspect = "iso", scale = list(draw = FALSE), xlab = NULL, ylab = NULL, par.settings = list(axis.line = list(col = "transparent")),
               panel = function(...) {
                 panel.polygon(c(A[1], B[1], C[1]), c(A[2], B[2], C[2]))
                 panel.polygon(c(a[1], b[1], c[1]), c(a[2], b[2], c[2]), col = grey(0.75))
               })
  invisible(aa)
}
