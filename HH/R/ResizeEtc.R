ResizeEtc <- function(c.list,
                      condlevelsName,
                      x.same, y.same,
                      layout,
                      strip=TRUE,
                      strip.left=TRUE,
                      strip.values, strip.left.values,
                      strip.par, strip.left.par,  ## only the second is effective when both are specified
                      resize.height, resize.width,
                      main,
                      ## sub,
                      ...) {

  resultA <- c.list
  if (!missing(condlevelsName))
    names(resultA$condlevels) <- condlevelsName ## this argument should be part of c.trellis

  resultB <- resultA
  if (!missing(x.same) && x.same) {
    resultB$x.limits <- range(resultB$x.limits)
    resultB$x.scales$relation <- "same"
  }

  resultC <- resultB
  if (!missing(y.same) && y.same) {
    resultC$y.limits <- range(resultC$y.limits)
    resultC$y.scales$relation <- "same"
  }

  resultD <-
    if (missing(layout))
      resultC
    else
      update(resultC, layout=layout)

  ## both strip and strip.left use the same par.strip.text
  ## I would like them to use par.strip.text and par.strip.left.text, where the
  ## new par.strip.left.text defaults to the current value of par.strip.text
  resultE <- resultD
  if (!missing(strip.values) && strip) {
    resultE <- update(resultE, strip=strip)
    resultE$condlevels[[condlevelsName]] <- strip.values
  }

  resultF <- resultE
  if (!missing(strip.left.values) && strip.left) {
    resultF <- update(resultF, strip.left=strip.left)
    resultF$condlevels[[condlevelsName]] <- strip.left.values
  }

  resultG <-
    if (missing(strip.par))
      resultF
    else
      update(resultF, par.strip.text=strip.par)

  resultH <-
    if (missing(strip.left.par))
      resultG
    else
      update(resultG, par.strip.text=strip.left.par)

  ##  w:width h:height
  wh <- paste(ifelse(missing(resize.width), '', "w"),
              ifelse(missing(resize.height), '', "h"),
              sep="")
  resultI <- switch(wh,
                    wh=resizePanels(resultH,
                      h=resize.height, w=resize.width),
                    w=resizePanels(resultH, w=resize.width+1),
                    h=resizePanels(resultH, h=resize.height+1),
                    resultH)
  resultI <- if (strip)
    resultI
  else
    update(resultI, strip=strip)

  resultJ <- resultI
    if (!missing(main) && !is.null(main)) {
      resultJ <- update(resultJ, main=main)
    }

##  resultK <- resultJ
##    if (!is.null(sub)) {
##      resultK <- update(resultK, sub=sub)
##    }
##
##  resultL <- resultK
  resultL <- resultJ
    if (length(list(...))) {
      resultL <- update(resultL, ...)
    }

  resultL
}

## source("c:/HOME/rmh/HH-R.package/HH/R/ResizeEtc.R")
## environment(ResizeEtc) <- environment(plot.likert)

if (FALSE) {
  ResizeEtc(A + as.layer(B))
  ResizeEtc(A + as.layer(B), x.same=FALSE)
  ResizeEtc(A + as.layer(B), y.same=FALSE)
  ResizeEtc(A + as.layer(B), layout=c(1,2))

  ResizeEtc(A + as.layer(B), strip.values=c("cccc","ddddd"))
  ResizeEtc(A + as.layer(B), strip.left.values=c("cccc","ddddd"))
  ResizeEtc(A + as.layer(B), resize.width=c(1,2))
  ResizeEtc(A + as.layer(B), main="abcd")
}
## environment(ResizeEtc) <- environment(plot.likert)
