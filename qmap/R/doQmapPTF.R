doQmapPTF <- function(x,fobj,...){
  if(!any(class(fobj)=="fitQmapPTF"))
    stop("class(fobj) should be fitQmapPTF")  
  UseMethod("doQmapPTF")
}

doQmapPTF.default <- function(x,fobj,...){
### transform 'x' using quantile mapping with
### transfere function (tfun) and and parameters (par)
### specifyed in 'fobj'
### based on the parameters specifyed in 'fobj'.
###
### x : a numeric vector
### fobj : output from fitQmap (class: "fitQmap")
### ... : not used but possibly needed in further development.
  ffun <- fobj$tfun
  funpar <- fobj$par
  wett <-  if(!is.null(fobj$wet.day)){
    x>=fobj$wet.day
  } else {
    rep(TRUE,length(x))
  }
  x[wett] <- do.call(ffun,c(list(x[wett]),as.list(funpar)))
  x[!wett] <- 0
  if(!is.null(fobj$wet.day))
    x[x<0] <- 0
  return(x)
}

doQmapPTF.matrix <- function (x, fobj, ...) {
### doQmap for matrix of time series to be
### transformed
###
### x: matrix with N columns
    if (ncol(x) != nrow(fobj$par))
        stop("'ncol(x)' and 'nrow(fobj$par)' should be eaqual\n")
    NN <- ncol(x)
    hind <- 1:NN
    names(hind) <- colnames(x)
    tfun <- fobj$tfun
    funpar <- fobj$par
    xx <- sapply(hind, function(i) {
        tr <- try({
            xh <- x[, i]
            ## start bug fix
            ## added: 27.04.2016
            wett <- if (!is.null(fobj$wet.day)) {
                xh >= fobj$wet.day[i]
            } else {
                rep(TRUE, length(xh))
            }
            ## wett <- xh >= fobj$wet.day[i]
            ## end bug fix
            xh[wett] <- do.call(tfun, c(list(xh[wett]), as.list(funpar[i,
                ])))
            xh[!wett] <- 0
            if (!is.null(fobj$wet.day))
                xh[xh < 0] <- 0
            xh
        }, silent = TRUE)
        if (class(tr) == "try-error") {
            warning("Quantile mapping for ", names(hind)[i],
                " failed NA's produced.")
            tr <- rep(NA, nrow(x))
        }
        return(tr)
    })
    rownames(xx) <- rownames(x)
    return(xx)
}


doQmapPTF.data.frame <-  function(x,fobj,...){
### doQmap for data.frame of time series to be
### transformed
###
### x: data.frame with N columns
  x <- as.matrix(x)
  x <- doQmapPTF.matrix(x,fobj,...)
  x <- as.data.frame(x)
  return(x)
}
