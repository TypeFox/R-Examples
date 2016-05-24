lmatRows <- function(x, focus) {
  UseMethod("lmatRows")
}

lmatRows.mmc.multicomp <-
  if.R(r=
       function(x, focus) {
         lmatRows(x$mca$glht$model, focus)
       }
       ,s=
       function(x, focus) {
         lmatRows(x$mca, focus)
       }
       )

     
lmatRows.multicomp <-
  if.R(r=
       function(x, focus) {
         lmatRows(x$glht$model, focus)
       }
       ,s=
       function(x, focus) {
         lmatRows(eval(x$lmcall), focus)
       }
       )   
  
lmatRows.glht <-
  if.R(r=
       function(x, focus) {
         lmatRows(x$model, focus)
       }
       ,s=
       lmatRows.glht <- function(x, focus) {
         stop("No glht class in S-Plus.")
       }
       )
       
lmatRows.lm <-
  if.R(r=
       function(x, focus) {
         if (missing(focus)) stop("'lmatRows' requires the 'focus' argument.")
         assign <- x$assign
         term.labels <- attr(x$terms, "term.labels")
         ## recover()
         which(term.labels[assign]==focus) + 1 ## allowing for "(Intercept)"
       }
       ,s=
       function(x, focus) {
         if (missing(focus)) stop("'lmatRows' requires the 'focus' argument.")
         terms <- attr(x$terms, "factors")
         levels <- sapply(x$contrasts, nrow)
         nl <- names(levels)
         tl <- terms
         tl[nl,] <- terms[nl,] * levels
         rowcounts <- apply(tl, 2, function(x) prod(x[x>0]))
         if (names(x$assign)[1] == "(Intercept)")
           rowcounts <- c("(Intercept)"=1, rowcounts)
         start <- c(1, cumsum(rev(rev(rowcounts)[-1]))+1)
         lmatRows <- apply(cbind(rowcounts, start=as.vector(start)), 1,
                           function(x) seq(from=x[2], length=x[1]))
         lmatRows[[focus]]
       }
       )

## lmatRows.lme <-
##   if.R(r=
##        function(x, focus) {
##          if (missing(focus)) stop("'lmatRows' requires the 'focus' argument.")
##          which(substring(names(x$coefficients$fixed), 1, nchar(focus)) == focus)
##        }
##        ,s=
##        function(x, focus) {
##          stop("'lmatRows.lme' not yet designed for S-Plus.")
##        }
##        )

## Used in both R and S-Plus
lmatContrast <- function(lmat.none, contrast.matrix) {
  levels <- dimnames(contrast.matrix)[[1]]
  lmat.none[,levels] %*% contrast.matrix
}

## source("~/HH-R.package/HH/R/lmatRows.R")
