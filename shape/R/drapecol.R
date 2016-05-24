
##==============================================================================
## drapecol: colors for draping persp (surface) plots
##==============================================================================

drapecol <- function(A, col=femmecol(100), NAcol = "white", lim = NULL) {

  nr <- nrow(A) ; nc <- ncol(A) ; ncol <- length(col)

 ## drape color matrix has one row and one column less than input matrix;
 ## take a weighted average
  AA <- 0.25 * (A[1:(nr-1),1:(nc-1)] +
                A[1:(nr-1),2:nc] +
                A[2:nr,1:(nc-1)] +
                A[2:nr,2:nc])

  if (! is.null(lim)) Ar <- lim
  else Ar <- range(AA, na.rm=TRUE)
  rn <- Ar[2] - Ar[1]

  ifelse (rn != 0, drape <- col[1+trunc((AA-Ar[1])/rn*(ncol-1))] ,
                   drape <- rep(col[1],ncol) )

  drape [is.na(drape)] <- NAcol

  return(drape)

}


