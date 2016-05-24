as.3dpoints <- function (...) 
{
    nv <- nargs()
    fargs <- list(...)
    if (nv == 3) {
        l1 <- length(fargs[[1]])
        l2 <- length(fargs[[2]])
        l3 <- length(fargs[[3]])
        if ((l1 == l2) & (l1 == l3)) {
            pts <- cbind(x=fargs[[1]], y=fargs[[2]], t=fargs[[3]])
        }
        else {
            stop("Cannot make points from different length vectors")
        }
    }
    else {
        if (nv == 1) {
            if (is.list(fargs[[1]])) {
                fargs <- fargs[[1]]
                if (any(names(fargs) == "x") & any(names(fargs) == 
                  "y") & any(names(fargs) == "t")) {
                  arx <- fargs$x
                  ary <- fargs$y
			art <- fargs$t
                  if ((length(arx) != length(ary)) | (length(arx) != length(art))) {
                    stop("Cannot make points from different length x, y and t list components!")
                  }
                  else {
                    pts <- cbind(x=arx, y=ary, t=art)
                  }
                }
                else {
                  stop("Cannot make points from list without x and y components.")
                }
            }
            else {
                if (is.3dpoints(fargs[[1]])) 
                  {
				pts <- fargs[[1]]
				colnames(pts)=c("x","y","t")	
			}
			
                else stop("Cannot make points from this object")
            }
        }
        else {
            stop("Cannot make object into points!")
        }
    }
    x=sort(pts[,3],index.return=TRUE)
    pts=pts[x$ix,]
    oldClass(pts)="stpp"
    pts
}
