ak2k <- function (obj, data, method = NULL, mpower = 2, verbose = 0)
{
    nm <- NULL
    for (i in obj$mstrataID) {
        cat(sprintf("\nmastrata=%d\n",i))
        idx <- which(obj$mstrata == i)
        tmp <- obj$groups[idx]
        tr <- idx[which(tmp == obj$g.names[1])]
        ct <- idx[which(tmp == obj$g.names[2])]
        n.tr <- length(tr)
        n.ct <- length(ct)
        m <- min(n.tr, n.ct)
        n <- n.tr + n.ct
        if (n.tr != n.ct) {
            idx2 <- c(tr, ct)
            if (is.null(method)) {
                mat <- matrix(runif(n * n), n, n)
                colnames(mat) <- rownames(data)[idx2]
                rownames(mat) <- rownames(data)[idx2]
            }
            else {
                mat <- as.matrix(dist(data[idx2, obj$vars], method = method,
                p = mpower))
            }
            m <- min(n.tr, n.ct)
            mat1 <- matrix(mat[1:n.tr, -(1:n.tr)], n.tr, n.ct)
            colnames(mat1) <- colnames(mat)[-(1:n.tr)]
            rownames(mat1) <- rownames(mat)[1:n.tr]
            if(n.tr>0 & n.ct>0){
             if (n.tr > n.ct) {
                for(k in 1:m){
                    print(k)
                    print(mat1)
                 mins <- apply(mat1, 2, function(x) which.min(x))[1]
                 coln <- colnames(mat1)[1]
                 rown <- rownames(mat1)[mins]
                 cat("\nrow-col\n")
                 print(c(length(coln),length(rown)))
                 nm <- c(nm, coln, rown)
                 mat1 <- mat1[-mins,]
                }
             }
             if(n.tr < n.ct) {
                for(k in 1:m){
                    print(k)
                    print(mat1)
                 mins <- apply(mat1, 1, function(x) which.min(x))[1]
                 rown <- rownames(mat1)[1]
                 coln <- colnames(mat1)[mins]
                 cat("\nrow-col2\n")
                 print(c(length(coln),length(rown)))
                 nm <- c(nm, coln, rown)
                 mat1 <- mat1[,mins]
                }
             }
         }
        } else { # n.tr != n.ct
            cat("\nhere\n")
            nm <- c(nm, rownames(obj$X)[c(ct, tr)])
        }
    }
    print(table(nm))
    idx <- match(nm, rownames(obj$X))
    print(idx)
    print(table(idx))
    idx <- idx[which(!is.na(idx))]
    idx <- unique(idx)
    print(table(idx))
    if (length(idx) > 0) {
        obj$matched[-idx] <- FALSE
        obj$mstrata[-idx] <- NA
        obj$w <- numeric(dim(data)[1])
        obj$w[idx] <- 1
        obj$k2k <- TRUE
        obj$tab <- cem.summary(obj = obj, verbose = verbose)
    }
    invisible(obj)
}

k2k <- function(obj, data, method=NULL, mpower=2, verbose=0){
		   nm <- NULL
        if(is.null(obj$X))
         stop("please first run cem() with option keep.all=TRUE")

        for(i in obj$mstrataID){
          idx <- which(obj$mstrata==i)
          #cat(sprintf("\nstrata=%d\n",i))
		  tmp <- obj$groups[idx] 
		  tr <- idx[which(tmp==obj$g.names[1])]
		  ct <- idx[which(tmp==obj$g.names[2])]
          n.tr <- length(tr)
          n.ct <- length(ct)
		  m <- min(n.tr,n.ct)
		  n <- n.tr+n.ct
		  if(n.tr != n.ct){
		   idx2 <- c(tr, ct)
		   if(is.null(method)){
		    mat <- matrix(runif(n*n), n,n)
			colnames(mat) <- rownames(data)[idx2]
			rownames(mat) <- rownames(data)[idx2]
			} else {
		     mat <- as.matrix(dist(data[idx2, obj$vars],method=method, p=mpower))
			}
		   m <- min(n.tr, n.ct)
		   mat1 <- matrix(mat[1:n.tr,-(1:n.tr)], n.tr, n.ct)

		   colnames(mat1) <- colnames(mat)[-(1:n.tr)]
		   rownames(mat1) <- rownames(mat)[1:n.tr]
		   if(n.tr > n.ct){
		    for(k in 1:m){
                #                print(mat1)
             mins <- apply(mat1, 2, function(x) min(x, na.rm=TRUE))
			 min.c <- min(mins, na.rm=TRUE)
			 col <- which(mins == min.c)[1]
			 row <- which(mat1[,col]==min.c)[1]
			 mat1[row, 1:n.ct] <- NA
			 mat1[1:n.tr ,col] <- NA
             #print(c(colnames(mat1)[col], rownames(mat1)[row]))
			 nm <- c(nm, colnames(mat1)[col], rownames(mat1)[row])
			}
		   } else {
		    for(k in 1:m){
                #  print(mat1)
             mins <- apply(mat1, 1, function(x) min(x, na.rm=TRUE))
			 min.r <- min(mins, na.rm=TRUE)
			 row <- which(mins == min.r)[1]
			 col <- which(mat1[row,]==min.r)[1]
			 mat1[row, ] <- NA
			 mat1[ ,col] <- NA
			 nm <- c(nm, colnames(mat1)[col], rownames(mat1)[row]) 
             #print(c(colnames(mat1)[col], rownames(mat1)[row]))
			}
		   }
		  } else {
		   nm <- c(nm, rownames(obj$X)[c(ct,tr)])
		  } 
		}


idx <- match(nm, rownames(obj$X))
idx <- idx[which(!is.na(idx))]
if(length(idx)>0){
 obj$matched[-idx] <- FALSE
 obj$mstrata[-idx] <- NA
 obj$w <- numeric(dim(data)[1])
 obj$w[idx] <- 1
 obj$k2k <- TRUE
 obj$tab <- cem.summary(obj=obj, verbose=verbose)
}
invisible(obj)
}


