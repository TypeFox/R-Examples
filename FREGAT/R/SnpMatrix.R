# Functions from package 'snpStats' v.1.20.0 (c) 2015 David Clayton

#setMethod("[", signature(x="SnpMatrix",i="ANY",j="ANY",drop="ANY"),
subset.SnpMatrix <-
          function(x, i, j, drop=FALSE) {
            if (drop!=FALSE)
              stop("dimensions cannot be dropped from a SnpMatrix object")
            
            if (missing(i)) i <- integer(0)
            else if (is.character(i)) {
              i <- match(i, rownames(x))
              if (any(is.na(i)))
                stop("No match for one or more row selections")
            }
            else if (is.numeric(i)) {
              if (any(i<0)) 
                  i <- (1:nrow(x))[i]
              if (min(i)<1 || max(i)>nrow(x)) 
                stop("One or more row selections out of range")
            }
            else if (is.logical(i)) {
              if (length(i)!=nrow(x))
                stop("logical row selection vector is wrong length")
              i <- (1:nrow(x))[i]
            }
            else
              stop("Illegal type for row selection, i")
            i <- as.integer(i)
            if (any(is.na(i)))
              stop("NAs in row selection vector")
            
            if (missing(j)) j <- integer(0)
            else if (is.character(j)) {
              j <- match(j, colnames(x))
              if (any(is.na(j)))
                stop("No match for one or more column selections")
            }
            else if (is.numeric(j)) {
              if (any(j<0))
                j <- (1:ncol(x))[j]
              if (min(j)<1 || max(j)>ncol(x)) 
                stop("One or more column selections out of range")
            }
            else if (is.logical(j)) {
              if (length(j)!=ncol(x))
                stop("logical column selection vector is wrong length")
              j <- (1:ncol(x))[j]
            }
            else
              stop("Illegal type for column selection, j")
            j <- as.integer(j)
            if (any(is.na(j)))
              stop("NAs in column selection vector")
            
            .Call("subset", x, i, j, PACKAGE="FREGAT")
          }

#setAs("SnpMatrix", "numeric",
SnpMatrix2numeric <-
      function(from) {
        .Call("asnum", from, PACKAGE="FREGAT")
      }
