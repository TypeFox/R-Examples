partial.order <-
function (x, type = c("strings", "galois"), labels = NULL) 
{
    if (match.arg(type) == "strings") {
        if (isTRUE(attr(x, "class") == "Strings") == TRUE) {
            x <- x$wt
            if (is.array(x) == FALSE) 
                stop("Data must be a stacked array of square matrices if a product of 'strings'.")
            if (is.na(dim(x)[3]) == FALSE) {
                tmp <- data.frame(matrix(ncol = (dim(x)[1] * 
                  dim(x)[2]), nrow = 0))
                for (i in 1:dim(x)[3]) {
                  tmp[i, ] <- as.vector(x[, , i])
                }
                rm(i)
                po <- as.matrix(array(0, dim = c(dim(x)[3], dim(x)[3])))
                for (j in 1:dim(x)[3]) {
                  for (i in 1:dim(x)[3]) {
                    if ((as.numeric(any(tmp[i, ] < tmp[j, ])) == 
                      1 && as.numeric(any(tmp[j, ] < tmp[i, ])) == 
                      0) | as.numeric(all(tmp[i, ] == tmp[j, 
                      ])) == 1) 
                      po[i, j] <- 1
                  }
                }
                rm(i, j)
                rownames(po) <- colnames(po) <- dimnames(x)[[3]]
            }
            else if (is.na(dim(x)[3]) == TRUE) {
                po <- 1
            }
        }
        else if (isTRUE(attr(x, "class") == "Strings") == FALSE) {
            stop("\"x\" should be an object of a \"Strings\" class.")
        }
    }
    if (match.arg(type) == "galois") {
        if (isTRUE(attr(x, "class")[1] == "Galois") == TRUE) {
            if (isTRUE(attr(x, "class")[2] == "full") == TRUE) {
                po <- matrix(0, nrow = length(x), ncol = length(x))
                for (j in 1:length(x)) {
                  for (i in 1:length(x)) {
                    ifelse(isTRUE(all(dhc(x[[i]]) %in% dhc(x[[j]]))) == 
                      TRUE, po[i, j] <- 1, NA)
                  }
                }
                rm(i, j)
            }
            else if (isTRUE(attr(x, "class")[2] == "reduced") == 
                TRUE) {
                po <- matrix(0, nrow = length(x$full), ncol = length(x$full))
                for (j in 1:length(x$full)) {
                  for (i in 1:length(x$full)) {
                    ifelse(isTRUE(all(dhc(x$full[[i]]) %in% dhc(x$full[[j]]))) == 
                      TRUE, po[i, j] <- 1, NA)
                  }
                }
                rm(i, j)
            }
            lb <- list()
            if (isTRUE(attr(x, "class")[2] == "full") == TRUE) {
                length(lb) <- length(x)
                for (i in 1:length(x)) {
                  if (isTRUE(is.na(attr(x, "names")[i])) == FALSE) {
                    lb[[i]] <- paste(paste("{", x[[i]], sep = ""), 
                      paste(attr(x, "names")[i], "}", sep = ""), 
                      sep = "} {")
                  }
                  else {
                    lb[[i]] <- paste(paste("{", x[[i]], sep = ""), 
                      paste(" ", "}", sep = ""), sep = "} {")
                  }
                }
                rm(i)
                colnames(po) <- rownames(po) <- lb
            }
            else if (isTRUE(attr(x, "class")[2] == "reduced") == 
                TRUE) {
                length(lb) <- length(x$reduc)
                for (i in 1:length(x$reduc)) {
                  if (isTRUE(is.na(attr(x$reduc, "names")[i])) == 
                    FALSE) {
                    lb[[i]] <- paste(paste("{", attr(x$reduc, 
                      "names")[i], sep = ""), paste(x$reduc[[i]], 
                      "}", sep = ""), sep = "} {")
                  }
                  else {
                    lb[[i]] <- paste(paste("{", " ", sep = ""), 
                      paste(x$reduc[[i]], "}", sep = ""), sep = "} {")
                  }
                }
                rm(i)
                colnames(po) <- rownames(po) <- lb
            }
            cp <- which(dimnames(po)[[1]] == "{} {}")
            dimnames(po)[[1]][cp] <- dimnames(po)[[2]][cp] <- cp
        }
        else if (isTRUE(attr(x, "class")[1] == "Galois") == FALSE) {
            stop("\"x\" should be an object of a \"Galois\" class.")
        }
    }
    if (is.null(labels) == FALSE) 
        dimnames(po)[[2]] <- dimnames(po)[[1]] <- labels
    po
}
