read.srt <-
function (file, header = TRUE, sep = "\t", attr = FALSE, toarray = TRUE, 
    dichot = FALSE, labels = NULL) 
{
    ifelse(isTRUE(is.data.frame(file) == FALSE) == TRUE, x <- utils::read.table(file, 
        header = header, sep = sep), x <- file)
    if (isTRUE(attr == TRUE) == TRUE) {
        xa <- x
        if (isTRUE(ncol(x) == 2L) == TRUE) {
            x <- as.data.frame(cbind(as.vector(x[, 1]), as.vector(x[, 
                1]), as.vector(x[, 2])))
        }
        else if (isTRUE(ncol(x) > 2L) == TRUE) {
            x <- cbind(x[, 1], x[, 1], x[, 2:ncol(x)])
        }
        attr(x, "names")[1:2] <- c("n", "n")
    }
    if (isTRUE(toarray == TRUE) == TRUE) {
        R <- (ncol(x) - 2L)
        if (R == 0L) 
            stop("You must specify at least one relation.")
        if (is.null(labels) == TRUE) {
            lbs <- unique(c(as.vector(x[, 1]), as.vector(x[, 
                2])))
            if (isTRUE(lbs == "") == TRUE) {
                warning("Node labels in the input are empty!")
                lbs <- 1:nrow(x)
            }
            else {
                NA
            }
            ifelse(isTRUE(R == 1L) == TRUE, MAT <- array(0L, 
                dim = c(length(lbs), length(lbs))), MAT <- array(0L, 
                dim = c(length(lbs), length(lbs), R)))
            dimnames(MAT)[[1]] <- dimnames(MAT)[[2]] <- lbs
        }
        else {
            ifelse(isTRUE(R == 1L) == TRUE, MAT <- array(0L, 
                dim = c(length(labels), length(labels))), MAT <- array(0L, 
                dim = c(length(labels), length(labels), R)))
            dimnames(MAT)[[1]] <- dimnames(MAT)[[2]] <- labels
        }
        ifelse(isTRUE(R == 1L) == FALSE, dimnames(MAT)[[3]] <- attr(x, 
            "names")[3:ncol(x)], NA)
        Dims <- attr(MAT, "dimnames")[[1]]
        if (isTRUE(ncol(x) > 3L) == TRUE) {
            for (r in 3:ncol(x)) {
                rel <- which(x[, r] != 0L)
                rrel <- x[rel, ]
                X <- integer(length(Dims))
                for (i in 1:length(Dims)) {
                  X[i] <- sum(as.numeric(rrel[, 1] == Dims[i]))
                }
                rm(i)
                attr(X, "names") <- Dims
                xx <- vector()
                for (i in 1:length(Dims)) {
                  ifelse(X[i] != 0L, xx[i] <- i, xx[i] <- NA)
                }
                rm(i)
                attr(xx, "names") <- Dims
                xx <- (stats::na.omit(xx))
                xx <- as.vector(attr(xx, "names"))
                nX <- X[which(X > 0L)]
                YY <- vector()
                if (isTRUE(length(xx) > 1L) == TRUE) {
                  for (i in 1:length(xx)) {
                    YY <- rrel[, 2][which(rrel[, 1] == attr(nX, 
                      "names")[i])]
                    if (isTRUE(length(YY) > 1L) == TRUE) {
                      for (j in 1:length(YY)) {
                        tmp <- MAT[(which((as.vector(rownames(MAT)) == 
                          xx[i]), arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                          YY[j]), arr.ind = TRUE)), (r - 2L)]
                        MAT[(which((as.vector(rownames(MAT)) == 
                          xx[i]), arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                          YY[j]), arr.ind = TRUE)), (r - 2L)] <- tmp + 
                          as.numeric(rrel[, r][which(rrel[, 1] == 
                            attr(nX, "names")[i])])[j]
                      }
                      rm(j)
                    }
                    else if (isTRUE(length(YY) == 1L) == TRUE) {
                      tmp <- MAT[(which((as.vector(rownames(MAT)) == 
                        xx[i]), arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                        YY), arr.ind = TRUE)), (r - 2L)]
                      MAT[(which((as.vector(rownames(MAT)) == 
                        xx[i]), arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                        YY), arr.ind = TRUE)), (r - 2L)] <- tmp + 
                        as.numeric(rrel[, r][which(rrel[, 1] == 
                          attr(nX, "names")[i])])
                    }
                  }
                  rm(i)
                  rm(xx, YY)
                }
                else if (isTRUE(length(xx) == 1L) == TRUE) {
                  YY <- rrel[, 2][which(rrel[, 1] == attr(nX, 
                    "names"))]
                  if (isTRUE(length(YY) > 1L) == TRUE) {
                    for (j in 1:length(YY)) {
                      MAT[(which((as.vector(rownames(MAT)) == 
                        xx), arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                        YY[j]), arr.ind = TRUE)), (r - 2L)] <- as.numeric(rrel[, 
                        r][which(rrel[, 1] == attr(nX, "names"))])[j]
                    }
                    rm(j)
                  }
                  else if (isTRUE(length(YY) == 1L) == TRUE) {
                    MAT[(which((as.vector(rownames(MAT)) == xx), 
                      arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                      YY), arr.ind = TRUE)), (r - 2L)] <- as.numeric(rrel[, 
                      r][which(rrel[, 1] == attr(nX, "names"))])
                  }
                  rm(xx, YY)
                }
            }
            rm(r)
        }
        else if (isTRUE(ncol(x) == 3L) == TRUE) {
            rel <- which(x[, 3] != 0L)
            rrel <- x[rel, ]
            X <- integer(length(Dims))
            for (i in 1:length(Dims)) {
                X[i] <- sum(as.numeric(rrel[, 1] == Dims[i]))
            }
            rm(i)
            attr(X, "names") <- Dims
            xx <- vector()
            for (i in 1:length(Dims)) {
                ifelse(X[i] != 0L, xx[i] <- i, xx[i] <- NA)
            }
            rm(i)
            attr(xx, "names") <- Dims
            xx <- (stats::na.omit(xx))
            xx <- as.vector(attr(xx, "names"))
            nX <- X[which(X > 0L)]
            YY <- vector()
            if (isTRUE(length(xx) > 1L) == TRUE) {
                for (i in 1:length(xx)) {
                  YY <- rrel[, 2][which(rrel[, 1] == attr(nX, 
                    "names")[i])]
                  if (isTRUE(length(YY) > 1L) == TRUE) {
                    for (j in 1:length(YY)) {
                      tmp <- MAT[(which((as.vector(rownames(MAT)) == 
                        xx[i]), arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                        YY[j]), arr.ind = TRUE))]
                      MAT[(which((as.vector(rownames(MAT)) == 
                        xx[i]), arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                        YY[j]), arr.ind = TRUE))] <- tmp + as.numeric(factor(rrel[, 
                        3][which(rrel[, 1] == attr(nX, "names")[i])]))[j]
                    }
                    rm(j)
                  }
                  else if (isTRUE(length(YY) == 1L) == TRUE) {
                    tmp <- MAT[(which((as.vector(rownames(MAT)) == 
                      xx[i]), arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                      YY), arr.ind = TRUE))]
                    MAT[(which((as.vector(rownames(MAT)) == xx[i]), 
                      arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                      YY), arr.ind = TRUE))] <- tmp + as.numeric(factor(rrel[, 
                      3][which(rrel[, 1] == attr(nX, "names")[i])]))
                  }
                }
                rm(i)
            }
            else if (isTRUE(length(xx) == 1L) == TRUE) {
                YY <- rrel[, 2][which(rrel[, 1] == attr(nX, "names"))]
                if (isTRUE(length(YY) > 1L) == TRUE) {
                  for (j in 1:length(YY)) {
                    MAT[(which((as.vector(rownames(MAT)) == xx), 
                      arr.ind = TRUE)), (which(as.vector(colnames(MAT) == 
                      YY[j]), arr.ind = TRUE)), 1] <- as.numeric(rrel[, 
                      3][which(rrel[, 1] == attr(nX, "names"))])[j]
                  }
                  rm(j)
                }
                else if (isTRUE(length(YY) == 1L) == TRUE) {
                  MAT[which(dimnames(MAT)[1][[1]] == YY), which(dimnames(MAT)[2][[1]] == 
                    YY)] <- as.numeric(as.vector(rrel[, 3]))
                }
                rm(xx, YY)
            }
        }
        if (isTRUE(dichot == TRUE) == TRUE) {
            MAT <- dichot(MAT)
        }
        return(MAT)
    }
    if (isTRUE(dichot == TRUE) == TRUE) {
        x[, 3:ncol(x)] <- dichot(x[, 3:ncol(x)])
    }
    if (isTRUE(attr == TRUE) == TRUE) {
        rownames(xa) <- xa[, 1]
        ifelse(isTRUE(ncol(xa) == 2L) == TRUE, NA, xa <- xa[, 
            2:ncol(xa)])
        return(xa)
    }
    else {
        return(x)
    }
}
