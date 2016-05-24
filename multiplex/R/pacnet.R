pacnet <-
function (file, prsep = ", ", toarray = FALSE, uniq = FALSE, 
    transp = FALSE) 
{
    arx <- scan(file, what = "character", nlines = -1, quiet = TRUE)
    if (isTRUE(length(grep("ALPHA", arx, fixed = TRUE)) == 0L) == 
        TRUE) {
        if (isTRUE(length(grep("ATOM", arx, fixed = TRUE)) == 
            0L) == TRUE) 
            stop("Input data does not have 'Generating and Induced Inclusions'")
    }
    ifelse(isTRUE(arx[1] == "FACTORIZATION") == TRUE, X <- arx[(utils::head(grep("Inclusions", 
        arx, fixed = TRUE), 1) + 1L):(utils::tail(grep("THERE", 
        arx, fixed = TRUE), 1) - 1L)], X <- arx[(utils::tail(grep("ALPHA", 
        arx, fixed = TRUE), 1) + 5L):(utils::tail(grep("THERE", 
        arx, fixed = TRUE), 1) - 1L)])
    x <- arx[utils::tail(grep("THERE", arx, fixed = TRUE), 1):length(arx)]
    v <- vector()
    for (i in 1:length(X)) {
        ifelse((isTRUE(as.logical(grep("(", X[i], fixed = TRUE))) == 
            FALSE & isTRUE(as.logical(grep(")", X[i], fixed = TRUE))) == 
            FALSE), v[length(v) + 1L] <- i, NA)
    }
    rm(i)
    m <- v[which(v%%2L == 0L)]
    n <- v[which(v%%2L != 0L)]
    n <- n[2:length(n)]
    rm(v)
    lt <- list()
    for (i in 1:length(m)) {
        if (is.na((n[i] - 1L)) == FALSE) 
            lt[[i]] <- X[(m[i] + 1):(n[i] - 1L)]
    }
    rm(i)
    rm(X)
    rm(n, m)
    ifelse(isTRUE(arx[1] == "FACTORIZATION") == TRUE, n <- as.numeric(arx[(utils::head(grep("Order", 
        arx, fixed = TRUE), 1) + 4L)]), n <- as.numeric(arx[(grep("(Semigroup", 
        arx, fixed = TRUE) - 1L)]))
    a <- grep("Inclusions:", x, fixed = TRUE)
    z <- grep("Atom", x, fixed = TRUE)
    lt2 <- list()
    for (i in 1:length(z)) {
        lt2[[i]] <- x[(a[i] + 1L):(z[i] - 1L)]
    }
    rm(i)
    for (i in 1:length(lt)) {
        lt[[i]] <- gsub("[[:punct:]]", "", lt[[i]])
    }
    rm(i)
    for (i in 1:length(lt2)) {
        lt2[[i]] <- gsub("[[:punct:]]", "", lt2[[i]])
    }
    rm(i)
    ii <- list()
    length(ii) <- length(lt)
    for (i in 1:length(lt)) {
        for (j in 1:(length(lt[[i]])/2L)) {
            ii[[i]] <- append(ii[[i]], paste(lt[[i]][seq(1, length(lt[[i]]), 
                by = 2)[j]], lt[[i]][seq(2, length(lt[[i]]), 
                by = 2)[j]], sep = prsep))
        }
        rm(j)
    }
    rm(i)
    ifelse(isTRUE(uniq == TRUE) == FALSE, NA, ii <- unique(ii))
    at <- list()
    length(at) <- length(lt2)
    for (i in 1:length(lt2)) {
        for (j in 1:(length(lt2[[i]])/2L)) {
            at[[i]] <- append(at[[i]], paste(lt2[[i]][seq(1, 
                length(lt2[[i]]), by = 2)[j]], lt2[[i]][seq(2, 
                length(lt2[[i]]), by = 2)[j]], sep = prsep))
        }
        rm(j)
    }
    rm(i)
    f <- utils::tail(grep("-|-", arx, fixed = TRUE), length(lt2))
    if (isTRUE(as.numeric(x[grep("-|-", x, fixed = TRUE)[1] - 
        1L]) != n) == TRUE) {
        g <- utils::tail(grep("-|-", x, fixed = TRUE))
        p <- as.numeric(x[grep("-|-", x, fixed = TRUE)[1] - 1L])
        q <- n - p
        arr1 <- array(NA, dim = c(n, p))
        arr2 <- array(NA, dim = c(n, q))
        for (i in 1:length(g)) {
            for (j in 1:n) {
                if (isTRUE((i%%2L) == 1L) == TRUE) 
                  arr1[j, ] <- as.numeric(x[(g[i] + 3L + ((j - 
                    1) * p + ((j - 1) * 2L))):(g[i] + (j * p) + 
                    (j * 2L))])
            }
            rm(j)
            for (j in 1:n) {
                if (isTRUE((i%%2L) == 0L) == TRUE) 
                  arr2[j, ] <- as.numeric(x[(g[i] + 3L + ((j - 
                    1) * q + ((j - 1) * 2L))):(g[i] + (j * q) + 
                    (j * 2L))])
            }
            rm(j)
            if (isTRUE((i%%2L) == 0L) == TRUE && isTRUE(i > 1) == 
                TRUE) {
                if (isTRUE(i == 2L) == TRUE) {
                  arr <- data.frame(arr1, arr2)
                }
                else if (isTRUE(i > 2) == TRUE) {
                  arr <- zbind(as.array(as.matrix(arr)), as.array(as.matrix(data.frame(arr1, 
                    arr2))))
                }
            }
        }
        rm(i)
    }
    else if (isTRUE(as.numeric(x[grep("-|-", x, fixed = TRUE)[1] - 
        1L]) == n) == TRUE) {
        arr <- array(NA, dim = c(n, n, length(f)))
        for (i in 1:length(f)) {
            for (j in 1:n) {
                arr[j, , i] <- as.numeric(arx[(f[i] + 3L + ((j - 
                  1) * n + ((j - 1) * 2L))):(f[i] + (j * n) + 
                  (j * 2L))])
            }
            rm(j)
        }
        rm(i)
    }
    if (isTRUE(toarray == FALSE) == TRUE) {
        if (transp) {
            for (i in 1:length(ii)) {
                ii[[i]] <- swp(ii[[i]])
            }
            rm(i)
            for (i in 1:length(at)) {
                at[[i]] <- swp(at[[i]])
            }
            rm(i)
            tmp <- arr
            if (is.na(dim(arr)[3]) == FALSE) {
                for (i in 1:dim(arr)[3]) {
                  arr[, , i] <- t(tmp[, , i])
                }
                rm(i)
            }
            else if (is.na(dim(arr)[3]) == TRUE) {
                arr <- t(tmp)
            }
        }
    }
    if (toarray) {
        iiarrs <- array(dim = c(n, n, length(ii)), dimnames = list(1:n, 
            1:n, 1:length(ii)))
        for (i in 1:length(ii)) {
            ifelse(isTRUE(transp == TRUE) == TRUE, iiarrs[, , 
                i] <- transf(t(ii[[i]]), type = "listmat", ord = n, 
                prsep = prsep), iiarrs[, , i] <- transf(ii[[i]], 
                type = "listmat", ord = n, prsep = prsep))
        }
        rm(i)
        tmp <- data.frame(matrix(ncol = dim(iiarrs)[1] * dim(iiarrs)[2], 
            nrow = 0L))
        for (i in 1:dim(iiarrs)[3]) {
            tmp[i, ] <- as.vector(iiarrs[, , i])
        }
        rm(i)
        iiarr <- array(dim = c(n, n, nrow(unique(tmp))))
        for (i in 1:nrow(unique(tmp))) {
            iiarr[, , i][1:(n * n)] <- as.numeric(unique(tmp)[i, 
                ])
        }
        rm(i)
        dimnames(iiarr)[[1]] <- dimnames(iiarr)[[2]] <- 1:n
        rm(tmp)
        atarrs <- array(dim = c(n, n, length(at)), dimnames = list(1:n, 
            1:n, 1:length(at)))
        for (i in 1:length(at)) {
            ifelse(isTRUE(transp == TRUE) == TRUE, atarrs[, , 
                i] <- transf(t(at[[i]]), type = "listmat", ord = n, 
                prsep = prsep), atarrs[, , i] <- transf(at[[i]], 
                type = "listmat", ord = n, prsep = prsep))
        }
        rm(i)
        tmp <- data.frame(matrix(ncol = dim(atarrs)[1] * dim(atarrs)[2], 
            nrow = 0))
        for (i in 1:dim(atarrs)[3]) {
            tmp[i, ] <- as.vector(atarrs[, , i])
        }
        rm(i)
        atarr <- array(dim = c(n, n, nrow(unique(tmp))))
        for (i in 1:nrow(unique(tmp))) {
            atarr[, , i][1:(n * n)] <- as.numeric(unique(tmp)[i, 
                ])
        }
        rm(i)
        dimnames(atarr)[[1]] <- dimnames(atarr)[[2]] <- 1:n
        rm(tmp)
        if (transp) {
            tmp <- arr
            if (is.na(dim(arr)[3]) == FALSE) {
                for (i in 1:dim(arr)[3]) {
                  arr[, , i] <- t(tmp[, , i])
                }
                rm(i)
            }
            else if (is.na(dim(arr)[3]) == TRUE) {
                arr <- t(tmp)
            }
            tmp <- iiarr
            for (i in 1:dim(iiarr)[3]) {
                iiarr[, , i] <- t(tmp[, , i])
            }
            rm(i)
            tmp <- iiarrs
            for (i in 1:dim(iiarrs)[3]) {
                iiarrs[, , i] <- t(tmp[, , i])
            }
            rm(i)
            tmp <- atarr
            for (i in 1:dim(atarr)[3]) {
                atarr[, , i] <- t(tmp[, , i])
            }
            rm(i)
            tmp <- atarrs
            for (i in 1:dim(atarrs)[3]) {
                atarrs[, , i] <- t(tmp[, , i])
            }
            rm(i)
        }
    }
    dimnames(arr)[[1]] <- dimnames(arr)[[2]] <- 1:n
    if (toarray) {
        ifelse(isTRUE(uniq == TRUE) == TRUE, lst <- list(ii = iiarr, 
            at = atarr, mc = arr), lst <- list(ii = iiarrs, at = atarrs, 
            mc = arr))
        class(lst) <- c("Pacnet", "toarray")
    }
    else {
        lst <- list(ii = ii, at = at, mc = arr)
        class(lst) <- "Pacnet"
    }
    if (uniq) 
        class(lst) <- append(class(lst), "uniq")
    if (transp) 
        class(lst) <- append(class(lst), "transp")
    return(lst)
}
