semigroup <-
function (x, type = c("numerical", "symbolic"), labels = NULL, 
    cmp = FALSE, smpl = FALSE) 
{
    if (is.array(x) == FALSE) 
        stop("Data must be a stacked array of square matrices.")
    if (is.na(dim(x)[3]) == FALSE) {
        if (is.null(labels) == FALSE) {
            dimnames(x)[[3]] <- labels
        }
        else {
            if (is.null(dimnames(x)[[3]]) == TRUE) 
                dimnames(x)[[3]] <- 1:dim(x)[3]
        }
    }
    x <- replace(x, x < 1, 0)
    x <- replace(x, x >= 1, 1)
    if (smpl == TRUE) {
        ifelse(is.null(dimnames(x)[[3]]) == FALSE, lbs <- dimnames(x)[[3]], 
            lbs <- 1:dim(x)[3])
        if (is.null(dimnames(x)[[3]]) == FALSE) {
            nlb <- list()
            for (i in 1:length(lbs)) {
                nlb[i] <- lbs[i]
            }
            rm(i)
            for (i in 1:length(nlb)) {
                lbs[i] <- (strsplit(nlb[[i]], "")[[1]][1])
            }
            dimnames(x)[[3]] <- lbs
        }
    }
    if (is.na(dim(x)[3]) == TRUE) {
        s0 <- data.frame(matrix(ncol = 1, nrow = 1))
        if (isTRUE(all.equal(replace(x %*% x, x %*% x >= 1, 1), 
            x) == TRUE)) 
            s0[1, 1] <- 1
        Bx <- array(dim = c(dim(x)[1], dim(x)[2], 2))
        Bx[, , 1] <- as.matrix(x)
        Bx[, , 2] <- replace(x %*% x, x %*% x >= 1, 1)
    }
    if (is.na(dim(x)[3]) == FALSE) {
        tmpo <- data.frame(matrix(ncol = (dim(x)[1] * dim(x)[2]), 
            nrow = 0))
        for (i in 1:dim(x)[3]) {
            ifelse(isTRUE(dim(x)[3] > 1) == TRUE, tmpo[i, ] <- as.vector(x[, 
                , i]), tmpo <- as.vector(x))
        }
        rm(i)
        if (isTRUE(is.character(dimnames(x)[[3]]) == TRUE) == 
            TRUE) 
            dimnames(x)[[3]][which(duplicated(dimnames(x)[[3]]))] <- 1:length(which(duplicated(dimnames(x)[[3]])))
        if (isTRUE(is.null(dim(tmpo)) == FALSE) == TRUE) 
            rownames(tmpo) <- dimnames(x)[[3]]
        tmpu <- unique(tmpo)
        if (isTRUE(dim(x)[3] < 2) == TRUE) 
            x <- array(tmpo, c(dim(x)[1], dim(x)[2]))
        if (isTRUE(dim(x)[3] > 1) == TRUE) {
            tmp <- array(dim = c(dim(x)[1], dim(x)[2], nrow(tmpu)))
            for (i in 1:nrow(tmpu)) {
                tmp[, , i][1:(dim(x)[1] * dim(x)[2])] <- as.numeric(tmpu[i, 
                  ])
            }
            rm(i)
            if (is.null(dimnames(tmp)[[1]]) == FALSE) 
                dimnames(tmp)[[3]] <- rownames(tmpu)
            if (is.null(dimnames(x)[[1]]) == FALSE) 
                dimnames(tmp)[[1]] <- dimnames(tmp)[[2]] <- dimnames(x)[[1]]
            x <- tmp
            dimnames(x)[[3]] <- as.list(rownames(tmpu))
        }
        rm(tmp)
        s0 <- data.frame(matrix(ncol = dim(x)[3], nrow = dim(x)[3]))
        for (k in 1:dim(x)[3]) {
            for (j in 1:dim(x)[3]) {
                tmp <- x[, , j] %*% x[, , k]
                tmp <- replace(tmp, tmp >= 1, 1)
                for (i in dim(x)[3]:1) {
                  if (isTRUE(all.equal(tmp, x[, , i]) == TRUE)) 
                    s0[j, k] <- i
                }
            }
        }
        rm(i, j, k)
        dimnames(s0)[[1]] <- 1:dim(x)[3]
        dimnames(s0)[[2]] <- 1:dim(x)[3]
        if (sum(as.numeric(is.na(s0))) == 0) 
            Bx <- x
        if (sum(as.numeric(is.na(s0))) > 0) {
            Bx <- array(dim = c(dim(x)[1], dim(x)[2], 0))
            for (i in 1:nrow(s0)) {
                for (j in 1:length(which(is.na(s0[i, ])))) {
                  if (length(which(is.na(s0[i, ]))) > 0) 
                    Bx <- zbnd(Bx, (replace(x[, , i] %*% x[, 
                      , which(is.na(s0[i, ]))[j]], x[, , i] %*% 
                      x[, , which(is.na(s0[i, ]))[j]] >= 1, 1)))
                }
            }
            rm(i, j)
            tmp <- data.frame(matrix(ncol = (dim(x)[1] * dim(x)[2]), 
                nrow = 0))
            for (i in 1:dim(Bx)[3]) {
                tmp[i, ] <- as.vector(Bx[, , i])
            }
            rm(i)
            xBx <- array(dim = c(dim(x)[1], dim(x)[2], nrow(unique(tmp))))
            for (i in 1:nrow(unique(tmp))) {
                xBx[, , i][1:(dim(Bx)[1] * dim(Bx)[2])] <- as.numeric(unique(tmp)[i, 
                  ])
            }
            rm(i)
            if (is.null(dimnames(xBx)) == FALSE) 
                dimnames(xBx)[[3]] <- (dim(x)[3] + 1):(dim(xBx)[3] + 
                  dim(x)[3])
            Bx <- zbnd(x, xBx)
            rm(xBx, tmp)
        }
    }
    while (sum(as.numeric(is.na(s0))) > 0) {
        BBx <- Bx
        for (i in 1:nrow(s0)) {
            for (j in 1:length(which(is.na(s0[i, ])))) {
                if (length(which(is.na(s0[i, ]))) > 0) 
                  BBx <- zbnd(BBx, (replace(Bx[, , i] %*% Bx[, 
                    , which(is.na(s0[i, ]))[j]], Bx[, , i] %*% 
                    Bx[, , which(is.na(s0[i, ]))[j]] >= 1, 1)))
            }
        }
        rm(i, j)
        tmp <- data.frame(matrix(ncol = (dim(Bx)[1] * dim(Bx)[2]), 
            nrow = 0))
        for (i in 1:dim(BBx)[3]) {
            tmp[i, ] <- as.vector(BBx[, , i])
        }
        rm(i)
        Bx <- array(dim = c(dim(x)[1], dim(x)[2], nrow(unique(tmp))))
        for (i in 1:nrow(unique(tmp))) {
            Bx[, , i][1:(dim(BBx)[1] * dim(BBx)[2])] <- as.numeric(unique(tmp)[i, 
                ])
        }
        rm(i)
        rm(tmp, BBx)
        if (is.na(dim(x)[3]) == TRUE) {
            s0 <- data.frame(matrix(ncol = 1, nrow = dim(Bx)[3]))
            for (j in 1:dim(Bx)[3]) {
                tmp <- Bx[, , j] %*% Bx[, , 1]
                tmp <- replace(tmp, tmp >= 1, 1)
                for (i in dim(Bx)[3]:1) {
                  if (isTRUE(all.equal(tmp, Bx[, , i]) == TRUE)) 
                    s0[j, 1] <- i
                }
            }
            rm(i, j)
        }
        if (is.na(dim(x)[3]) == FALSE) {
            s0 <- data.frame(matrix(ncol = dim(x)[3], nrow = dim(Bx)[3]))
            for (k in 1:dim(x)[3]) {
                for (j in 1:dim(Bx)[3]) {
                  tmp <- Bx[, , j] %*% Bx[, , k]
                  tmp <- replace(tmp, tmp >= 1, 1)
                  for (i in dim(Bx)[3]:1) {
                    if (isTRUE(all.equal(tmp, Bx[, , i]) == TRUE)) 
                      s0[j, k] <- i
                  }
                }
            }
            rm(i, j, k)
        }
    }
    ifelse(isTRUE(is.na(dim(x)[3])) == TRUE, dimnames(s0)[[2]] <- 1, 
        dimnames(s0)[[2]] <- 1:dim(x)[3])
    tmpO <- data.frame(matrix(ncol = (dim(Bx)[1] * dim(Bx)[2]), 
        nrow = 0))
    for (i in 1:dim(Bx)[3]) {
        tmpO[i, ] <- as.vector(Bx[, , i])
    }
    rm(i)
    rownames(tmpO) <- dimnames(Bx)[[3]]
    tmpU <- unique(tmpO)
    tmp <- array(dim = c(dim(Bx)[1], dim(Bx)[2], nrow(tmpU)))
    for (i in 1:nrow(tmpU)) {
        tmp[, , i][1:(dim(Bx)[1] * dim(Bx)[2])] <- as.numeric(tmpU[i, 
            ])
    }
    rm(i)
    Bx <- tmp
    dimnames(Bx)[[3]] <- as.list(rownames(tmpU))
    E <- s0
    rm(s0)
    if (dim(Bx)[3] == ncol(E)) {
        W <- rbind(cbind(1:ncol(E), NA, NA))
        colnames(W) <- c("", "n", "g")
    }
    if (dim(Bx)[3] > ncol(E)) {
        tmp <- (ncol(E) + 1):dim(Bx)[3]
        z <- vector()
        for (i in 1:length(tmp)) {
            z[i] <- which(t(E) == tmp[i])[1]
        }
        rm(i)
        g <- vector()
        for (i in 1:length(tmp)) {
            ifelse(z[i]%%ncol(E) == 0, g[i] <- ncol(E), g[i] <- z[i]%%ncol(E))
        }
        rm(i)
        n <- vector()
        for (i in 1:length(tmp)) {
            ifelse(z[i]%%ncol(E) == 0, n[i] <- (z[i]%/%ncol(E)), 
                n[i] <- ((z[i]%/%ncol(E)) + 1))
        }
        rm(i)
        W <- rbind(cbind(1:ncol(E), NA, NA), cbind(((ncol(E) + 
            1):nrow(E)), n, g))
        rm(z, n, g)
    }
    rm(tmp)
    if (is.na(dim(x)[3]) == TRUE) {
        ifelse(is.null(dimnames(Bx)[[3]]) == TRUE, lbl <- 1:dim(Bx)[3], 
            lbl <- dimnames(Bx)[[3]])
    }
    if (is.na(dim(x)[3]) == FALSE) {
        if (is.null(dimnames(x)[[3]]) == TRUE) 
            lbl <- 1:dim(Bx)[3]
        if (is.null(dimnames(x)[[3]]) == FALSE) {
            if (isTRUE(dim(Bx)[3] == dim(x)[3]) == TRUE) 
                lbl <- dimnames(x)[[3]]
            if (isTRUE(dim(Bx)[3] < dim(x)[3]) == TRUE) 
                lbl <- rownames(tmpO)
            if (isTRUE(dim(Bx)[3] > dim(x)[3]) == TRUE) 
                lbl <- c(dimnames(x)[[3]], (dim(x)[3] + 1):dim(Bx)[3])
        }
        dimnames(Bx)[[3]] <- lbl
        for (i in which(W[, 2] < (ncol(E) + 1))) {
            lbl[(i)] <- paste(dimnames(Bx)[[3]][W[, 2][i]], dimnames(Bx)[[3]][W[, 
                3][i]], sep = "")
        }
        rm(i)
        for (i in which(W[, 2] < (which(W[, 2] > ncol(E))[1]))[(length(which(W[, 
            2] < (ncol(E) + 1))) + 1):length(which(W[, 2] < (which(W[, 
            2] > ncol(E))[1])))]) {
            lbl[(i)] <- paste(dimnames(Bx)[[3]][W[W[, 2][i], 
                ][2]], dimnames(Bx)[[3]][W[W[, 2][i], ][3]], 
                dimnames(Bx)[[3]][W[, 3][i]], sep = "")
        }
        rm(i)
        for (i in which(W[, 2] < (which(W[, 2] >= (which(W[, 
            2] > ncol(E)))[1])[1]))[(length(which(W[, 2] < (which(W[, 
            2] > ncol(E))[1]))) + 1):length(which(W[, 2] < (which(W[, 
            2] >= (which(W[, 2] > ncol(E)))[1])[1])))]) {
            lbl[(i)] <- paste(dimnames(Bx)[[3]][W[W[W[W[, 2][i], 
                ], ][1, ][2], ][2]], dimnames(Bx)[[3]][W[W[W[W[, 
                2][i], ], ][1, ][2], ][3]], dimnames(Bx)[[3]][W[W[, 
                2][i], ][3]], dimnames(Bx)[[3]][W[, 3][i]], sep = "")
        }
        rm(i)
        for (i in which(W[, 2] < (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] > ncol(E))[1])))[1])[1]))[(length(which(W[, 
            2] < (which(W[, 2] >= (which(W[, 2] > ncol(E)))[1])[1]))) + 
            1):length(which(W[, 2] < (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] > ncol(E))[1])))[1])[1])))]) {
            lbl[(i)] <- paste(dimnames(Bx)[[3]][W[W[W[W[W[, 2][i], 
                ], ][1, ][2], ][2], ][2]], dimnames(Bx)[[3]][W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[, 
                2][i], ], ][1, ][2], ][3]], dimnames(Bx)[[3]][W[W[, 
                2][i], ][3]], dimnames(Bx)[[3]][W[, 3][i]], sep = "")
        }
        rm(i)
        for (i in which(W[, 2] < (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] > ncol(E))[1]))[1])))[1])[1]))[(length(which(W[, 
            2] < (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] > ncol(E))[1])))[1])[1]))) + 1):length(which(W[, 
            2] < (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] > ncol(E))[1]))[1])))[1])[1])))]) {
            lbl[(i)] <- paste(dimnames(Bx)[[3]][W[W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][2], ][2]], dimnames(Bx)[[3]][W[W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[, 
                2][i], ], ][1, ][2], ][3]], dimnames(Bx)[[3]][W[W[, 
                2][i], ][3]], dimnames(Bx)[[3]][W[, 3][i]], sep = "")
        }
        rm(i)
        for (i in which(W[, 2] < (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] > ncol(E))[1]))[1]))[1])))[1])[1]))[(length(which(W[, 
            2] < (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] > ncol(E))[1]))[1])))[1])[1]))) + 
            1):length(which(W[, 2] < (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] > ncol(E))[1]))[1]))[1])))[1])[1])))]) {
            lbl[(i)] <- paste(dimnames(Bx)[[3]][W[W[W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][2], ][2], ][2]], 
                dimnames(Bx)[[3]][W[W[W[W[W[W[W[, 2][i], ], ][1, 
                  ][2], ][2], ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[W[W[, 
                  2][i], ], ][1, ][2], ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[W[, 
                  2][i], ], ][1, ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[, 
                  2][i], ], ][1, ][2], ][3]], dimnames(Bx)[[3]][W[W[, 
                  2][i], ][3]], dimnames(Bx)[[3]][W[, 3][i]], 
                sep = "")
        }
        rm(i)
        for (i in which(W[, 2] < (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] > ncol(E))[1]))[1]))[1]))[1])))[1])[1]))[(length(which(W[, 
            2] < (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] > ncol(E))[1]))[1]))[1])))[1])[1]))) + 
            1):length(which(W[, 2] < (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] > ncol(E))[1]))[1]))[1]))[1])))[1])[1])))]) {
            lbl[(i)] <- paste(dimnames(Bx)[[3]][W[W[W[W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][2], ][2], ][2], 
                ][2]], dimnames(Bx)[[3]][W[W[W[W[W[W[W[W[, 2][i], 
                ], ][1, ][2], ][2], ][2], ][2], ][2], ][3]], 
                dimnames(Bx)[[3]][W[W[W[W[W[W[W[, 2][i], ], ][1, 
                  ][2], ][2], ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[W[W[, 
                  2][i], ], ][1, ][2], ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[W[, 
                  2][i], ], ][1, ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[, 
                  2][i], ], ][1, ][2], ][3]], dimnames(Bx)[[3]][W[W[, 
                  2][i], ][3]], dimnames(Bx)[[3]][W[, 3][i]], 
                sep = "")
        }
        rm(i)
        for (i in which(W[, 2] < (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] > ncol(E))[1]))[1]))[1]))[1]))[1])))[1])[1]))[(length(which(W[, 
            2] < (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] > ncol(E))[1]))[1]))[1]))[1])))[1])[1]))) + 1):length(which(W[, 
            2] < (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] >= (which(W[, 2] >= (which(W[, 
            2] >= (which(W[, 2] > ncol(E))[1]))[1]))[1]))[1]))[1])))[1])[1])))]) {
            lbl[(i)] <- paste(dimnames(Bx)[[3]][W[W[W[W[W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][2], ][2], ][2], 
                ][2], ][2]], dimnames(Bx)[[3]][W[W[W[W[W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][2], ][2], ][2], 
                ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][2], ][2], ][2], 
                ][3]], dimnames(Bx)[[3]][W[W[W[W[W[W[W[, 2][i], 
                ], ][1, ][2], ][2], ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[W[, 
                2][i], ], ][1, ][2], ][2], ][3]], dimnames(Bx)[[3]][W[W[W[W[, 
                2][i], ], ][1, ][2], ][3]], dimnames(Bx)[[3]][W[W[, 
                2][i], ][3]], dimnames(Bx)[[3]][W[, 3][i]], sep = "")
        }
        rm(i)
        dimnames(Bx)[[3]] <- lbl
    }
    if (is.null(dimnames(x)[[1]]) == FALSE) 
        dimnames(Bx)[[1]] <- dimnames(Bx)[[2]] <- dimnames(x)[[1]]
    switch(match.arg(type), numerical = S <- cbind(E, data.frame(matrix(ncol = (nrow(E) - 
        ncol(E)), nrow = nrow(E)))), symbolic = {
        e <- E
        for (i in 1:nrow(E)) {
            for (j in 1:ncol(E)) {
                e[i, j] <- lbl[as.matrix(E[i, j])]
            }
            rm(j)
        }
        rm(i)
        S <- cbind(as.data.frame(e), data.frame(matrix(ncol = (nrow(E) - 
            ncol(E)), nrow = nrow(E))))
        rm(e)
    })
    if (ncol(S) != ncol(E)) {
        for (j in (ncol(E) + 1):nrow(E)) {
            for (i in 1:nrow(E)) {
                tmp <- replace(Bx[, , i] %*% (Bx[, , as.numeric(which(E == 
                  j, arr.ind = TRUE)[1, ][1])] %*% Bx[, , as.numeric(which(E == 
                  j, arr.ind = TRUE)[1, ][2])]), Bx[, , i] %*% 
                  (Bx[, , as.numeric(which(E == j, arr.ind = TRUE)[1, 
                    ][1])] %*% Bx[, , as.numeric(which(E == j, 
                    arr.ind = TRUE)[1, ][2])]) >= 1, 1)
                for (k in 1:dim(Bx)[3]) {
                  if (isTRUE(all.equal(Bx[, , k], tmp)) == TRUE) {
                    switch(match.arg(type), numerical = S[i, 
                      j] <- k, symbolic = S[i, j] <- lbl[k])
                    break
                  }
                }
                rm(k)
            }
            rm(j)
        }
        rm(i)
        rm(tmp)
    }
    else {
        NA
    }
    switch(match.arg(type), numerical = dimnames(S)[[1]] <- dimnames(S)[[2]] <- 1:dim(Bx)[3], 
        symbolic = dimnames(S)[[1]] <- dimnames(S)[[2]] <- lbl)
    switch(match.arg(type), numerical = S <- as.matrix(S), symbolic = NA)
    if (isTRUE(dim(x)[3] < dim(Bx)[3]) == TRUE) 
        Bx <- Bx[, , (dim(x)[3] + 1):dim(Bx)[3]]
    if (is.na(dim(x)[3]) == FALSE) {
        if (isTRUE(nrow(tmpo) == nrow(tmpu)) == TRUE) {
            ifelse(cmp == TRUE, lst <- list(dim = dim(x)[1], 
                gens = x, cmps = Bx, ord = nrow(S), st = lbl, 
                S = S), lst <- list(dim = dim(x)[1], gens = x, 
                ord = nrow(S), st = lbl, S = S))
        }
        else {
            ifelse(cmp == TRUE, lst <- list(dim = dim(x)[1], 
                gens = x, cmps = Bx, ord = nrow(S), Note = paste("Relation ", 
                  dimnames(x)[[3]][which(duplicated(tmpo) == 
                    TRUE)], " is repeated in the input data and has been equated", 
                  sep = "'"), st = lbl, S = S), lst <- list(dim = dim(x)[1], 
                gens = x, ord = nrow(S), Note = paste("Relation ", 
                  dimnames(x)[[3]][which(duplicated(tmpo) == 
                    TRUE)], " is repeated in the input data and has been equated", 
                  sep = "'"), st = lbl, S = S))
        }
    }
    else if (is.na(dim(x)[3]) == TRUE) {
        ifelse(cmp == TRUE, lst <- list(dim = dim(x)[1], gens = x, 
            cmps = Bx[, , 2:dim(Bx)[3]], ord = nrow(S), st = lbl, 
            S = S), lst <- list(dim = dim(x)[1], gens = x, ord = nrow(S), 
            st = lbl, S = S))
    }
    class(lst) <- c("Semigroup", match.arg(type))
    return(lst)
}
