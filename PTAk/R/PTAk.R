#PTAk  is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. (see file LICENSE)
#    see also <http://www.gnu.org/licenses/>. 
#   Dr Didier G.Leibovici CC 2001-20010-2012

"howtoPTAk" <-
function()
{
	print(packageDescription("PTAk"))
		cat( "\n")
		cat("********************", "\n",
	"          Copyright GPL >=2 2000, 2007, 2010, 2012 Didier Leibovici" , "\n",
    "            see the citation file","\n",
    "            and for a good introduction ","\n",
    "           Leibovici, D.G. (2010) JSS:34(10), www.jstatsoft.org/v34/i10/","\n",
    "          contact me at c3s2i@free.fr","\n")
cat("         see some examples on http://c3s2i.free.fr","\n")
   
}
"CANDPARA" <-
function (X, dim = 3, test = 1e-08, Maxiter = 1000, smoothing = FALSE,
    smoo = list(NA), verbose = getOption("verbose"), file = NULL, 
    modesnam = NULL, addedcomment = "")
{
   datanam <- substitute(X)
    sym <- NULL
    if (is.list(X)) {
        if (is.list(X$met)) 
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        diX <- length(dim(X$data))
        for (d in 1:diX) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]], 
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:diX)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:diX)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:diX, laperm))
            }
            else X$data <- X$data * sqrt(X$met[[d]])
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
    }
    debtime <- proc.time()
    pass. <- function(a, r) {
        pasta <- a
        for (i in 2:r) pasta <- paste(pasta, a, sep = "")
        return(pasta)
    }
    if (verbose) {
        cat("\n", "       ----------+++++++++++------------", 
            "\n", ifelse(smoothing, paste("Smoothed ", "\n"), 
                ""), "              PARAFAC/CANDECOMP ", "\n", 
            file = ifelse(is.null(file), "", file), append = TRUE)
        cat("       ----------+++++++++++------------", "\n", 
            file = ifelse(is.null(file), "", file), append = TRUE)
        cat(" Data is ... ", deparse(datanam), "...", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  .... Tensor of order ", length(dim(X)), file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  ....  with dimensions: ", dim(X), "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        if (!is.null(modesnam)) 
            cat("modes are ", modesnam, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
        if (metrics) 
            cat("---Analysis with non-Identity metrics  ------", 
                "\n", file = ifelse(is.null(file), "", file), 
                append = TRUE)
        if (!addedcomment == "") 
            cat("\n", addedcomment, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
    }
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    ord <- length(dim(X))
    if (ord <=2) {
        stop(paste("--- X must be an array  of k > 2 entries! ---"))
    }
    if (is.null(modesnam)) {
        modesnam <- paste(rep("mo", ord), 1:ord)
    }
    if (smoothing) {
        if (length(smoo) < ord) 
            smoo <- rep(list(smoo[[1]]), ord)
    }
    else smoo <- list(NULL)
    sval0 <- INITIA(X, modesnam = modesnam, method = "svd", dim = dim)
    test0 <- 1
    atest <- 0
    sval <- sval0
    iter <- 0
    if (smoothing) {
        for (j in 1:ord) if (!is.list(smoo[[j]])) 
            smoo[[j]] <- list(smoo[[j]])
        for (a in 2:dim) {
            for (j in 1:ord) if (length(smoo[[j]]) < a) 
                smoo[[j]][[a]] <- smoo[[j]][[a - 1]]
        }
    }
    while (test0 > test) {
        iter <- iter + 1
        if (verbose & iter%%100 == 1) 
            cat("\n", " ----------- iteration-", iter, file = ifelse(is.null(file), 
                "", file), append = TRUE)
        for (i in 1:ord) {
            if (iter == 1) {
                if (verbose) 
                  cat("\n", i, "^", sval0[[i]]$d, file = ifelse(is.null(file), 
                    "", file), append = TRUE)
                sval[[i]]$d <- NULL
            }
        }
        if (i == 1) {
            tzz <- 1
            Z <- 1
        }
        else {
            ifelse(dim == 1, tzz <- 1, tzz <- sval[[1]]$v %*% 
                t(sval[[1]]$v))
            Z <- t(sval[[1]]$v)
        }
        for (j in 2:ord) {
            if (!j == i) {
                if (dim > 1) 
                  tzz <- tzz * (sval[[j]]$v %*% t(sval[[j]]$v))
                Z <- RaoProd(t(sval[[j]]$v), Z)
            }
        }
        sval[[i]]$v <- t(matrix(aperm(X, c(i, (1:length(dim(X)))[-i])), 
            nrow = dim(X)[i]) %*% Z %*% Ginv(tzz))
        if (smoothing) {
            for (a in 1:dim) if (is.function(smoo[[i]][[a]])) 
                sval[[i]]$v[a, ] <- smoo[[i]][[a]](sval[[i]]$v[a, 
                  ])
        }
        sval[[i]]$d <- sqrt(diag(sval[[i]]$v %*% t(sval[[i]]$v)))
        sval[[i]]$v <- sval[[i]]$v/sval[[i]]$d
        atest <- atest + sum((sval[[i]]$v - sval0[[i]]$v)^2)
        if (!is.null(sym)) {
            for (i in ord:1) {
                if (!i == sym[i]) 
                  sval[[sym[i]]] <- sval[[i]]
            }
        }
        sval0 <- sval
        test0 <- sqrt(atest)
        atest <- 0
        if (verbose & (iter%%100) == 1) 
            cat("\n", "----------- test =         ", test0, "\n", 
                file = ifelse(is.null(file), "", file), append = TRUE)
        if (iter > (Maxiter - 1) & (iter - Maxiter)%%Maxiter == 0) {
            cat("\n \n \n \n \n ", " WARNING ****** Iteration already =  ", 
                iter, "\n")
            cat(" ** type anything to STOP ** just RETURN to carry on **", 
                "\n")
            conti <- scan("", what = "", n = 1, quiet = TRUE, 
                flush = TRUE, )
            if (length(conti) > 0) 
                stop(paste(" ---- Aborted by request ---- "))
        }
    }
    pourRR2 <- function() {
        tens <- t(sval[[1]]$v) %*% diag(sval[[ord]]$d)
        for (r in 2:ord) {
            tens <- RaoProd(t(sval[[r]]$v), tens)
        }
        return(summary(lm(as.vector(X) ~ tens - 1)))
    }
    pass. <- function(a, r) {
        pasta <- a
        for (i in 2:r) pasta <- paste(pasta, a, sep = "")
        return(pasta)
    }
    PCnam <- paste("v", pass.(1:dim, ord), sep = "")
    ssX <- sum(X^2)
    sstens <- (sval[[i]]$d^2)
    PCT <- 100 * sstens/ssX
    sval[[i]]$lm <- pourRR2()
    if (verbose) {
        cat(" --------optimisation  done ", "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat(" --------Final iteration----", iter, "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat(" ----------- test =         ", test0, "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("\n", " --Norms-- ", sval[[i]]$d, "\n", " --Percent-- ", 
            PCT)
        cat("\n", " ---Total R2 ", sval[[i]]$lm$r.squared * 100, 
            "%", "\n")
    }
    cat("-----Execution Time-----", (proc.time() - debtime)[3], 
        "\n")
    if (metrics) {
        for (d in 1:length(sval)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  sval[[d]]$v <- sval[[d]]$v %*% Powmat(met[[d]], 
                    -1/2)
                }
                else {
                  sval[[d]]$v <- t(1/sqrt(met[[d]]) * t(sval[[d]]$v))
                }
            }
            else sval[[d]]$v <- sval[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    sval[[i]]$pct <- as.vector(PCT)
    sval[[i]]$ssX <- as.vector(ssX)
    sval[[i]]$vsnam <- PCnam
    sval[[i]]$datanam <- datanam
    sval[[i]]$method <- match.call()
    class(sval) <- c("CANDPARA", "PTAk")
    invisible(return(sval))
}
"PCAn" <-
function (X, dim = c(2, 2, 2, 3), test = 1e-12, Maxiter = 400, 
    smoothing = FALSE, smoo = list(NA), verbose = getOption("verbose"), 
    file = NULL, modesnam = NULL, addedcomment = "") 
{
    datanam <- substitute(X)
    sym <- NULL
    if (is.list(X)) {
        if (is.list(X$met)) 
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        diX <- length(dim(X$data))
        for (d in 1:diX) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]], 
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:diX)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:diX)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:diX, laperm))
            }
            else X$data <- X$data * sqrt(X$met[[d]])
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
    }
    debtime <- proc.time()
    pass. <- function(a, r) {
        pasta <- a
        for (i in 2:r) pasta <- paste(pasta, a, sep = "")
        return(pasta)
    }
    if (verbose) {
        cat("----------+++++++++++------------", "\n", ifelse(smoothing, 
            paste("Smoothed ", "\n"), ""), " PCA-n modes  ", 
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
        cat(" Data is ... ", deparse(datanam), "...", "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  .... Tensor of order ", length(dim(X)), file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  ....  with dimensions: ", dim(X), "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        if (!is.null(modesnam)) 
            cat("modes are ", modesnam, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
        if (metrics) 
            cat("---Analysis with non-Identity metrics  ------", 
                "\n", file = ifelse(is.null(file), "", file), 
                append = TRUE)
        if (!addedcomment == "") 
            cat("\n", addedcomment, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
    }
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    solutions <- NULL
    ord <- length(dim(X))
    if (ord <=2) {
        stop(paste("--- X must be an array  of k > 2 entries! ---"))
    }
    if (smoothing) {
        if (length(smoo) < ord) 
            smoo <- rep(list(smoo[[1]]), ord)
    }
    else smoo <- list(NULL)
    if (is.null(modesnam)) {
        modesnam <- paste(rep("mo", ord), 1:ord)
    }
    if (!length(dim) == ord) 
        stop(" Wrong length for dim argument (= Rank-spaces !)")
    for (j in 1:ord) if (dim[j] > dim(X)[j]) 
        stop(" (dim argument) some Rank-spaces are too big!")
    sval0 <- INITIA(X, modesnam = modesnam, method = "svd", dim = dim)
    test0 <- 1
    atest <- 0
    sval <- sval0
    iter <- 0
    if (smoothing) {
        for (j in 1:ord) {
            if (!is.list(smoo[[j]])) 
                smoo[[j]] <- list(smoo[[j]])
            for (a in 2:dim[j]) if (length(smoo[[j]]) < a) 
                smoo[[j]][[a]] <- smoo[[j]][[a - 1]]
        }
    }
    while (test0 > test) {
        iter <- iter + 1
        if (verbose & iter%%100 == 1) {
            cat(" ----------- iteration-", iter, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
        }
        for (i in 1:ord) {
            if (iter == 1) {
                if (verbose) {
                  cat(" ", i, "^", sval0[[i]]$d, file = ifelse(is.null(file), 
                    "", file), append = TRUE)
                }
            }
            sval[[i]]$d <- NULL
        }
        corematv <- X
        for (j in 1:ord) {
            if (j < i) 
                corematv <- CONTRACTION(corematv, sval[[j]]$v, 
                  Xwiz = 1, zwiX = 2)
            if (j > i) 
                corematv <- CONTRACTION(corematv, sval[[j]]$v, 
                  Xwiz = 2, zwiX = 2)
        }
        corematv <- matrix(corematv, nrow = dim(X)[i])
        if (smoothing) 
            svdcormatv <- svdsmooth(corematv, nomb = dim[i], 
                smooth = list(NA, smoo[[i]]))
        else svdcormatv <- svd(corematv)
        sval[[i]]$dopt <- svdcormatv$d[1:dim[i]]
        sval[[i]]$v <- t(svdcormatv$u[, 1:dim[i]])
        if (all(svdcormatv$u[, 1] < 0)) 
            sval[[i]]$v <- -sval[[i]]$v
        coremat <- array(t(corematv) %*% t(sval[[i]]$v), c(dim[-i], 
            dim[i]))
        atest <- atest + sum((sval[[i]]$v - sval0[[i]]$v)^2)
        if (!is.null(sym)) {
            for (i in ord:1) {
                if (!i == sym[i]) 
                  sval[[sym[i]]] <- sval[[i]]
            }
        }
        sval0 <- sval
        if (verbose & iter%%100 == 0) {
            cat(" --", coremat, file = ifelse(is.null(file), 
                "", file), append = TRUE)
        }
        test0 <- sqrt(atest)
        atest <- 0
        if (verbose & (iter%%100) == 1) {
            cat("\n", "----------- test =         ", test0, "\n", 
                file = ifelse(is.null(file), "", file), append = TRUE)
        }
        if (iter > (Maxiter - 1) & (iter - Maxiter)%%Maxiter == 0) {
            cat("\n \n \n \n \n ", " WARNING ****** Iteration already =  ", 
                iter, "\n")
            cat(" ** type anything to STOP ** just RETURN to carry on **", 
                "\n")
            conti <- scan("", what = "", n = 1, quiet = TRUE, 
                flush = TRUE, )
            if (length(conti) > 0) 
                stop(paste(" ---- Aborted by request ---- "))
        }
    }
    PCnam <- outer(outer(1:dim[1], 1:dim[2], FUN = "paste", sep = ""), 
        1:dim[3], FUN = "paste", sep = "")
    if (ord > 3) {
        for (t in 4:ord) {
            PCnam <- outer(PCnam, 1:dim[q], FUN = "paste", sep = "")
        }
    }
    PCnam <- paste("v", PCnam, sep = "")
    ssX <- sum(X^2)
    sstens <- sum(coremat^2)
    totPCT <- 100 * sstens/ssX
    if (verbose) {
        cat(" --------optimisation  done ", "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat(" --------Final iteration----", iter, "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat(" ----------- test =         ", test0, "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("\n", " --Core Matrix-- ", coremat, "\n", " --  Percent -- ", 
            totPCT, "%", "\n")
    }
    if (metrics) {
        for (d in 1:length(sval)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  sval[[d]]$v <- sval[[d]]$v %*% Powmat(met[[d]], 
                    -1/2)
                }
                else {
                  sval[[d]]$v <- t(1/sqrt(met[[d]]) * t(sval[[d]]$v))
                }
            }
            else sval[[d]]$v <- sval[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    cat("-----Execution Time-----", (proc.time() - debtime)[3], 
        "\n")
    sval[[i]]$d <- as.vector(coremat)
    sval[[i]]$coremat <- coremat
    sval[[i]]$pct <- as.vector(totPCT)
    sval[[i]]$ssX <- as.vector(ssX)
    sval[[i]]$vsnam <- PCnam
    sval[[i]]$datanam <- datanam
    sval[[i]]$method <- match.call()
    class(sval) <- c("PCAn", "PTAk")
    invisible(return(sval))
}
"REBUILDPCAn" <-
function (solu) 
{
    lo <- length(solu)
    recon <- t(solu[[1]]$v)
    for (k in 2:lo) {
        recon <- recon %o% t(solu[[k]]$v)
    }
    reconf <- CONTRACTION(recon, solu[[lo]]$coremat, Xwiz = (1:(lo - 
        1)) * 2, zwiX = 1:(lo - 1))
    cat("\n", "--- RMSerror---", sqrt(mean((eval(solu[[lo]]$datanam) - 
        reconf)^2)))
    invisible(return(reconf))
}
"CONTRACTION" <-
function (X, z, Xwiz = NULL, zwiX = NULL, rezwiX = FALSE, usetensor = TRUE)
{
    if (usetensor) {
        if (is.null(zwiX)) {
            if (is.vector(z))
                zwiX <- 1
            else zwiX <- 1:length(dim(z))
        }
        if (is.null(Xwiz)) {
            if (is.vector(X))
                Xwiz <- 1
            else Xwiz <- 1:length(dim(X))
        }
        return(tensor(X, z, Xwiz, zwiX))
    }
    else {
        non <- function(A, awi) {
            (1:length(dim(A)))[!(1:length(dim(A))) %in% awi]
        }
        zbX <- FALSE
        if (length(dim(as.array(X))) < length(dim(as.array(z)))) {
            zbX <- TRUE
            temp <- X
            X <- z
            z <- temp
            temp <- Xwiz
            Xwiz <- zwiX
            zwiX <- temp
        }
        if (is.vector(z)) {
            zwiX <- 1
            zz <- z
            lacolz <- NULL
            lacolaz <- NULL
            if (is.null(Xwiz))
                Xwiz <- (1:length(dim(X)))[dim(X) %in% length(z)]
        }
        else {
            if (is.null(zwiX)) {
                zwiX <- 1:length(dim(z))
            }
            if (is.null(Xwiz))
                Xwiz <- match(dim(z)[zwiX], dim(X))
            Xwiz <- Xwiz[!is.na(Xwiz)]
            if (rezwiX)
                zwiX <- match(dim(X)[Xwiz], dim(z))
            zwiX <- zwiX[!is.na(zwiX)]
            if (!all(dim(X)[Xwiz] == dim(z)[zwiX]))
                stop(paste(" @@@@@ WRONG matching for contraction!@@@@@"))
            if (all(dim(z) %in% dim(X)[Xwiz]) & length(dim(z)) <
                length(dim(X)[Xwiz])) {
                zz <- as.vector(aperm(z, zwiX))
                lacolz <- NULL
                lacolaz <- NULL
            }
            else {
                czwiX <- non(z, zwiX)
                lacolaz <- (1:length(dim(z)))[czwiX]
                lacolz <- (dim(z))[lacolaz]
                zz <- matrix(aperm(z, c(zwiX, lacolaz)), ncol = prod(lacolz))
            }
        }
        lacola <- (1:length(dim(X)))[non(X, Xwiz)]
        laperm <- c(Xwiz, lacola)
        lacol <- (dim(X))[lacola]
        toconz <- matrix(aperm(X, laperm), ncol = prod(lacol))
        Xz <- t(toconz) %*% zz
        dinam <- function(A) {
            namA <- rep(paste(1:length(dim(A))), dim(A))
            dinamA <- list(namA[1:dim(A)[1]])
            for (e in 2:length(dim(A))) {
                dinamA <- c(dinamA, list(namA[(sum(dim(A)[1:(e -
                  1)]) + 1):sum(dim(A)[1:e])]))
            }
            return(dinamA)
        }
        if (!is.null(dimnames(X)))
            dimnamX <- dimnames(X)
        else dimnamX <- dinam(X)
        if (!is.vector(z)) {
            if (!is.null(dimnames(z)))
                dimnamz <- dimnames(z)
            else dimnamz <- dinam(z)
        }
        if (is.null(lacolaz))
            ladim <- dimnamX[lacola]
        else ladim <- c(dimnamX[lacola], dimnamz[lacolaz])
        Xz <- array(Xz, c(lacol, lacolz), dimnames = ladim)
        if (zbX)
            Xz <- aperm(Xz, c((1:length(dim(Xz)))[-lacola], lacola))
        return(Xz)
    }
}
"CONTRACTION.list" <-
function (X, zlist, moins = 1, zwiX = NULL, usetensor = TRUE,
    withapply = FALSE)
{
    mplu <- 1
    lz <- length(zlist)
    for (tu in 1:lz) {
        if (!tu %in% moins) {
            if (withapply)
                X <- apply(X, (1:length(dim(X)))[-mplu], FUN = function(x) as.vector(x) %*%
                  zlist[[tu]]$v)
            else X <- CONTRACTION(X, zlist[[tu]]$v, Xwiz = mplu,
                zwiX = zwiX[tu], usetensor = usetensor)
        }
        else mplu <- mplu + 1
    }
    return(X)
}
"INITIA" <-
function (X, modesnam = NULL, method = "Presvd", dim = 1, ...)
{
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    VV <- list(NULL)
    if (is.null(modesnam))
        modesnam <- paste(rep("m", length(dim(X))), 1:length(dim(X)))
    if (!is.function(method) && method == "Presvd")
        dim <- 1
    if (length(dim) == 1)
        dim <- rep(dim, length(dim(X)))
    for (i in 1:length(dim(X))) {
        cci <- (1:length(dim(X)))[-i]
        if (is.function(method))
            VV[[i]] <- method(matrix(aperm(X, c(cci, i)), ncol = dim(X)[i]),
                ...)
        else {
            if (method == "Presvd")
                VV[[i]] <- PPMA(matrix(aperm(X, c(cci, i)), ncol = dim(X)[i]),
                  pena = list(NULL, NULL))
            if (method == "svd")
                VV[[i]] <- svd(matrix(aperm(X, c(cci, i)), ncol = dim(X)[i]))
        }
        if (dim[i] > dim(X)[i])
            dimi <- dim(X)[i]
        else dimi <- dim[i]
        VV[[i]]$d <- VV[[i]]$d[1:dimi]
        VV[[i]]$modesnam <- modesnam[[i]]
        VV[[i]]$n <- dimnames(X)[[i]]
        VV[[i]]$v <- t(VV[[i]]$v[, 1:dimi])
        if (dimi == 1)
            VV[[i]]$v <- as.vector(VV[[i]]$v)
        VV[[i]]$u <- NULL
    }
    return(VV)
}
"PROJOT" <-
function (X, solu, numo = 1, bortho = TRUE, Ortho = TRUE, metrics = NULL)
{
    txDy <- function(x, D, y) {
        if (!is.null(D)) {
            if (is.vector(D))
                y <- D * y
            if (is.matrix(D))
                y <- D %*% y
        }
        if (is.vector(x) & is.vector(y))
            return(sum(x * y))
        else return(t(x) %*% y)
    }
    projmat <- function(Y, x, D = NULL, bortho = TRUE) {
        if (is.vector(x))
            Y <- x %*% (1/txDy(x, D, x) * txDy(x, D, Y))
        if (is.matrix(x)) {
            if (!bortho)
                Y <- x %*% Powmat(txDy(x, D, x), -1) %*% txDy(x,
                  D, Y)
            else Y <- x %*% ((1/diag(txDy(x, D, x))) * txDy(x,
                D, Y))
        }
        return(Y)
    }
    ldx <- length(dim(X))
    if (!is.list(numo))
        numo <- rep(list(numo), ldx)
    if (!is.list(bortho))
        bortho <- rep(list(bortho), ldx)
    if (!is.list(Ortho))
        Ortho <- rep(list(Ortho), ldx)
    if (!is.list(metrics))
        metrics <- rep(list(metrics), ldx)
    for (i in 1:ldx) {
        if (!is.null(numo[[i]])) {
            z <- solu[[i]]$v
            if (is.matrix(z)) {
                if (!dim(X)[i] == dim(z)[2])
                  stop("----WRONG DIMENSIONS----")
                else z <- z[numo[[i]], ]
            }
            else {
                if (!dim(X)[i] == length(z))
                  stop("----WRONG DIMENSIONS----")
            }
            lacola <- (1:length(dim(X)))[-i]
            laperm <- c(i, lacola)
            lacol <- (dim(X))[lacola]
            toconz <- matrix(aperm(X, laperm), ncol = prod(lacol))
            if (!is.vector(z))
                z <- t(z)
            PXz <- projmat(toconz, z, D = metrics[[i]], bortho = bortho[[i]])
            PXz <- array(PXz, c(dim(X)[i], lacol))
            if (i == ldx)
                PXz <- aperm(PXz, c(2:ldx, 1))
            if (!i == 1 & !i == ldx)
                PXz <- aperm(PXz, c(2:(i), 1, (i + 1):ldx))
            dimnames(PXz) <- dimnames(X)
            if (Ortho[[i]])
                X <- X - PXz
            else X <- PXz
        }
    }
    return(X)
}
"REBUILD" <-
function (solutions, nTens = 1:2, testvar = 1, redundancy = FALSE)
{
    if (!is.list(solutions)) {
        stop(" should be a solutions object see PTA3")
    }
    ord <- length(solutions)
    if (as.character(solutions[[ord]]$method)[1] == "PCA")
        REBUILDPCAn(solutions)
    else {
        tensfin <- 0
        deja <- NULL
        dejaTP <- NULL
        testpass <- length(nTens)
        for (cp in nTens) {
            if (100 * (solutions[[ord]]$d[cp]^2/solutions[[ord]]$ssX[1]) >
                testvar) {
                if (!solutions[[ord]]$d[cp] %in% deja || (redundancy &
                  substr(solutions[[ord]]$vsnam[cp], 2, 1) ==
                    "t")) {
                  if (!substr(solutions[[ord]]$vsnam[cp], 1,
                    1) == "*") {
                    deja <- c(deja, solutions[[ord]]$d[cp])
                    dejaTP <- c(dejaTP, cp)
                  }
                  tens <- solutions[[1]]$v[cp, ] * solutions[[ord]]$d[cp]
                  names(tens) <- solutions[[1]]$n
                  for (d in 2:ord) {
                    atens <- solutions[[d]]$v[cp, ]
                    names(atens) <- solutions[[d]]$n
                    tens <- tens %o% atens
                  }
                  tensfin <- tensfin + tens
                }
            }
        }
        pcre <- 100 * sum(deja^2)/solutions[[ord]]$ssX[1]
        cat("-- Variance Percent rebuilt", solutions[[ord]]$datanam,
            " at ", pcre, "% ", "\n")
        cat("-- MSE ", mean((eval(solutions[[ord]]$datanam) -
            tensfin)^2), "\n")
        cat("-- with ", length(deja), " Principal Tensors out of ",
            length(nTens), " given", "\n")
        if (pcre > 100) {
            cat("\n", "--WARNING !--- redundancy in choice of solutions to rebuild !!!",
                "\n")
            print(pcre, digits = 20)
        }
        comp <- 100 - 100 * (sum(dim(eval(solutions[[ord]]$datanam))) +
            1) * length(deja)/prod(dim(eval(solutions[[ord]]$datanam)))
        cat("-- compression    ", comp, " %", "\n")
        if (comp < 0) {
            cat("******no compression ....", "\n")
            dejadedans <- cbind(dejaTP, deja)
            rownames(dejadedans) <- solutions[[ord]]$vsnam[dejaTP]
            print(dejadedans)
        }
        return(tensfin)
    }
}

"RESUM" <-
function (solb, sola = NULL, numass = NULL, verbose = getOption("verbose"),
    file = NULL, summary = FALSE, testvar = 0.1, with=TRUE)
{
    if (!is.null(sola)) {
        numlast <- length(sola[[length(sola)]]$d)
        if (is.null(numass))
            num <- numlast
        if (!is.null(numass))
            num <- numass
        for (i in 1:length(solb)) {
            for (j in 1:length(sola)) {
                if (as.vector(solb[[i]]$modesnam) == as.vector(sola[[j]]$modesnam)) {
                  sola[[j]]$v <- rbind(sola[[j]]$v, solb[[i]]$v)
                  if ("iter" %in% names(solb[[i]]))
                    sola[[j]]$iter <- c(sola[[j]]$iter, solb[[i]]$iter)
                  if ("test" %in% names(solb[[i]]))
                    sola[[j]]$test <- c(sola[[j]]$test, solb[[i]]$test)
                }
            }
        }
        for (k in 1:length(sola)) {
            if (is.matrix(sola[[k]]$v)) {
                if ((dim(sola[[k]]$v)[[1]]) == numlast) {
                  sola[[k]]$v <- rbind(sola[[k]]$v, rep(1, length(solb[[length(solb)]]$d)) %x%
                    t(sola[[k]]$v[num, ]))
                }
            }
            if (!is.matrix(sola[[k]]$v)) {
                sola[[k]]$v <- rbind(sola[[k]]$v, rep(1, length(solb[[length(solb)]]$d)) %x%
                  t(sola[[k]]$v))
            }
        }
        sola[[k]]$d <- c(sola[[k]]$d, solb[[i]]$d)
        sola[[k]]$pct <- c(sola[[k]]$pct, solb[[i]]$pct)
        sola[[k]]$ssX <- c(sola[[k]]$ssX, solb[[i]]$ssX)
        if ("smoocheck" %in% names(solb[[i]]))
            sola[[k]]$smoocheck <- cbind(sola[[k]]$smoocheck,
                solb[[i]]$smoocheck)
        for (m in 1:length(solb[[i]]$vsnam)) {
            for (n in 1:length(sola[[k]]$vsnam)) {
                if ((round(sola[[k]]$d[n], digits = 10) == round(solb[[i]]$d[m],
                  digits = 10)) & (!substr(solb[[i]]$vsnam[m],
                  1, 1) == "*")) {
                  solb[[i]]$vsnam[m] <- paste("*t", solb[[i]]$vsnam[m],
                    sep = "")
                }
            }
        }
        sola[[k]]$vsnam <- c(sola[[k]]$vsnam, solb[[i]]$vsnam)
    }
    else {
        sola <- solb
        k <- length(sola)
    }
    pctota <- (100 * (sola[[k]]$d)^2)/sola[[k]]$ssX[1]
    if (verbose & !summary) {
        cat("                ------Percent Rebuilt from Selected ----",
            sum(pctota[!substr(sola[[k]]$vsnam, 1, 1) == "*" &
                pctota > testvar]), "%", "\n")
    }
    if (!is.null(file)) {
        cat("               ------Percent Rebuilt from Selected ----",
            sum(pctota[!substr(sola[[k]]$vsnam, 1, 1) == "*" &
                pctota > testvar]), "%", "\n", file = file, append = TRUE)
        if (verbose) {
            sink(file = file, append = TRUE)
            summ <- as.matrix(cbind(1:length(sola[[k]]$d), sola[[k]]$d,
                sola[[k]]$ssX, sola[[k]]$pct, pctota))
            dimnames(summ) <- list(sola[[k]]$vsnam, c("-no-",
                "--Sing Val--", "--ssX--", "--local Pct--", "--Global Pct--"))
            summ <- summ[pctota > testvar, ]
            print(summ, digits = 5)
            sink()
        }
    }
    if (summary) {
         if(class(sola)[1]=="PCAn")  cat("\n", "++++ PCA- ", k, "modes ++++ ", "\n","summary function not available yet using summary.PTAk!","\n","Core tensor taken like Sing Val!")   
         else {
         if(class(sola)[1]=="CANDPARA")  cat("\n", "++++ CANDECOMP/PARAFAC- ", k, "modes ++++ ", "\n","\n")
           else cat("\n", "++++ PTA- ", k, "modes ++++ ", "\n")
           }
        di <- NULL
        for (r in 1:length(sola)) di <- c(di, length(sola[[r]]$v[1,
            ]))
        nostar <- !substr(sola[[k]]$vsnam, 1, 1) == "*"
        cat("               data= ", deparse(sola[[k]]$datanam),
            " ", di, "\n")
        if(sola[[k]]$addedcomment!="")cat("   ", sola[[k]]$addedcomment, "\n")
        cat("                ------Percent Rebuilt----", sum(pctota[nostar]),
            "%", "\n")
        summ <- matrix(cbind(1:length(sola[[k]]$d), sola[[k]]$d,
            sola[[k]]$ssX, sola[[k]]$pct, pctota), ncol = 5)
        summ <- summ[(pctota > testvar) & with, ]
        summ <- matrix(summ, ncol = 5)
        dimnames(summ) <- list(sola[[k]]$vsnam[pctota > testvar & with], c("-no-", "--Sing Val--", "--ssX--", "--local Pct--",
            "--Global Pct--"))
        sumex <- sum(summ[!substr(dimnames(summ)[[1]], 1, 1) ==
            "*", 5])
        cat("                ------Percent Rebuilt from Selected ----",
            sumex, "%", "\n")
        print(summ, digits = 5)
        cat("\n", "++++               ++++", "\n")
        if (!is.null(testvar) & !testvar == 0)
            cat(" Shown are selected over ", length(sola[[k]]$vsnam[nostar]),
                " PT  with var>", testvar, "% total", "\n")
        else cat(" over ", length(sola[[k]]$vsnam[nostar]), " PT ",
            "\n")
       invisible(summ)
    }
    else invisible(sola)
}
"TENSELE" <-
function (T, moins = NULL, asarray = TRUE, order = NULL, id = NULL)
{
    dim <- NULL
    if (is.null(order))
        order <- length(T):1
    if (is.list(id))
        asarray <- TRUE
    vu <- 0
    for (i in order) {
        if (!(i %in% moins)) {
            if (is.null(id[[i]]))
                Tv <- T[[i]]$v
            else Tv <- T[[i]]$v[id[[i]], ]
            if (vu == 0)
                tensel <- Tv
            if (!vu == 0) {
                if (asarray) {
                  tensel <- Tv %o% tensel
                }
                if (!asarray) {
                  tensel <- as.vector(tensel %x% Tv)
                  dim <- c(length(Tv), dim)
                }
            }
            vu <- 1
        }
    }
    return(tensel)
}
"summary.FCAk" <-
function (object, testvar = 0.5, dontshow = "*",...)
{        solution <-object
      #dobjectontshow NULL == "*"
      nostar <- (!substr(object[[length(object)]]$vsnam, 1, 1) == "*")    
      show <-nostar
       
    	if (is.character(dontshow)){
    		if(dontshow=="*") dontshow="[*]"
        	  show <- !(grepl(dontshow, object[[length(object)]]$vsnam)) 
        	  }
   	 	if (is.logical(dontshow))show <-!dontshow 
    
   
    k <- length(object)
    ismodel <- "E=" %in% class(object)
    wha="complete independence"
    if(ismodel)wha=" model(E=) "
    pctota <- (100 * (object[[k]]$d)^2)/object[[k]]$ssX[1]
    pctotafc <-pctota
    if(!ismodel)pctotafc <- (100 * (object[[k]]$d)^2)/(object[[k]]$ssX[1] - 1)
    cat("\n", "+++ FCA- ",wha," ++ ", k, "modes+++ ", "\n")
    di <- NULL
    for (r in 1:length(object)) di <- c(di, length(object[[r]]$v[1,
        ]))
    cat("     ++ Contingency Table ", deparse(object[[k]]$datanam),
        " ", di, " ++", "\n")
    if(object[[k]]$addedcomment!="")cat("   ", object[[k]]$addedcomment, "\n")
    cat("     -----Total Percent Rebuilt----", sum(pctota[nostar]),
        "%", "\n")
    cat("     ++ Percent of lack of ",wha, " rebuilt  ++ ",
        ifelse(ismodel,sum(pctotafc[show]),sum(pctotafc[show][-1])), "%", "\n")
    cat("                    selected pctoafc > ",
        testvar, "%  total= ",  ifelse(ismodel,sum(pctotafc[show & pctotafc >
            testvar]),sum(pctotafc[show & pctotafc >
            testvar][-1])), "\n")
    summ <- matrix(cbind(1:length(object[[k]]$d), object[[k]]$d,
        object[[k]]$ssX, pctota, pctotafc), ncol = 5)
    summ <- summ[pctotafc > testvar & show, ]
    summ <- matrix(summ, ncol = 5)
    if(!ismodel)summ[1, 5] <- NA
    dimnames(summ) <- list(object[[k]]$vsnam[pctotafc > testvar &
        show], c("-no-", "--Sing Val--", "--ssX--", "--Global Pct--",
        "--FCA--"))
    print(summ, digits = 5)
    cat("\n", "++++               ++++", "\n")
    if (!is.null(testvar) & !testvar == 0)
        cat(" Shown are selected  over ", length(object[[k]]$vsnam[show]) -
            1, " PT  with pct FCA >", testvar, "% ", "\n")
    else cat(" over ", length(object[[k]]$vsnam), " PT (with*)",
        "\n")
  invisible(summ)
}
"summary.PTAk" <-
function (object, testvar = 1, dontshow = "*",...)
{
           if (is.character(dontshow)){
    		if(dontshow=="*") dontshow="[*]"
        	  show <- !(grepl(dontshow, object[[length(object)]]$vsnam)) 
        	  }
   	 	if (is.logical(dontshow))show <-!dontshow 
    RESUM(object, summary = TRUE, testvar = testvar, with = show)
}
"APSOLU3" <-
function (X, solu, pt3 = NULL, nbPT2 = 1, smoothing = FALSE,
    smoo = list(NA), verbose = getOption("verbose"), file = NULL, ...)
{
    if (is.list(X)) {
        if (is.list(X$met))
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        for (d in 1:3) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]],
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:3)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:3)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:3, laperm))
            }
            else X$data <- X$data * sqrt(X$met[[d]])
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
        for (d in 1:length(solu)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solu[[d]]$v <- solu[[d]]$v %*% Powmat(met[[d]],
                    1/2)
                }
                else {
                  solu[[d]]$v <- t(sqrt(met[[d]]) * t(solu[[d]]$v))
                }
            }
            else solu[[d]]$v <- solu[[d]]$v * sqrt(met[[d]])
        }
    }
    if (!is.null(pt3))
        numsol <- pt3
    else numsol <- length(solu[[length(solu)]]$d)
    Zsol <- list(NULL, NULL)
    if (smoothing & !length(smoo) == 3)
        stop(paste("--- Smoothing list must be of length 3! ---"))
    for (i in 1:3) {
        if (verbose) {
            cat(" ----------APSOLU3------------", file = ifelse(is.null(file),
                "", file), append = TRUE)
            cat(" ---- Associated solution to entry ---", i,
                file = ifelse(is.null(file), "", file), append = TRUE)
            cat("  ....  of dimension: ", dim(X)[i], "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
        }
        tracei <- i
        Z <- CONTRACTION(X, if (is.matrix(solu[[i]]$v))
            solu[[i]]$v[numsol, ]
        else solu[[i]]$v, Xwiz = i, zwiX = ifelse(length(numsol) ==
            1, 1, 2))
        i <- tracei
        if (nbPT2 == 1)
            nomb <- min(dim(Z))
        else nomb <- min(dim(Z), nbPT2)
        if (smoothing == TRUE)
            solq <- svdsmooth(Z, nomb = nomb, smooth = smoo[-i],...)
        else solq <- svd(Z)
        Zsol[[1]]$modesnam <- solu[[((1:3)[-i])[1]]]$modesnam
        nomb <- min(nomb, length(solq$d))
        Zsol[[1]]$v <- t(solq$u[, 1:nomb])
        Zsol[[2]]$modesnam <- solu[[((1:3)[-i])[2]]]$modesnam
        Zsol[[2]]$v <- t(solq$v[, 1:nomb])
        if (smoothing == TRUE) {
            Zsol[[2]]$smoocheck <- array(NA, c(3, nomb))
            Zsol[[2]]$smoocheck[(1:3)[-i], ] <- solq$smoocheck
        }
        if (!dim(Z)[2] > length(solq$d))
            ssX <- sum(solq$d^2)
        else ssX <- sum(as.vector(Z)^2)
        Zsol[[2]]$d <- solq$d[1:nomb]
        Zsol[[2]]$pct <- (100 * (solq$d^2)/ssX)[1:nomb]
        Zsol[[2]]$ssX <- rep(ssX, nomb)
        Zsol[[2]]$vsnam <- c(paste("*", dim(X)[i], solu[[length(solu)]]$vsnam[numsol],
            dim(X)[-i][1], dim(X)[-i][2], sep = ""), rep(paste(dim(X)[i],
            solu[[length(solu)]]$vsnam[numsol], dim(X)[-i][1],
            dim(X)[-i][2]), nomb - 1))
        solu <- RESUM(Zsol, solu, numass = numsol, verbose = verbose,
            file = file)
    }
    if (metrics) {
        for (d in 1:length(solu)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solu[[d]]$v <- solu[[d]]$v %*% Powmat(met[[d]],
                    -1/2)
                }
                else {
                  solu[[d]]$v <- t(1/sqrt(met[[d]]) * t(solu[[d]]$v))
                }
            }
            else solu[[d]]$v <- solu[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    return(solu)
}
"APSOLUk" <-
function (X, solu, nbPT, nbPT2 = 1, smoothing = FALSE, smoo = list(NA),
    minpct = 0.1, ptk = NULL, verbose = getOption("verbose"),
    file = NULL, modesnam = NULL, ...)
{
    if (is.list(X)) {
        if (is.list(X$met))
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        diX <- length(dim(X$data))
        for (d in 1:diX) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]],
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:diX)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:diX)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:diX, laperm))
            }
            else X$data <- X$data * sqrt(X$met[[d]])
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
        for (d in 1:length(solu)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solu[[d]]$v <- solu[[d]]$v %*% Powmat(met[[d]],
                    1/2)
                }
                else {
                  solu[[d]]$v <- t(sqrt(met[[d]]) * t(solu[[d]]$v))
                }
            }
            else solu[[d]]$v <- solu[[d]]$v * sqrt(met[[d]])
        }
    }
    numsol <- length(solu[[length(solu)]]$d)
    if (!is.null(ptk))
        numsol <- ptk
    Zsol <- list(NULL, NULL)
    kor <- length(dim(X))
    if (smoothing & !length(smoo) == kor)
        stop(paste("--- Smoothing list must be of length ", kor,
            "! ---"))
    for (i in 1:kor) {
        if (verbose) {
            cat("\n", "\n", "            ++++++++++++++++ --APSOLUk-- ",
                file = ifelse(is.null(file), "", file), append = TRUE)
            cat(solu[[kor]]$vsnam[numsol], " Associated solution to entry ---",
                i, file = ifelse(is.null(file), "", file), append = TRUE)
            cat("  ....  of dimension: ", dim(X)[i], "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
        }
        tracei <- i
        Z <- CONTRACTION(X, matrix(solu[[i]]$v, ncol = dim(X)[i])[numsol,
            ], Xwiz = i)
        i <- tracei
        if (length(dim(Z)) == 3) {
            solZ <- PTA3(Z, nbPT = nbPT[1], nbPT2 = nbPT2, smoothing = smoothing,
                smoo = smoo[-i], minpct = minpct, verbose = verbose,
                file = file, modesnam = modesnam[-i], ...)
        }
        if (length(dim(Z)) > 3) {
            solZ <- PTAk(Z, nbPT = nbPT, nbPT2 = nbPT2, smoothing = smoothing,
                smoo = smoo[-i], minpct = minpct, verbose = verbose,
                file = file, modesnam = modesnam[-i], ...)
        }
        nno <- length(solZ[[length(solZ)]]$vsnam)
        for (n in 1:nno) {
            if (!substr(solZ[[length(solZ)]]$vsnam[n], 1, 1) ==
                "*") {
                solZ[[length(solZ)]]$vsnam[n] <- paste(dim(X)[i],
                  solZ[[length(solZ)]]$vsnam[n], sep = "-")
            }
            else {
                solZ[[length(solZ)]]$vsnam[n] <- paste("*", paste(dim(X)[i],
                  substr(solZ[[length(solZ)]]$vsnam[n], 2, 100),
                  sep = "-"), sep = "")
            }
        }
        solZ[[length(solZ)]]$vsnam[1] <- paste("*", solZ[[length(solZ)]]$vsnam[1],
            sep = "")
        if (smoothing == TRUE) {
            smooche <- solZ[[length(solZ)]]$smoocheck
            solZ[[length(solZ)]]$smoocheck <- array(NA, c(kor,
                nno))
            solZ[[length(solZ)]]$smoocheck[(1:kor)[-i], ] <- smooche
        }
        solu <- RESUM(solZ, solu, numass = numsol, verbose = verbose,
            file = file)
    }
    if (metrics) {
        for (d in 1:length(solu)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solu[[d]]$v <- solu[[d]]$v %*% Powmat(met[[d]],
                    -1/2)
                }
                else {
                  solu[[d]]$v <- t(1/sqrt(met[[d]]) * t(solu[[d]]$v))
                }
            }
            else solu[[d]]$v <- solu[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    return(solu)
}
"FCAk" <-
function (X, nbPT = 3, nbPT2 = 1, minpct = 0.01, smoothing = FALSE,
    smoo = rep(list(function(u) ksmooth(1:length(u), u, kernel = "normal",
        bandwidth = 3, x.points = (1:length(u)))$y), length(dim(X))),
    verbose = getOption("verbose"), file = NULL, modesnam = NULL,
    addedcomment = "", chi2 = TRUE, E = NULL, ...)
{
    ldx <- length(dim(X))
    if (ldx <=2) {
        stop(paste("--- X must be an array  of k > 2 entries! ---"))
    }
    if (verbose) {
        cat("\n", "       ----------+++++++++++------------",
            "\n", ifelse(smoothing, paste("Penalised ", "\n"),
                ""), "       Correspondence Analysis on ", ldx,
            " modes ", "\n", file = ifelse(is.null(file), "",
                file), append = TRUE)
        if (!is.null(modesnam))
            cat("modes are ", modesnam, "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
        cat("       ----------+++++++++++------------", "\n",
            fil = ifelse(is.null(file), "", file), append = TRUE)
        cat("Data   = complete independence    + lack of independence ...",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
        cat(" lack of independence = partial independence + lack of independence ... etc ...",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
    }
    Y <- FCAmet(X, chi2 = chi2, E = E)
    if(ldx==3){
     solutions <- PTA3(Y, nbPT = nbPT, nbPT2 = nbPT2, smoothing = smoothing,
        smoo = smoo, minpct = minpct, verbose = verbose, file = file,
        modesnam = modesnam, addedcomment = addedcomment, ...)}
    else{
    solutions <- PTAk(Y, nbPT = nbPT, nbPT2 = nbPT2, smoothing = smoothing,
        smoo = smoo, minpct = minpct, verbose = verbose, file = file,
        modesnam = modesnam, addedcomment = addedcomment, ...)
    }
    solutions[[ldx]]$datanam <- substitute(X)
    solutions[[ldx]]$method <- match.call()
    solutions[[ldx]]$addedcomment <- addedcomment
    class(solutions) <- c("FCAk","PTAk")
    if(!is.null(E)) class(solutions) <- c("FCAk","PTAk","E=")
    invisible(solutions)
}
"FCA2" <-
function (X, nbdim =NULL, minpct = 0.01, smoothing = FALSE,
    smoo = rep(list(function(u) ksmooth(1:length(u), u, kernel = "normal",
        bandwidth = 3, x.points = (1:length(u)))$y), length(dim(X))),
    verbose = getOption("verbose"), file = NULL, modesnam = NULL,
    addedcomment = "", chi2 = FALSE, E = NULL, ...)
{
    ldx <- length(dim(X))
    if (ldx !=2) {
        stop(paste("--- X must be an array  of k = 2 entries! ---"))
    }
    if (verbose) {
        cat("\n", "       ----------+++++++++++------------",
            "\n", ifelse(smoothing, paste("Penalised ", "\n"),
                ""), "       Correspondence Analysis on ", ldx,
            " modes ", "\n", file = ifelse(is.null(file), "",
                file), append = TRUE)
        if (!is.null(modesnam))
            cat("modes are ", modesnam, "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
        cat("       ----------+++++++++++------------", "\n",
            fil = ifelse(is.null(file), "", file), append = TRUE)
        cat("Data   = complete independence    + lack of independence ...",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
        cat(" lack of independence = partial independence + lack of independence ... etc ...",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
    }
    Y <- FCAmet(X, chi2 = chi2, E = E)
     solutions <- SVDgen(Y,smoothing = smoothing, nomb=nbdim, smoo = smoo)
    
    solutions[[ldx]]$datanam <- substitute(X)
    solutions[[ldx]]$method <- match.call()
    solutions[[ldx]]$addedcomment <- addedcomment
    class(solutions) <- c("FCAk", "FCA2","PTAk")
    if(!is.null(E)) class(solutions) <- c("FCAk","FCA2","PTAk","E=")
    invisible(solutions)
}
"PPMA" <-
function (X, test = 1e-10, pena = list(function(u) ksmooth(1:length(u),
    u, kernel = "normal", bandwidth = 3, x.points = (1:length(u)))$y,
    NA), ini = mean, vsmin = 1e-20, Maxiter = 2000,...)
{
    v0 <- apply(X, 2, FUN = ini)
    if (all(v0 < 1e-04)) {
        v0 <- (X[sample(1:dim(X)[1], 1), ])
        if (max(abs(X)) < test * 1e-08) {
            cat(" Sum of squares veryyyy smallll  .......", "\n")
            return(list(u = as.matrix(rep(0, dim(X)[1])), v = as.matrix(rep(0, dim(X)[2])),
                d = 0, iter = 0, test = NA))
        }
    }
	PTnam="vsassocie"
	 if(exists("PTnam",envir=.GlobalEnv))PTnam=get("PTnam",envir=.GlobalEnv)
    test0 <- 1
    iter <- 1
    while (test0 > test) {
        u <- as.vector(X %*% v0)
        if (is.function(pena[[1]]))
            u <- pena[[1]](u)
        d <- sqrt(u %*% u)
        if (!d == 0)
            u <- u/d
        v <- as.vector(u %*% X)
        if (is.function(pena[[2]]))
            v <- pena[[2]](v)
        d <- sqrt(v %*% v)
        if (!d == 0)
            v <- v/d
        if (test0 == 1)
            u0 <- u
        if (!d < vsmin)
            test0 <- sum((u - u0)^2) + sum((v - v0)^2)
        else test0 <- 0
        v0 <- v
        u0 <- u
        iter <- iter + 1
        if (iter > (Maxiter - 1) && (iter - Maxiter)%%Maxiter == 0) {
            cat("\n \n \n \n \n ", " WARNING ****** Iteration already =  ",
                iter, "test= ", test0, "\n")
            cat(" ** type  999  to STOP ** just RETURN to carry on **",
                "\n")
            cat(" or type a new test value initial was", test,
                "\n")
            conti <- scan("", n = 1, quiet = TRUE, flush = TRUE)
            if (length(conti) > 0) {
                if (conti == 999)
                  stop(paste(" ---- Aborted by request ---- "))
                if (is.numeric(conti))
                  test <- conti
            }
        }
    }
    return(list("u" = as.matrix(u), "v" = as.matrix(v), "d" = as.vector(d),
        "iter" = iter, "test" = test0))
}
"PTA3" <-
function (X, nbPT = 2, nbPT2 = 1, smoothing = FALSE, smoo = list(function(u) ksmooth(1:length(u), 
    u, kernel = "normal", bandwidth = 4, x.points = (1:length(u)))$y, 
    function(u) smooth.spline(u, df = 3)$y, NA), minpct = 0.1, 
    verbose = getOption("verbose"), file = NULL, modesnam = NULL, 
    addedcomment = "", ...) 
{
    datanam <- substitute(X)
    if (is.list(X)) {
        if (is.list(X$met)) 
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        for (d in 1:3) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]], 
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:3)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:3)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:3, laperm))
            }
            else X$data <- X$data * sqrt(X$met[[d]])
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
    }
    debtime <- proc.time()
    pass. <- function(a, r) {
        pasta <- a
        for (i in 2:r) pasta <- paste(pasta, a, sep = "")
        return(pasta)
    }
    if (verbose) {
        cat("\n", "       ----------+++++++++++------------", 
            "\n", ifelse(smoothing, paste("Smoothed ", "\n"), 
                ""), "              PTA 3modes ", "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
        cat("       ----------+++++++++++------------", "\n", 
            file = ifelse(is.null(file), "", file), append = TRUE)
        cat(" Data is ... ", deparse(datanam), "...", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  .... Tensor of order ", length(dim(X)), file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  ....  with dimensions: ", dim(X), "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        if (!is.null(modesnam)) 
            cat("modes are ", modesnam, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
        if (metrics) 
            cat("---Analysis with non-Identity metrics  ------", 
                "\n", file = ifelse(is.null(file), "", file), 
                append = TRUE)
        if (!addedcomment == "") 
            cat("\n", addedcomment, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
    }
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    if (length(dim(X)) <=2) {
        stop(paste("--- X must be an array  of k > 2 entries! ---"))
    }
    solutions <- NULL
    if (smoothing) {
        if (smoothing & !length(smoo) == 3) 
            stop(paste("--- Smoothing list must be of length 3! ---"))
        for (j in 1:3) if (!is.list(smoo[[j]])) 
            smoo[[j]] <- list(smoo[[j]])
    }
    for (t in 1:nbPT) {
gc()
        if (verbose) {
            cat("----- Principal Tensor ---- ", paste("vs", pass.(t, 
                3), sep = ""), file = ifelse(is.null(file), "", 
                file), append = TRUE)
        }
        if (smoothing) {
            if (t > 1) {
                for (j in 1:3) if (length(smoo[[j]]) == t - 1) 
                  smoo[[j]][[t]] <- smoo[[j]][[t - 1]]
            }
            tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]], smoo[[3]][[t]])
        }
        else tosmoo <- list(NA)
 gc()
        solut <- SINGVA(X, verbose = verbose, file = file, PTnam = paste("vs", 
            pass.(t, 3), sep = ""), smoothing = smoothing, smoo = tosmoo, 
            modesnam = modesnam, ...)
        if (is.null(solutions) & verbose) 
            cat(" --- GLobal Percent --- ", (100 * solut[[3]]$d^2)/solut[[3]]$ssX[1], 
                "%", "\n", file = ifelse(is.null(file), "", file), 
                append = TRUE)
        if (verbose & !is.null(solutions)) {
            cat("                 -- GLobal Percent -- ", (100 * 
                solut[[3]]$d^2)/solutions[[3]]$ssX[1], "%", "\n", 
                file = ifelse(is.null(file), "", file), append = TRUE)
        }
        if (!is.null(solutions)) {
            if (100 * solut[[length(solut)]]$d^2/solutions[[length(solutions)]]$ssX[1] < 
                minpct) {
                cat("\n", "\n", " ++ Last 3-modes vs < ", minpct, 
                  "% stopping this level and under ++", "\n")
                solutions <- RESUM(solut, solutions, verbose = verbose, 
                  file = file)
                break
            }
        }
        if (nbPT2 >= 1) 
            solut <- APSOLU3(X, solut, pt3 = NULL, nbPT2 = nbPT2, 
                smoothing = smoothing, smoo = tosmoo, verbose = verbose, 
                file = file,...)
        if (verbose) 
            cat("\n", "+++ PTA 3modes  ------After ---", paste("vs", 
                pass.(t, 3), sep = ""), file = ifelse(is.null(file), 
                "", file), append = TRUE)
        solutions <- RESUM(solut, solutions, verbose = verbose, 
            file = file)
gc()
        if (t < nbPT) 
            X <- PROJOT(X, solut)
    }
    if (metrics) {
        for (d in 1:3) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solutions[[d]]$v <- solutions[[d]]$v %*% Powmat(met[[d]], 
                    -1/2)
                }
                else {
                  solutions[[d]]$v <- t(1/sqrt(met[[d]]) * t(solutions[[d]]$v))
                }
            }
            else solutions[[d]]$v <- solutions[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    solutions[[3]]$method <- match.call()
    solutions[[3]]$addedcomment <- addedcomment
    solutions[[length(solutions)]]$datanam <- datanam
    cat("\n", "-----Execution Time-----", (proc.time() - debtime)[3], 
        "\n")
    class(solutions) <- c("PTAk")
    invisible(solutions)
}
"PTAk" <-
function (X, nbPT = 2, nbPT2 = 1, minpct = 0.1, smoothing = FALSE,
    smoo = list(NA), verbose = getOption("verbose"), file = NULL,
    modesnam = NULL, addedcomment = "", ...)
{
    datanam <- substitute(X)
    if (is.list(X)) {
        if (is.list(X$met))
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        diX <- length(dim(X$data))
        for (d in 1:diX) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]],
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:diX)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:diX)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:diX, laperm))
            }
            else X$data <- X$data * sqrt(X$met[[d]])
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
    }
    debtime <- proc.time()
    pass. <- function(a, r) {
        pasta <- a
        for (i in 2:r) pasta <- paste(pasta, a, sep = "")
        return(pasta)
    }
    if (verbose) {
        cat("----------+++++++++++------------", "\n", ifelse(smoothing,
            paste("Penalised ", "\n"), ""), " Principal Tensor Analysis on k modes ",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
        cat(" Data is ... ", deparse(datanam), "...", "\n", file = ifelse(is.null(file),
            "", file), append = TRUE)
        cat("  .... Tensor of order ", length(dim(X)), file = ifelse(is.null(file),
            "", file), append = TRUE)
        cat("  ....  with dimensions: ", dim(X), "\n", file = ifelse(is.null(file),
            "", file), append = TRUE)
        if (!is.null(modesnam))
            cat("modes are ", modesnam, "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
        if (metrics)
            cat("---Analysis with non-Identity metrics  ------",
                "\n", file = ifelse(is.null(file), "", file),
                append = TRUE)
        if (!addedcomment == "")
            cat("\n", addedcomment, "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
    }
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    kor <- length(dim(X))
    if (kor <=2) {
        stop(paste("--- X must be an array  of k > 2 entries! ---"))
    }
    if (length(nbPT) < kor - 2) {
        nbPT <- rep(nbPT[1], kor - 2)
    }
    if (is.null(modesnam)) {
        modesnam <- paste(rep("mo", kor), 1:kor)
    }
    solutions <- NULL
    if (smoothing) {
        if (smoothing & !length(smoo) == kor)
            stop(paste("--- Smoothing list must be of length ",
                kor, "! ---"))
        for (j in 1:kor) if (!is.list(smoo[[j]]))
            smoo[[j]] <- list(smoo[[j]])
    }
    for (t in 1:(nbPT[kor - 2])) {
    gc()
        if (verbose)
            cat("\n", "\n", "                ++++++  k-modes Solutions  ---- k=",
                kor, paste(", vs", pass.(t, kor), sep = ""),
                "++++++", "\n", "\n", file = ifelse(is.null(file),
                  "", file), append = TRUE)
        tosmoo <- list(NA)
        if (smoothing) {
            if (t > 1) {
                for (j in 1:kor) if (length(smoo[[j]]) == t -
                  1)
                  smoo[[j]][[t]] <- smoo[[j]][[t - 1]]
            }
            if (kor == 3)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]])
            if (kor == 4)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]])
            if (kor == 5)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]], smoo[[5]][[t]])
            if (kor == 6)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]], smoo[[5]][[t]],
                  smoo[[6]][[t]])
            if (kor == 7)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]], smoo[[5]][[t]],
                  smoo[[6]][[t]], smoo[[7]][[t]])
            if (kor == 8)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]], smoo[[5]][[t]],
                  smoo[[6]][[t]], smoo[[7]][[t]], smoo[[8]][[t]])
        }
  gc()
        solut <- SINGVA(X, verbose = verbose, file = file,
            PTnam = paste("vs", pass.(t, kor), sep = ""), 
            smoothing = smoothing, smoo = tosmoo, modesnam = modesnam, ...)
        if (is.null(solutions) & verbose)
            cat("                 -- GLobal Percent -- ", solut[[kor]]$pct,
                "%", "\n", file = ifelse(is.null(file), "", file),
                append = TRUE)
        if (verbose & !is.null(solutions)) {
            cat("                 -- GLobal Percent -- ", (100 *
                solut[[length(solut)]]$d^2)/solutions[[length(solutions)]]$ssX[1],
                "%", "\n", file = ifelse(is.null(file), "", file),
                append = TRUE)
        }
        if (!is.null(solutions)) {
            if (100 * solut[[length(solut)]]$d^2/solutions[[length(solutions)]]$ssX[1] <
                minpct) {
                cat("\n", "\n", " ++ Last ", kor, "-modes vs < ",
                  minpct, "% stopping this level and under ++",
                  "\n", file = ifelse(is.null(file), "", file),
                  append = TRUE)
                solutions <- RESUM(solut, solutions, verbose = verbose,
                  file = file)
                break
            }
        }
        if (kor - 3 > 0) {
            if (!nbPT[kor - 3] == 0) {
            gc()
                solut <- APSOLUk(X, solut, nbPT = nbPT, nbPT2 = nbPT2,
                  smoothing = smoothing, smoo = tosmoo, minpct = minpct,
                  ptk = NULL, verbose = verbose, file = file,
                  modesnam = modesnam, ...)
            }
        }
        if (kor == 3 & nbPT2 >= 1) {
        gc()
            ptk <- NULL
            solut <- APSOLU3(X, solut, pt3 = ptk, nbPT2 = nbPT2,
                smoothing = smoothing, smoo = tosmoo, verbose = verbose,
                file = file, ...)
        }
        solutions <- RESUM(solut, solutions, verbose = verbose,
            file = file)
        if (is.null(solutions[[length(solutions)]]$datanam))
            solutions[[length(solutions)]]$datanam <- datanam
        gc()
        if (t < nbPT[kor - 2])
            X <- PROJOT(X, solut)
    }
    if (metrics) {
        for (d in 1:length(solutions)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solutions[[d]]$v <- solutions[[d]]$v %*% Powmat(met[[d]],
                    -1/2)
                }
                else {
                  solutions[[d]]$v <- t(1/sqrt(met[[d]]) * t(solutions[[d]]$v))
                }
            }
            else solutions[[d]]$v <- solutions[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    solutions[[kor]]$method <- match.call()
    solutions[[kor]]$addedcomment <- addedcomment
    cat("-----Execution Time-----", (proc.time() - debtime)[3],
        "\n")
    class(solutions) <- c("PTAk")
    invisible(solutions)
}
"SINGVA" <-
function (X, test = 1e-12, PTnam = "vs111", Maxiter = 2000, verbose = getOption("verbose"),
    file = NULL, smoothing = FALSE, smoo = list(NA), modesnam = NULL,
    Ini = "Presvd", sym = NULL)
{
    datanam <- substitute(X)
    if (is.list(X)) {
        if (is.list(X$met))
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        diX <- length(dim(X$data))
        for (d in 1:diX) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]],
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:diX)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:diX)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:diX, laperm))
            }
            else X$data <- X$data * sqrt(X$met[[d]])
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
    }
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    ord <- length(dim(X))
    if (verbose) {
        cat("\n", "       ----------+++++++++++------------ RPVSCC algorithm ",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
        cat("                             ------------ Singular Value  ",
            PTnam, "\n", file = ifelse(is.null(file), "", file),
            append = TRUE)
        cat("                                       ----  dimensions:  ",
            dim(X), "\n", file = ifelse(is.null(file), "", file),
            append = TRUE)
    }
    if(class(Ini)=="PTAk") 
          sval0 <- Ini 
    else{
      sval0 <- INITIA(X, modesnam = modesnam, method = Ini)
    }
    if (!is.null(sym)) {
        if (!ord == length(sym))
            stop(paste("--- Wrong length for parameter sym ! ---"))
        for (i in 1:ord) {
            if (!i == sym[i])
                sval0[[i]] <- sval0[[sym[i]]]
        }
    }
    if (verbose) {
        cat(" ----------------------", "\n", "Initialisation  done",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
    }
    sval <- sval0
    if (smoothing) {
        sval[[ord]]$smoocheck <- array(FALSE, c(ord, 1))
        if (!length(smoo) == ord)
            stop(paste("--- Smoothing list must be of length ",
                ord, "! ---"))
        for (i in 1:ord) smoo[[i]] <- toplist(smoo[[i]])
    }
    test0 <- 1
    atest <- 0
    iter <- 0
    while (test0 > test) {
        iter <- iter + 1
        if (verbose & iter%%100 == 1) {
            cat(" ----------- iteration-", iter, "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
        }
        for (i in 1:ord) {
            if (iter == 1) {
                if (verbose) {
                  cat(" ", i, "^", sval0[[i]]$d, file = ifelse(is.null(file),
                    "", file), append = TRUE)
                }
                sval[[i]]$d <- NULL
            }
            tracei <- i
            v <- CONTRACTION.list(X, sval0, moins = i)
            i <- tracei
            if (smoothing) {
                if (is.function(smoo[[i]])) {
                  v <- smoo[[i]](as.vector(v))
                  sval[[ord]]$smoocheck[i, 1] <- TRUE
                }
            }
            sval[[i]]$v <- as.vector(v)
            sval0[[i]]$v <- as.vector(sval0[[i]]$v)
            sd <- sqrt(sval[[i]]$v %*% sval[[i]]$v)
            if (verbose & iter%%100 == 0) {
                cat(" --", sd, file = ifelse(is.null(file), "",
                  file), append = TRUE)
            }
            if (!sd == 0)
                sval[[i]]$v <- (sval[[i]]$v)/sd
            atest <- atest + (sval[[i]]$v - sval0[[i]]$v) %*%
                (sval[[i]]$v - sval0[[i]]$v)
            if (!is.null(sym)) {
                for (i in ord:1) {
                  if (!i == sym[i])
                    sval[[sym[i]]] <- sval[[i]]
                }
            }
            sval0 <- sval
        }
        test0 <- sqrt(atest)
        atest <- 0
        if (verbose & (iter%%100) == 1) {
            cat("\n", "----------- test =         ", test0, "\n",
                file = ifelse(is.null(file), "", file), append = TRUE)
        }
        if (iter > (Maxiter - 1) && (iter - Maxiter)%%Maxiter == 0) {
            cat("\n \n \n \n \n ", " WARNING ****** Iteration already =  ",
                iter, "test= ", test0, "\n")
            cat(" ** type  999  to STOP ** just RETURN to carry on **",
                "\n")
            cat(" or type a new test value initial was", test,
                "\n")
            conti <- scan("", n = 1, quiet = TRUE, flush = TRUE)
            if (length(conti) > 0) {
                if (conti == 999)
                  stop(paste(" ---- Aborted by request ---- "))
                if (is.numeric(conti))
                  test <- conti
            }
        }
        sval0 <- sval
    }
    ssX <- sum(X^2)
    sstens <- sd^2
    totPCT <- 100 * sstens/ssX
    if (!verbose) {
        cat(" ---Final iteration--- ", iter, "\n")
        cat(" --Singular Value-- ", sd, " -- Local Percent -- ",
            totPCT, "%", "\n")
    }
    else {
        cat(" --------Final iteration----", iter, "\n", file = ifelse(is.null(file),
            "", file), append = TRUE)
        cat(" ----------- test =         ", test0, "\n", file = ifelse(is.null(file),
            "", file), append = TRUE)
        cat("\n", " --Singular Value-- ", sd, " -- Local Percent -- ",
            totPCT, "%", "\n", file = ifelse(is.null(file), "",
                file), append = TRUE)
    }
    sval[[i]]$iter <- iter
    sval[[i]]$test <- test
    sval[[i]]$d <- as.vector(sd)
    sval[[i]]$pct <- as.vector(totPCT)
    sval[[i]]$ssX <- as.vector(ssX)
    sval[[i]]$vsnam <- PTnam
    if (metrics) {
        for (d in 1:length(sval)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  sval[[d]]$v <- as.vector(sval[[d]]$v %*% Powmat(met[[d]],
                    -1/2))
                }
                else {
                  sval[[d]]$v <- 1/sqrt(met[[d]]) * sval[[d]]$v
                }
            }
            else sval[[d]]$v <- sval[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    class(sval) <- c("PTAk")
    return(sval)
}
"SVDgen" <-
function (Y, D2 = 1, D1 = 1, smoothing = FALSE, nomb = NULL,
    smoo = list(function(u) ksmooth(1:length(u), u, kernel = "normal",
        bandwidth = 3, x.points = (1:length(u)))$y)){
    
    datanam <- substitute(Y)   
 if(is.list(Y)) {
    D1=Y$met[[1]]
    D2=Y$met[[2]]
    Y=Y$data
  }
    nomb <- min(nomb, dim(Y))
    
    dinam <- dimnames(Y)
    if (length(D1) == dim(Y)[1]^2) {
        Y <- Powmat(D1, 1/2) %*% Y
    }
    else if (length(D1) == dim(Y)[1]) {
        Y <- sqrt(D1) * Y
    }
    else if (length(D1) == 1) {
        Y <- sqrt(D1) * Y
    }
    else stop(paste(" ----- Wrong DIMENSION for Metric of the First Entry ------ !!!!@@@@"))
    if (length(D2) == dim(Y)[2]^2) {
        Y <- Y %*% Powmat(D2, 1/2)
    }
    else if (length(D2) == dim(Y)[2]) {
        Y <- t(sqrt(D2) * t(Y))
    }
    else if (length(D2) == 1) {
        Y <- sqrt(D2) * Y
    }
    else stop(paste(" ----- Wrong DIMENSION for Metric of the First Entry ------ !!!!@@@@"))
    dimnames(Y) <- dinam
    if (smoothing == TRUE)
        result <- svdsmooth(Y, nomb = nomb, smooth = smoo)
    else result <- svd(Y)
    if (length(D1) == dim(Y)[1]^2) {
        result$u <- Powmat(D1, -1/2) %*% result$u
    }
    else if (length(D1) == dim(Y)[1]) {
        result$u <- 1/sqrt(D1) * result$u
    }
    else if (length(D1) == 1) {
        result$u <- (1/sqrt(D1)) * result$u
    }
    if (length(D2) == dim(Y)[2]^2) {
        result$v <- Powmat(D2, -1/2) %*% result$v
    }
    else if (length(D2) == dim(Y)[2]) {
        result$v <- 1/sqrt(D2) * result$v
    }
    else if (length(D2) == 1) {
        result$v <- 1/sqrt(D2) * result$v
    }
    if (all(result$u[, 1] < 0)) {
        result$u <- result$u * (-1)
        result$v <- result$v * (-1)
    }
    solutions <- list(NULL, NULL)
    solutions[[1]]$v <- t(result$u)[1:nomb, ]
    solutions[[1]]$modesnam <- "lignes"
    solutions[[1]]$n <- dimnames(Y)[[1]]
    solutions[[2]]$v <- t(result$v)[1:nomb, ]
    solutions[[2]]$modesnam <- "colonnes"
    solutions[[2]]$n <- dimnames(Y)[[2]]
    solutions[[2]]$d <- result$d[1:nomb]
    if (smoothing | nomb < min(dim(Y)))
        ssX <- sum(as.vector(Y)^2)
    else ssX <- sum(result$d^2)
    solutions[[2]]$pct <- (100 * result$d^2/ssX)[1:nomb]
    solutions[[2]]$ssX <- rep(ssX, nomb)
    solutions[[2]]$vsnam <- paste("vs", 1:nomb, sep = "")
    solutions[[2]]$datanam <- datanam
    solutions[[2]]$addedcomment <- ""
    solutions[[2]]$method <- match.call()
    if (smoothing)
        solutions[[2]]$smoocheck <- result$smoocheck
    class(solutions) <- c("PTAk")
    return(solutions)
}
"svdsmooth" <-
function (X, nomb = min(dim(X)), smooth = list(function(u) ksmooth(1:length(u),
    u, kernel = "normal", bandwidth = 3, x.points = (1:length(u)))$y),
    vsmin = 1e-16,...)
{
    if (!is.list(smooth[[1]]))
        smooth[[1]] <- list(smooth[[1]])
    if (length(smooth) < 2)
        smooth <- rep(list(smooth[[1]]), 2)
    if (!is.list(smooth[[2]]))
        smooth[[2]] <- list(smooth[[2]])
    solu <- list(NULL, NULL)
    solu[[1]]$v <- array(0, c(nomb, dim(X)[1]))
    solu[[2]]$v <- array(0, c(nomb, dim(X)[2]))
    solu[[2]]$d <- rep(0, nomb)
    solu[[2]]$smoocheck <- array(NA, c(2, nomb))
    solu[[2]]$smoocheck[1:2, 1] <- c(is.function(smooth[[1]][[1]]),
        is.function(smooth[[2]][[1]]))
    fi <- PPMA(X, pena = list(smooth[[1]][[1]], smooth[[2]][[1]]),...)
    solu[[1]]$v[1, ] <- fi$u
    solu[[2]]$v[1, ] <- fi$v
    solu[[2]]$d[1] <- fi$d
    for (qi in 2:nomb) {
        X <- PROJOT(X, solu, numo = (qi - 1))
        if (length(smooth[[1]]) == qi - 1)
            smooth[[1]][[qi]] <- smooth[[1]][[qi - 1]]
        if (length(smooth[[2]]) == qi - 1)
            smooth[[2]][[qi]] <- smooth[[2]][[qi - 1]]
        tempi <- list(toplist(smooth[[1]][[qi]]), toplist(smooth[[2]][[qi]]))
        fi <- PPMA(X, pena = tempi,...)
        solu[[1]]$v[qi, ] <- fi$u
        solu[[2]]$v[qi, ] <- fi$v
        solu[[2]]$d[qi] <- fi$d
        solu[[2]]$smoocheck[1:2, qi] <- c(is.function(smooth[[1]][[qi]]),
            is.function(smooth[[2]][[qi]]))
        if (fi$d < vsmin)
            break
    }
    return(list(u = t(solu[[1]]$v), d = solu[[2]]$d, v = t(solu[[2]]$v),
        smoocheck = solu[[2]]$smoocheck))
}
"toplist" <-
function (li)
{
    while (is.list(li)) {
        li <- li[[1]]
    }
    return(li)
}
"CauRuimet" <-
function (Z, ker = 1, m0 = 1, withingroup = TRUE, loc = substitute(apply(Z,
    2, mean, trim = 0.1)), matrixmethod = TRUE, Nrandom=3000)
{

    debtime <- proc.time()
   if(Nrandom < dim(Z)[1]) {
         sN=sample(1:(dim(Z)[1]),Nrandom)
         Z=Z[sN,]
       if(is.matrix(m0)){m0=m0[sN,]}
   }
  if(is.character(m0)){
    if (m0 == "tridiag") {
        m0 <- array(as.integer(0), c(dim(Z)[1], dim(Z)[1]))
        m0[1:2, 1] <- c(1, 1)
        m0[(dim(Z)[1] - 1):dim(Z)[1], dim(Z)[1]] <- c(1, 1)
        for (j in 2:(dim(Z)[1] - 1)) {
            m0[j - 1, j] <- 1
            m0[j, j] <- 1
            m0[j + 1, j] <- 1
        }
    }
  }
    mz <- eval(loc)
    Sz <- sweep(Z, 2, mz)
    Sz <- t(Sz) %*% Sz/(dim(Z)[1] - 1)
    norm2S <- function(u, S = Powmat(Sz, (-1))) {
        return(t(u) %*% S %*% u)
    }
    if (is.numeric(ker)) {
        g <- ker
        ker <- function(t) {
            return(exp(-(g * t)))
        }
    }
    if (withingroup) {
        if (matrixmethod) {
            distZiZj <- norm2S(t(Z))
            diadis <- diag(distZiZj)/2
            distZiZj <- 2 * sweep(sweep(-distZiZj, 2, -diadis),
                1, -diadis)
            M <- m0 * ker(distZiZj)
            sumM <- (sum(as.vector(M)) - dim(Z)[1])/2
            M <- diag(apply(M, 2, sum)) - M
            W <- norm2S(Z, M)/sumM
        }
        else {
            W <- matrix(0, nrow = dim(Z)[2], ncol = dim(Z)[2])
            totad <- 0
            for (i in 1:(dim(Z)[1] - 1)) for (j in (i + 1):dim(Z)[1]) {
                ad <- as.double(ker(norm2S(Z[i, ] - Z[j, ])))
                if (is.matrix(m0))
                  ad <- ad * m0[i, j]
                W <- W + ad * ((Z[i, ] - Z[j, ]) %o% (Z[i, ] -
                  Z[j, ]))
                totad <- totad + ad
            }
            totad <- totad
            W <- W/totad
        }
        cat("--- W local variance like ---","\n")
    }
    else {
    if (matrixmethod | !matrixmethod) {
       
        W <- matrix(0, nrow = dim(Z)[2], ncol = dim(Z)[2])
        totad <- 0
         if (is.matrix(m0)){
         for (i in 1:(dim(Z)[1] - 1)) for (j in (i + 1):dim(Z)[1]) {
                ad <- as.double(ker(norm2S(Z[i, ] - Z[j, ])))
                  ad <- ad * m0[i, j]
                W <- W + ad * ((Z[i, ] - mz) %o% (Z[j, ] -
                  mz))
                totad <- totad + ad
            }
            totad <- totad
            W <- W/totad
             cat("--- W global variance like  ---","\n")   
         }
         else {
        for (i in 1:(dim(Z)[1])) {
            ad <- as.double(ker(norm2S(Z[i, ] - mz)))
            W <- W + ad * ((Z[i, ] - mz) %o% (Z[i, ] - mz))
            totad <- totad + ad
        }
        totad <- totad #* dim(Z)[1]^2
        W <- W/totad
         cat("--- W robust total  variance like ---","\n")
        }
       }
     }
    cat("-----Execution Time-----", (proc.time() - debtime)[3],
        "\n")
    return(W)
}
"Detren" <-
function (dat, Mm = c(1, 3), rsd = TRUE, tren = function(x) smooth.spline(as.vector(x),
    df = 5)$y)
{
    tre <- apply(dat, Mm, FUN = tren)
    dimi <- c(dim(dat)[-Mm], dim(dat)[Mm])
    tre <- aperm(array(tre, dimi), match(dimi, dim(dat)))
    if (rsd)
        return(dat - tre)
    else return(tre)
}
"FCAmet" <-
function (X, chi2 = FALSE, E = NULL,No0margins=TRUE)
{
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    
    datanam <- substitute(X)
    ord <- length(dim(X)) 
     N <- sum(X)
      metafc <- rep(list(NULL), ord)
       
    if(No0margins){
      pass. <- function(a, r) {
        pasta <- a
        if(r==0)return("")
        if(r==1)return(pasta)
        for (i in 2:r) pasta <- paste(pasta, a, sep = "")
        return(pasta)
    }
    evalCh.f<-function(st){
      #st is a expression quoted e.g."x=2"
      tmp <- tempfile()
       writeLines(st, tmp)
       return(eval.parent(parse(tmp)))      
         } #end of evalCh.f          
     #library(tensorA)
      dnam=dimnames(X)
      #X=to.tensor(as.vector(X),dim(X))
       for(t in 1:ord){
          metafc[[t]] <- apply(X, t, sum)
          if (any(metafc[[t]]==0)){
             cat("-------Zeros for margin:",t,"\n")
             cat("-------zero values replaced by min margin on dim:",t,"\n")
             cat("--------old $count N:",N,"\n")
               the0=(1:length(metafc[[t]]))[metafc[[t]]==0]
               cat(the0,"\n")
             amin=min(metafc[[t]][metafc[[t]]!=0])/prod(dim(X)[-t])
            
          # evalCh.f(paste("X[[",names(X)[t],"=the0]]=amin",sep=""))
             # t th position X[,,,theo,,]=amin
           evalCh.f(paste("X[",  pass.(",",(t-1)),"the0", pass.(",",(ord-t)), "]=rep(amin,length(the0))",sep=""))
            
          }
       }
       #X=array(as.vector(X),dim(X))
       #dimnames(X)=dnam 
     } # X rebuilt if zeos margins
      
      N <- sum(X)  
     metafc[[1]] <- apply(X, 1, sum)/N
     Indep <-  metafc[[1]]
    for (t in 2:ord) {
        metafc[[t]] <- apply(X, t, sum)/N
        Indep <- Indep %o% metafc[[t]]
    }
    if (chi2) {
        Indep <- array(Indep, dim(X))
        Chi2 <- N * sum((X/N - Indep)^2/Indep)
        cat("\n", " --")
        cat("\n", "++++ Data is            ", deparse(datanam),
            "        +++++++")
        cat("\n", "-------------- Multiple contingency Table of dimensions ",
            dim(X), "  ----", "\n")
        cat("\n", "-------------- Chi2 = ", Chi2, " with ddl = ",
            prod(dim(X) - 1))
        cat("\n", " ------------- p(>Chi2)= ", pchisq(Chi2, df = prod(dim(X) -
            1), lower.tail = FALSE), "\n")
        cat("\n", " --", "\n")
    }
    
    if (!is.null(E))
        invisible(list(data = (X/N - E)/Indep, met = metafc,
            count = N))
    else invisible(list(data = (X/N)/Indep, met = metafc, count = N))
}
"Ginv" <-
function (A)
{
    Powmat(A, -1)
}
"IterMV" <-
function (n = 10, dat, Mm = c(1, 3), Vm = c(2, 3),
    fFUN = mean, usetren = FALSE, tren = function(x) smooth.spline(as.vector(x),
        df = 5)$y, rsd = TRUE)
{
    sdev <- function(x) {
        sd(as.vector(x))
    }
    for (i in 1:n) {
        if (usetren) {
            dat <- Detren(dat = dat, Mm = Mm, tren = tren, rsd = rsd)
        }
        else {
            mean.dat <- apply(dat, Mm, FUN = fFUN)
            dat <- sweep(dat, Mm, mean.dat)
        }
        sd.dat <- apply(dat, Vm, sdev)
        if (sd.dat == 1)
            warning("zero variances were replaced by 1")
        sd.dat <- ifelse(sd.dat == 0, 1, sd.dat)
        dat <- sweep(dat, Vm, sd.dat, FUN = "/")
    }
    cat("\n", "---Max of the means: ", max(apply(dat, Mm, mean),
        nam = TRUE), "\n")
    return(dat)
}
"Multcent" <-
function (dat, bi = c(1, 2), by = 3, centre = mean,
    centrebyBA = c(TRUE, FALSE), scalebyBA = c(TRUE, FALSE))
{
    if (centrebyBA[1]) {
        me <- apply(dat, by, FUN = centre)
        dat <- sweep(dat, by, me)
    }
    sdev <- function(x) {
        sd(as.numeric(x))
    }
    if (scalebyBA[1]) {
        sca <- apply(dat, by, sdev)
        if (sca == 1)
            warning("zero variances were replaced by 1")
        sca <- ifelse(sca == 0, 1, sca)
        dat <- sweep(dat, by, sca, FUN = "/")
    }
    if (!is.null(bi)) {
        for (g in 1:length(bi)) {
            me <- apply(dat, c(bi[g], by), FUN = centre)
            dat <- sweep(dat, c(bi[g], by), me)
        }
    }
    if (centrebyBA[2]) {
        me <- apply(dat, by, FUN = centre)
        dat <- sweep(dat, by, me)
    }
    if (scalebyBA[2]) {
        sca <- apply(dat, by, sdev)
        if (sca == 1)
            warning("zero variances were replaced by 1")
        sca <- ifelse(sca == 0, 1, sca)
        dat <- sweep(dat, by, sca, FUN = "/")
    }
    return(dat)
}
"Powmat" <-
function (A, pw, eltw = FALSE)
{
    A <- as.matrix(A)
    if (eltw) {
        dimA <- dim(A)
        A <- as.vector(A)
        RR <- A^pw
        RR[abs(RR) == Inf] <- A[abs(RR) == Inf]
        if (dimA[2] > 1)
            RR <- matrix(RR, ncol = dimA[2])
    }
    else {
        valsin <- svd(A)
        diago <- valsin$d[valsin$d > 1e-06]
        diago <- diago^pw
        if (length(diago) == 0) {
            RR <- matrix(0, ncol(A), nrow(A))
            return(RR)
        }
        if (length(diago) == 1)
            RR <- t(as.matrix(valsin$v[, 1]) %*% t(as.matrix(valsin$u[,
                1]))) * diago
        else RR <- valsin$u[, 1:length(diago)] %*% diag(diago) %*%
            t(valsin$v[, 1:length(diago)])
        RR <- as.matrix(RR)
        if (pw < 0 & (!min(dim(RR)) == 1))
            RR <- t(RR)
        if (length(RR) == 1)
            RR <- as.numeric(RR)
        else if (dim(RR)[1] == 1)
            RR <- as.vector(RR)
    }
    return(RR)
}
"RaoProd" <-
function (A, B)
{
    A <- as.matrix(A)
    B <- as.matrix(B)
    if (min(dim(A)) == 1 & min(dim(B)) == 1)
        return(as.vector(A) %x% as.vector(B))
    else {
        if (length(A) == 1 || length(B) == 1) {
            ifelse(length(B) == 1, return(A * as.vector(B)),
                return(as.vector(A) * B))
        }
        else {
            if (!dim(A)[2] == dim(B)[2])
                stop("Wrong number of columns")
            f <- dim(A)[2]
            re <- array(0, c(dim(A)[1] * dim(B)[1], f))
            for (w in 1:f) {
                re[, w] <- A[, w] %x% B[, w]
            }
            return(re)
        }
    }
}
"RiskJackplot" <-
function (x, nbvs = 1:20, mod = NULL, max = NULL, rescaled = TRUE, 
    ...)
{     solution <- x
    qchoix <- nbvs
    ord <- length(solution)
    if (is.null(max))
        max <- length(solution[[ord]]$d)
    if(max(nbvs) >= 0.8*max)warning("nbvs could be too high!")    
    if (is.null(mod))
        mod <- 1:length(solution)
    val <- solution[[ord]]$d[!substr(solution[[ord]]$vsnam, 1,
        1) == "*"]
    if (ord > 2)
        iden <- order(val)[length(val):1]
    else iden <- 1:max
    covalid <- function() {
        mindiff <- min((val[iden][-length(solution[[ord]]$d)]^2 -
            val[iden][-1]^2))
        for (mode in mod) {
            if (length(solution[[mode]]$v[1, ]) < solution[[ord]]$ssX[1]^2/mindiff/length(solution[[mode]]$v[1,
                ]^2)) {
                cat(" WARNING ..mode ", mode, " ..n= ", length(solution[[mode]]$v[1,
                  ]), " validity condition  >", solution[[ord]]$ssX[1]^2/mindiff/length(solution[[mode]]$v[1,
                  ]^2), "\n")
            }
        }
    }
    RJack <- matrix(rep(0, max(mod) * length(qchoix)), c(max(mod),
        length(qchoix)))
    for (m in mod) {
        for (q in qchoix) {
            tl <- 0
            if (q > (max - 1))
                q <- max - 1
            for (k in 1:q) {
                for (j in (q + 1):max(iden)) {
                  l1 <- solution[[ord]]$d[iden[j]]^2
                  l2 <- solution[[ord]]$d[iden[k]]^2
                  tjk <- mean(solution[[m]]$v[iden[j], ]^2 *
                    solution[[m]]$v[iden[k], ]^2) * l1 * l2
                  diff <- (l1 - l2)^2
                  tl <- tl + tjk/diff
                }
            }
            RJack[m, match(q, qchoix)] <- tl * 1/(length(solution[[m]]$v[j,
                ]) - 1)
            if (q == (max - 1))
                q <- max(qchoix)
        }
    }
    for (u in mod) {
        if (rescaled)
            RJack[u, ] <- (RJack[u, ] - min(RJack[u, ]))/(max(RJack[u,
                ]) - min(RJack[u, ]))
        plot(qchoix, RJack[u, ], xlab = "Nb of dimensions", ylab = "Risk's approx",
            lty = u, col = u, type = "b", ...)
        par(new = TRUE)
    }
    legend(2, max(RJack)/2, paste("Risk-mode",
        mod), col = mod, lty = mod, bty = "n", cex = 0.7)
    invisible(par(new = FALSE))
}
"Susan1D" <-
function (y, x = NULL, sigmak = NULL, sigmat = NULL, ker = list(function(u) return(exp(-0.5 *
    u^2))))
{
    if (is.null(x))
        x <- 1:length(y)
    else {
        if (!length(x) == length(y))
            stop("Wrong length for x")
        y <- y[order(x)]
        x <- sort(x)
    }
    if (is.null(sigmat))
        sigmat <- 8 * (length(y)^(-1/5))
    if (is.null(sigmak))
        sigmak <- 1/2 * (range(y)[2] - range(y)[1])
    if (length(ker) < 2)
        ker <- list(t = ker[[1]], k = ker[[1]])
    knei <- max(1, round(2 * sigmat))
    resul <- y
    for (t in 1:length(y)) {
        xt <- 0
        wjt <- 0
        for (j in max(1, t - knei):min(length(y), t + knei)) {
            wj <- ker$t((x[j] - x[t])/sigmat) * ker$k((y[j] -
                y[t])/sigmak)
            xt <- xt + wj * y[j]
            wjt <- wjt + wj
        }
        resul[t] <- xt/wjt
    }
    return(resul)
}
"COS2" <-function(solu,mod=1,solnbs=2:4){
	# as (t(phi_s) D phi_s)i/(t(vec_phi) D_-phi vec_phi ) if normed to lambda_s otherwise has to be times lambda_s
	# can be added on the upto for accumulated "rendering" 
	
    X <-eval(solu[[length(solu)]]$datanam)
                   
		if("FCAk" %in% class(solu)){
				if("E" %in% names(solu[[length(solu)]]$method)){
					leE=eval((solu[[length(solu)]]$method)$E)
				 X <-FCAmet(X,E=leE)				
				}
				else{				
						X <-FCAmet(X)
						X$data <- X$data -1
				} 
		}					 
	 # norm of the mod vv
	
	 if(is.list(X)) {
	 	# XD1/2 except mod
        nam <- dimnames(X$data)
        diX <- length(dim(X$data))
        for (d in 1:diX) {
        	if (d==mod){X$met[[d]]=1}
           
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]],
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:diX)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:diX)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:diX, laperm))
            }
            else X$data <- X$data * sqrt(X$met[[d]])
        }
       # metAll <- X$met
        X <- X$data
        dimnames(X) <- nam
    }
    
     diX <-length(dim(X))
       lacola <- (1:diX)[-mod]
           laperm <- c(mod, lacola)
              lacol <- (dim(X))[lacola]
       Xmod <- matrix(aperm(X, laperm), ncol = prod(lacol))	 	 
	 
	   
	  cos2 <- (solu[[mod]]$v[solnbs,]**2)*rep(solu[[length(solu)]]$d[solnbs]**2,dim(solu[[mod]]$v)[2])
	  vv <- apply(Xmod**2,1,sum)
	cos2 <- t(cos2)/vv
	
	 colnames(cos2) <- paste(rep("cos2",length(solnbs)),paste(":",solnbs,":",sep=""),solu[[length(solu)]]$vsnam[solnbs],sep="_")
	 rownames(cos2) <- solu[[mod]]$n
	return(round(cos2*1000))
}
##############
"CTR"  <-function(solu, mod=1,solnbs=1:4){
	# as ()t(phi_s) D phi_s)i/lambda_s  if normed to lambda_s otherwise not divided by lambda
                    
	if(is.list(eval(solu[[length(solu)]]$datanam))) {
	   	met <- eval(solu[[length(solu)]]$datanam)$met[[mod]]}
	else{
		if("FCAk" %in% class(solu)){
			met <-FCAmet(eval(solu[[length(solu)]]$datanam))$met[[mod]]
			}
		else{
	   		met <-1
	   	}
	 } 
	if(is.vector(met)){
		ctr <- t(solu[[mod]]$v[solnbs,]**2)*met 
		}
	else{
		ctr <- t( (solu[[mod]]$v[solnbs,]%*%met)*solu[[mod]]$v[solnbs,]) 
	}
	colnames(ctr) <- paste(rep("ctr",length(solnbs)),paste(":",solnbs,":",sep=""),solu[[length(solu)]]$vsnam[solnbs],sep="_")
	rownames(ctr) <- solu[[mod]]$n
	return(round(ctr*1000))
}
"plot.PTAk" <-
function (x, labels = TRUE, mod = 1, nb1 = 1, nb2 = NULL,
    coefi = list(NULL, NULL), xylab = TRUE, ppch = (1:length(solution)),
    lengthlabels = 2, scree = FALSE, ordered = TRUE,
    nbvs = 40, RiskJack = NULL, method = "",ZoomInOut=NULL, Zlabels=NULL, Zcol=NULL,...)
{      solution <- x
	awaybor=1.04
if(class(solution)[1]=="PCAn" | class(solution)[1]=="CANDPARA" )cat("\n","Plot function not available yet using Plot.PTAk!","\n")
    if (is.null(coefi[[1]]))
        coefi[[1]] <- rep(1, length(solution))
    if (is.null(coefi[[2]]))
        coefi[[2]] <- rep(1, length(solution))
    if (is.null(lengthlabels))
        lengthlabels <- rep(10, length(solution))
    if (length(lengthlabels) == 1)
        lengthlabels <- rep(lengthlabels, length(solution))
    ord <- length(solution)
    if ("FCAk" %in% class(x) && !("E=" %in% class(x)) ) {
        divv <- solution[[ord]]$ssX[1] - 1
        perclab <- "% FCA"
        if (length(nbvs) == 1)
            nbvs <- 2:nbvs
        else if (1 %in% nbvs)
            nbvs <- nbvs[-match(1, nbvs)]
    }
    else {
        divv <- solution[[ord]]$ssX[1]
        perclab <- "% global"
    }
    di <- NULL
    for (r in 1:length(solution)) di <- c(di, length(solution[[r]]$v[1,
        ]))
    if (!scree) {
        xlab <- ""
        ylab <- ""
        ylim <- NULL
        xlim <- NULL
        if (is.null(nb2)) {
            xlim <- c(1, max(di[mod]) + 1)
            if (xylab)
                ylab <- paste(solution[[ord]]$vsnam[nb1], " local",
                  round(solution[[ord]]$pct[nb1], 2), "% ", round((100 *
                    (solution[[ord]]$d[nb1])^2)/divv, 2), perclab)
        }
        else {
            if (xylab)
                xlab <- paste(solution[[ord]]$vsnam[nb1], " local",
                  round(solution[[ord]]$pct[nb1], 2), "% ", round((100 *
                    (solution[[ord]]$d[nb1])^2)/divv, 2), perclab)
            if (xylab)
                ylab <- paste(solution[[ord]]$vsnam[nb2], " local",
                  round(solution[[ord]]$pct[nb2], 2), "% ", round((100 *
                    (solution[[ord]]$d[nb2])^2)/divv, 2), perclab)
        }
        for (u in mod) {  
            if (!is.null(nb2)) {
                xyn <- t(solution[[u]]$v[c(nb1, nb2), ]) %*% diag(c(coefi[[1]][u], coefi[[2]][u]))                 
                xaxt <- "s"
                ylim <- range(xyn[,2],ylim)
                xlim <- range(xyn[,1],xlim)
                
            }
            else {
                xyn <- solution[[u]]$v[nb1, ] * coefi[[1]][u]
               if (!"xaxt" %in% names(list(...)))
                  xaxt <- "n" 
				ylim <- range(xyn,ylim)               
            }
		}
		ylim <- awaybor *ylim
		if (!is.null(nb2))xlim <- awaybor *xlim
        if(!is.null(ZoomInOut[[2]])) ylim <- ZoomInOut[[2]]
          if(!is.null(ZoomInOut[[1]]) &&  !is.null(nb2)) xlim <-ZoomInOut[[1]]

        for (u in mod) { 
        	 if (!is.null(nb2)) {
                 xy <- t(solution[[u]]$v[c(nb1, nb2), ]) %*% diag(c(coefi[[1]][u],coefi[[2]][u]))
             }
             else {
                 xy <- solution[[u]]$v[nb1, ] * coefi[[1]][u]
             }         	
            if (labels) {
                if ("xlab" %in% names(list(...))) {
                  if ("ylab" %in% names(list(...)))
                    plot(xy, xlim = xlim, ylim = ylim, pch = ppch[u],
                      xaxt = xaxt,...)
                  else plot(xy, xlim = xlim, ylim = ylim, ylab = ylab,
                    pch = ppch[u], xaxt = xaxt, ...)
                  if (is.null(nb2))
                    axis(1, 1:length(xy))
                }
                else {
                  if ("ylab" %in% names(list(...)))
                    plot(xy, xlim = xlim, ylim = ylim, pch = ppch[u],
                      xaxt = xaxt, xlab = xlab, ...)
                  else plot(xy, xlim = xlim, ylim = ylim, ylab = ylab,
                    pch = ppch[u], xaxt = xaxt, xlab = xlab,
                    ...)
                }
              if(!is.null(Zlabels[[u]])){
              	ZlabUN <- Zlabels[[u]] 
              	       }
              else {
              	ZlabUN <- solution[[u]]$n
              }  
                if (!is.null(ZlabUN)) {
                	ZcolUN= u
              	if(!is.null(Zcol[[u]]))ZcolUN=Zcol[[u]]

                  if (is.factor(ZlabUN)) {
                    if ("cex" %in% names(list(...)))
                      cex <- list(...)$cex
                    else cex <- par("cex")
                    text(xy, labels = substr(levels(ZlabUN),
                      1, lengthlabels[u]), col = ZcolUN,
                      pos = 4,...)
                    if (is.null(nb2)) {
                      par(new = TRUE)
                      plot(xy ~ ZlabUN, xlab = "", ylab = "",
                        ylim = ylim, cex = cex)
                      par(new = FALSE)
                    }
                  }
                  else { 
                  	ZcolUN= u
              	if(!is.null(Zcol[[u]]))ZcolUN=Zcol[[u]]
                  	text(xy, labels = substr(ZlabUN,
                    1, lengthlabels[u]), col = ZcolUN, pos = 4,...)
                    }
                }
            } 
            else 
            if ("xlab" %in% names(list(...))) {
                if ("ylab" %in% names(list(...)))
                  plot(xy, xlim = xlim, ylim = ylim, pch = ppch[u],
                    col = u, xaxt = xaxt, ...)
                else plot(xy, xlim = xlim, ylim = ylim, ylab = ylab,
                  pch = ppch[u], col = u, xaxt = xaxt, ...)
                if (is.null(nb2))
                  axis(1, 1:length(xy))
            }
            else {
                if ("ylab" %in% names(list(...)))
                  plot(xy, xlim = xlim, ylim = ylim, pch = ppch[u],
                    col = u, xaxt = xaxt, xlab = xlab, ...)
                else plot(xy, xlim = xlim, ylim = ylim, ylab = ylab,
                  pch = ppch[u], col = u, xaxt = xaxt, xlab = xlab,
                  ...)
            }
            if(!is.null(nb2)) {
            abline(h = 0, col = "green", lty = 2)
            abline(v = 0, col = "green", lty = 2)
            }
            par(new = TRUE)
        }
        invisible(par(new = FALSE))
    }
    else {
        if (!is.null(ordered)) {
            if (ordered == TRUE) {
                ld <- length(solution[[ord]]$d[!substr(solution[[ord]]$vsnam,
                  1, 1) == "*"])
                if (length(nbvs) == 1) {
                  nbvs <- min(max(3, nbvs), ld)
                  nbvs <- 1:nbvs
                }
                scre <- 100 * ((solution[[ord]]$d[!substr(solution[[ord]]$vsnam,
                  1, 1) == "*"])^2)/divv 
                scre <- (sort(scre[nbvs]))            
                scre <- scre[length(scre):1]
                if (!is.null(RiskJack)) scre <- scre[1:min(max(RiskJack,2),max(nbvs-2))]
                nbvs <- nbvs[1:length(scre)]
                plot(nbvs, scre, xlab = "Ordered ", ylab = "Squared Singular Values (%)",
                  xaxt = "n", ...)
                axis(1, at = nbvs)
                par(new = TRUE)
                plot(nbvs, ylim = c(0, 100), cumsum(scre), axes = FALSE,
                  lwd = 2, lty = 1, type = "b", pch = "c", col = 3,
                  xlab = "", ylab = "")
                axis(4, at = atpc <- seq(0, 100, 10), labels = formatC(atpc,
                  format = "fg"), col.axis = 3)
                par(new = TRUE)
                if (!is.null(RiskJack))
                  RiskJackplot(solution, nbvs = nbvs, mod = NULL,
                    max = NULL, rescaled = TRUE,
                    axes = FALSE, ann = FALSE, pch = "r")
                par(new = FALSE)
            } 
            
            if (ordered == FALSE) {
                ld <- length(solution[[ord]]$d[!substr(solution[[ord]]$vsnam,
                  1, 1) == "*"])
                if (length(nbvs) == 1) {
                  nbvs <- min(max(5, nbvs), ld)
                  nbvs <- 1:nbvs
                }
                scre <- ((solution[[ord]]$d)^2)[nbvs]
                plot(nbvs, scre, xlab = "Unordered with redundancy",
                  ylab = "Squared Singular Values", ...)
            }
        }
    }
    invisible(par(new = FALSE))
}
