rel.sys <-
function (x, type = c("matlist", "listmat"), bonds = c("entire", 
    "strong", "weak"), sel = NULL, prsep = ", ", loops = FALSE, 
    attr = NULL) 
{
    if (isTRUE(attr == 0L) == TRUE) {
        attr <- NULL
    }
    else {
        NA
    }
    if (match.arg(type) == "matlist") {
        if (is.null(attr) == FALSE) {
            if (is.na(dim(x)[3]) == FALSE) {
                if (isTRUE(max(attr) > dim(x)[3]) == TRUE) 
                  stop("Value of 'attr' greater than dim(x)[3]")
            }
            else if (is.na(dim(x)[3]) == TRUE) {
                if (isTRUE(max(attr) > 1L) == TRUE) 
                  stop("Value of 'attr' greater than dim(x)[3]")
            }
            ats <- bundles(x, collapse = FALSE, loops = TRUE, 
                prsep = prsep)[[7]][attr]
        }
        else if (is.null(attr) == TRUE) {
            ats <- logical(0)
        }
        if (is.na(dim(x)[3]) == FALSE) {
            if (isTRUE(all(seq(dim(x)[3]) %in% attr)) == FALSE) {
                bd <- bundles(x[, , which(!(seq(dim(x)[3]) %in% 
                  attr))], collapse = FALSE, loops = loops, prsep = prsep)
            }
            else if (isTRUE(all(seq(dim(x)[3]) %in% attr)) == 
                TRUE) {
                bd <- NULL
            }
        }
        else {
            bd <- bundles(x, collapse = FALSE, loops = loops, 
                prsep = prsep)
        }
        switch(match.arg(bonds), entire = lbd <- bd, strong = lbd <- list(bd$recp, 
            bd$txch, bd$mixd, bd$full), weak = lbd <- list(bd$asym, 
            bd$tent))
        if (is.null(lbd) == FALSE) {
            if (is.na(dim(x)[3]) == FALSE && isTRUE((dim(x)[3] - 
                length(attr)) == 0L) == FALSE) {
                stb <- list()
                for (k in 1:(dim(x)[3] - length(attr))) {
                  tmp <- vector()
                  for (i in 1:length(lbd)) {
                    if (isTRUE(length(lbd[[i]]) > 0L) == TRUE) {
                      ifelse(is.na(dim(x[, , which(!(seq(dim(x)[3]) %in% 
                        attr))])[3]) == TRUE, tmp <- append(tmp, 
                        lbd[[i]]), tmp <- append(tmp, lbd[[i]][k]))
                    }
                  }
                  rm(i)
                  stb[[k]] <- as.vector(unlist(tmp))
                }
                rm(k)
            }
            else {
                stb <- vector()
                for (i in 1:length(lbd)) {
                  stb <- append(stb, lbd[[i]])
                }
                rm(i)
            }
        }
        else {
            stb <- lbd
        }
        if (is.null(sel) == FALSE) {
            ntsel <- list()
            for (k in 1:length(stb)) {
                tss <- which(dhc(stb[[k]]) %in% sel)
                if (isTRUE(length(tss) > 0) == TRUE) {
                  tmpsel <- vector()
                  for (i in 1:length(tss)) {
                    if (isTRUE((tss[i]%%2L) == 1L) == TRUE) {
                      tmpsel <- append(tmpsel, stb[[k]][ceiling(tss[i]/2L)])
                    }
                    else {
                      tmpsel <- append(tmpsel, stb[[k]][floor(tss[i]/2L)])
                    }
                  }
                  rm(i)
                  ntsel[[k]] <- as.vector(unlist(tmpsel))
                }
            }
            rm(k)
            rm(tmpsel, tss)
            stb <- ntsel
        }
        else {
            NA
        }
        if (length(stb) > 0L) {
            ties <- vector()
            for (k in 1:length(stb)) {
                for (i in 1:length(stb[[k]])) {
                  if (isTRUE(length(stb[[k]]) > 0L) == TRUE) {
                    ties <- append(ties, dhc(stb[[k]][i], prsep = prsep))
                  }
                }
                rm(i)
            }
            rm(k)
        }
        else {
            ties <- stb <- character(0)
        }
        if (is.null(attr) == TRUE) {
            if (is.na(dim(x)[3]) == FALSE) {
                ifelse(is.null(dimnames(x)[[3]]) == TRUE, attr(stb, 
                  "names") <- 1:(dim(x)[3] - length(attr)), attr(stb, 
                  "names") <- dimnames(x)[[3]])
            }
        }
        else {
            ifelse(is.null(dimnames(x)[[3]]) == TRUE, attr(stb, 
                "names") <- which(!(seq(dim(x)[3]) %in% attr)), 
                attr(stb, "names") <- dimnames(x)[[3]][which(!(seq(dim(x)[3]) %in% 
                  attr))])
        }
        if (is.null(dimnames(x)[[1]]) == TRUE) {
            note <- "Input labels in 'x' are NULL."
            lbs <- 1:dim(x)[1]
        }
        else {
            note <- NULL
            lbs <- dimnames(x)[[1]]
        }
        if (isTRUE(length(ats) > 0L) == TRUE) {
            ifelse(length(note) > 0L, RS <- (list(ord = dim(x)[1], 
                nodes = lbs, sel = sel, sys.ord = nlevels(factor(ties)), 
                incl = lbs[which(lbs %in% levels(factor(ties)))], 
                excl = lbs[which(!(lbs %in% levels(factor(ties))))], 
                bond.type = bonds, size = length(unlist(stb)), 
                Note = note, prsep = prsep, Ties = stb, Attrs.ord = length(unlist(ats)), 
                Attrs = jnt(dhc(ats, prsep = prsep), prsep = prsep))), 
                RS <- (list(ord = dim(x)[1], nodes = lbs, sel = sel, 
                  sys.ord = nlevels(factor(ties)), incl = lbs[which(lbs %in% 
                    levels(factor(ties)))], excl = lbs[which(!(lbs %in% 
                    levels(factor(ties))))], bond.type = bonds, 
                  size = length(unlist(stb)), prsep = prsep, 
                  Ties = stb, Attrs.ord = length(unlist(ats)), 
                  Attrs = jnt(dhc(ats, prsep = prsep), prsep = prsep))))
        }
        else {
            ifelse(isTRUE(length(note) > 0L) == TRUE, RS <- (list(ord = dim(x)[1], 
                nodes = lbs, sel = sel, sys.ord = nlevels(factor(ties)), 
                incl = lbs[which(lbs %in% levels(factor(ties)))], 
                excl = lbs[which(!(lbs %in% levels(factor(ties))))], 
                bond.type = bonds, size = length(unlist(stb)), 
                Note = note, prsep = prsep, Ties = stb)), RS <- (list(ord = dim(x)[1], 
                nodes = lbs, sel = sel, sys.ord = nlevels(factor(ties)), 
                incl = lbs[which(lbs %in% levels(factor(ties)))], 
                excl = lbs[which(!(lbs %in% levels(factor(ties))))], 
                bond.type = bonds, size = length(unlist(stb)), 
                prsep = prsep, Ties = stb)))
        }
        class(RS) <- "Rel.System"
        return(RS)
    }
    else if (match.arg(type) == "listmat") {
        if (isTRUE(attr(x, "class") == "Rel.System") == FALSE) 
            stop("Relational system must be a \"Rel.System\" class.")
        if (isTRUE(x$sys.ord == 0L) == TRUE) 
            stop("Relational system chosen is empty!")
        lbst <- attr(x$Ties, "names")
        if (is.null(sel) == FALSE) {
            if (isTRUE(sel == 1L) == TRUE) {
                sel <- x$nodes[which(x$nodes %in% unlist(dhc(x$Attrs)))]
            }
            else if (isTRUE(sel == 0L) == TRUE) {
                sel <- x$nodes[which(!(x$nodes %in% unlist(dhc(x$Attrs))))]
            }
            else if (isTRUE(any(sel %in% x$nodes)) == FALSE) {
                return(sel)
            }
            else {
                NA
            }
            ntsel <- list()
            for (k in 1:length(x$Ties)) {
                tss <- which(dhc(x$Ties[[k]]) %in% sel)
                if (isTRUE(length(tss) > 0) == TRUE) {
                  tmpsel <- vector()
                  for (i in 1:length(tss)) {
                    if (isTRUE((tss[i]%%2L) == 1L) == TRUE) {
                      tmpsel <- append(tmpsel, x$Ties[[k]][ceiling(tss[i]/2L)])
                    }
                    else {
                      tmpsel <- append(tmpsel, x$Ties[[k]][floor(tss[i]/2L)])
                    }
                  }
                  rm(i)
                  ntsel[[k]] <- as.vector(unlist(tmpsel))
                }
                else {
                  lbst <- lbst[which(!(lbst %in% lbst[k]))]
                }
            }
            rm(k)
            rm(tmpsel, tss)
            ntsel <- ntsel[unlist(lapply(ntsel, length) != 0)]
            x$Ties <- ntsel
            lbs <- unique(dhc(unlist(ntsel)))
            n <- length(lbs)
            r <- length(lbst)
        }
        else if (is.null(sel) == TRUE) {
            n <- x$sys.ord
            r <- length(x$Ties)
            lbs <- x$incl
        }
        arr <- array(0, dim = c(n, n, r))
        for (i in 1:r) {
            if (isTRUE(length(x$Ties[[i]]) > 0) == TRUE) {
                arr[, , i] <- transf(x$Ties[[i]], type = "listmat", 
                  ord = n, labels = lbs)
            }
            else {
                NA
            }
        }
        rm(i)
        dimnames(arr)[[1]] <- dimnames(arr)[[2]] <- lbs
        dimnames(arr)[[3]] <- lbst
        if (is.null(x$Attrs) == FALSE) {
            arra <- array(0, dim = c(n, n, length(x$Attrs)))
            dimnames(arra)[[1]] <- dimnames(arra)[[2]] <- lbs
            dimnames(arra)[[3]] <- attr(x$Attrs, "names")
            for (i in 1:length(x$Attrs)) {
                act <- dhc(x$Attrs[[i]], prsep = prsep)
                if (isTRUE(length(act) > 0) == TRUE) {
                  diag(arra[, , i])[which(lbs %in% dhc(x$Attrs[[i]]))] <- 1L
                }
            }
            rm(i)
            attrs <- dim(arr)[3]
            arr <- zbind(arr, arra)
            class(arr) <- c("Rel.System", paste("Attrs.", paste(attrs + 
                1L, dim(arr)[3], sep = "="), sep = ": "))
        }
        return(arr)
    }
}
