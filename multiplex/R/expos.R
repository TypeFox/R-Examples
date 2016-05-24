expos <-
function (rs, classes = FALSE, allClasses = FALSE, allNodes = TRUE) 
{
    if (isTRUE(attr(rs, "class") == "Rel.System") == FALSE) 
        stop("Relational system must be a \"Rel.System\" class.")
    if (isTRUE(rs$sys.ord == 0L) == TRUE) 
        stop("Relational system chosen is empty!")
    if (isTRUE(rs$Attrs.ord == 0L) == TRUE | is.null(rs$Attrs.ord) == 
        TRUE) 
        return("There are no attributes in the relational system.")
    if (isTRUE(any(duplicated(attr(rs$Attrs, "names")))) == TRUE) {
        rs$attrs <- rs$Attrs[which(!(duplicated(attr(rs$Attrs, 
            "names"))))]
        dpl <- which(duplicated(attr(rs$Attrs, "names")))
        for (i in 1:length(dpl)) {
            rs$attrs[which(attr(rs$attrs, "names") == attr(rs$Attrs, 
                "names")[dpl[i]])] <- jnt(dhc(unlist(rs$attrs[which(attr(rs$attrs, 
                "names") == attr(rs$Attrs, "names")[dpl[i]])], 
                rs$Attrs[dpl[i]]), prsep = rs$prsep), prsep = rs$prsep)
        }
        rm(i)
    }
    else {
        rs$attrs <- rs$Attrs
    }
    rs$attrs2 <- list()
    attrs2 <- vector()
    for (i in 1:length(rs$attrs)) {
        if (isTRUE(is.null(rs$attrs[[i]])) == FALSE) {
            rs$attrs2[[length(rs$attrs2) + 1L]] <- rs$attrs[[i]]
            attrs2 <- append(attrs2, attr(rs$attrs, "names")[i])
        }
        else {
            NA
        }
    }
    rm(i)
    attr(rs$attrs2, "names") <- attrs2
    if (isTRUE(length(rs$attrs2) == 0L) == TRUE) 
        stop("There are no attributes in the relational system.")
    natt <- length(unique(attr(rs$attrs2, "names")))
    At <- data.frame(matrix(0L, ncol = rs$ord, nrow = natt))
    colnames(At) <- rs$nodes
    rownames(At) <- unique(attr(rs$attrs2, "names"))
    for (i in 1:natt) {
        At[i, which(colnames(At) %in% dhc(rs$attrs2[[i]], prsep = rs$prsep))] <- 1L
    }
    rm(i)
    if (is.null(rs$incl) == FALSE) {
        ifelse(isTRUE(is.null(nrow(At)) == TRUE) == TRUE, at <- At[which(attr(At, 
            "names") %in% rs$incl)], at <- as.data.frame(At)[, 
            which(colnames(At) %in% rs$incl)])
    }
    if (isTRUE(nrow(at) > 0L) == TRUE) {
        Adpt <- list()
        Adpt[[1]] <- colnames(at)[which(apply(at, 2, sum) == 
            0L)]
        Adpt[[2]] <- colnames(at)[which(apply(at, 2, sum) == 
            nrow(at))]
        if (isTRUE(nrow(at) > 1L) == TRUE) {
            rst <- which(!(colnames(at) %in% unlist(Adpt)))
            if (isTRUE(length(rst) >= 1L) == TRUE) {
                aat <- at[, rst]
                bnch <- as.data.frame(t(expand.grid(rep(list(0L:1L), 
                  natt))))
                rownames(bnch) <- unique(attr(rs$attrs2, "names"))
                colnames(bnch)[1] <- "Null"
                for (i in 2:ncol(bnch)) colnames(bnch)[i] <- jnt(rownames(bnch)[which(bnch[, 
                  i] == 1L)], prsep = ", ")
                rm(i)
                bn <- bnch[, 2:(ncol(bnch) - 1L)]
                tmp <- vector()
                for (k in 1:ncol(bn)) {
                  if (isTRUE(is.null(ncol(aat)) == FALSE) == 
                    TRUE) {
                    for (i in 1:ncol(aat)) {
                      if (isTRUE(all(bn[, k] == aat[, i])) == 
                        TRUE) {
                        tmp <- append(tmp, colnames(aat)[i])
                      }
                      else {
                        NA
                      }
                    }
                    rm(i)
                  }
                  else {
                    ifelse(isTRUE(all(as.vector(bn[, k]) == aat)) == 
                      TRUE, tmp <- append(tmp, colnames(aat)[rst]), 
                      NA)
                  }
                  Adpt[[k + 2L]] <- as.character(tmp)
                  tmp <- vector()
                }
                rm(k)
                attr(Adpt, "names") <- c("NONE", "ALL", colnames(bn))
            }
            else {
                NA
            }
            adpt <- Adpt[c(1, 3:length(Adpt))]
        }
        else if (isTRUE(nrow(at) == 1L) == TRUE) {
            attr(Adpt, "names") <- c("NONE", "ALL")
            adpt <- Adpt[1]
        }
    }
    else if (isTRUE(nrow(at) > 0L) == FALSE) {
        Adpt <- NULL
    }
    if (isTRUE(length(unlist(adpt)) > 0L) == TRUE) {
        ladpt <- list()
        for (l in 1:length(adpt)) {
            if (isTRUE(length(adpt[[l]]) != 0L) == TRUE) {
                ex <- list()
                for (i in 1:length(adpt[[l]])) {
                  ngbsadpt <- ngbs(adpt[[l]][i], rs, type = "und")
                  if (isTRUE(length(ngbsadpt) == 1L) == TRUE) {
                    ex <- at[, which(colnames(at) %in% ngbsadpt)]
                    ex <- dichot(ex - at[, which(adpt[[l]][i] == 
                      colnames(at))])
                    attr(ex, "names") <- rownames(at)
                  }
                  else {
                    tmp <- at[, which(colnames(at) %in% ngbsadpt)]
                    tmp <- dichot(tmp - at[, which(adpt[[l]][i] == 
                      colnames(at))])
                    ex <- apply(tmp, 1, sum)/ncol(tmp)
                  }
                  ladpt[[length(ladpt) + 1L]] <- ex
                }
                rm(i)
            }
            else {
                NA
            }
        }
        rm(l)
        attr(ladpt, "names") <- as.vector(unlist(adpt))
        uladpt <- unlist(ladpt)
        slct <- (1:length(uladpt))%%nrow(at)
        slct <- replace(slct, slct == 0L, nrow(at))
        Exps <- list()
        for (i in 1:nrow(at)) {
            temp <- uladpt[slct == i]
            Exps[[i]] <- round(temp[which(temp > 0L)], 2)
        }
        rm(i)
        exx <- list()
        exx2 <- list()
        for (i in 1:length(Exps)) {
            if (isTRUE(length(Exps[[i]]) > 0L) == TRUE) {
                exx[[i]] <- dhc(names(Exps[[i]]), "[.]")[which(1:(length(Exps[[i]]) * 
                  2L)%%2L == 1L)]
                exx2[[i]] <- noquote(as.vector(Exps[[i]]))
            }
            else {
                exx2[[i]] <- exx[[i]] <- NULL
            }
        }
        rm(i)
        exps <- list()
        for (i in 1:length(Exps)) {
            if (isTRUE(length(Exps[[i]]) > 0L) == TRUE) {
                exps[[i]] <- paste(exx[[i]], paste0(exx2[[i]] * 
                  100L, "%"), sep = "=")
            }
        }
        rm(i)
        attr(exps, "names") <- paste0("to_", rownames(at))
        exps <- noquote(exps)
    }
    else {
        exps <- NULL
    }
    if (classes) {
        if (allNodes) {
            Adpt[[1]] <- colnames(At)[which(apply(At, 2, sum) == 
                0L)]
            Adpt[[2]] <- colnames(At)[which(apply(At, 2, sum) == 
                nrow(At))]
            rst <- which(!(colnames(At) %in% unlist(Adpt)))
            if (isTRUE(length(rst) >= 1L) == TRUE) {
                Aat <- At[, rst]
                tmp <- vector()
                for (k in 1:ncol(bn)) {
                  if (isTRUE(is.null(ncol(Aat)) == FALSE) == 
                    TRUE) {
                    for (i in 1:ncol(Aat)) {
                      if (isTRUE(all(bn[, k] == Aat[, i])) == 
                        TRUE) {
                        tmp <- append(tmp, colnames(Aat)[i])
                      }
                      else {
                        NA
                      }
                    }
                    rm(i)
                  }
                  else {
                    ifelse(isTRUE(all(as.vector(bn[, k]) == Aat)) == 
                      TRUE, tmp <- append(tmp, colnames(At)[rst]), 
                      NA)
                  }
                  Adpt[[k + 2L]] <- as.vector(c(Adpt[[k + 2L]], 
                    tmp))
                  tmp <- vector()
                }
                rm(k)
            }
            else {
                NA
            }
        }
        if (allClasses) {
            clss <- Adpt
        }
        else {
            clss <- list()
            ncls <- vector()
            for (i in 1:length(Adpt)) {
                if (isTRUE(length(Adpt[[i]]) > 0L) == TRUE) {
                  clss[[length(clss) + 1L]] <- Adpt[[i]]
                  ncls <- append(ncls, attr(Adpt, "names")[i])
                }
                else {
                  NA
                }
            }
            rm(i)
            attr(clss, "names") <- ncls
        }
        return(list(Classes = noquote(clss), Bonds=rs$bond.type, 
          Exposure = exps))
    }
    else {
        return(list(Bonds=rs$bond.type, Exposure = exps))
    }
}
