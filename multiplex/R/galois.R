galois <-
function (x, labeling = c("full", "reduced")) 
{
    ifelse(isTRUE(is.data.frame(x)) == TRUE, NA, x <- as.data.frame(x))
    eq0 <- list()
    for (i in which(duplicated(x))) {
        tmp <- vector()
        for (j in 1:nrow(x)) {
            if (isTRUE(all(x[i, ] == x[j, ]) == TRUE) == TRUE) {
                tmp <- append(tmp, rownames(x)[j])
            }
            else {
                NA
            }
        }
        rm(j)
        eq0[[length(eq0) + 1L]] <- paste(rownames(x)[i], jnt(tmp), 
            sep = ", ")
    }
    rm(i)
    if (isTRUE(length(eq0) > 0L) == TRUE) {
        for (i in 1:length(eq0)) {
            eq0[[i]] <- jnt(dhc(eq0[[i]]))
        }
        rm(i)
        eq0 <- unique(eq0)
        X <- x
        x <- unique(x)
    }
    else {
        NA
    }
    conj <- list()
    length(conj) <- ncol(x)
    for (i in 1:ncol(x)) {
        conj[[i]] <- jnt(rownames(x)[which(x[, i] != 0L)])
    }
    rm(i)
    attr(conj, "names") <- colnames(x)
    conj1 <- list()
    length(conj1) <- nrow(x)
    for (i in 1:nrow(x)) {
        conj1[[i]] <- jnt(colnames(x)[which(x[i, ] != 0L)])
    }
    rm(i)
    attr(conj1, "names") <- rownames(x)
    conj2 <- list()
    conj2n <- vector()
    if (isTRUE(length(conj) > 1L) == TRUE) {
        for (k in 2:length(conj)) {
            for (i in k:length(conj)) {
                if (isTRUE(length(intersect(dhc(conj[[k - 1]]), 
                  dhc(conj[[i]]))) > 1L) == TRUE) {
                  conj2[[length(conj2) + 1L]] <- jnt(intersect(dhc(conj[[k - 
                    1]]), dhc(conj[[i]])))
                }
                else if (isTRUE(length(intersect(dhc(conj[[k - 
                  1L]]), dhc(conj[[i]]))) > 1L) == FALSE) {
                  conj2[[length(conj2) + 1L]] <- intersect(dhc(conj[[k - 
                    1L]]), dhc(conj[[i]]))
                }
                conj2n[length(conj2n) + 1L] <- paste(attr(conj, 
                  "names")[k - 1L], attr(conj, "names")[i], sep = ", ")
            }
            rm(i)
        }
        rm(k)
        attr(conj2, "names") <- conj2n
    }
    else {
        conj2 <- conj
    }
    if (isTRUE(length(which(conj2 %in% conj == FALSE)) != 0L) == 
        TRUE) {
        if (isTRUE(jnt(rownames(x)) %in% c(conj, conj2)) == TRUE) {
            conj3 <- list(conj, unique(conj2[which(conj2 %in% 
                conj == FALSE)]))
        }
        else if (isTRUE(jnt(rownames(x)) %in% c(conj, conj2)) == 
            FALSE) {
            conj3 <- list(conj, unique(conj2[which(conj2 %in% 
                conj == FALSE)]), jnt(rownames(x)))
        }
    }
    else if (isTRUE(length(which(conj2 %in% conj == FALSE)) != 
        0L) == FALSE) {
        conj3 <- list(conj, jnt(rownames(x)))
    }
    for (i in 1:length(conj)) {
        if (isTRUE(length(which(conj2 %in% conj[[i]] == TRUE)) != 
            0L) == TRUE) {
            attr(conj3[[1]], "names")[i] <- jnt(attr(conj2, "names")[which(conj2 %in% 
                conj[[i]])])
        }
        else {
            NA
        }
    }
    rm(i)
    for (i in 1:length(conj3[[2]])) {
        if (isTRUE(length(attr(conj2, "names")[which(conj3[[2]][[i]] == 
            conj2)]) != 0L) == TRUE) {
            if (isTRUE(length(dhc(conj3[[2]][[i]])) == 1L) == 
                TRUE) {
                attr(conj3[[2]], "names")[i] <- conj1[[which(conj3[[2]][[i]] == 
                  attr(conj1, "names"))]]
            }
            else if (isTRUE(length(dhc(conj3[[2]][[i]])) == 1L) == 
                FALSE) {
                attr(conj3[[2]], "names")[i] <- jnt(attr(conj2, 
                  "names")[which(conj3[[2]][[i]] == conj2)])
            }
        }
        else {
            attr(conj3[[2]], "names")[i] <- jnt(attr(conj2, "names")[which(conj2 %in% 
                levels(factor(unlist(conj2))))])
        }
    }
    rm(i)
    if (isTRUE(jnt(rownames(x)) %in% c(conj, conj2)) == FALSE) {
        ifelse(isTRUE(length(conj3) > 2L) == TRUE, con <- c(conj3[[1]], 
            conj3[[2]], conj3[[3]]), con <- c(conj3[[1]], conj3[[2]]))
    }
    else {
        con <- c(conj3[[1]], conj3[[2]])
    }
    der <- unique(con)
    attr(der, "names") <- unique(attr(con, "names"))
    po <- matrix(0L, nrow = length(der), ncol = length(der))
    for (j in 1:length(der)) {
        for (i in 1:length(der)) {
            ifelse(isTRUE(all(dhc(der[[i]]) %in% dhc(der[[j]]))) == 
                TRUE, po[i, j] <- 1L, NA)
        }
    }
    rm(i, j)
    ints <- attr(der, "names")
    exts <- der
    rownames(po) <- colnames(po) <- ints
    for (i in (length(conj3[[1]]) + 1L):length(der)) {
        if (isTRUE(length(der) > length(conj3[[1]])) == TRUE) {
            ifelse(isTRUE(length(flt(i, po, rclos = FALSE)) > 
                0L) == TRUE, ints[i] <- jnt(dhc(ints[flt(i, po, 
                rclos = FALSE)])), NA)
        }
        else {
            NA
        }
    }
    rm(i)
    attr(der, "names") <- attr(exts, "names") <- ints
    if (isTRUE(length(der) > 2L) == TRUE) {
        for (k in 2:length(der)) {
            for (i in k:length(der)) {
                if (isTRUE((k - 1L) == i) == FALSE) {
                  if (isTRUE(any(isTRUE(all(po[, i] - po[, (k - 
                    1L)] != -1)) == TRUE | isTRUE(all(po[, (k - 
                    1L)] - po[, i] != -1)) == TRUE)) == TRUE) {
                    if (isTRUE(all(po[, i] - po[, (k - 1L)] != 
                      -1)) == TRUE) {
                      if (isTRUE(ints[(k - 1L)] == "") == FALSE) {
                        ints[(k - 1L)] <- jnt(dhc(ints[(k - 1L)])[which(!(dhc(ints[(k - 
                          1L)]) %in% dhc(ints[i])))])
                      }
                      else {
                        NA
                      }
                      if (isTRUE(length(exts[[i]]) == 0L) == 
                        FALSE) {
                        exts[[i]] <- jnt(dhc(exts[[i]])[which(!(dhc(exts[[i]]) %in% 
                          dhc(exts[[(k - 1L)]])))])
                      }
                      else {
                        NA
                      }
                    }
                    else if (isTRUE(all(po[, i] - po[, (k - 1L)] != 
                      -1)) == FALSE) {
                      if (isTRUE(all.equal(dhc(ints[i]), dhc(ints[(k - 
                        1L)]))) == TRUE) {
                        ints[i] <- ""
                      }
                      else if (isTRUE(all.equal(dhc(ints[i]), 
                        dhc(ints[(k - 1L)]))) == FALSE) {
                        ifelse(isTRUE(length(jnt(dhc(ints[i])[which(!(dhc(ints[i]) %in% 
                          dhc(ints[(k - 1L)])))])) == 0L) == 
                          TRUE, NA, ints[i] <- jnt(dhc(ints[i])[which(!(dhc(ints[i]) %in% 
                          dhc(ints[(k - 1L)])))]))
                      }
                      if (isTRUE(all.equal(dhc(exts[[(k - 1L)]]), 
                        dhc(exts[[i]]))) == TRUE) {
                        exts[[(k - 1L)]] <- ""
                      }
                      else if (isTRUE(all.equal(dhc(exts[[(k - 
                        1L)]]), dhc(exts[[i]]))) == FALSE) {
                        exts[[(k - 1L)]] <- jnt(dhc(exts[[(k - 
                          1L)]])[which(!(dhc(exts[[(k - 1L)]]) %in% 
                          dhc(exts[[i]])))])
                      }
                    }
                  }
                  else {
                    NA
                  }
                }
                else {
                  NA
                }
            }
            rm(i)
        }
        rm(k)
    }
    else if (isTRUE(length(der) == 2L) == TRUE) {
        if (isTRUE(any(isTRUE(all(po[, 2] - po[, 1] != -1)) == 
            TRUE | isTRUE(all(po[, 1] - po[, 2] != -1)) == TRUE)) == 
            TRUE) {
            if (isTRUE(all(po[, 2] - po[, 1] != -1)) == TRUE) {
                if (isTRUE(is.na(ints[2])) == FALSE) {
                  ints[2] <- jnt(dhc(ints[1])[which(!(dhc(ints[1]) %in% 
                    dhc(ints[2])))])
                }
                else {
                  NA
                }
                if (isTRUE(length(exts[[2]]) == 0L) == FALSE) {
                  exts[[2]] <- jnt(dhc(exts[[2]])[which(!(dhc(exts[[2]]) %in% 
                    dhc(exts[[1]])))])
                }
                else {
                  NA
                }
            }
            else if (isTRUE(all(po[, 2] - po[, 1] != -1)) == 
                FALSE) {
                if (isTRUE(all.equal(dhc(ints[2]), dhc(ints[1]))) == 
                  TRUE) {
                  ints[2] <- ""
                }
                else if (isTRUE(all.equal(dhc(ints[2]), dhc(ints[1]))) == 
                  FALSE) {
                  ifelse(isTRUE(length(jnt(dhc(ints[2])[which(!(dhc(ints[2]) %in% 
                    dhc(ints[1])))])) == 0L) == TRUE, NA, ints[2] <- jnt(dhc(ints[2])[which(!(dhc(ints[2]) %in% 
                    dhc(ints[1])))]))
                }
                if (isTRUE(all.equal(dhc(exts[[1]]), dhc(exts[[2]]))) == 
                  TRUE) {
                  exts[[1]] <- ""
                }
                else if (isTRUE(all.equal(dhc(exts[[1]]), dhc(exts[[2]]))) == 
                  FALSE) {
                  exts[[1]] <- jnt(dhc(exts[[1]])[which(!(dhc(exts[[1]]) %in% 
                    dhc(exts[[2]])))])
                }
            }
        }
        else {
            NA
        }
    }
    else {
        NA
    }
    dupl <- levels(factor(unlist(dhc(exts))[which(duplicated(unlist(dhc(exts))) == 
        TRUE)]))
    dder <- der
    if (isTRUE(length(dupl) > 0L) == TRUE) {
        for (i in 1:length(dupl)) dder[[(length(dder) + 1L)]] <- dupl[i]
        for (i in 1:length(dupl)) {
            for (j in 1:nrow(x)) {
                if (isTRUE(any(isTRUE(all(x[j, ] - x[which(rownames(x) == 
                  dupl[i]), ] != -1)) == TRUE & isTRUE(all(x[which(rownames(x) == 
                  dupl[i]), ] - x[j, ] != -1)) == FALSE)) == 
                  TRUE) {
                  dder[[(length(der) + i)]] <- jnt(append(dder[[(length(der) + 
                    i)]], rownames(x)[j]))
                }
                else {
                  NA
                }
            }
            rm(j)
        }
        rm(i)
        dpo <- matrix(0L, nrow = length(dder), ncol = length(dder))
        for (j in 1:length(dder)) {
            for (i in 1:length(dder)) {
                ifelse(isTRUE(all(dhc(dder[[i]]) %in% dhc(dder[[j]]))) == 
                  TRUE, dpo[i, j] <- 1L, NA)
            }
        }
        rm(i, j)
        attr(dder, "names")[which(is.na(attr(dder, "names")) == 
            TRUE)] <- ""
        ints <- attr(dder, "names")
        exts <- dder
        rownames(dpo) <- colnames(dpo) <- ints
        PO <- dpo
        for (i in (length(der) + 1L):length(dder)) {
            if (isTRUE(length(dder) > length(der)) == TRUE) {
                ifelse(isTRUE(length(flt(i, dpo, rclos = FALSE)) > 
                  0L) == TRUE, ints[i] <- jnt(dhc(ints[flt(i, 
                  dpo, rclos = FALSE)])), NA)
            }
            else {
                NA
            }
        }
        rm(i)
        attr(dder, "names") <- ints
        for (k in 2:length(dder)) {
            for (i in k:length(dder)) {
                if (isTRUE((k - 1L) == i) == FALSE) {
                  if (isTRUE(any(isTRUE(all(dpo[, i] - dpo[, 
                    (k - 1L)] != -1)) == TRUE | isTRUE(all(dpo[, 
                    (k - 1L)] - dpo[, i] != -1)) == TRUE)) == 
                    TRUE) {
                    if (isTRUE(all(dpo[, i] - dpo[, (k - 1L)] != 
                      -1)) == TRUE) {
                      if (isTRUE(ints[(k - 1L)] == "") == FALSE) {
                        ints[(k - 1L)] <- jnt(dhc(ints[(k - 1L)])[which(!(dhc(ints[(k - 
                          1L)]) %in% dhc(ints[i])))])
                      }
                      else {
                        NA
                      }
                      if (isTRUE(length(exts[[i]]) == 0L) == 
                        FALSE) {
                        exts[[i]] <- jnt(dhc(exts[[i]])[which(!(dhc(exts[[i]]) %in% 
                          dhc(exts[[(k - 1L)]])))])
                      }
                      else {
                        NA
                      }
                    }
                    else if (isTRUE(all(dpo[, i] - dpo[, (k - 
                      1L)] != -1)) == FALSE) {
                      if (isTRUE(all.equal(dhc(ints[i]), dhc(ints[(k - 
                        1L)]))) == TRUE) {
                        ints[i] <- ""
                      }
                      else if (isTRUE(all.equal(dhc(ints[i]), 
                        dhc(ints[(k - 1L)]))) == FALSE) {
                        ifelse(isTRUE(length(jnt(dhc(ints[i])[which(!(dhc(ints[i]) %in% 
                          dhc(ints[(k - 1L)])))])) == 0L) == 
                          TRUE, NA, ints[i] <- jnt(dhc(ints[i])[which(!(dhc(ints[i]) %in% 
                          dhc(ints[(k - 1L)])))]))
                      }
                      if (isTRUE(all.equal(dhc(exts[[(k - 1L)]]), 
                        dhc(exts[[i]]))) == TRUE) {
                        exts[[(k - 1L)]] <- ""
                      }
                      else if (isTRUE(all.equal(dhc(exts[[(k - 
                        1L)]]), dhc(exts[[i]]))) == FALSE) {
                        exts[[(k - 1L)]] <- jnt(dhc(exts[[(k - 
                          1L)]])[which(!(dhc(exts[[(k - 1L)]]) %in% 
                          dhc(exts[[i]])))])
                      }
                    }
                  }
                  else {
                    NA
                  }
                }
                else {
                  NA
                }
            }
            rm(i)
        }
        rm(k)
    }
    else if (isTRUE(length(dupl) > 0L) == FALSE) {
        PO <- po
    }
    if (isTRUE(length(eq0) > 0L) == TRUE) {
        for (k in 1:length(eq0)) {
            cmb <- which(exts %in% dhc(eq0)[[k]])
            exts[[cmb]] <- eq0[[k]]
            for (i in 1:length(dder)) {
                ifelse(isTRUE(any(dhc(dder)[[i]] %in% dhc(eq0)[[k]]) == 
                  TRUE) == TRUE, dder[[i]] <- jnt(c(dhc(dder)[[i]], 
                  dhc(eq0)[[k]]), prsep = ", "), NA)
            }
            rm(i)
        }
        rm(k, cmb)
    }
    mi <- NULL
    for (i in 1:dim(PO)[1]) {
        flt <- flt(i, PO, rclos = TRUE)
        ifelse(isTRUE(length(flt) == dim(PO)[1]) == TRUE, mi <- i, 
            NA)
    }
    rm(i)
    if (isTRUE(jnt(colnames(x)) %in% attr(dder, "names")) == 
        FALSE) {
        if (isTRUE(length(mi) == 0L) == TRUE) {
            dder[[length(dder) + 1L]] <- exts[[length(exts) + 
                1L]] <- ints[length(ints) + 1L] <- ""
            attr(dder, "names")[length(dder)] <- jnt(colnames(x))
        }
        else {
            attr(dder, "names")[mi] <- jnt(colnames(x))
        }
    }
    else {
        NA
    }
    attr(dder, "names")[which(is.na(attr(dder, "names")) == TRUE)] <- ""
    class(dder) <- c("Galois", "full")
    derr <- exts
    attr(derr, "names") <- ints
    attr(derr, "names")[which(is.na(attr(derr, "names")) == TRUE)] <- ""
    lst <- (red = redl(dder, derr))
    class(lst) <- c("Galois", "reduced")
    switch(match.arg(labeling), full = dder, reduced = lst)
}
