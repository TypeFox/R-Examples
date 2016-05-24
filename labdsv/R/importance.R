importance <- function (taxa, clustering, minval = 0, digits = 2, show = minval, 
    sort = FALSE, typical=TRUE, spcord) 
{
    clustering <- clustify(clustering) 

    res <- matrix(0, nrow = ncol(taxa), ncol = length(levels(clustering)))
    x <- apply(taxa,2,function(x){tapply(x,clustering,sum)})
    if (typical) {
        y <- apply(taxa,2,function(x){tapply(x>0,clustering,sum)})
    } else {
        y <- apply(taxa,2,function(x){tapply(x>=0,clustering,sum)})
    }
    y[x==0] <- 1
    res <- x/y
    keep <- as.logical(apply(res, 2, max) >= minval)
    res <- res[, keep]
    tmp <- as.data.frame(t(res))
    row.names(tmp) <- names(taxa)[keep]
    if (!missing(spcord)) {
        tmp <- tmp[rev(order(spcord[keep])), ]
    }

    tmp <- format(round(tmp, digits = digits))
    tmp[tmp < show] <- substring(" .  ", 1, nchar(tmp[1, 1]))
    if (sort) {
        cat("\nConstancy Table\n\n")
        print(tmp)
        repeat {
            plots <- readline(" enter the species: ")
            if (plots == "") {
                break
            } else {
                pnt <- readline(" in front of        : ")
            }
            for (i in (strsplit(plots, ",")[[1]])) {
                ord <- 1:nrow(tmp)
                x <- match(i, row.names(tmp))
                if (!is.na(x)) {
                  ord <- ord[-x]
                  y <- match(pnt, row.names(tmp[ord, ]))
                  if (!is.na(y)) {
                    if (y == 1) {
                      ord <- c(x, ord)
                    } else {
                      first <- ord[1:(y - 1)]
                      last <- ord[y:length(ord)]
                      ord <- c(first, x, last)
                    }
                    tmp <- tmp[ord, ]
                    print(tmp)
                  } else {
                    print(paste("species", pnt, "does not exist"))
                  }
                } else {
                  print(paste("species", i, "does not exist"))
                }
            }
        }
        attr(tmp,'call') <- match.call()
        attr(tmp,'taxa') <- deparse(substitute(taxa))
        return(tmp)
    } else {
        attr(tmp,'call') <- match.call()
        attr(tmp,'taxa') <- deparse(substitute(taxa))
        return(tmp)
    }
}
