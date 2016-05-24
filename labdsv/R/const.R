const <- function (taxa, clustering, minval = 0, show = minval, digits = 2, 
    sort = FALSE, spcord = NULL) 
{
    clustering <- clustify(clustering)

    res <- matrix(0, nrow = ncol(taxa), ncol = length(levels(clustering)))
    x <- apply(taxa,2,function(x){tapply(x>0,clustering,sum)})
    y <- as.numeric(table(clustering))
    res <- x/y
    keep <- as.logical(apply(res, 2, max) >= minval)
    res <- res[, keep]
    tmp <- as.data.frame(t(res))
    row.names(tmp) <- names(taxa)[keep]
    if (!is.null(spcord)) {
        tmp <- tmp[rev(order(spcord[keep])), ]
    }
    tmp <- format(round(tmp, digits = digits))
    tmp[tmp < show] <- substring(" .  ", 1, digits + 2)
    names(tmp) <- attr(clustering,'levels')

    if (sort) {
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
                  z <- ord[x]
                  ord <- ord[-x]
                  y <- match(pnt, row.names(tmp))
                  if (!is.na(y)) {
                    if (y > 1) {
                      first <- ord[1:(y - 1)]
                      last <- ord[y:length(ord)]
                      ord <- c(first, z, last)
                    } else {
                      last <- ord[y:length(ord)]
                      ord <- c(z, last)
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
        return(tmp)
    } else {
        return(tmp)
    }
}
