"stats.bplot" <-

function (x, by, style = "tukey", outlier = TRUE) 



{

    if (!missing(by)) {

        x <- cat.to.list(c(x), by)

    }

    if (!is.list(x) & !is.matrix(x)) 

        x <- matrix(x, ncol = 1)

    if (is.list(x)) {

        ncol <- length(x)

        out <- as.list(1:ncol)

        names(out) <- names(x)

        for (j in (1:ncol)) {

            if (is.numeric(x[[j]])) {

                out[[j]] <- describe.bplot(x[[j]], style = style, 

                  outlier = outlier)

            }

        }

        return(out)

    }

    if (is.matrix(x)) {

        nc <- ncol(x)

        out <- as.list(1:nc)

        names(out) <- dimnames(x)[[2]]

        for (j in (1:nc)) {

            out[[j]] <- describe.bplot(x[, j], style = style, 

                outlier = outlier)

        }

        return(out)

    }



}

