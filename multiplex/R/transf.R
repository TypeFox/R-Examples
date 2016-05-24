transf <-
function (x, type = c("matlist", "listmat"), lb2lb = FALSE, labels = NULL, 
    prsep = ", ", ord = NULL) 
{
    if (isTRUE(is.character(x) == TRUE) == TRUE) 
        type <- "listmat"
    if (match.arg(type) == "matlist") {
        if (is.na(dim(x)[3]) == FALSE) 
            stop("Use the \"rel.sys\" function for 3D arrays.")
        if (isTRUE(is.matrix(x) == TRUE) == FALSE) 
            x <- as.matrix(x)
        if (isTRUE(lb2lb == TRUE) == TRUE) {
            if (isTRUE(is.null(labels) == TRUE) == TRUE) {
                if (isTRUE(is.null(dimnames(x)[[1]]) == TRUE | 
                  is.null(dimnames(x)[[2]]) == TRUE) == TRUE) 
                  stop("To use the \"lb2lb\" option you need to specify the labels.")
                labelsr <- dimnames(x)[[1]]
                labelsc <- dimnames(x)[[2]]
            }
            else {
                labelsr <- labelsc <- labels
            }
        }
        else {
            NA
        }
        if (isTRUE(sum(x) > 0) == TRUE) {
            inc <- list()
            rws <- vector()
            cls <- vector()
            for (k in 1:max(x)) {
                X <- dichot(x, c = k)
                for (i in 1:length(which((X) == 1))) {
                  cls[i] <- (ceiling(which((X) == 1)/dim(x)[1]))[i]
                  ifelse((which((X) == 1)%%dim(x)[1])[i] == 0, 
                    rws[i] <- (which((X) == 1)%%dim(x)[1])[i] + 
                      dim(x)[1], rws[i] <- (which((X) == 1)%%dim(x)[1])[i])
                  ifelse(isTRUE(lb2lb == TRUE) == TRUE, inc[[length(inc) + 
                    1]] <- paste(labelsr[rws[i]], labelsc[cls[i]], 
                    sep = prsep), inc[[length(inc) + 1]] <- paste(rws[i], 
                    cls[i], sep = prsep))
                }
                rm(i)
            }
            rm(k)
            return(sort(unlist(inc)))
        }
        else {
            return(paste(0, 0, sep = prsep))
        }
    }
    if (match.arg(type) == "listmat") {
        if (is.character(x) == FALSE) 
            stop("Input data has to be a list with character format.")
        if (isTRUE(lb2lb == TRUE) == TRUE) {
        }
        if (isTRUE(is.null(labels) == TRUE) == TRUE) {
            if (isTRUE(is.null(ord) == FALSE) == TRUE) {
                mat <- matrix(0, ncol = ord, nrow = ord)
            }
            else if (isTRUE(is.null(ord) == TRUE) == TRUE) {
                mat <- matrix(0, nrow = length(dhc(jnt(x))), 
                  ncol = length(dhc(jnt(x))))
            }
            tmp <- as.list(unlist(x))
            nx <- vector()
            for (i in 1:length(tmp)) {
                if (isTRUE(length(tmp[[i]]) > 0) == TRUE) {
                  for (j in 1:length(tmp[[i]])) {
                    nx <- append(nx, strsplit(tmp[[i]][j], prsep)[[1]][1])
                    nx <- append(nx, strsplit(tmp[[i]][j], prsep)[[1]][2])
                  }
                }
            }
            rm(i)
            rm(tmp)
            lbs <- levels(factor(nx))
        }
        else {
            mat <- matrix(0, ncol = length(labels), nrow = length(labels))
            lbs <- labels
        }
        for (i in 1:length(lbs)) if (isTRUE(is.numeric(lbs[i]) == 
            FALSE) == TRUE) {
            lb2lb <- TRUE
            break
        }
        else {
            NA
        }
        for (i in 1:length(x)) {
            ifelse(isTRUE(lb2lb == TRUE) == TRUE, mat[which(strsplit(x[i], 
                prsep)[[1]][1] == lbs), which(strsplit(x[i], 
                prsep)[[1]][2] == lbs)] <- (mat[which(strsplit(x[i], 
                prsep)[[1]][1] == lbs), which(strsplit(x[i], 
                prsep)[[1]][2] == lbs)] + 1), mat[as.numeric(strsplit(x[i], 
                prsep)[[1]])[1], as.numeric(strsplit(x[i], prsep)[[1]])[2]] <- (mat[as.numeric(strsplit(x[i], 
                prsep)[[1]])[1], as.numeric(strsplit(x[i], prsep)[[1]])[2]] + 
                1))
        }
        rm(i)
        if (isTRUE(length(lbs) == nrow(mat)) == TRUE) 
            rownames(mat) <- colnames(mat) <- lbs
        return(mat)
    }
}
