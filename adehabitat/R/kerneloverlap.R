kerneloverlap <- function(xy, id = NULL,
                          method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"),
                          lev=95, conditional=FALSE, ...)
{
    ## Verifications
    method <- match.arg(method)

    ## UD estimation
    x <- kernelUD(xy, id, same4all=TRUE, ...)
    vol <- getvolumeUD(x)

    ## Matrix of results
    res <- matrix(0, ncol=length(x), nrow=length(x))

    ## loop for each animal
    for (i in 1:length(x)) {
        for (j in 1:i) {

            if (method=="HR") {
                vi <- vol[[i]]$UD
                vj <- vol[[j]]$UD
                vi[vi<=lev] <- 1
                vi[vi>lev] <- 0
                vj[vj<=lev] <- 1
                vj[vj>lev] <- 0
                vk <- vi*vj
                res[i,j] <- sum(vk)/sum(vi)
                res[j,i] <- sum(vk)/sum(vj)
            }

            if (method=="PHR") {
                vi <- x[[i]]$UD
                vj <- x[[j]]$UD
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[j,i] <- sum(vi*aj)*(attr(vi,"cellsize")^2)
                    res[i,j] <- sum(vj*ai)*(attr(vi,"cellsize")^2)
                } else {
                    res[j,i] <- sum(vi*aj)*(attr(vi,"cellsize")^2)
                    res[i,j] <- sum(vj*ai)*(attr(vi,"cellsize")^2)
                }
            }



            if (method=="VI") {
                vi <- c(x[[i]]$UD)
                vj <- c(x[[j]]$UD)
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[i,j] <- res[j,i] <- sum(pmin(vi, vj))*(attr(x[[i]]$UD,"cellsize")^2)
                } else {
                    res[i,j] <- res[j,i] <- sum(pmin(vi, vj))*(attr(x[[i]]$UD,"cellsize")^2)
                }
            }

            if (method=="BA") {
                vi <- x[[i]]$UD
                vj <- x[[j]]$UD
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[j,i] <- res[i,j] <- sum(sqrt(vi)*sqrt(vj))*(attr(vi,"cellsize")^2)
                } else {
                    res[j,i] <- res[i,j] <- sum(sqrt(vi)*sqrt(vj))*(attr(vi,"cellsize")^2)
                }
            }

            if (method=="UDOI") {
                vi <- x[[i]]$UD
                vj <- x[[j]]$UD
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    ak <- sum(ai*aj)*(attr(vi,"cellsize")^2)
                    res[j,i] <- res[i,j] <- ak * sum(vi*vj)*(attr(vi,"cellsize")^2)
                } else {
                    ak <- sum(ai*aj)*(attr(vi,"cellsize")^2)
                    res[j,i] <- res[i,j] <- ak * sum(vi*vj)*(attr(vi,"cellsize")^2)
                }
            }

            if (method=="HD") {
                vi <- x[[i]]$UD
                vj <- x[[j]]$UD
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[j,i] <- res[i,j] <- sqrt(sum((sqrt(vi) - sqrt(vj))^2*(attr(vi,"cellsize")^2)))
                } else {
                    res[j,i] <- res[i,j] <- sqrt(sum((sqrt(vi) - sqrt(vj))^2*(attr(vi,"cellsize")^2)))
                }
            }
        }
    }
    rownames(res) <- names(x)
    colnames(res) <- names(x)
    return(res)
}




kerneloverlaphr <- function(x,
                            method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"),
                            lev=95, conditional=FALSE)
{
    ## Verifications
    method <- match.arg(method)

    if ((!inherits(x, "khrud"))&(!inherits(x, "kbbhrud")))
        stop("x should be of class khrud or kbbhrud")

    liii <- lapply(x, function(x) x$UD)
    names(liii) <- letters[1:length(liii)]
    verif <- as.kasc(liii)

    ## UD estimation
    vol <- getvolumeUD(x)

    ## Matrix of results
    res <- matrix(0, ncol=length(x), nrow=length(x))

    ## loop for each animal
    for (i in 1:length(x)) {
        for (j in 1:i) {

            if (method=="HR") {
                vi <- vol[[i]]$UD
                vj <- vol[[j]]$UD
                vi[vi<=lev] <- 1
                vi[vi>lev] <- 0
                vj[vj<=lev] <- 1
                vj[vj>lev] <- 0
                vk <- vi*vj
                res[i,j] <- sum(vk)/sum(vi)
                res[j,i] <- sum(vk)/sum(vj)
            }

            if (method=="PHR") {
                vi <- x[[i]]$UD
                vj <- x[[j]]$UD
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[j,i] <- sum(vi*aj)*(attr(vi,"cellsize")^2)
                    res[i,j] <- sum(vj*ai)*(attr(vi,"cellsize")^2)
                } else {
                    res[j,i] <- sum(vi*aj)*(attr(vi,"cellsize")^2)
                    res[i,j] <- sum(vj*ai)*(attr(vi,"cellsize")^2)
                }
            }



            if (method=="VI") {
                vi <- c(x[[i]]$UD)
                vj <- c(x[[j]]$UD)
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[i,j] <- res[j,i] <- sum(pmin(vi, vj))*(attr(x[[i]]$UD,"cellsize")^2)
                } else {
                    res[i,j] <- res[j,i] <- sum(pmin(vi, vj))*(attr(x[[i]]$UD,"cellsize")^2)
                }
            }

            if (method=="BA") {
                vi <- x[[i]]$UD
                vj <- x[[j]]$UD
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[j,i] <- res[i,j] <- sum(sqrt(vi)*sqrt(vj))*(attr(vi,"cellsize")^2)
                } else {
                    res[j,i] <- res[i,j] <- sum(sqrt(vi)*sqrt(vj))*(attr(vi,"cellsize")^2)
                }
            }

            if (method=="UDOI") {
                vi <- x[[i]]$UD
                vj <- x[[j]]$UD
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    ak <- sum(ai*aj)*(attr(vi,"cellsize")^2)
                    res[j,i] <- res[i,j] <- ak * sum(vi*vj)*(attr(vi,"cellsize")^2)
                } else {
                    ak <- sum(ai*aj)*(attr(vi,"cellsize")^2)
                    res[j,i] <- res[i,j] <- ak * sum(vi*vj)*(attr(vi,"cellsize")^2)
                }
            }

            if (method=="HD") {
                vi <- x[[i]]$UD
                vj <- x[[j]]$UD
                ai <- vol[[i]]$UD
                aj <- vol[[j]]$UD
                ai[ai<=lev] <- 1
                ai[ai>lev] <- 0
                aj[aj<=lev] <- 1
                aj[aj>lev] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[j,i] <- res[i,j] <- sqrt(sum((sqrt(vi) - sqrt(vj))^2*(attr(vi,"cellsize")^2)))
                } else {
                    res[j,i] <- res[i,j] <- sqrt(sum((sqrt(vi) - sqrt(vj))^2*(attr(vi,"cellsize")^2)))
                }
            }
        }
    }
    rownames(res) <- names(x)
    colnames(res) <- names(x)
    return(res)
}
