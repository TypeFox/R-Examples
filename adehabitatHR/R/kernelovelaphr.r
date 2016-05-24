kerneloverlaphr <- function (x, method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"),
                           percent = 95, conditional = FALSE, ...)
{
    method <- match.arg(method)
    if (!inherits(x, "estUDm"))
        stop("x should be of class estUDm")
    if (length(x)==1)
        stop("several animals are needed for this function")
    if (slot(x[[1]],"vol"))
        stop("x should not be a volume under UD")
    vol <- getvolumeUD(x)
    x <- lapply(x, function(y) {
        coo <- coordinates(y)
        y[order(coo[,1], coo[,2]),]
    })
    vol <- lapply(vol, function(y) {
        coo <- coordinates(y)
        y[order(coo[,1], coo[,2]),]
    })

    gp <- gridparameters(vol[[1]])
    res <- matrix(0, ncol = length(x), nrow = length(x))
    for (i in 1:length(x)) {
        for (j in 1:i) {
            if (method == "HR") {
                vi <- vol[[i]][[1]]
                vj <- vol[[j]][[1]]
                vi[vi <= percent] <- 1
                vi[vi > percent] <- 0
                vj[vj <= percent] <- 1
                vj[vj > percent] <- 0
                vk <- vi * vj
                res[i, j] <- sum(vk)/sum(vi)
                res[j, i] <- sum(vk)/sum(vj)
            }
            if (method == "PHR") {
                vi <- x[[i]][[1]]
                vj <- x[[j]][[1]]
                ai <- vol[[i]][[1]]
                aj <- vol[[j]][[1]]
                ai[ai <= percent] <- 1
                ai[ai > percent] <- 0
                aj[aj <= percent] <- 1
                aj[aj > percent] <- 0
                if (conditional) {
                  vi <- vi * ai
                  vj <- vj * aj
                  res[j, i] <- sum(vi * aj) * (gp[1, 2]^2)
                  res[i, j] <- sum(vj * ai) * (gp[1, 2]^2)
                }
                else {
                  res[j, i] <- sum(vi * aj) * (gp[1, 2]^2)
                  res[i, j] <- sum(vj * ai) * (gp[1, 2]^2)
                }
            }
            if (method == "VI") {
                vi <- x[[i]][[1]]
                vj <- x[[j]][[1]]
                ai <- vol[[i]][[1]]
                aj <- vol[[j]][[1]]
                ai[ai <= percent] <- 1
                ai[ai > percent] <- 0
                aj[aj <= percent] <- 1
                aj[aj > percent] <- 0
                if (conditional) {
                  vi <- vi * ai
                  vj <- vj * aj
                  res[i, j] <- res[j, i] <- sum(pmin(vi, vj)) *
                    (gp[1, 2]^2)
                }
                else {
                  res[i, j] <- res[j, i] <- sum(pmin(vi, vj)) *
                    (gp[1, 2]^2)
                }
            }
            if (method == "BA") {
                vi <- x[[i]][[1]]
                vj <- x[[j]][[1]]
                ai <- vol[[i]][[1]]
                aj <- vol[[j]][[1]]
                ai[ai <= percent] <- 1
                ai[ai > percent] <- 0
                aj[aj <= percent] <- 1
                aj[aj > percent] <- 0
                if (conditional) {
                  vi <- vi * ai
                  vj <- vj * aj
                  res[j, i] <- res[i, j] <- sum(sqrt(vi) * sqrt(vj)) *
                    (gp[1, 2]^2)
                }
                else {
                  res[j, i] <- res[i, j] <- sum(sqrt(vi) * sqrt(vj)) *
                    (gp[1, 2]^2)
                }
            }
            if (method == "UDOI") {
                vi <- x[[i]][[1]]
                vj <- x[[j]][[1]]
                ai <- vol[[i]][[1]]
                aj <- vol[[j]][[1]]
                ai[ai <= percent] <- 1
                ai[ai > percent] <- 0
                aj[aj <= percent] <- 1
                aj[aj > percent] <- 0
                if (conditional) {
                  vi <- vi * ai
                  vj <- vj * aj
                  ak <- sum(ai * aj) * (gp[1, 2]^2)
                  res[j, i] <- res[i, j] <- ak * sum(vi * vj) *
                    (gp[1, 2]^2)
                }
                else {
                  ak <- sum(ai * aj) * (gp[1, 2]^2)
                  res[j, i] <- res[i, j] <- ak * sum(vi * vj) *
                    (gp[1, 2]^2)
                }
            }
            if (method == "HD") {
                vi <- x[[i]][[1]]
                vj <- x[[j]][[1]]
                ai <- vol[[i]][[1]]
                aj <- vol[[j]][[1]]
                ai[ai <= percent] <- 1
                ai[ai > percent] <- 0
                aj[aj <= percent] <- 1
                aj[aj > percent] <- 0
                if (conditional) {
                  vi <- vi * ai
                  vj <- vj * aj
                  res[j, i] <- res[i, j] <- sqrt(sum((sqrt(vi) -
                    sqrt(vj))^2 * (gp[1, 2]^2)))
                }
                else {
                  res[j, i] <- res[i, j] <- sqrt(sum((sqrt(vi) -
                    sqrt(vj))^2 * (gp[1, 2]^2)))
                }
            }
        }
    }
    rownames(res) <- names(x)
    colnames(res) <- names(x)
    return(res)
}
