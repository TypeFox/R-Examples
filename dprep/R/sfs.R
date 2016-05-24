sfs <-
function (data, method = c("lda", "knn", "rpart"), kvec = 5, 
    repet = 10) 
{
#    require("MASS")
#    require("class")
#    require("rpart")
    if (sum(is.na(data))> 0) 
        stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    if (!(method %in% c("lda", "knn", "rpart"))) {
        cat("The classifier entered is not supported by this function.\n")
        return(method)
    }
    n = dim(data)[1]
    p = dim(data)[2]
    numbersel = rep(0, repet)
    fsel = rep(0, p - 1)
    for (i in 1:repet) {
        indic <- rep(0, p - 1)
        output <- indic
        varia <- p
        for (k in 1:(p - 1)) {
            correct <- rep(0, p - 1)
            if (k > 1) {
                varia <- c(where, varia)
            }
            for (m in 1:(p - 1)) {
                if (indic[m] == 0) {
                  which <- c(m, varia)
                  if (method == "lda") 
                    correct[m] <- cv10lda2(data[, which])
                  else if (method == "knn") 
                    {correct[m] <- cv10knn2(data[, which], kvec)
    }              else correct[m] <- cv10rpart2(data[, which])
                }
            }
            prov <- correct + runif(p - 1)
            where <- sum((1:(p - 1)) * as.numeric(max(prov) == 
                prov))
            output[k] <- correct[where]/n
            indic[where] <- 1
            if (k > 1) {
                if (output[k] <= output[k - 1]) {
                  indic <- rep(1, p - 1)
                }
            }
        }
        which <- rev(which)
        which <- which[-1]
        which1 <- which[1:(length(which) - 1)]
        numbersel[i] = length(which1)
        fsel[which1] = fsel[which1] + 1
    }
    bestsize = round(mean(numbersel))
    rev(order(fsel))
    bestsubset = rev(order(fsel))[1:bestsize]
    cat("The best subset of features is:")
    cat("\n")
    bestsubset
}
