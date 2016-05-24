fun.theo.mv.gld <-
function (L1, L2, L3, L4, param, normalise = "N") 
{
    if (length(L1) > 1) {
        L4 <- L1[4]
        L3 <- L1[3]
        L2 <- L1[2]
        L1 <- L1[1]
    }
    result <- rep(NA, 4)
    if (tolower(param) == "rs") {
        v1 <- fun.rsb(L3, L4, 1)
        v2 <- fun.rsb(L3, L4, 2)
        v3 <- fun.rsb(L3, L4, 3)
        v4 <- fun.rsb(L3, L4, 4)
        result[1] <- L1 + 1/L2 * ((L3 + 1)^-1 - (L4 + 1)^-1)
        result[2] <- (v2 - v1^2)/L2^2
        result[3] <- (v3 - 3 * v1 * v2 + 2 * (v1)^3)/((L2^3) * 
            result[2]^1.5)
        result[4] <- (v4 - 4 * v1 * v3 + 6 * (v1^2) * v2 - 3 * 
            v1^4)/((L2^4) * result[2]^2)
        result[!as.logical((L3 > -1/1:4) * (L4 > -1/1:4))] <- NA
    }
    if (tolower(param) == "fmkl" | tolower(param) == "fkml") {
        if (L3 == 0 & L4 == 0) {
            v1 <- 0
            v2 <- fun.fmkl0(2)
            v3 <- 0
            v4 <- fun.fmkl0(4)
            result <- rep(NA, 4)
            result[1] <- v1/L2 + L1
            result[2] <- (v2 - v1^2)/L2^2
            result[3] <- (v3 - 3 * v1 * v2 + 2 * (v1)^3)/((L2^3) * 
                result[2]^1.5)
            result[4] <- (v4 - 4 * v1 * v3 + 6 * (v1^2) * v2 - 
                3 * v1^4)/((L2^4) * result[2]^2)
            if (normalise == "Y") {
                result[4] <- result[4] - 3
            }
            names(result) <- c("mean", "variance", "skewness", 
                "kurtosis")
            return(result)
        }
        if (L3 == 0 & L4 != 0) {
            v1 <- fun.fmkl.L30(1, L4)
            v2 <- fun.fmkl.L30(2, L4)
            v3 <- fun.fmkl.L30(3, L4)
            v4 <- fun.fmkl.L30(4, L4)
            result <- rep(NA, 4)
            result[1] <- v1/L2 + L1
            result[2] <- (v2 - v1^2)/L2^2
            result[3] <- (v3 - 3 * v1 * v2 + 2 * (v1)^3)/((L2^3) * 
                result[2]^1.5)
            result[4] <- (v4 - 4 * v1 * v3 + 6 * (v1^2) * v2 - 
                3 * v1^4)/((L2^4) * result[2]^2)
            if (normalise == "Y") {
                result[4] <- result[4] - 3
            }
            names(result) <- c("mean", "variance", "skewness", 
                "kurtosis")
            return(result)
        }
        if (L3 != 0 & L4 == 0) {
            v1 <- fun.fmkl.L40(1, L3)
            v2 <- fun.fmkl.L40(2, L3)
            v3 <- fun.fmkl.L40(3, L3)
            v4 <- fun.fmkl.L40(4, L3)
            result <- rep(NA, 4)
            result[1] <- v1/L2 + L1
            result[2] <- (v2 - v1^2)/L2^2
            result[3] <- (v3 - 3 * v1 * v2 + 2 * (v1)^3)/((L2^3) * 
                result[2]^1.5)
            result[4] <- (v4 - 4 * v1 * v3 + 6 * (v1^2) * v2 - 
                3 * v1^4)/((L2^4) * result[2]^2)
            if (normalise == "Y") {
                result[4] <- result[4] - 3
            }
            names(result) <- c("mean", "variance", "skewness", 
                "kurtosis")
            return(result)
        }
        v1 <- fun.fmklb(L3, L4, 1)
        v2 <- fun.fmklb(L3, L4, 2)
        v3 <- fun.fmklb(L3, L4, 3)
        v4 <- fun.fmklb(L3, L4, 4)
        result[1] <- L1 - 1/L2 * ((L3 + 1)^-1 - (L4 + 1)^-1)
        result[2] <- (v2 - v1^2)/L2^2
        result[3] <- (v3 - 3 * v1 * v2 + 2 * (v1)^3)/((L2^3) * 
            result[2]^1.5)
        result[4] <- (v4 - 4 * v1 * v3 + 6 * (v1^2) * v2 - 3 * 
            v1^4)/((L2^4) * result[2]^2)
        result[as.logical(min(L3, L4) < (-1/1:4))] <- NA
    }
    if (normalise == "Y") {
        result[4] <- result[4] - 3
    }
    names(result) <- c("mean", "variance", "skewness", "kurtosis")
    return(result)
}
