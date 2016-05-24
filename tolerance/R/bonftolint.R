bonftol.int <- function (fn, P1 = 0.005, P2 = 0.005, alpha = 0.05, ...) 
{
    P <- 1 - (P1 + P2)
    if (P < 0) {
        stop(paste("Invalid values for P1 and P2.", 
            "\n"))
    }
    Lower <- suppressMessages(fn(P = 1 - P1, side = 1, alpha = alpha, ...))
    Upper <- suppressMessages(fn(P = 1 - P2, side = 1, alpha = alpha, ...))
    message("These are ", (1 - alpha) * 100, "%/", P * 100, "% 2-sided tolerance limits controlling ", P1 * 100, "% in the lower tail and ", P2 * 100, "% in the upper tail.")
    if(class(Lower) != "list") {
        out <- cbind(Lower[, -which(colnames(Lower) == "1-sided.upper")],
            Upper[, which(colnames(Upper) == "1-sided.upper")])
        colnames(out)[(ncol(out) - 1):ncol(out)] <- c("2-sided.lower", "2-sided.upper")
        out[, 2] <- P
    }   else {
            for(i in 1:length(out)) {
                out[[i]] <- cbind(Lower[[i]][, -which(colnames(Lower[[i]]) == "1-sided.upper")],
                    Upper[[i]][, which(colnames(Upper[[i]]) == "1-sided.upper")])
                colnames(out)[(ncol(out[[i]]) - 1):ncol(out[[i]])] <- c("2-sided.lower", "2-sided.upper")
            }
        }
    out
}

