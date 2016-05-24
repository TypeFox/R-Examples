samp.postopt <- function (opt.cov.d, opt.mu.d, d.keep, prior, B = 400, B0 = 8000, 
    d = 10) 
{	
	lo <- apply(prior, 2, min)
	hi <- apply(prior, 2, max)
    H.k <- prior
    for (i in 1:d) {
        if (!is.na(opt.mu.d[i, 1])) {
            H.new <- mvrnorm(n = B, mu = opt.mu.d[i, ], Sigma = opt.cov.d[, 
                , i])
            for (j in 1:B) {
                lt0 <- table(H.new[j, ] <= 0)
                while (as.numeric(lt0[1]) != ncol(prior) | H.new[j, 1] < 
                  lo[1] | H.new[j, 1] > hi[1] | H.new[j, 2] < 
                  lo[2] | H.new[j, 2] > hi[2] | H.new[j, 3] < 
                  lo[3] | H.new[j, 3] > hi[3] | H.new[j, 4] < 
                  lo[4] | H.new[j, 4] > hi[4] | H.new[j, 5] < 
                  lo[5] | H.new[j, 5] > hi[5] | H.new[j, 6] < 
                  lo[6] | H.new[j, 6] > hi[6] | H.new[j, 7] < 
                  lo[7] | H.new[j, 7] > hi[7] | H.new[j, 8] < 
                  lo[8] | H.new[j, 8] > hi[8]) {
                  H.new[j, ] <- mvrnorm(1, mu = opt.mu.d[i, ], 
                    Sigma = opt.cov.d[, , i])
                  lt0 <- table(H.new[j, ] < 0)
                }
            }
            H.k <- rbind(H.k, H.new)
        }
    }
    #H.k <- rbind(prior, H.new0)
    H.k <- cbind(H.k[, 1:8])
    H.new <- H.k[(B0 + 1):nrow(H.k), ]
    B1 <- d.keep * B
    return(list(H.k = H.k, H.new = H.new, B1 = B1))
}



