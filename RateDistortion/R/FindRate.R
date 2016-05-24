FindRate <-
function(x, px, y, rho.fn, D,
                     tol = 1e-5, eps = 0.001,
                     max.iters = Inf,
                     slope.range = c(-100, 0),
                     search.bracket = TRUE,
                     max.search.steps = 10,
                     rho.scale = TRUE, verbose = FALSE, ...) {
    ## This function will find the channel that achieves the minimum
    ## information rate necessary to achieve the level of performance
    ## (distortion) given by D.

    if(verbose) {
        cat(sprintf("% 10s", "s"), ": ")
        cat(sprintf("% 10s", "R"), "\t")
        cat(sprintf("% 10s", "D"), "\n")
        cat(paste(rep("=", 43), collapse = ""), "\n")
    }
    
    obj.fn <- function(s) {
        if(verbose) cat(sprintf("% 10.4g", s), ": ")
        channel <- BlahutAlgorithm(x, px, y, rho.fn, s, eps = eps,
                                   max.iters = max.iters,
                                   rho.scale = rho.scale, ...)
        if(verbose) cat(
            sprintf("% 10.4g", channel$R), "\t",
            sprintf("% 10.4g", channel$D), "\n")
        (channel$D - D)
    }

    sign.lower <- sign(obj.fn(slope.range[1]))
    sign.upper <- sign(obj.fn(slope.range[2]))

    s.opt <- NA
    search.count <- 0
    while(sign.lower == sign.upper) {
        search.count <- search.count + 1
        if(search.count > max.search.steps) {
            if(sign.lower > 0) {
                s.opt <- slope.range[1]
            } else {
                s.opt <- slope.range[2]
            }
            msg <- paste("Exceeded maximum number of iterations in bracket search.")
            warning(msg, call. = FALSE)
            break
        }
        if(search.bracket) {
            if(sign.lower > 0) {
                slope.range[1] <- -1 * 1.5 * abs(slope.range[1])
                sign.lower <- sign(obj.fn(slope.range[1]))
            } else {
                if(slope.range[2] < 0) {
                    slope.range[2] <- 0
                    sign.upper <- sign(obj.fn(slope.range[2]))
                } else {
                    s.opt <- slope.range[2]
                    msg <- paste("Could not bracket slope range for D = ", D, ".", sep = "")
                    warning(msg, call. = FALSE)
                    break;
                }
            }
        } else {
            break
        }
    }

    if(is.na(s.opt)) {
        s.opt <- uniroot(obj.fn, slope.range, tol = tol)$root
    }
    if(verbose) cat("  s.opt = ", s.opt, "\n")
    channel <- BlahutAlgorithm(x, px, y, rho.fn, s.opt, eps = eps,
                               max.iters = max.iters, ...)
    return(channel)
}
