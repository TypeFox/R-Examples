FindOptimalChannel <-
function(x, px, y, rho.fn, R,
                               tol = 1e-5, eps = 0.001,
                               slope.range = c(-100, 0),
                               max.iters = Inf,
                               rho.scale = 1,
                               verbose = FALSE,
                               search.bracket = TRUE,
                               max.search.steps = 10, ...) {
    ## Search for an optimal capacity limited channel,
    ## i.e., a channel that minimizes distortion subject to the constraint
    ## that the information rate does not exceed a specified value R over
    ## the information source described by x.
    
    if(verbose) {
        cat(sprintf("\n% 10s", "s"), ": ")
        cat(sprintf("% 10s", "R"), "\t")
        cat(sprintf("% 10s", "D"), "\n")
        cat(paste(rep("=", 43), collapse = ""), "\n")
    }
    
    obj.fn <- function(s) {
        if(s > 0) {
            stop("FindOptimalChannel: Slope parameter must be negative.")
        }
        if(verbose) cat(sprintf("% 10.4g", s), ": ")
        channel <- BlahutAlgorithm(x, px, y, rho.fn, s, eps = eps,
                                   max.iters = max.iters,
                                   rho.scale = rho.scale, ...)
        if(verbose) cat(
            sprintf("% 10.4g", channel$R), "\t",
            sprintf("% 10.4g", channel$D), "\n")
        (channel$R - R)
    }

    sign.lower <- sign(obj.fn(slope.range[1]))
    sign.upper <- sign(obj.fn(slope.range[2]))

    s.opt <- NA
    search.count <- 0
    while(sign.lower == sign.upper) {
        search.count <- search.count + 1
        if(search.count > max.search.steps) {
            if(sign.lower < 0) {
                s.opt <- slope.range[1]
            } else {
                s.opt <- slope.range[2]
            }
            msg <- paste("Exceeded maximum number of iterations in bracket search.")
            warning(msg, call. = FALSE)
            break
        }
        if(search.bracket) {
            if(sign.lower < 0) {
                slope.range[1] <- -1 * 1.5 * abs(slope.range[1])
                sign.lower <- sign(obj.fn(slope.range[1]))
            } else {
                if(slope.range[2] < 0) {
                    slope.range[2] <- 0
                    sign.upper <- sign(obj.fn(slope.range[2]))
                } else {
                    s.opt <- slope.range[2]
                    msg <- paste("Could not bracket slope range for R = ", R, ".", sep = "")
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
                               max.iters = max.iters,
                               rho.scale = rho.scale, ...)
    return(channel)
}
