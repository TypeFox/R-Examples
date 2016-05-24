setMethod("profile",
    signature(fitted = "mlePP"),
    function (fitted, which = 1L:p, maxsteps = 100, 
        alpha = 0.01, zmax = sqrt(qchisq(1 - alpha, 1L)), del = zmax/5, 
        trace = FALSE, ...) 
    {

        onestep <- function(step) {
            bi <- B0[i] + sgn * step * del * std.err[i]
            fix <- list(bi)
            names(fix) <- pi
            call$fixed <- c(fix, fix0)
            pfit <- tryCatch(eval.parent(call, 2L), error = identity)
            if (inherits(pfit, "error")) 
                return(NA)
            else {
                zz <- 2 * (pfit@min - fitted@min)
                ri <- pv0
                ri[, names(pfit@coef)] <- pfit@coef
                ri[, pi] <- bi
                if (zz > -0.001) 
                  zz <- max(zz, 0)
                else stop("profiling has found a better solution, so original fit had not converged")
                z <- sgn * sqrt(zz)
                pvi <<- rbind(pvi, ri)
                zi <<- c(zi, z)
            }
            if (trace) 
                cat(bi, z, "\n")
            z
        }
        summ <- summary(fitted)
        std.err <- summ@coef[, "Std. Error"]
        Pnames <- names(B0 <- fitted@coef)
        pv0 <- t(as.matrix(B0))
        p <- length(Pnames)
        prof <- vector("list", length = length(which))
        names(prof) <- Pnames[which]
        call <- fitted@call
        call$modSim <- TRUE
        call$modCI <- FALSE
        call$dplot <- FALSE
        ndeps <- eval.parent(call$control$ndeps)
        parscale <- eval.parent(call$control$parscale)
        fix0 <- eval.parent(call$fixed)
        for (i in which) {
            zi <- 0
            pvi <- pv0
            pi <- Pnames[i]
            if (!is.null(ndeps)) 
                call$control$ndeps <- ndeps[-i]
            if (!is.null(parscale)) 
                call$control$parscale <- parscale[-i]
            for (sgn in c(-1, 1)) {
                if (trace) 
                  cat("\nParameter:", pi, c("down", "up")[(sgn + 
                    1)/2 + 1], "\n")
                step <- 0
                z <- 0
                lastz <- 0
                while ((step <- step + 1) < maxsteps && abs(z) < 
                  zmax) {
                  z <- onestep(step)
                  if (is.na(z)) 
                    break
                  lastz <- z
                }
                if (abs(lastz) < zmax) {
                  for (dstep in c(0.2, 0.4, 0.6, 0.8, 0.9)) {
                    z <- onestep(step - 1 + dstep)
                    if (is.na(z) || abs(z) > zmax) 
                      break
                  }
                }
                else if (length(zi) < 5L) {
                  mxstep <- step - 1L
                  step <- 0.5
                  while ((step <- step + 1L) < mxstep) onestep(step)
                }
            }
            si <- order(pvi[, i])
            prof[[pi]] <- data.frame(z = zi[si])
            prof[[pi]]$par.vals <- pvi[si, , drop = FALSE]
        }
        new("profile.mle", profile = prof, summary = summ)
    }
)
