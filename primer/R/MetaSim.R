`MetaSim` <-
function (Time = 50, NSims = 1, method = "hanski", ci = 0.25, 
    e = 0.25, phi = 0.75, p0 = 0.5, D = 0.5) 
{
    require(deSolve)
    out.stochastic <- matrix(NA, nrow = Time + 1, ncol = NSims)
    for (run in 1:NSims) out.stochastic[, run] <- {
        dpdt <- switch(pmatch(method, c("hanski", "levins", "gotelli", 
            "lande"), nomatch = ""), `1` = function(t, y, parms) {
            p <- y[1]
            with(as.list(parms), {
                dp <- cir * p * (1 - p) - er * p * (1 - p)
                return(list(dp))
            })
        }, `2` = function(t, y, parms) {
            p <- y[1]
            with(as.list(parms), {
                dp <- cir * p * (1 - p) - er * p
                return(list(dp))
            })
        }, `3` = function(t, y, parms) {
            p <- y[1]
            with(as.list(parms), {
                dp <- cir * (1 - p) - er * p
                return(list(dp))
            })
        }, `4` = function(t, y, parms) {
            p <- y[1]
            with(as.list(parms), {
                dp <- cir * p * (1 - D - p) - er * p
                return(list(dp))
            })
        })
        nout <- numeric(Time + 1)
        nout[1] <- p0
        for (j in 1:Time) nout[j + 1] <- {
            cir <- runif(1, ci * (1 - phi), ci * (1 + phi))
            er <- runif(1, e * (1 - phi), e * (1 + phi))
            params <- c(cir = cir, er = er, D = D)
            n <- ode(nout[j], 0:1, dpdt, params)[2, 2]
            if (n < 0) 
                n = 0
            if (n > 1) 
                n = 1
            n
        }
        nout
    }
    list(method = method, t = 0:(Time), Ns = as.matrix(out.stochastic), 
        Parameters = c(NSims = NSims, ci = ci, e = e, phi = phi, 
            p0 = p0, D = D))
}
