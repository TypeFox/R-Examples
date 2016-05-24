
fastesa <- function (g0, sigma, noccasions, traps, mask) {
    dk2 <- outer (traps$x, mask$x, '-')^2 + outer (traps$y, mask$y, '-')^2
    pk <- g0 * exp(-dk2/2/sigma^2)
    temp2 <- 1-apply (1-pk, 2, prod, na.rm = TRUE) ^ noccasions
    attr(mask, 'area') * sum(temp2, na.rm = TRUE)
}

## Expected number of captures (Efford Dawson and Borchers 2009 Eq. B.2)
Ecap <- function (g0, sigma, noccasions, traps, mask) {
    dk2 <- outer (traps$x, mask$x, '-')^2 + outer (traps$y, mask$y, '-')^2
    pkm <- g0 * exp(-dk2/2/sigma^2)
    n.m <- noccasions * apply (pkm, 2, sum, na.rm=T) # over traps
    p.m <- 1-apply (1-pkm, 2, prod, na.rm=T)^noccasions # over traps
    sum(n.m) / sum(p.m)
}

scenariosFromStatistics <- function (sigma, noccasions, traps, mask, nval, rval, g0.int = c(0.001, 0.999)) {
    caproot <- function (g0, cn) {
        Ecap(g0, sigma, noccasions, traps, mask) - cn }
    fn <- function (nr) {
        n <- nr[1]; r <- nr[2]; noccasions = nr[3]; sigma = nr[4]
        g0 <- uniroot (caproot, interval = g0.int, cn = (n+r)/n)$root
        a <- fastesa (g0, sigma, noccasions, traps, mask)
        D <- n/a
        c (n, r, noccasions, sigma, D, g0, a)
    }
    output <- as.data.frame(t(apply (expand.grid(nval, rval, noccasions, sigma), 1, fn)))
    names(output) <- c('n', 'r', 'noccasions', 'sigma', 'D', 'g0', 'a')
    nrw <- nrow(output)
    extra <- data.frame(trapsindex = 1, nrepeats = 1, detectfn = 0, recapfactor = 1,
                   popindex = 1, detindex = 1, fitindex = 1, group = 1)
    out <- cbind(output, extra[rep(1,nrw),])
    out$scenario <- 1:nrw
    scenariofields <- c('scenario','trapsindex', 'noccasions', 'nrepeats', 'D', 'g0', 'sigma', 'detectfn', 'recapfactor',
                           'popindex', 'detindex', 'fitindex')
    out <- out[,scenariofields]
    rownames(out) <- paste('n=', output$n, ',r=', output$r, ',nocc=', output$noccasions, ',sigma=', output$sigma, sep='')
    out
}

# library(secr)
# grid36 <- make.grid(nx = 6, ny = 6, spacing = 200)
# mask <- make.mask(grid36, buffer = 2000)
# tmp <- scenariosFromStatistics (sigma = 300, noccasions = 44, traps = grid36, mask = mask, nval = 14, rval = 34)
# run.scenarios(tmp, nrepl=5, traps=grid36, mask=mask)


