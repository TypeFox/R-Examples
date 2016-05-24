TangleInference <-
function (gen, D_G, D_E, nit, thinning, theta.max, theta.init, 
    run, ud, set) 
{
    ud_save = ud
    nsite = dim(gen)[1]
    nloc = dim(gen)[2]
    nalM = dim(gen)[3]
    mdim = c(nsite, nloc, nalM)
    zcur = array(dim = mdim, rnorm(n = prod(mdim)))
    nal <- rep(2, nloc)
    for (i in 1:nloc) nal[i] = max(which(apply(gen[, i, ], 2, 
        function(x) !all(x %in% 0)) == T))
    tvt = array(data = NA, dim = c(5, nit/thinning, sum(run)))
    fvt = array(data = NA, dim = c(nit/thinning, nsite, nloc, 
        nalM, sum(run)))
    fcur = fprop = xcur = xprop = ycur = yprop = zprop = array(dim = mdim, 
        data = -999)
    zsav = fsav = array(dim = c(nit/thinning, nsite, nloc, nalM), 
        data = -999)
    thma = zeta = matrix(nrow = 5, ncol = nit/thinning, data = -999)
    thetaprop = rep(-999, 5)
    Sigma = Uchol = diag(nsite)
    if (set > 0) {
        na = sample(nsite * nloc, set)
        coo = matrix(ncol = 2, nrow = nsite * nloc, data = NA)
        coo[, 1] = 1:nsite
        for (i in 1:nloc) coo[((i - 1) * nsite + 1):(i * nsite), 
            2] = i
        sot = matrix(ncol = 2, data = coo[na, ])
        gen_tr = gen
        for (i in 1:set) gen_tr[sot[i, 1], sot[i, 2], ] = NA
        gen_tr[array(is.na(gen_tr), length(gen_tr))] = 0
    }
    else if (set == 0) 
        gen_tr = gen
    salva = ls(all.names = FALSE)
    j = 1
    for (i in 1:length(run)) {
        if (run[i] == 1) {
            swi <- i
            ud <- ud_save
            switch(swi, {
                cat("Will check now model 'geog+envi'\n")
            }, {
                cat("Will check now model 'geog'\n")
                ud[3] = 0
            }, {
                cat("Will check now model 'envi'\n\n")
                ud[2] = 0
            })
            res <- .Fortran("tangle", as.integer(nsite), as.integer(nloc), 
                as.integer(nal), as.integer(nalM), as.integer(gen_tr), 
                as.double(D_G), as.double(D_E), as.double(thma), 
                as.double(theta.init), as.double(thetaprop), 
                as.integer(nit), as.integer(thinning), as.double(fcur), 
                as.double(fprop), as.double(xcur), as.double(xprop), 
                as.double(ycur), as.double(yprop), as.double(zcur), 
                as.double(zprop), as.double(Sigma), as.double(Sigma), 
                as.double(Uchol), as.double(theta.max[1]), as.double(theta.max[2]), 
                as.double(theta.max[3]), as.double(theta.max[4]), 
                as.double(theta.max[5]), as.double(zeta), as.double(fsav), 
                as.integer(ud), as.double(zsav), as.integer(swi))
            tvt[, , j] <- matrix(nrow = 5, ncol = nit/thinning, 
                res[[8]])
            fvt[, , , , j] <- array(dim = c(nit/thinning, nsite, 
                nloc, nalM), res[[30]])
            rm(list = setdiff(ls(all.names = TRUE), c(salva, 
                "salva", "tvt", "fvt", "j")))
            j = j + 1
        }
    }
    nna = sum(nal[sot[, 2]])
    f.me = array(dim = c(nna, sum(run)))
    g.na = rep(NA, nna)
    kon = 0
    for (i in 1:set) {
        for (ial in 1:nal[sot[i, 2]]) {
            kon = kon + 1
            g.na[kon] = gen[sot[i, 1], sot[i, 2], ial]
            k = 1
            for (j in 1:length(run)) {
                if (run[j] == 1) {
                  f.me[kon, k] = mean(fvt[, sot[i, 1], sot[i, 
                    2], ial, k])
                  k = k + 1
                }
            }
        }
    }
    mod.lik = array(data = -Inf, dim = length(run))
    k = 1
    for (j in 1:length(run)) {
        if (run[j] == 1) {
            mod.lik[j] = 0
            for (i in 1:nna) {
                mod.lik[j] = mod.lik[j] + g.na[i] * log(f.me[i, 
                  k])
            }
            k = k + 1
        }
    }
    names(mod.lik) = c("G+E", "G", "E")
    output <- list(theta = tvt, frequencies = fvt, mod.lik = mod.lik)
    return(output)
}
