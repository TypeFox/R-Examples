MCMCCV <-
function (gen, D_G, D_E, nit, thinning, theta.max = c(10, 100 * 
    max(D_G), 100 * max(D_E), 1, 0.9), theta.init = c(1, max(D_G), 
    max(D_E), 0.7, 0.1), run = c(TRUE, FALSE, FALSE), ud = c(TRUE, 
    TRUE, TRUE, TRUE, TRUE), n.validation.set = 0, print.pct = TRUE) 
{
    max_D_G = max(D_G)
    max_D_E = max(D_E)
    D_G_scaled = D_G/max_D_G
    D_E_scaled = D_E/max_D_E
    theta.max[2] = theta.max[2]/max_D_G
    theta.max[3] = theta.max[3]/max_D_E
    theta.init[2] = theta.init[2]/max_D_G
    theta.init[3] = theta.init[3]/max_D_E
    run <- as.numeric(run)
    ud <- as.numeric(ud)
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
    if (n.validation.set > 0) {
        na = sample(nsite * nloc, n.validation.set)
        coo = matrix(ncol = 2, nrow = nsite * nloc, data = NA)
        coo[, 1] = 1:nsite
        for (i in 1:nloc) coo[((i - 1) * nsite + 1):(i * nsite), 
            2] = i
        sot = matrix(ncol = 2, data = coo[na, ])
        gen_tr = gen
        for (i in 1:n.validation.set) gen_tr[sot[i, 1], sot[i, 
            2], ] = NA
        gen_tr[array(is.na(gen_tr), length(gen_tr))] = 0
    }
    else if (n.validation.set == 0) 
        gen_tr = gen
    salva = ls(all.names = FALSE)
    j = 1
    for (i in 1:length(run)) {
        if (run[i] == 1) {
            swi <- i
            ud <- ud_save
            switch(swi, {
                cat("Computations for model 'geog+envi'\n")
            }, {
                cat("Computations for model 'geog'\n")
                ud[3] = 0
            }, {
                cat("Computations for model 'envi'\n\n")
                ud[2] = 0
            })
            res <- .Fortran("tangle", PACKAGE = "Sunder", as.integer(nsite), 
                as.integer(nloc), as.integer(nal), as.integer(nalM), 
                as.integer(gen_tr), as.double(D_G_scaled), as.double(D_E_scaled), 
                as.double(thma), as.double(theta.init), as.double(thetaprop), 
                as.integer(nit), as.integer(thinning), as.double(fcur), 
                as.double(fprop), as.double(xcur), as.double(xprop), 
                as.double(ycur), as.double(yprop), as.double(zcur), 
                as.double(zprop), as.double(Sigma), as.double(Sigma), 
                as.double(Uchol), as.double(theta.max[1]), as.double(theta.max[2]), 
                as.double(theta.max[3]), as.double(theta.max[4]), 
                as.double(theta.max[5]), as.double(zeta), as.double(fsav), 
                as.integer(ud), as.double(zsav), as.integer(swi), 
                as.integer(print.pct))
            tvt[, , j] <- matrix(nrow = 5, ncol = nit/thinning, 
                res[[8]])
            fvt[, , , , j] <- array(dim = c(nit/thinning, nsite, 
                nloc, nalM), res[[30]])
            rm(list = setdiff(ls(all.names = TRUE), c(salva, 
                "salva", "tvt", "fvt", "j")))
            j = j + 1
        }
    }
    mod.lik = array(data = NA, dim = length(run))
    if (n.validation.set > 0) {
        nna = sum(nal[sot[, 2]])
        f.me = array(dim = c(nna, sum(run)))
        g.na = rep(NA, nna)
        kon = 0
        for (i in 1:n.validation.set) {
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
    }
    tvt[2, , ] = tvt[2, , ] * max_D_G
    tvt[3, , ] = tvt[3, , ] * max_D_E
    theta.max[2] = theta.max[2] * max_D_G
    theta.max[3] = theta.max[3] * max_D_E
    theta.init[2] = theta.init[2] * max_D_G
    theta.init[3] = theta.init[3] * max_D_E
    dimnames(tvt) = list(theta_components = c("alpha", "beta_G", 
        "beta_E", "gamma", "delta"), time = paste(1:(nit/thinning)), 
        model = c("G+E", "G", "E")[as.logical(run)])
    names(mod.lik) = c("G+E", "G", "E")
    output <- list(theta = tvt, frequencies = fvt, mod.lik = mod.lik, 
        theta.max = theta.max, theta.init = theta.init, run = as.logical(run), 
        ud = as.logical(ud), nit = nit, thinning = thinning, 
        n.validation.set = n.validation.set)
    return(output)
}
