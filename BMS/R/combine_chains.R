combine_chains <-
function (...) 
{
    combine_topmods <- function(topmodobj1, topmodobj2) {
        nregs1 = topmodobj1$nregs
        nregs2 = topmodobj2$nregs
        if (nregs1 != nregs2) {
            stop("The number of regressors in both BMA chains has to be the same!")
        }
        k1 = length(topmodobj1$ncount())
        k2 = length(topmodobj2$ncount())
        nbmodels1 = topmodobj1$nbmodels
        nbmodels2 = topmodobj2$nbmodels
        ncount1 = topmodobj1$ncount()
        ncount2 = topmodobj2$ncount()
        lik1 = topmodobj1$lik()
        lik2 = topmodobj2$lik()
        bool1 = topmodobj1$bool()
        bool2 = topmodobj2$bool()
        betas1 = topmodobj1$betas()
        betas2 = topmodobj2$betas()
        betas2_1 = topmodobj1$betas2()
        betas2_2 = topmodobj2$betas2()
        fv1 = topmodobj1$fixed_vector()
        fv2 = topmodobj2$fixed_vector()
        if (all(betas1 == 0) | all(betas2 == 0)) {
            dobetas = FALSE
        }
        else {
            dobetas = TRUE
        }
        if (all(betas2_1 == 0) | all(betas2_2 == 0)) {
            dobetas2 = FALSE
        }
        else {
            dobetas2 = TRUE
        }
        idxin2_boolof1in2 = match(bool1, bool2)
        idxin1_boolof1in2 = which(!is.na(idxin2_boolof1in2))
        idxin2_boolof1in2 = idxin2_boolof1in2[!is.na(idxin2_boolof1in2)]
        ncount2[idxin2_boolof1in2] = ncount2[idxin2_boolof1in2] + 
            ncount1[idxin1_boolof1in2]
        if (any(idxin1_boolof1in2)) {
            ncount1 = ncount1[-idxin1_boolof1in2]
            lik1 = lik1[-idxin1_boolof1in2]
            bool1 = bool1[-idxin1_boolof1in2]
        }
        lika = c(lik2, lik1)
        orderlika = order(lika, decreasing = TRUE)
        lika = lika[orderlika]
        ncounta = c(ncount2, ncount1)[orderlika]
        boola = c(bool2, bool1)[orderlika]
        if (dobetas) {
            if (any(idxin1_boolof1in2)) 
                betas1 = betas1[, -idxin1_boolof1in2]
            betasa = cbind(betas2, betas1)[, orderlika]
            betasa_not0 = betasa != 0
            vecka = colSums(betasa_not0)
            vbetaa = as.vector(betasa[as.vector(betasa_not0)])
        }
        else {
            vecka = 0
            vbetaa = numeric(0)
        }
        if (dobetas2) {
            if (any(idxin1_boolof1in2)) 
                betas2_1 = betas2_1[, -idxin1_boolof1in2]
            betasa2 = cbind(betas2_2, betas2_1)[, orderlika]
            vbetaa2 = as.vector(betasa2[as.vector(betasa_not0)])
        }
        else {
            vbetaa2 = numeric(0)
        }
        fva = numeric(0)
        lfv = 0
        if ((nrow(fv1) == nrow(fv2)) & ((nrow(fv1) > 0) & (nrow(fv2) > 
            0))) {
            if (any(idxin1_boolof1in2)) 
                fv1 = fv1[, -idxin1_boolof1in2]
            if (!is.matrix(fv1)) 
                fv1 = t(fv1)
            fva = as.vector(cbind(fv2, fv1)[, orderlika])
            lfv = nrow(fv1)
        }
        return(.top10(nmaxregressors = nregs1, nbmodels = length(lika), 
            bbeta = dobetas, lengthfixedvec = lfv, bbeta2 = dobetas2, 
            inivec_lik = lika, inivec_bool = boola, inivec_count = ncounta, 
            inivec_vbeta = vbetaa, inivec_vbeta2 = vbetaa2, inivec_veck = vecka, 
            inivec_fixvec = fva))
    }
    combine_2chains <- function(flso1, flso2) {
        topmod.combi = combine_topmods(flso1$topmod, flso2$topmod)
        gpi <- flso1$gprior.info
        gpi$shrinkage.moments = numeric(length(gpi$shrinkage.moments))
        io1 = flso1$info
        io2 = flso2$info
        obj.combi = .post.calc(gprior.info = gpi, add.otherstats = io1$add.otherstats + 
            io2$add.otherstats, k.vec = (io1$k.vec[-1] + io2$k.vec[-1]), 
            null.count = (io1$k.vec[1] + io2$k.vec[1]), flso1$X.data, 
            topmods = topmod.combi, b1mo = io1$b1mo + io2$b1mo, 
            b2mo = io1$b2mo + io2$b2mo, iter = io1$iter + io2$iter, 
            burn = io1$burn + io2$burn, inccount = io1$inccount + 
                io2$inccount, models.visited = io1$models.visited + 
                io2$models.visited, K = io1$K, N = io1$N, msize = io1$msize + 
                io2$msize, timed = io1$timed + io2$timed, cumsumweights = io1$cumsumweights + 
                io2$cumsumweights, mcmc = flso1$arguments$mcmc, 
            possign = io1$pos.sign + io2$pos.sign)
        stpos1 = as.matrix(flso1$start.pos)
        stpos2 = as.matrix(flso2$start.pos)
        startpos.combi = cbind(rbind(stpos1, matrix(0, max(0, 
            nrow(stpos2) - nrow(stpos1)), ncol(stpos1))), rbind(stpos2, 
            matrix(0, max(0, nrow(stpos1) - nrow(stpos2)), ncol(stpos2))))
        call.combi = c(flso1$bms.call, flso2$bms.call)
        args.combi = flso1$arguments
        args2 = flso2$arguments
        args.combi$burn = args.combi$burn + args2$burn
        args.combi$iter = args.combi$iter + args2$iter
        if ((length(args.combi$mprior.size) == 1) | (length(args.combi$mprior.size) == 
            1)) {
            args.combi$mprior.size = mean(c(args.combi$mprior.size, 
                args2$mprior.size))
        }
        args.combi$nmodel = topmod.combi$nbmodels
        args.combi$user.int = (args.combi$user.int & args2$user.int)
        args.combi$g.stats = (args.combi$g.stats & args2$g.stats)
        mp1 = flso1$mprior.info
        mp2 = flso2$mprior.info
        if (mp1$mp.mode != mp2$mp.mode) {
            mpall = list()
        }
        else {
            mpall = mp1
            mpall$mp.msize = 0.5 * mp1$mp.msize + 0.5 * mp2$mp.msize
            mpall$origargs$mpparam = 0.5 * mp1$origargs$mpparam + 
                0.5 * mp2$origargs$mpparam
            mpall$mp.Kdist = 0.5 * mp1$mp.Kdist + 0.5 * mp2$mp.Kdist
        }
        result = list(info = obj.combi$info, arguments = args.combi, 
            topmod = topmod.combi, start.pos = startpos.combi, 
            gprior.info = obj.combi$gprior.info, mprior.info = mpall, 
            X.data = flso1$arguments$X.data, reg.names = obj.combi$reg.names, 
            bms.call = call.combi)
        class(result) = "bma"
        return(result)
    }
    arglist = list(...)
    if (!all(unlist(lapply(arglist, is.bma)))) 
        stop("All of the input arguments must be BMA objects!")
    if (!all(lapply(arglist, function(xx) xx$info$K) == arglist[[1]]$info$K)) 
        stop("All of the input BMA objects must have an equal number of max regressors (i.e. equal (X.data))!")
    if (!all(lapply(arglist, function(xx) xx$info$N) == arglist[[1]]$info$N)) 
        stop("All of the input BMA objects must have equal X.data!")
    if (!all(lapply(arglist, function(xx) xx$gprior.info$gtype) == 
        arglist[[1]]$gprior.info$gtype)) 
        stop("All of the input BMA objects must have the same type of g-prior (bms-argument g)")
    if (length(arglist) == 1) 
        return(arglist[[1]])
    combined_output <- combine_2chains(arglist[[1]], arglist[[2]])
    if (nargs() > 2) {
        for (inarg in 3:nargs()) {
            combined_output <- combine_2chains(arglist[[inarg]], 
                combined_output)
        }
    }
    return(combined_output)
}
