hfdenoise.wav <-
function (x, binsize, transform = "binhf", meth = "u", van = 1, 
    fam = "DaubExPhase", min.level = 3, coarse = FALSE) 
{
    if (transform == "binhf") {
        hfx1 <- binhf.wd(x, binsize)
        hfx <- hfx1$transformed
    }
    if (transform == "ansc") {
        hfx <- ansc(x, binsize)
    }
    tmp.wd <- wd(hfx, filter.number = van, family = fam)
    if (meth == "u") {
        thr.wd <- threshold.wd(tmp.wd, levels = (min.level):(tmp.wd$nlevels - 
            1), policy = "universal", type = "hard")
    }
    if (meth == "c") {
        thr.wd <- threshold.wd(tmp.wd, levels = (min.level):(tmp.wd$nlevels - 
            1), policy = "cv", type = "hard")
    }
    if (coarse == TRUE) {
    	# Manual "nullevels" for DNA denoise and new version of wavethresh:

       levelstonull = (thr.wd$nlevels - 11):(thr.wd$nlevels - 1)
	
	for (lev in levelstonull) {
            d <- accessD(wd, level = lev)
            d <- rep(0, length(d))
            thr.wd <- putD(wd, level = lev, v = d)
        }
    }
    
    tmp <- wr(thr.wd)
    if (transform == "binhf") {
        tmp <- list(transformed = tmp, cnew = hfx1$cnew)
        fhat <- invbinhf.wd(tmp, binsize)
    }
    if (transform == "ansc") {
        fhat <- invansc(tmp, binsize)
    }
    fhat
}

