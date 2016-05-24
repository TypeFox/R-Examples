topmod <-
function (nbmodels, nmaxregressors = NA, bbeta = FALSE, lengthfixedvec = 0, 
    liks = numeric(0), ncounts = numeric(0), modelbinaries = matrix(0, 
        0, 0), betas = matrix(0, 0, 0), betas2 = matrix(0, 0, 
        0), fixed_vector = matrix(0, 0, 0)) 
{
    if (!is.numeric(nbmodels)) 
        stop("argument 'nbmodels' needs to be an integer>0")
    nbmodels = floor(nbmodels[[1]])
    if (nbmodels[[1]] < 0) 
        stop("argument 'nbmodels' needs to be an integer>0")
    if (bbeta > 0) 
        bbeta2 = TRUE
    else bbeta2 = FALSE
    bbeta = as.logical(bbeta)
    if (!bbeta & (length(betas) > 0)) 
        bbeta = TRUE
    if (!bbeta2 & (length(betas2) > 0)) 
        bbeta2 = TRUE
    if (is.na(lengthfixedvec[1])) 
        lengthfixedvec = 0
    if ((lengthfixedvec == 0) & length(fixed_vector) > 0) {
        lengthfixedvec = nrow(fixed_vector)
    }
    if (length(liks) > nbmodels) 
        stop("liks longer than nbmodels allows")
    if (length(ncounts) > nbmodels) 
        stop("ncounts longer than nbmodels allows")
    if ((length(modelbinaries) == 0) & (length(betas) > 0)) {
        modelbinaries = as.logical(betas)
        dim(modelbinaries) = dim(betas)
    }
    if (ncol(modelbinaries) > nbmodels) 
        stop("modelbinaries has more columns than than nbmodels allows")
    bindim = dim(modelbinaries)
    modelbinaries = as.logical(modelbinaries)
    dim(modelbinaries) = bindim
    if ((is.na(nmaxregressors[1])) & (length(modelbinaries) > 
        0)) {
        nmaxregressors = nrow(modelbinaries)
    }
    if (is.na(nmaxregressors)) 
        stop("argument 'nmaxregressors' is missing")
    nmaxregressors = floor(nmaxregressors[[1]])
    if (nmaxregressors <= 0) 
        stop("argument 'nmaxregressors' needs to be a positive integer")
    if ((length(ncounts) == 0) & (length(liks) > 0)) {
        ncounts = rep(1, length(liks))
    }
    if (length(modelbinaries) > 0) 
        if (nmaxregressors != nrow(modelbinaries)) 
            stop("nrow() of modelbinaries does not match nmaxregressors")
    if (bbeta & (length(betas) > 0)) 
        if (nmaxregressors != nrow(betas)) 
            stop("nrow() of betas does not match nmaxregressors")
    if (bbeta2 & (length(betas2) > 0)) 
        if (nmaxregressors != nrow(betas2)) 
            stop("nrow() of betas2 does not match nmaxregressors")
    N = length(liks)
    if (length(ncounts) != length(liks)) 
        stop("lengths of arguments 'liks' and 'ncounts' do not conform")
    if (ncol(modelbinaries) != length(liks)) 
        stop("nrow of argument 'modelbinaries' does not conform to length of argument 'liks'")
    if (bbeta) 
        if (ncol(betas) != length(liks)) 
            stop("nrow of argument 'betas' does not conform to length of argument 'liks'")
    if (bbeta2) 
        if (ncol(betas2) != length(liks)) 
            stop("nrow of argument 'betas2' does not conform to length of argument 'liks'")
    if (lengthfixedvec) 
        if (ncol(fixed_vector) != length(liks)) 
            stop("nrow of argument 'fixed_vector' does not conform to length of argument 'liks'")
    morder = order(liks, decreasing = TRUE)
    liks = liks[morder]
    modelbinaries = modelbinaries[, morder, drop = FALSE]
    ncounts = ncounts[morder]
    if (bbeta) {
        betas = betas[, morder, drop = FALSE]
    }
    if (bbeta2) {
        betas2 = betas2[, morder, drop = FALSE]
    }
    if (lengthfixedvec) {
        fixed_vector = fixed_vector[, morder, drop = FALSE]
    }
    hexobj = .hexcode.binvec.convert(nmaxregressors)
    bools = as.vector(sapply(as.list(as.data.frame(modelbinaries)), 
        hexobj$as.hexcode))
    if (length(bools) == 0) {
        bools = character(0)
    }
    veck = numeric(0)
    betas_raw = numeric(0)
    betas2_raw = numeric(0)
    if (bbeta & (length(bbeta) > 0)) {
        veck = colSums(modelbinaries)
        betas_raw = as.vector(betas)[as.vector(modelbinaries)]
    }
    if (bbeta2 & (length(bbeta2) > 0)) {
        veck = colSums(modelbinaries)
        betas2_raw = as.vector(betas2)[as.vector(modelbinaries)]
    }
    fixedvec = as.vector(fixed_vector)
    .top10(nmaxregressors = nmaxregressors, nbmodels = nbmodels, 
        bbeta = bbeta, lengthfixedvec = lengthfixedvec, bbeta2 = bbeta2, 
        inivec_lik = liks, inivec_bool = bools, inivec_count = ncounts, 
        inivec_vbeta = betas_raw, inivec_vbeta2 = betas2_raw, 
        inivec_veck = veck, inivec_fixvec = fixedvec)
}
