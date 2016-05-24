read.plotSEMM_wACOV <- function(read){
    
    #apparently I need to do read in the results.dat file for now
    tmp <- read$input$savedata
    if(is.null(tmp))
        stop('ACOV matrix could not be read in')
    tmp <- gsub(' ', '', tmp)
    tmp <- gsub(';', '', tmp)
    tmp1 <- strsplit(tmp, '=')
    tmp2 <- strsplit(tmp, 'is')
    if(any(sapply(tmp1, length) == 2L)) tmp <- tmp1
    else tmp <- tmp2
    results_f <- NULL
    for(i in 1L:length(tmp))
        if(grepl('results', tmp[[i]][1L])) results_f <- tmp[[i]][2L]
    if(is.null(results_f))
        stop('Mplus results file not found. Please save using \'results = filename.dat \' ')
    dat <- t(read.table(results_f, sep=' ', header=FALSE, fill=TRUE))
    dat <- dat[!is.na(dat)]
    read$means <- matrix(dat)
    
    ovars <- strsplit(toupper(read$input$variable$usevariables), split = ' ')[[1L]]
    pi <- read$class_counts$modelEstimated$proportion
    pars <- read$parameters[[1L]]
    pars <- pars[!(pars$param %in% ovars), ] #latents only
    tmp <- min(which(grepl("*\\.ON$", pars$paramHeader)))
    ON <- pars$paramHeader[tmp]
    DV <- strsplit(ON, '.ON')[[1L]]
    IV <- pars$param[tmp]
    acov <- read$tech3$paramCov.savedata
    if(is.null(acov))
        stop('acov matrix must also be saved from Mplus')
    tech1 <- read$tech1$parameterSpecification
    nclass <- length(pi)
    
    #pars
    alphaarray <- rbind(pars[pars[,'paramHeader'] == 'Means' & pars[,'param'] == IV,'est'],
                        pars[pars[,'paramHeader'] == 'Intercepts' & pars[,'param'] == DV,'est'])
    psiarray <- rbind(pars[pars[,'paramHeader'] == 'Variances' & pars[,'param'] == IV,'est'],
                      pars[pars[,'paramHeader'] == 'Residual.Variances' & pars[,'param'] == DV,'est'])
    betavec <- pars[pars[,'paramHeader'] == ON,'est']
    ci_v <- c(pars[pars[,'paramHeader'] == 'Means' & pars[,'param'] != IV,'est'], 0)
    retpars <- list(alphaarray=alphaarray, psiarray=psiarray, betavec=betavec, ci_v=ci_v)
    
    #locations in acov
    c_loc <- tech1$LATENT.CLASS.REGRESSION.MODEL.PART$alpha.c
    if(length(c_loc) > 1L){
        c_loc <- c_loc[-length(c_loc)]
    } else c_loc <- 0
    alpha_loc <- beta_loc <- psi_loc <- c()
    for(i in 1L:nclass){
        tmp <- tech1[[i]]
        alpha_loc <- c(alpha_loc, tmp$alpha)
        beta_loc <- c(beta_loc, tmp$beta[2,1])
        psi_loc <- c(psi_loc, diag(tmp$psi))        
    }    
    loc <- list(c_loc=c_loc, alpha_loc=alpha_loc, beta_loc=beta_loc, psi_loc=psi_loc)   
    
    ret <- list(nclass=nclass, nparm=ncol(acov), acov=acov, loc=loc, pars=retpars, 
                means=read$means)
    ret
}
