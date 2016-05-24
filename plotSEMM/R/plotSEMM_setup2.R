plotSEMM_setup2 <- function(setup, alpha = .025, boot = NULL, boot.CE=FALSE, boot.CI=TRUE,
                            points = 50, fixed_value=NA){
    
    #only supports 2 or more classes 
    #requires Mplus read in file as input; setup <- read.plotSEMM_wACOV(read)
    if(!is.na(fixed_value)) points <- points + 1

    nclass <- classes <- setup$nclass
    nparam <- setup$nparm; acov <- setup$acov; loc <- setup$loc
    c_loc <- loc$c_loc; alpha_loc <- loc$alpha_loc; beta_loc <- loc$beta_loc; 
    psi_loc <- loc$psi_loc
    locations <- c(c_loc,alpha_loc,beta_loc,psi_loc[seq(from=1, to=length(psi_loc), by=2)])
    p <- nclass*(4)+(nclass-1)  #number of variables to be used later in generating curves
    p1 <- length(unique(locations))
    pars <- setup$pars
    alphaarray <- pars$alphaarray; psiarray <- pars$psiarray; 
    gamma <- betavec <- pars$betavec; ci_v <- pars$ci_v
    means <- setup$means
    
    if(any(psiarray < 0))
        stop('Negative variances supplied. Please fix.')
        
    sum_expi <- sum(exp(ci_v))
    pi_v <- exp(ci_v) / sum_expi
    
    # ----------------------------------------------------------------------------------------------------------------------
    # code for computing I(z) taken from older plotting code
    # ----------------------------------------------------------------------------------------------------------------------
    
    # these need to be in 2x2 matrices/ 2x1 vectors for computations
    
    alphaarray2 <- array(data = NA, c(2, 1, classes))
    betaarray <- array(data = NA, c(2, 2, classes))
    psiarray2 <- array(data = NA, c(2, 2, classes))
    
    for (i in 1:classes) {
        alphaarray2[, , i] <- matrix(c(alphaarray[1, i], alphaarray[2, i]), 2, 1, byrow = TRUE)
        betaarray[, , i] <- matrix(c(0, 0, betavec[i], 0), 2, 2, byrow = TRUE)
        psiarray2[, , i] <- matrix(c(psiarray[1, i], 0, 0, psiarray[2, i]), 2, 2, byrow = TRUE)
    }
    
    IMPCOV <- array(data = NA, c(2, 2, classes))
    IMPMEAN <- array(data = NA, c(2, 2, classes))
    
    for (i in 1:classes) {
        IMPCOV[, , i] <- solve(diag(x = 1, nrow = 2, ncol = 2) - betaarray[, , i]) %*% (psiarray2[, , i]) %*% t(solve(diag(x = 1, 
            nrow = 2, ncol = 2) - betaarray[, , i]))
        IMPMEAN[, , i] <- solve(diag(x = 1, nrow = 2, ncol = 2) - betaarray[, , i]) %*% (alphaarray2[, , i])
    }
    
    
    MuEta_1 <- vector(mode = "numeric", length = classes)
    MuEta_2 <- vector(mode = "numeric", length = classes)
    VEta_1 <- vector(mode = "numeric", length = classes)
    VEta_2 <- vector(mode = "numeric", length = classes)
    COVKSIETA <- vector(mode = "numeric", length = classes)
    
    for (i in 1:classes) {
        MuEta_1[i] = IMPMEAN[1, 1, i]
        MuEta_2[i] = IMPMEAN[2, 2, i]
        VEta_1[i] = IMPCOV[1, 1, i]
        VEta_2[i] = IMPCOV[2, 2, i]
        COVKSIETA[i] = IMPCOV[1, 2, i]
    }
    # upper and lower bounds for plots
    
    muEta1 <- 0
    muEta2 <- 0
    
    for (i in 1:classes) {
        
        muEta1 <- muEta1 + pi_v[i] * MuEta_1[i]
        muEta2 <- muEta2 + pi_v[i] * MuEta_2[i]
    }
    
    vEta1 <- 0
    vEta2 <- 0
    for (i in 1:classes) {
        for (j in 1:classes) {
            if (i < j) {
                vEta1 <- vEta1 + pi_v[i] * pi_v[j] * (MuEta_1[i] - MuEta_1[j])^2
                vEta2 <- vEta2 + pi_v[i] * pi_v[j] * (MuEta_2[i] - MuEta_2[j])^2
            }
        }
    }
    for (i in 1:classes) {
        
        vEta1 <- vEta1 + pi_v[i] * VEta_1[i]
        vEta2 <- vEta2 + pi_v[i] * VEta_2[i]
    }
    
    LEta1 = muEta1 - 3 * sqrt(vEta1)
    UEta1 = muEta1 + 3 * sqrt(vEta1)
    LEta2 = muEta2 - 3 * sqrt(vEta2)
    UEta2 = muEta2 + 3 * sqrt(vEta2)
    
    LB = min(LEta1, LEta2)
    UB = max(UEta1, UEta2)
    
    
    # computations for contour plot
    
    Eta1 <- seq(LEta1, UEta1, length = points - !is.na(fixed_value))
    Eta2 <- seq(LEta2, UEta2, length = points - !is.na(fixed_value))
    if(!is.na(fixed_value)){
        Eta1 <- c(Eta1, fixed_value)
        Eta2 <- c(Eta2, fixed_value)
    }
    
    
    r <- vector(mode = "numeric", length = classes)
    
    for (i in 1:classes) {
        r[i] <- COVKSIETA[i]/sqrt(VEta_1[i] * VEta_2[i])
    }
    
    denKE <- function(Eta1, Eta2) {
        
        placeholder <- 0
        denKE_ <- matrix(data = 0, nrow = length(Eta1), ncol = classes)
        
        for (i in 1:classes) {
            z <- ((Eta1 - MuEta_1[i])^2)/VEta_1[i] + ((Eta2 - MuEta_2[i])^2)/VEta_2[i] - 2 * r[i] * (Eta1 - MuEta_1[i]) * (Eta2 - 
                MuEta_2[i])/sqrt(VEta_1[i] * VEta_2[i])
            denKE_[, i] <- (1/(2 * 22/7 * sqrt(VEta_1[i]) * sqrt(VEta_2[i]) * sqrt(1 - r[i]^2))) * exp(-z/(2 * (1 - r[i]^2)))
        }
        
        for (i in 1:classes) {
            placeholder <- placeholder + pi_v[i] * denKE_[, i]
        }
        
        denKE <- placeholder
    }
    
    z <- outer(Eta1, Eta2, denKE)
    
    # ---------------------------------------------------------------------------------------------------------------------- End
    # of code for computing I(z)
    # ----------------------------------------------------------------------------------------------------------------------
    
    x <- Eta1
    x2 <- Eta2
    
    phi <- array(data = 0, c(points, classes))
    for (i in 1:classes) {
        phi[, i] <- dnorm(x, mean = alphaarray[1, i], sd = sqrt(psiarray[1, i]))
    }
    
    
    a_pi <- array(data = 0, c(points, classes))
    a_pi2 <- array(data = 0, c(points, classes))
    for (i in 1:classes) {
        a_pi[, i] <- pi_v[i] * dnorm(x, mean = MuEta_1[i], sd = sqrt(VEta_1[i]))
        a_pi2[, i] <- pi_v[i] * dnorm(x2, mean = MuEta_2[i], sd = sqrt(VEta_2[i]))
    }
    
    sumpi <- array(data = 0, c(points, 1))
    sumpi2 <- array(data = 0, c(points, 1))
    for (i in 1:classes) {
        sumpi[, 1] <- sumpi[, 1] + a_pi[, i]
        sumpi2[, 1] <- sumpi2[, 1] + a_pi2[, i]
    }
    
    pi <- array(data = 0, c(points, classes))
    for (i in 1:classes) {
        pi[, i] <- a_pi[, i]/sumpi[, 1]
    }
    y <- 0
    for (i in 1:classes) {
        y <- y + pi[, i] * (alphaarray[2, i] + gamma[i] * x)
    }
    
    # Derivatives for delta method CIs
    D <- 0
    for (i in 1:classes) {
        D = D + exp(ci_v[i]) * phi[, i]
    }
    
    
    dalpha <- array(data = 0, c(points, classes))
    dphi <- array(data = 0, c(points, classes))
    dc <- array(data = 0, c(points, classes - 1))
    for (i in 1:classes) {
        for (j in 1:classes) {
            if (i != j) {
                dalpha[, i] <- dalpha[, i] + ((exp(ci_v[j]) * phi[, j] * ((alphaarray[2, i] - alphaarray[2, j]) + (gamma[i] - 
                  gamma[j]) * x)))
                dphi[, i] <- dphi[, i] + ((exp(ci_v[j]) * phi[, j] * ((alphaarray[2, i] - alphaarray[2, j]) + (gamma[i] - gamma[j]) * 
                  x)))
                
            }
        }
        dalpha[, i] <- (dalpha[, i] * exp(ci_v[i]) * phi[, i] * ((x - alphaarray[1, i])/psiarray[1, i])) * (1/D^2)
        dphi[, i] <- dphi[, i] * exp(ci_v[i]) * phi[, i] * (((x - alphaarray[1, i])^2 - 1)/psiarray[1, i]) * (1/(2 * psiarray[1, 
            i])) * (1/D^2)
    }
    
    if(classes > 1){
        for (i in 1:(classes - 1)) {
            for (j in 1:(classes)) {
                if (i != j) {
                    dc[, i] <- dc[, i] + (exp(ci_v[j]) * phi[, j] * ((alphaarray[2, i] - alphaarray[2, j]) + (gamma[i] - gamma[j]) * 
                      x))
                    
                    
                }
            }
            dc[, i] <- dc[, i] * exp(ci_v[i]) * phi[, i] * (1/D^2)
        }
    }
    
    dkappa <- array(data = 0, c(points, classes))
    dgamma <- array(data = 0, c(points, classes))
    for (i in 1:classes) {
        dkappa[, i] <- exp(ci_v[i]) * phi[, i]/D
        dgamma[, i] <- exp(ci_v[i]) * phi[, i] * x/D
    }
    
    
    ct <- 0
    varordered <- c()
    for (i in 1:classes) {
        varordered <- c(varordered, alpha_loc[i + ct], alpha_loc[i + 1 + ct], beta_loc[i])
        ct <- ct + 1
    }
    
    varordered <- c(varordered, c_loc)
    
    ct <- 0
    for (i in 1:classes) {
        varordered <- c(varordered, psi_loc[i + ct])
        ct <- ct + 1
    }
    
    acovd <- acov[varordered, varordered]
    
    
    deriv <- c()
    for (i in 1:classes) {
        deriv <- cbind(deriv, dalpha[, i], dkappa[, i], dgamma[, i])
    }
    
    deriv <- c(deriv, dc, dphi)
    
    deriv <- matrix(deriv, nrow = points, ncol = p)
    
    se <- sqrt(diag(deriv %*% acovd %*% t(deriv)))
    q <- abs(qnorm(alpha/2, mean = 0, sd = 1))
    sq <- sqrt(qchisq(1 - alpha, p))
    
    # delta method nonsimultaneous confidence intervals
    lo <- y - q * se
    hi <- y + q * se
    
    # delta method simultaneous confidence intervals
    slo <- y - sq * se
    shi <- y + sq * se
                
    #boostrap CI's
    if(boot.CI){
        draws <- 1000
        bs <- bs.CE(boot, x=x, alpha=alpha, boot=TRUE)
        yall <- bs$bs.yall
        yall <- apply(yall, 2, sort)
        LCL = (alpha/2) * draws
        UCL = (1 - (alpha/2)) * draws    
        LCLall = yall[LCL, ]
        UCLall = yall[UCL, ]
    } else {
        LCLall <- UCLall <- numeric(length(lo))
    }
    
    #old prep coding
    Ksi <- Eta1
    Eta <- Eta2
    denKsi <- sumpi
    denEta <- sumpi2
    post <- pi
    pKsi <- a_pi
    pEta <- a_pi2
    
    etahmat <- matrix(data = 0, nrow = length(Ksi), ncol = classes)
    for (i in 1:classes) {
        etahmat[, i] <- alphaarray[2, i] + gamma[i] * Ksi
    }
    
    etah_ <- y
    lo_ <- lo
    hi_ <- hi
    slo_ <- slo
    shi_ <- shi
    LCLall_ <- LCLall
    UCLall_ <- UCLall   
    
    bs_lo <- bs_high <- 0
    if(boot.CE){
        bs <- bs.CE(boot, x=x, alpha=alpha, boot=FALSE)
        bs_lo <- bs$lb
        bs_high <- bs$ub
    }
    
    SEMLIdatapks <- data.frame(Eta1=Ksi, Eta2=Eta, agg_denEta1=denKsi, agg_denEta2=denEta, 
                               agg_pred=etah_, class_pred=I(etahmat), 
                               contour=I(z), classes=classes, class_prob=I(post), class_denEta1=I(pKsi), 
                               class_denEta2=I(pEta), bs_CIlo=LCLall_,
                               bs_CIhi=UCLall_, delta_CIlo=lo_, delta_CIhi=hi_, delta_CElo=slo_, 
                               delta_CEhi=shi_, x, alpha=alpha, setup2=TRUE,
                               boot=boot.CE, bs_CElo=bs_lo, bs_CEhi=bs_high)
    SEMLIdatapks
}
