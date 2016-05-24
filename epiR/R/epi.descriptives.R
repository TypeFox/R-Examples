epi.descriptives <- function(dat, conf.level = 0.95){
    tcrit <- qt((1 - conf.level) / 2, length(dat), lower.tail = FALSE)

    an <- length(dat)
    amean <- mean(dat, na.rm = TRUE)
    amedian <- as.vector(quantile(dat, probs = 0.50, na.rm = TRUE))
    avar <- var(dat, na.rm = TRUE)
    asd <- sd(dat, na.rm = TRUE)
    ase <- asd / sqrt(an)

    ana <- is.na(dat); ana <- sum(as.numeric(ana))
    
    aq25 <- as.vector(quantile(dat, probs = 0.25, na.rm = TRUE))
    aq75 <- as.vector(quantile(dat, probs = 0.75, na.rm = TRUE))
    alcl <- mean(dat, na.rm = TRUE) - (tcrit * ase)
    aucl <- mean(dat, na.rm = TRUE) + (tcrit * ase)
    
    amin <- min(dat, na.rm = TRUE)
    amax <- max(dat, na.rm = TRUE)
    
    # Geometric mean. Make sure all values positive:
    tdat <- dat[dat > 0]
    
    gn <- length(tdat)
    gmean <- exp(mean(log(tdat), na.rm = TRUE)) 
    gsd <- sd(log(tdat), na.rm = TRUE)
    gse <- gsd / sqrt(gn)
    
    gna <- is.na(tdat); gna <- sum(as.numeric(gna))
    
    gq25 <- as.vector(exp(quantile(log(tdat), probs = 0.25, na.rm = TRUE)))
    gq75 <- as.vector(exp(quantile(log(tdat), probs = 0.75, na.rm = TRUE)))

    glcl <- exp(mean(log(tdat), na.rm = TRUE) - (tcrit * gse))
    gucl <- exp(mean(log(tdat), na.rm = TRUE) + (tcrit * gse))
    
    gmin <- as.vector(exp(min(log(tdat), na.rm = TRUE)))
    gmax <- as.vector(exp(max(log(tdat), na.rm = TRUE)))
    
    # Skewness:
    x <- dat - mean(dat, na.rm = TRUE)
    skew <- sqrt(an) * sum(x^3, na.rm = TRUE) / (sum(x^2, na.rm = TRUE)^(3/2))
    
    # Kurtosis:
    x <- dat - mean(dat, na.rm = TRUE)
    r <- an * sum(x^4, na.rm = TRUE) / (sum(x^2, na.rm = TRUE)^2)
    kurt <- ((an + 1) * (r - 3) + 6) * (an - 1)/((an - 2) * (an - 3))
    
    rval <- list(
        arithmetic = data.frame(n = an, mean = amean, sd = asd, q25 = aq25, q75 = aq75, lower = alcl, upper = aucl, 
            min = amin, max = amax, na = ana),
        geometric = data.frame(n = gn, mean = gmean, sd = gsd,  q25 = gq25, q75 = gq75, lower = glcl, upper = gucl, 
            min = gmin, max = gmax, na = gna),
        symmetry = data.frame(skewness = skew, kurtosis = kurt)
    )
        
   return(rval)
 }
 
 
