"epi.smd" <- function(mean.trt, sd.trt, n.trt, mean.ctrl, sd.ctrl, n.ctrl, names, method = "cohens", conf.level = 0.95)
{
    # Declarations:
    N <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N, mean = 0, sd = 1)

    k <- length(names)
    N.i <- n.trt + n.ctrl

    # Pooled standard deviation of the two groups:
    s.i <- sqrt((((n.trt - 1) * sd.trt^2) + ((n.ctrl - 1) * sd.ctrl^2)) / (N.i - 2))

    if(method == "cohens") {
    # Standardised mean difference method using Cohen's d:
    MD.i <- (mean.trt - mean.ctrl) / s.i
    SE.MD.i <- sqrt((N.i / (n.trt * n.ctrl)) + (MD.i^2 / (2 * (N.i - 2))))
    lower.MD.i <- MD.i - (z * SE.MD.i)
    upper.MD.i <- MD.i + (z * SE.MD.i)
    }

    if(method == "hedges") {
    # Standardised mean difference method using Hedge's adjusted g:
    MD.i <- ((mean.trt - mean.ctrl) / s.i) * (1 - (3/ ((4 * N.i) - 9)))
    SE.MD.i <- sqrt((N.i / ((n.trt * n.ctrl)) + (MD.i^2 / (2 * (N.i - 3.94)))))
    lower.MD.i <- MD.i - (z * SE.MD.i)
    upper.MD.i <- MD.i + (z * SE.MD.i)
    }
    
    else if(method == "glass") {
    # Standardised mean difference method using Glass's delta:
    MD.i <- (mean.trt - mean.ctrl) / sd.ctrl
    SE.MD.i <- sqrt((N.i / ((n.trt * n.ctrl)) + (MD.i^2 / (2 * (n.ctrl - 1)))))
    lower.MD.i <- MD.i - (z * SE.MD.i)
    upper.MD.i <- MD.i + (z * SE.MD.i)
    }

    # IV pooled standardised mean difference:
    w.i <- 1 / (SE.MD.i)^2
    MD.iv <- sum(w.i * MD.i) / sum(w.i)
    SE.MD.iv <- 1/sqrt((sum(w.i)))
    lower.MD.iv <- MD.iv - (z * SE.MD.iv)
    upper.MD.iv <- MD.iv + (z * SE.MD.iv)
    
    # Heterogeneity statistic:
    Q <- sum(w.i * (MD.i - MD.iv)^2)
    df <- k - 1
    p.heterogeneity <- 1 - pchisq(Q, df)
    
    tau.sq <- (Q - (k - 1)) / (sum(w.i) - (sum((w.i)^2) / sum(w.i)))

    # If Q is less than (k - 1) tau.sq equals zero:
    tau.sq <- ifelse(Q < (k - 1), 0, tau.sq)
    w.dsl.i <- 1 / (SE.MD.i^2 + tau.sq)
    MD.dsl <- sum(w.dsl.i * MD.i) / sum(w.dsl.i)
    SE.MD.dsl <- 1 / sqrt(sum(w.dsl.i))
    lower.MD.dsl <- MD.dsl - (z * SE.MD.dsl)
    upper.MD.dsl <- MD.dsl + (z * SE.MD.dsl)
    
    # Results:
    md <- as.data.frame(cbind(MD.i, SE.MD.i, lower.MD.i, upper.MD.i))
    names(md) <- c("est", "se", "lower", "upper")
    md.invar <- as.data.frame(cbind(MD.iv, SE.MD.iv, lower.MD.iv, upper.MD.iv))
    names(md.invar) <- c("est", "se", "lower", "upper")
    md.dsl <- as.data.frame(cbind(MD.dsl, lower.MD.dsl, upper.MD.dsl))
    names(md.dsl) <- c("est", "lower", "upper")
    heterogeneity = c(Q = Q, df = df, p.value = p.heterogeneity)
    
    rval <- list(md = md, md.invar = md.invar, md.dsl = md.dsl, heterogeneity = heterogeneity)
    return(rval)
}
