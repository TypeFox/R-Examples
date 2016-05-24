"epi.mh" <- function(ev.trt, n.trt, ev.ctrl, n.ctrl, names, method = "odds.ratio", alternative = c("two.sided", "less", "greater"), conf.level = 0.95)
    {
        # Declarations:
        k <- length(names)
        a.i <- ev.trt
        b.i <- n.trt - ev.trt
        c.i <- ev.ctrl
        d.i <- n.ctrl - ev.ctrl

        N. <- 1 - ((1 - conf.level) / 2)
        z <- qnorm(N., mean = 0, sd = 1)
                
        # Test each strata for zero values. Add 0.5 to all cells if any cell has a zero value:
        for(i in 1:k){
        if(a.i[i] < 1 | b.i[i] < 1 | c.i[i] < 1 | d.i[i] < 1){
           a.i[i] <- a.i[i] + 0.5; b.i[i] <- b.i[i] + 0.5; c.i[i] <- c.i[i] + 0.5; d.i[i] <- d.i[i] + 0.5
           }
        }

        n.1i <- a.i + b.i
        n.2i <- c.i + d.i
        N.i <- a.i + b.i + c.i + d.i
 
        # For summary odds ratio:
        R <-  sum((a.i * d.i) / N.i)
        S <-  sum((b.i * c.i) / N.i)
        E <-  sum(((a.i + d.i) * a.i * d.i) / N.i^2)
        F. <- sum(((a.i + d.i) * b.i * c.i) / N.i^2)
        G <-  sum(((b.i + c.i) * a.i * d.i) / N.i^2) 
        H <-  sum(((b.i + c.i) * b.i * c.i) / N.i^2)

        # For summary risk ratio:
        P <- sum(((n.1i * n.2i * (a.i + c.i)) - (a.i * c.i * N.i)) / N.i^2)
        R. <- sum((a.i * n.2i) / N.i)
        S. <- sum((c.i * n.1i) / N.i)
        
        if(method == "odds.ratio"){  
            # Individual study odds ratios:
            OR.i <- (a.i * d.i) / (b.i * c.i)   
            lnOR.i <- log(OR.i)
            SE.lnOR.i <- sqrt(1/a.i + 1/b.i + 1/c.i + 1/d.i)
            SE.OR.i <- exp(SE.lnOR.i)
            lower.lnOR.i <- lnOR.i - (z * SE.lnOR.i)
            upper.lnOR.i <- lnOR.i + (z * SE.lnOR.i)
            lower.OR.i <- exp(lower.lnOR.i)
            upper.OR.i <- exp(upper.lnOR.i)
                
            # Weights:
            w.i <- (b.i * c.i) / N.i
            # w.i <- 1 / (1/a.i + 1/b.i + 1/c.i + 1/d.i)
            w.iv.i <- 1 / (SE.lnOR.i)^2
        
            # MH pooled odds ratios (relative effect measures combined in their natural scale):
            OR.mh <- sum(w.i * OR.i) / sum(w.i)                
            lnOR.mh <- sum(w.i * log(OR.i)) / sum(w.i)
            
            # Same method for calculating confidence intervals around pooled OR as epi.2by2 so results differ from page 304 of Egger, Smith and Altman:            
            G <- a.i * d.i / N.i
            H <- b.i * c.i / N.i
            P <- (a.i + d.i) / N.i
            Q <- (b.i + c.i) / N.i
            GQ.HP <- G * Q + H * P
            sumG <- sum(G)
            sumH <- sum(H)
            sumGP <- sum(G * P)
            sumGH <- sum(G * H)
            sumHQ <- sum(H * Q)
            sumGQ <- sum(G * Q)
            sumGQ.HP <- sum(GQ.HP)
            var.lnOR.mh <- sumGP / (2 * sumG^2) + sumGQ.HP/(2 * sumGH) + sumHQ/(2 * sumH^2)
            SE.lnOR.mh <- sqrt(var.lnOR.mh)
            SE.OR.mh <- exp(SE.lnOR.mh)
            lower.OR.mh <- exp(lnOR.mh - z * SE.lnOR.mh)
            upper.OR.mh <- exp(lnOR.mh + z * SE.lnOR.mh)
               
            # Test of heterogeneity (based on inverse variance weights):
            Q <- sum(w.iv.i * (lnOR.i - lnOR.mh)^2)
            df <- k - 1
            p.heterogeneity <- 1 - pchisq(Q, df)
    
            # Higgins and Thompson (2002) H^2 and I^2 statistic:
            Hsq <- Q / (k - 1)
            lnHsq <- log(Hsq)
            if(Q > k) {
              lnHsq.se <- (1 * log(Q) - log(k - 1)) / (2 * sqrt(2 * Q) - sqrt((2 * (k - 3))))
              }
            if(Q <= k) {
              lnHsq.se <- sqrt((1/(2 * (k - 2))) * (1 - (1 / (3 * (k - 2)^2))))
              }
            lnHsq.l <- lnHsq - (z * lnHsq.se)
            lnHsq.u <- lnHsq + (z * lnHsq.se)
            Hsq.l <- exp(lnHsq.l)
            Hsq.u <- exp(lnHsq.u)
            Isq <- ((Hsq - 1) / Hsq) * 100
            Isq.l <- ((Hsq.l - 1) / Hsq.l) * 100
            Isq.u <- ((Hsq.u - 1) / Hsq.u) * 100      
            
            # Test of effect. Code for p-value taken from z.test function in TeachingDemos package:
            effect.z <- lnOR.mh / SE.lnOR.mh
            alternative <- match.arg(alternative)
            p.effect <- switch(alternative, two.sided = 2 * pnorm(abs(effect.z), lower.tail = FALSE), less = pnorm(effect.z), greater = pnorm(effect.z, lower.tail = FALSE))
                
            # Results:
            OR <- as.data.frame(cbind(OR.i, SE.OR.i, lower.OR.i, upper.OR.i))
            names(OR) <- c("est", "se", "lower", "upper")
               
            OR.summary <- as.data.frame(cbind(OR.mh, SE.OR.mh, lower.OR.mh, upper.OR.mh))
            names(OR.summary) <- c("est", "se", "lower", "upper")
                                
            weights <- as.data.frame(cbind(w.i, w.iv.i))
            names(weights) <- c("raw", "inv.var")
                
            Hsq <- as.data.frame(cbind(Hsq, Hsq.l, Hsq.u))
            names(Hsq) <- c("est", "lower", "upper")
        
            Isq <- as.data.frame(cbind(Isq, Isq.l, Isq.u))
            names(Isq) <- c("est", "lower", "upper")
        
            rval <- list(OR = OR, OR.summary = OR.summary, weights = weights,
            heterogeneity = c(Q = Q, df = df, p.value = p.heterogeneity),
            Hsq = Hsq,
            Isq = Isq,
            effect = c(z = effect.z, p.value = p.effect))
            }
            
        else
        if(method == "risk.ratio"){
            # Individual study risk ratios:
            RR.i <- (a.i / n.1i) / (c.i / n.2i) 
            lnRR.i <- log(RR.i)
            SE.lnRR.i <- sqrt(1/a.i + 1/c.i - 1/n.1i - 1/n.2i)
            SE.RR.i <- exp(SE.lnRR.i)
            lower.lnRR.i <- lnRR.i - (z * SE.lnRR.i)
            upper.lnRR.i <- lnRR.i + (z * SE.lnRR.i)
            lower.RR.i <- exp(lower.lnRR.i)
            upper.RR.i <- exp(upper.lnRR.i)

            # Weights:
            w.i <- (c.i * n.1i) / N.i
            w.iv.i <- 1 / (SE.lnRR.i)^2
        
            # MH pooled odds ratios (relative effect measures combined in their natural scale):
            RR.mh <- sum(w.i * RR.i) / sum(w.i)
            lnRR.mh <- log(RR.mh)
            SE.lnRR.mh <- sqrt(P / (R. * S.))
            SE.RR.mh <- exp(SE.lnRR.mh)
            lower.lnRR.mh <- log(RR.mh) - (z * SE.lnRR.mh)
            upper.lnRR.mh <- log(RR.mh) + (z * SE.lnRR.mh)
            lower.RR.mh <- exp(lower.lnRR.mh)
            upper.RR.mh <- exp(upper.lnRR.mh)

            # Test of heterogeneity (based on inverse variance weights):
            Q <- sum(w.iv.i * (lnRR.i - lnRR.mh)^2)
            df <- k - 1
            p.heterogeneity <- 1 - pchisq(Q, df)
                
            # Higgins and Thompson (2002) H^2 and I^2 statistic:
            Hsq <- Q / (k - 1)
            lnHsq <- log(Hsq)
            if(Q > k) {
               lnHsq.se <- (1 * log(Q) - log(k - 1)) / (2 * sqrt(2 * Q) - sqrt((2 * (k - 3))))
               }
            if(Q <= k) {
               lnHsq.se <- sqrt((1/(2 * (k - 2))) * (1 - (1 / (3 * (k - 2)^2))))
               }
            lnHsq.l <- lnHsq - (z * lnHsq.se)
            lnHsq.u <- lnHsq + (z * lnHsq.se)
            Hsq.l <- exp(lnHsq.l)
            Hsq.u <- exp(lnHsq.u)
            Isq <- ((Hsq - 1) / Hsq) * 100
            Isq.l <- ((Hsq.l - 1) / Hsq.l) * 100
            Isq.u <- ((Hsq.u - 1) / Hsq.u) * 100
                
            # Test of effect. Code for p-value taken from z.test function in TeachingDemos package:
            effect.z <- log(RR.mh) / SE.lnRR.mh
            alternative <- match.arg(alternative)
            p.effect <- switch(alternative, two.sided = 2 * pnorm(abs(effect.z), lower.tail = FALSE), less = pnorm(effect.z), greater = pnorm(effect.z, lower.tail = FALSE))
               
            # Results:
            RR <- as.data.frame(cbind(RR.i, SE.RR.i, lower.RR.i, upper.RR.i))
            names(RR) <- c("est", "se", "lower", "upper")
                
            RR.summary <- as.data.frame(cbind(RR.mh, SE.RR.mh, lower.RR.mh, upper.RR.mh))
            names(RR.summary) <- c("est", "se", "lower", "upper")
                
            weights <- as.data.frame(cbind(w.i, w.iv.i))
            names(weights) <- c("raw", "inv.var")
               
            Hsq <- as.data.frame(cbind(Hsq, Hsq.l, Hsq.u))
            names(Hsq) <- c("est", "lower", "upper")
        
            Isq <- as.data.frame(cbind(Isq, Isq.l, Isq.u))
            names(Isq) <- c("est", "lower", "upper")
        
            rval <- list(RR = RR, RR.summary = RR.summary, weights = weights,
            heterogeneity = c(Q = Q, df = df, p.value = p.heterogeneity),
            Hsq = Hsq,
            Isq = Isq,
            effect = c(z = effect.z, p.value = p.effect))
            }
    return(rval)
}
