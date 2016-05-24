buck.model <- 
    glm(dead~conc,data=buckthorn, family=binomial)
summary(buck.model)
# calculate the percentage dead at each concentration used.
tbl <- aggregate(buckthorn$dead,
        by=list(buckthorn$conc),
        FUN=function(x){sum(x)/length(x)}) 
names(tbl) <- c('conc','propDead'); tbl
concentrations = seq(0,0.5,by=0.02)
fits <- predict(buck.model,new=data.frame(conc=concentrations), 
            type="response")
buck.plot <- xyplot(fits~concentrations,type='l',
    ylab="predicted death rate",xlab="concentration of glyphosate",
    panel=function(x,y,...) {
        panel.xyplot(tbl$conc,tbl$propDead)
        panel.xyplot(x,y,...)
    })

observed <- xtabs(~dead+conc,data=buckthorn); observed
colTotals <- apply(observed,2,sum); colTotals
fits <- predict(buck.model,new=data.frame(conc=unique(buckthorn$conc)), 
            type="response")
fits
expected <- rbind( (1-fits)*colTotals, fits*colTotals ); expected
lrt <- 2 * sum( observed *  log (observed/expected) ); lrt
pearson <- sum( ( observed - expected )^2 / expected ); pearson
# pvals
1 - pchisq(pearson, df=2)
1 - pchisq(lrt, df=2)
