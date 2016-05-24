pamr.plotstrata <-
function (fit, survival.time, censoring.status)
{
    group <-apply(fit$proby,1,which.is.max)
#    require(survival)
    n.class <- length(unique(group))
    junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
    junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))
 
  pv <- 1-pchisq(2*(junk2$loglik[2]-junk2$loglik[1]),df=n.class-1)

if(!is.null(fit$cutoffs.survival)){
    labels <- rep(NULL,n.class)
    labels[1] <- paste("(1)   ","<= ", round(fit$cutoffs.survival[1],2),sep="")
    if(n.class>2){
        for(i in 2:(n.class-1)){
          labels[i] <- paste("(",as.character(i),")  ", " > ",
        round(fit$cutoffs.survival[i-1],2), "  & <= ", 
        round(fit$cutoffs.survival[i],2), sep="")
     }}
    labels[n.class] <-  paste("(",as.character(n.class),")  ", " > ",round(fit$cutoffs.survival[n.class-1],2),sep="")
  }

else{labels <- as.character(1:n.class)}

#    win.metafile()
    plot(junk, col = 2:(2 + n.class - 1), xlab = "Time", ylab = "Probability of survival", main="Survival Strata Plot")
 #   legend(0.7 * max(fit$survival.time), 0.9, col = 2:(2 + n.class -
     legend(.01* max(fit$survival.time), 0.2, col = 2:(2 + n.class -
        1), lty = rep(1, n.class), legend = labels)
     text(0.1 * max(fit$survival.time), .25, paste("pvalue=",as.character(round(pv,4))))

#   dev.off()
#   return(TRUE)
  }
