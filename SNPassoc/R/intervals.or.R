`intervals.or` <-
function(o, level, x.b, var, ...)
{
 x<-o
 z<-abs(qnorm((1-level)/2))
 mat.coef <- merge(x$coef, summary(x)$coef, by=0, all.x=TRUE, sort=FALSE)
 nom.pos <- data.frame(names(x$coef), ordre=1:length(x$coef))
 mat.ordre <- merge(nom.pos, mat.coef, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
 mat.ordre <- mat.ordre[order(mat.ordre$ordre),]
 
 a <- as.matrix(mat.ordre[,c("Estimate")])
 se <- as.matrix(mat.ordre[,c("Std. Error")])
    
 or <- exp(a)
 li <- exp(a - z * se)
 ls <- exp(a + z * se)
 if(missing(var))
 {
  focus <- dim(a)[1.]
  or.ic <- round(cbind(or[focus,  ], li[focus,  ], ls[focus,  ]), 2.)
  or.ic[or.ic > 999.] <- NA
  t1 <- anova(x, x.b, test = "Chi")
  p.as <- t1[2, grep("^P.*Chi",names(t1))]
  or.ic <- cbind(or.ic, p.as)
  dimnames(or.ic) <- NULL
 }
    else
 {
  focus <- dim(a)[1.] - length(levels(var)) + 2:length(levels(var))
  or.ic <- round(cbind(or[focus,  ], li[focus,  ], ls[focus,  ]), 2)
  or.ic[or.ic > 999.] <- NA
  or.ic <- round(rbind(c(1, NA, NA), or.ic), 2)
  t1 <- anova(x, x.b, test = "Chi")
  p.as <- t1[2, grep("^P.*Chi",names(t1))]
  or.ic <- cbind(or.ic, c(p.as, rep(NA, times = length(levels(var)) - 1)))
  dimnames(or.ic) <- list(levels(var), c("   OR ", "lower", "upper", "p-value"))
 }

 list(or.ic = or.ic)
 
}

