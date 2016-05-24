`intervals.dif` <-
function(o, level, x.b, var, pval=TRUE, ...)
{

 x<-o
 z<-abs(qnorm((1-level)/2))

 mat.coef <- merge(x$coef, summary(x)$coef, by=0, all.x=TRUE, sort=FALSE)
 nom.pos <- data.frame(names(x$coef), ordre=1:length(x$coef))
 mat.ordre <- merge(nom.pos, mat.coef, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
 mat.ordre <- mat.ordre[order(mat.ordre$ordre),]
 
 a <- as.matrix(mat.ordre[,c("Estimate")])
 se <- as.matrix(mat.ordre[,c("Std. Error")])
 
 li <- a - (z * se)
 ls <- a + (z * se)
 
 if (missing(var))
 {
  focus <- nrow(a);
  m <- cbind(a[focus,],li[focus,],ls[focus,])
 
  if (pval)
  {
   p.as <- anova(x, x.b, test = "F")$"Pr(>F)"[2]
   m <- cbind(m, p.as)
   colnames(m) <- c("dif","lo","up","pval")
  }
  else
  {
   colnames(m) <- c("dif","lower","uppper")
  }
 }
 else
 {
  focus <- nrow(a) - length(levels(var)) + 2:length(levels(var))
  m <- cbind(a[focus,],li[focus,],ls[focus,])
  m <- rbind(c(0,NA,NA),m)
 
  if (pval)
  {
   p.as <- anova(x, x.b, test = "F")$"Pr(>F)"[2]
   m <- cbind(m, c(p.as,rep(NA,times=length(levels(var))-1)))
   colnames(m) <- c("dif","lower","upper","pval")
  }
  else
  {
   colnames(m) <- c("dif","lower","upper")
  }
 }
 
 list(m=m);
}

