mtdt2 <- function(x,verbose=TRUE,n.sim=NULL,...)
{
  for(p in c("BradleyTerry2")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!require(p, quietly = TRUE, character.only=TRUE))
        warning(paste("mtdt2 needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  dims <- dim(x)[1]
  colnames(x) <- paste(1:dims,sep="")
  rownames(x) <- paste(1:dims,sep="")
  c2b <- as.data.frame(BradleyTerry2::countsToBinomial(x))
  names(c2b) <- c("allele1","allele2","transmitted","nontransmitted")
  allele1 <- with(c2b, allele1)
  allele2 <- with(c2b, allele2)
  transmitted <- with(c2b, transmitted)
  nontransmitted <- with(c2b, nontransmitted)
  btx <- BradleyTerry2::BTm(cbind(transmitted,nontransmitted), allele1, allele2, ~ allele, id="allele", data=c2b, ...)
  t1 <- btx$null.deviance
  t2 <- btx$deviance
  t3 <- t1-t2
  f1 <- btx$df.residual
  f2 <- btx$df.null
  f3 <- f2-f1
  ctest <- c(t3,t1,t2)
  df <- c(f3,f2,f1)
  pexp <- pchisq(ctest,df,lower.tail=FALSE,log.p=TRUE)/log(10)
  p <- 10^pexp
  transmissions <- with(c2b,transmitted+nontransmitted)
  if(!missing(n.sim))
  {
    transmissionlen <- length(transmissions)
    transmittedn <- with(c2b,transmitted)
    nontransmittedn <- with(c2b,nontransmitted)
    pn <- rep(0,3)
    for(i in 1:n.sim)
    {
       for(j in 1:transmissionlen) transmittedn[j] <- rbinom(1,transmissions[j],0.5)
       nontransmittedn <- transmissions-transmittedn
       c2bn <- data.frame(c2b,transmittedn,nontransmittedn)
       btn <- BradleyTerry2::BTm(cbind(transmittedn,nontransmittedn), allele1, allele2, ~ allele, id="allele", data=c2bn, ...)
       t1 <- btn$null.deviance
       t2 <- btn$deviance
       t3 <- t1-t2
       ctestn <- c(t3,t1,t2)
       if(ctestn[1]>=ctest[1]) pn[1] <- pn[1]+1
       if(ctestn[2]>=ctest[2]) pn[2] <- pn[2]+1
       if(ctestn[3]>=ctest[3]) pn[3] <- pn[3]+1
    }
    pn <- (pn+1)/(n.sim+1)
  }
  if(verbose)
  {
    print(cbind(c2b,transmissions,fitted=btx$fitted.values*transmissions))
    print(summary(btx,corr=TRUE))
    cat("Chi-square for allele-wise TDT =",ctest[1],"df =",df[1],"p =",p[1],"\n")
    cat("Chi-square for genotype-wise TDT =",ctest[2],"df =",df[2],"p =",p[2],"\n")
    cat("Chi-square for goodness-of-fit of allele-wise model =",ctest[3],"df =",df[3],"p =",p[3],"\n")
    if(!missing(n.sim)) cat("The corresponding Monte Carlo p values are",pn,"\n")
  }
  if(missing(n.sim)) invisible(list(c2b=c2b,BTm=btx,X2=ctest,df=df,p=p))
  else invisible(list(c2b=c2b,BTM=btx,X2=ctest,df=df,p=p,pn=pn))
}

# created on 20-4-2010
# last updated on 22-4-2010
# As BTm(data=) is not working, we need get around with attachment or global assignment.
