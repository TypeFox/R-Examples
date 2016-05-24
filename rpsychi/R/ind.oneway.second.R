ind.oneway.second <- function(m, sd, n, unbiased=TRUE, contr=NULL, sig.level=.05, digits=3){

##if orthogonal
if(length(n)==1){
     n <- rep(n, length(m))
}

#if biased standard deviation
if(unbiased==FALSE){
         sd <- ssd2sd(n,sd)
}


##(a) anova table
k  <- length(m)           #number of groups
Xg  <- sum(n*m)/sum(n)

dfb <- k - 1              #degree of freedom
dfw   <- sum(n) - k       #degree of freedom

MSb <- sum(n * (m - Xg)^2)/(k-1)  #MS between
MSw <-  sum((n-1)*sd^2)/dfw       #MS within
SSb <- dfb * MSb
SSw <- dfw * MSw
SSt <- SSb + SSw

f.value <- MSb/MSw                #f value


anova.table  <- data.frame(matrix(NA,ncol=4, nrow=3))
rownames(anova.table) <- c("Between (A)", "Within", "Total")
colnames(anova.table) <- c("SS", "df", "MS", "F")
anova.table$SS <- c(SSb, SSw,SSt)
anova.table$df <- c(dfb, dfw, dfb+dfw)
anova.table$MS <- c(MSb,MSw,NA)
anova.table$F  <- c(f.value, NA,NA)
class(anova.table) <- c("anova", "data.frame")
anova.table <- round(anova.table, digits)



##(b) omnibus effect size eta and omega squared

#eta square
etasq <- SSb / SSt
delta.lower <- delta.upper <- numeric(length(etasq))
delta.lower <- try(FNONCT(f.value, dfb, dfw, prob=1-sig.level/2), silent=TRUE)
delta.upper <- try(FNONCT(f.value, dfb, dfw, prob=sig.level/2), silent=TRUE)
if(is.character(delta.lower)){
  delta.lower <- 0
}

etasq.lower <- delta.lower / (delta.lower + dfb + dfw + 1)
etasq.upper <- delta.upper / (delta.upper + dfb + dfw + 1)


#omega square
omegasq <- (SSb - dfb * MSw)/(SSt + MSw)
sosb_L  <- SSt * etasq.lower
msw_L   <- (SSt - sosb_L)/dfw
omegasq.lower <- (sosb_L - (dfb*msw_L))/(SSt+msw_L)

sosb_U  <- SSt * etasq.upper
msw_U   <- (SSt - sosb_U)/dfw
omegasq.upper <- (sosb_U - (dfb*msw_U))/(SSt+msw_U)

omnibus.es <- round(c(etasq=etasq, etasq.lower=etasq.lower, etasq.upper=etasq.upper),
                      digits)


##(c) raw contrasts
temp  <- combinations(k,2)
cont1 <- matrix(0, nrow=nrow(temp),ncol=k)
cont1.lab <- rep(0,nrow(temp))

#in case did not specify contrasts
for(i in 1:nrow(temp)){
     cont1[i, temp[i,1]] <- 1
     cont1[i, temp[i,2]] <- -1
     cont1.lab[i] <- paste(temp[i,1],"-",temp[i,2], sep="")
     rownames(cont1) <- cont1.lab
}


#in case specify contrasts
if(!is.null(contr)){
      if(is.vector(contr)){
                cont1 <- t(as.matrix(contr))
      }else{
                cont1 <- contr
      }
}


#F test for contrasts
psi <-  colSums(t(cont1)  * as.vector(m))                #raw contrasts
SSpsi <- (psi^2)/colSums(t(cont1^2) / as.vector(n))

nmat <- matrix(n, nrow=nrow(cont1), ncol=length(n), byrow=TRUE)
psi.std <- sqrt(MSw * rowSums(cont1 ^ 2/nmat))

psi.lower <- psi + psi.std * qt(sig.level/2, dfw)
psi.upper <- psi + psi.std * qt(sig.level/2, dfw, lower.tail=FALSE)

raw.contrasts <- round(data.frame(mean.diff=psi, lower=psi.lower, upper=psi.upper, std=psi.std), digits)
rownames(raw.contrasts) <- rownames(cont1)



##(d) standardized contrasts
gpsi <- psi/sqrt(MSw)       #effect size
gpsi.std <- sqrt(rowSums(cont1 ^ 2/nmat))
gpsi.lower <- gpsi + gpsi.std * qt(sig.level/2, dfw)
gpsi.upper <- gpsi + gpsi.std * qt(sig.level/2, dfw, lower.tail=FALSE)
standardized.contrasts <- round(data.frame(es=gpsi, lower=gpsi.lower, upper=gpsi.upper, std=gpsi.std), digits)
rownames(standardized.contrasts) <- rownames(cont1)



##(e) statistical power
c.delta <- c(.10, .25, .4)
criterion.power <- round(power.f(sig.level=sig.level, u=dfb, n=sum(n)/k,delta=c.delta), digits)
names(criterion.power) <- c("small", "medium", "large")



##(e) output
output <- list(anova.table=anova.table, omnibus.es=omnibus.es, raw.contrasts=raw.contrasts, standardized.contrasts = standardized.contrasts, power=criterion.power)
return(output)
}

