dep.oneway.second <- function(m, sd, n, corr, unbiased=TRUE, contr=NULL, sig.level=.05, digits=3){

##if orthogonal
if(length(n)==1){
     n <- rep(n, length(m))
}

#if biased standard deviation
if(unbiased==FALSE){
      sd <- ssd2sd(n,sd)
}


##(a) anova table
cov <- r2cov(sd=sd, R=corr)
Mcov <- mean(cov[lower.tri(cov)])
a   <- length(m)
Mt  <- mean(m)
N   <- sum(n)
      
#between effect
dfa <- a-1
MSa <- sum(n*(m-Mt)^2)/(a-1)
SSa <- dfa * MSa
          
#within effect
dfw <- N-a
MSw <- sum((n-1)*sd^2)/sum(n-1)
SSw <- dfw * MSw
      

#a*S error
MSas <- MSw-Mcov
dfas <- (a-1)*(n[[1]]-1)
SSas <- dfas * MSas
      
#Subject effect
dfs <- n[[1]]-1
SSs <- SSw - SSas
MSs <- SSs/dfs
      
f.value <- MSa/MSas
p.value <- pf(f.value, dfa, dfas, lower.tail=FALSE)
           
#total df/total SS
dft <- dfa + dfw
SSt <- SSa + SSw


#
anova.table  <- data.frame(matrix(NA,ncol=5, nrow=5))
rownames(anova.table) <- c("Between (A)", "Within", "Subjects (S)", "A * S (Error)", "Total")
colnames(anova.table) <- c("SS", "df", "MS", "F", "p.value")
anova.table$SS <- c(SSa, SSw, SSs, SSas, SSt)
anova.table$df <- c(dfa, dfw, dfs, dfas, dft)
anova.table$MS <- c(MSa, MSw,  MSs, MSas ,NA)
anova.table$F  <- c(f.value, NA, NA, NA, NA)
anova.table$p.value   <- c(p.value, NA,NA, NA,NA)
class(anova.table) <- c("anova", "data.frame")
anova.table <- round(anova.table, digits)



##(b) omnibus effect size eta and omega squared
etasq <- SSa / (SSa + SSas)
omegasq <- (dfa * (MSa - MSas)) / (dfa * (MSa - MSas) + a* n[1]*MSas)
   
omnibus.es <- round(c(partial.etasq=etasq), digits)


##(c) raw contrasts
temp  <- combinations(a, 2)
cont1 <- matrix(0, nrow=nrow(temp),ncol=a)
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
var.d.psi <- numeric(nrow(cont1))
for(i in 1:length(var.d.psi)){
  sect1 <- sum(cont1[i,] ^ 2 * sd^2)
  sect2 <- cont1[i,] %*% t(cont1[i,]) * cov
  var.d.psi[i] <- sect1 + sum(sect2[lower.tri(sect2)]) * 2 
}

psi.std <- sqrt(var.d.psi/n[1])
psi.lower <- psi + psi.std * qt(sig.level/2, dfs)
psi.upper <- psi + psi.std * qt(sig.level/2, dfs, lower.tail=FALSE)

raw.contrasts <- round(data.frame(mean.diff=psi, lower=psi.lower, upper=psi.upper, std=psi.std), digits)
rownames(raw.contrasts) <- rownames(cont1)



##(d) standardized contrasts
gpsi <- psi/sqrt(MSw)
gpsi.lower <- psi.lower/sqrt(MSw)
gpsi.upper <- psi.upper/sqrt(MSw)
standardized.contrasts <- round(data.frame(es=gpsi, lower=gpsi.lower, upper=gpsi.upper), digits)
rownames(standardized.contrasts) <- rownames(cont1)



##(e) output
output <- list(anova.table=anova.table, omnibus.es=omnibus.es, raw.contrasts=raw.contrasts, standardized.contrasts = standardized.contrasts)
return(output)
}
