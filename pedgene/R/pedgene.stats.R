
## Title: pedgene.stats.R
## Purpose: worker function to calculate ped-based kernel and burden tests for rare variants
##          on one gene at a time. Wrapper function over multiple genes is pedgene()
## Author: Jason Sinnwell
## Created: 8/2013
## Updated: 6/4/2015

pedgene.stats <- function(geno, c.factor, chrom, male.dose, sex, resid,
                          weights=NULL, weights.mb=FALSE, weights.beta=c(1,25),
                          method="davies", acc.davies=1e-6) {

    x.chrom <- as.character(chrom)=="X"
    
    # remove monomorphic markers
    # for x.chrom, genotype codes may be different so make sure not monomorphic in both sexes
    if(x.chrom) {
      v.males <- apply(geno[sex==1,,drop=FALSE], 2, var, na.rm = TRUE)
      v.females <- apply(geno[sex==2,,drop=FALSE], 2, var, na.rm = TRUE)
      v <- pmax(v.males, v.females)
    } else {
      v <- apply(geno, 2, var, na.rm = TRUE)
    }
  
    nvariant0 <- ncol(geno)
    geno <- geno[, v > 0,drop=FALSE]
    geno <- as.matrix(geno)
    nvariant <- ncol(geno)
    nvariant.noninform <- nvariant0 - nvariant
    if(nvariant==0) {
        return(list(stat.kernel = NA, pval.kernel = NA,  
                    stat.burden = NA, pval.burden = NA,
                    nvariant = nvariant, noninform=nvariant.noninform))
      }
    
    # Account for missing genotypes. Could exclude subjects with
    # any missing, but  here just replace missing with mean values.
    # Would be best for user to decide how to handle missing prior
    # to calling this function.

    col.miss <- apply(is.na(geno), 2, any)
   
    if(any(col.miss == TRUE))
      {
        for(j in 1:ncol(geno))
          {
            if(col.miss[j])
              {
                if(x.chrom)
                  {
                    mn.male <- mean(geno[sex==1,j], na.rm=TRUE)
                    geno[is.na(geno[,j])&sex==1, j] <- mn.male
                    mn.female <- mean(geno[sex==2,j], na.rm=TRUE)
                    geno[is.na(geno[,j])&sex==2, j] <- mn.female          
                  } else
                {
                  mn <- mean(geno[,j], na.rm=TRUE)
                  geno[is.na(geno[,j]), j] <- mn
                }
              }
          }
      }
   
    if(x.chrom) {
      count.male   <- apply(geno[sex==1,,drop=FALSE], 2, sum)
      count.female <- apply(geno[sex==2,,drop=FALSE], 2, sum)
      n.male   <- sum(sex==1)
      n.female <- sum(sex==2)
      maf <- (count.male + count.female)/(n.male + 2*n.female)
    } else {
      maf <- apply(geno/2, 2, mean)
    }   
    if(is.null(weights)) {    
      if(weights.mb==TRUE) {        
        wt <- 1/sqrt(maf * (1-maf))
      } else {
        ## weights=NULL and weights.MB=FALSE, all default to do weights.beta
        wt <- dbeta(maf, weights.beta[1], weights.beta[2])
      }
    } else {
      wt <- weights[ v > 0 ]
    }

    # estimate cor among markers
    r.mat <- cor(geno)
   
    ## compute f from top of p412.  results in f=1 in M-B weights,
    ## still needed for beta and user-specified weights, b/c it is binomial variance
    
    f <- wt * sqrt(maf*(1-maf))
    
    fRmat <- (f %o% f) * r.mat

    var.z <- fRmat * c.factor    
    
    # score males according to male.dose
    if(x.chrom & (male.dose !=1) )  {
      geno.score <- geno     
      geno.score[sex==1,] <- geno[sex==1,]*male.dose       
      ## kernel stat info
      kmat <- geno.score %*% diag(wt^2,nrow=length(wt),ncol=length(wt)) %*% t(geno.score)
      ## Burden stat info
      burden.score.subject <- as.vector(geno.score %*% wt)        
    } else {     
      ## kernel stat info
      kmat <- geno %*% diag(wt^2,nrow=length(wt),ncol=length(wt)) %*% t(geno)
      ## Burden stat info
      burden.score.subject <- as.vector(geno %*% wt)
    }

    ########### Burden stat (2-sided)     
    stat.num <- (resid %*% burden.score.subject)
    factor.sum <- sum(fRmat)
    stat.denom <- sqrt(factor.sum * c.factor)    
    burden.stat <- stat.num / stat.denom
    burden.pval <- pchisq(burden.stat^2, 1,lower.tail=FALSE)
    
    ########## Quadratic kernel stat
    ## if only 1 variant, reduces to burden stat, o/w do kernel test
    if(nvariant>1) {
      stat.kernel <- as.vector(resid %*% kmat %*% resid)
      eig <-eigen(var.z, symmetric=T, only.values=T)
      evals <-eig$values[eig$values>1e-6*eig$values[1]]
      if(method=="davies") {
        pval.kernel <- davies(stat.kernel , evals, acc=acc.davies)$Qq        
        ## Davies' method sometimes instable, returns out-of-range p-value
        ## Help file suggests lerge acc.davies of 1e-4 or so to fix
      }
      if(method=="kounen") {
        ## in main method, require(survery) for kounen
        pval.kernel <- pchisqsum(stat.kernel, rep(1, length(evals)),
                           evals, lower.tail=FALSE, method="saddle")
      }

      ## davies and kounen sometimes return out of range p-values.
      ## fix to 1 or 0 on either side.
      pval.kernel <- min(pval.kernel, 1)
      pval.kernel <- max(pval.kernel,0)
     
    } else {
      stat.kernel <- burden.stat^2  
      pval.kernel <- burden.pval
    }     
    
    lst <- list(stat.kernel = stat.kernel, 
                     pval.kernel = pval.kernel,
                     stat.burden = burden.stat,
                     pval.burden = burden.pval,
                     nvariant = nvariant,
                     noninform=nvariant.noninform)

    return(lst)
    
}
