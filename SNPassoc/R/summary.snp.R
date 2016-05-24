`summary.snp` <-
function (object, ...) 
{
    n <- length(object)
    nas <- is.na(object)
    n.typed <- n - sum(nas)
    ll <- levels(object)
    tbl <- table(object)
    tt <- c(tbl)
    names(tt) <- dimnames(tbl)[[1]]
    if (any(nas)) 
      { 
        tt.g <- c(tt, "NA's" = sum(nas))
        missing.allele<-sum(nas)/(sum(tt)+sum(nas))
      }
    else
      {
       tt.g <- tt
       missing.allele<-0
      }
    tt.g.prop <- prop.table(tbl)
    if (any(nas)) 
        tt.g.prop <- c(tt.g.prop, NA)
    ans.g <- cbind(frequency = tt.g, percentage = tt.g.prop * 100)

    alle <- attr(object, "allele.names")
    alle1 <- length(grep(paste(alle[1], "/", sep = ""), as.character(object))) + 
             length(grep(paste("/", alle[1], sep = ""), as.character(object)))
    if (length(alle) > 1) {

        alle2 <- length(grep(paste(alle[2], "/", sep = ""), as.character(object))) + 
            length(grep(paste("/", alle[2], sep = ""), as.character(object)))
        tt.a <- c(alle1, alle2)
        tt.a.prop <- prop.table(tt.a)
        ans.a <- cbind(frequency = tt.a, percentage = tt.a.prop * 100)
        pvalueHWE <- SNPHWE(c(tbl,0,0)[1:3])  # VM make sure 3 genotypes are sent i			
        dimnames(ans.a)[[1]] <- alle
     }
     else {
        tt.a <- alle1
        tt.a.prop <- prop.table(tt.a)
        ans.a <- t(c(frequency = tt.a, percentage = tt.a.prop * 100))
        rownames(ans.a)<-alle
        pvalueHWE <- NA
     }
     if (any(nas)) 
         ans.a <- rbind(ans.a, "NA's" = c(2 * sum(nas), NA))
    
    ans <- list(allele.names = alle, allele.freq = ans.a, genotype.freq = ans.g, 
        n = n, n.typed = n.typed, HWE = pvalueHWE, missing.allele=missing.allele)
    class(ans) <- "summary.snp"
    ans
}

