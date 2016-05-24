`expandsetupSNP` <-
function (o) 
{
	if(!inherits(o,"logical")){ # all missings
    x <- summary(o)
    control<-!is.na(x$allele.freq[,2]) & x$allele.freq[,2]!=0
    o<-order(x$allele.freq[control,2],decreasing=TRUE)

    alleles <- rbind(x$allele.names)[o]
    if (length(alleles) > 1) alleles <- paste(alleles, collapse = "/")
    aux<-ifelse(any(!is.na(x$allele.freq[, 2])),
             max(x$allele.freq[, 2], na.rm = TRUE),NA)
    out <- data.frame(alleles = alleles, 
         major.allele.freq = ifelse(is.na(aux),NA,round(aux, 1)), 
         HWE = round(x$HWE, 6), missing = round(x$missing * 100, 1))
	} else {
    out <- data.frame(alleles = NA, major.allele.freq = NA, HWE = NA, missing = 100)
}
    out
}

