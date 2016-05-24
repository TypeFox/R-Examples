`getSignificantSNPs` <-
function (x, chromosome, model, sig = 1e-15) 
{
    if (!inherits(x, "WGassociation")) 
        stop("x must be an object of class 'WGassociation'")

    if(is.null(attr(x, "gen.info")))
     {
        pvalues <- attr(x, "pvalues")
        mm<-charmatch(model,dimnames(pvalues)[[2]])
        if (is.na(mm))
          stop("model selected is not correct")  
        sel2<-pvalues[,mm]<=sig
        SNPs.sel <- pvalues[sel2, ]
        pos.sel <- attr(x, "colSNPs")[sel2]
        out <- list(names = dimnames(SNPs.sel)[[1]], column = pos.sel)
     }

    else
     {
      if (!chromosome %in% c(1:22) & chromosome != "X") 
        stop("chromosome should be either a number between 1 and 22 or X")
      if (chromosome == "X") 
        chromosome <- 23
      gen.info <- attr(x, "gen.info")
      pvalues <- attr(x, "pvalues")
      chrs <- gen.info[, 2]
      chr.l <- unique(chrs)
      chr <- chr.l[orderChromosome(chr.l)]
      sel <- chr[chromosome]
      SNPs <- gen.info[gen.info[, 2] %in% sel, 1]
      sel2 <- dimnames(pvalues)[[1]] %in% SNPs & !is.na(pvalues[, 
        2]) & pvalues[, 2] <= sig
      SNPs.sel <- pvalues[sel2, ]
      pos.sel <- attr(x, "colSNPs")[sel2]
      out <- list(names = dimnames(SNPs.sel)[[1]], column = pos.sel)
     }
    out
}

