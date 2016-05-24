`[.setupSNP` <-
function(x,i,j,...){

    out<-NextMethod("[")

    if (!is.null(dim(out))){
      k<- match(attr(x, "label.SNPs"), names(out)) # nuevas columnas con snps
      k<-k[!is.na(k)]
      ik<- match(names(out), attr(x, "label.SNPs")) # nuevas columnas con snps
      ik<-ik[!is.na(ik)]

      snps<- attr(x, "label.SNPs")[ik]
#      for (l in snps) attr(out[,l],"allele.names")<- attr(x[,l],"allele.names")
            
      attr(out, "colSNPs")        <- sort(k)
      attr(out, "label.SNPs")     <- attr(x, "label.SNPs")[ik]
      attr(out, "gen.info")       <- attr(x, "gen.info")[ik,]
    }
    out }

