summary.RSNPset<-function(object, verbose=TRUE, ...) {
  if(verbose){
    cat("- Efficient score statistics based on", attr(object,'n'), "samples.\n")
    cat("- SNP sets range in size from ", min(attr(object, 'mAna')), " to ", max(attr(object, 'mAna')), ".\n" ,sep="")
    cat("-", attr(object, "KSub")-attr(object, "KAna"),"SNP sets were not included in the analysis\n")
    s<-attr(object,'mSub')
    a<-attr(object, 'mAna')
    cat("-",sum( s[names(a)] != a ), "SNP sets contained SNPs that were not included in the analysis.\n")
    B<-attr(object, 'B')
    if(B>0) {
      cat("-", attr(object, 'B'), attr(object, 'r.method'), "replicates were computed.\n")
      if(attr(object, 'r.method') == "permutation") {
        cat("- ret.rank = ", attr(object, 'ret.rank'), " : The ranks of the permutation variance matrices were ", ifelse( attr(object, 'ret.rank'), "","not "), "returned.\n",sep="")
        cat("- v.permute = ",attr(object, 'v.permute')," : Variance was ", ifelse( attr(object, 'v.permute'), "","not "), "recomputed for each permutation replicate.\n",sep="")
      }
    }
    if(any(!is.na(attr(object, 'pinv.check')))) cat("- pinv.tol =",attr(object, 'pinv.tol'),"\n") 
  }
  return(attr(object, 'pinv.check'))
}
