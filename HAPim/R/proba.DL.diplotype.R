`proba.DL.diplotype` <-
function(DL.chrom1,DL.chrom2){


 DL.QQ    = DL.chrom1*DL.chrom2
 DL.qq    = 1 - DL.chrom1 - DL.chrom2 + DL.QQ
 DL.Qq    =  DL.chrom1 - DL.QQ
 DL.qQ    = DL.chrom2 - DL.QQ
 DL.d     = cbind(DL.QQ,DL.qq,DL.Qq,DL.qQ)

 DL.d


}

