myReadVCF <- function(vcffile){

vcf       <- .Call("myReadVCFC",vcffile)
matrix    <- vcf[[1]]
positions <- vcf[[2]]
rownames(matrix) <- vcf[[3]]

rm(vcf)
gc()

o_b_j_sub    <- list(matrix=matrix,reference=NaN,positions=positions)

return(o_b_j_sub)
}
