cluster.strains <-
function(variants) {

  # If the quality columns are still in the SNPs, remove them.
  variants = strip.quality.columns(variants)

  # Separate the SNP header from the numeric SNP data.
  col.names = colnames(variants)
  snp.hdr = variants[,1:5]
  variants = as.matrix(variants[,-1:-5])
  
  # If we don't have numeric SNPs, alert the user.
  if(!is.numeric(variants)) {
    stop(paste("The SNPs in cluster.variants must be numeric.",
         "Run convert.variants.to.numeric() first"))
  } # if(!is.numeric(variants))

  # Create a matrix that indictes if two SNPs are equal.
  sim = matrix(0, ncol(variants), ncol(variants), dimnames = list(colnames(variants), 
               colnames(variants)))
  for(i in 1:ncol(variants)) {
    sim[i,] = colMeans(variants[,i] == variants, na.rm = T)
  } # for(i)

  # Cluster.
  cl = hclust(as.dist(1.0 - sim), method = "average")
  variants = variants[,cl$order]
  col.names[-1:-5] = col.names[-1:-5][cl$order]

  variants = data.frame(snp.hdr, variants, stringsAsFactors = F)
  colnames(variants) = col.names
  return(variants)
}
