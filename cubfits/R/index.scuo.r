# This is originally coded by Drew Schmidt and adjusted by WCC.

# Calculate the Synonymous Codon Usage Order (SCUO) index for each gene.
# Used as a substitute for Phi in the without phi case
# See: http://www.tandfonline.com/doi/abs/10.1080/03081070500502967
calc_scuo_values <- function(codon.counts)
{
  # Get the gene ids
  gene.names <- unique(codon.counts$ORF)
  codons <- codon.counts[, 3:8]
 
  # Compute SCUO
  scuo <- unlist(
    lapply(gene.names,
      function(i) # i'th gene
      {
        gene <- which(codon.counts$ORF == i)
        codons <- codons[gene, ]
       
        # Sum and Prob of codons within an AA within a gene
        sums <- rowSums(codons, na.rm=T)
        p_ij <- sweep(codons, 1, sums, "/")
       
        # Shannon's Entropy for each AA within each gene H_i
        H_i <- -rowSums(p_ij * log(p_ij), na.rm=T)
       
        Hmax_i <- -log(1/rowSums(!is.na(codons)))
       
        # SCUO for each amino acid
        O_i <- (Hmax_i - H_i)/Hmax_i
       
        # Composition ratio of i'th amino acid
        denom <- sum(sums)
        F_i <- rowSums(codons/denom, na.rm=T)
       
        # Average SCUO for each gene
        O <- sum(F_i*O_i)
       
        return( O )
      }
    )
  )
 
 
  # Return
  output.df <- as.data.frame(gene.names, stringsAsFactors = FALSE)
 
  output.df$scuo <- scuo
  names(output.df) <- c("ID", "SCUO")
 
  return( output.df )
}

