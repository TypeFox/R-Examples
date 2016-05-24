### Convert colors.

get.color <- function(u.codon, color = .CF.PT$color){
  if(length(u.codon) <= 4){
    ### Cases of less than or equal to 4 synonymous codon are all fine, but
    ### the colors should be matched the position. e.g. C = c("TGC", "TGT")
    ### should be colors 2 and 4 rather than colors 1 and 2.
    color <- color[c("A", "C", "G", "T") %in% gsub("..(.)", "\\1", u.codon)]
  } else{
    if(all(u.codon %in% .CF.GV$synonymous.codon$L)){
      ### L = c("CTA", "CTC", "CTG", "CTT", "TTA", "TTG")
      ### Do nothing.
    } else if(all(u.codon %in% .CF.GV$synonymous.codon$R)){
      ### R = c("AGA", "AGG", "CGA", "CGC", "CGG", "CGT")
      color <- color[c(5:6, 1:4)]
    } else if(all(u.codon %in% .CF.GV$synonymous.codon$S)){
      ### S = c("AGC", "AGT", "TCA", "TCC", "TCG", "TCT")
      ### This is for un-split case. The split case is ok in 4 codons case.
      color <- color[c(5:6, 1:4)]
    } else{
      stop("Not in right order nor right unique codon.")
    }
  }

  color
} # End of get.color().

