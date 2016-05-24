### For adjusting focal codons.
###
### Case 1: Only 2 codons, differ in one codon.
###   - A: C1 C2 C3 C4
###   - B: C1 C4 C3 C2    where C2 is the new focal codon.
###   In this case, B needs to be swapped C4 and C2.
###
### Case 2: Differ in 2 codons
###   - A: C1 C2 C3 C6 C5 C4    where C4 is the new focal codon.
###   - B: C1 C6 C3 C4 C5 C2    where C2 is the new focal codon.
###   In this case, C6 is the common codon, B needs to be swapped to as A.
###   First, swap C6 and C2, then swap C4 and C6.

adjust.focal.codon <- function(y, label, label.focal, y.ci = NULL){
  label.aa <- gsub("(.*)\\.(.*)", "\\1", label)
  label.codon <- gsub("(.*)\\.(.*)", "\\2", label)

  focal.aa <- gsub("(.*)\\.(.*)", "\\1", label.focal)
  focal.codon <- gsub("(.*)\\.(.*)", "\\2", label.focal)

  aa.names <- unique(label.aa)
  for(i.aa in aa.names){
    id <- label.codon[label.aa == i.aa] != focal.codon[focal.aa == i.aa]

    if(sum(id) == 1){    ### Differ in one codon.
      ### Swap y.
      tmp.y <- y[label.aa == i.aa]
      max.y <- tmp.y[which(id)]
      tmp.y <- tmp.y - max.y 
      tmp.y[which(id)] <- -max.y
      y[label.aa == i.aa] <- tmp.y

      ### Swap y.ci if any.
      if(!is.null(y.ci)){
        tmp.y.ci <- matrix(y.ci[label.aa == i.aa,], ncol = 2)
        max.y.ci <- tmp.y.ci[which(id),]
        tmp.y.ci <- tmp.y.ci - max.y
        tmp.y.ci[which(id),] <- -max.y.ci
        y.ci[label.aa == i.aa,] <- tmp.y.ci
      }
    } else if(sum(id) == 2){    ### Differ in two codons.
      ### Find intersection first then find common and swapping codons.
      diff.codon <- label.codon[label.aa == i.aa][id]
      common.codon <- diff.codon[diff.codon %in% focal.codon[focal.aa == i.aa]]
      id.common <- which(label.codon[label.aa == i.aa] == common.codon)
      id.swap <- which(label.codon[label.aa == i.aa] ==
                       diff.codon[diff.codon != common.codon])

      ### Swap y.
      tmp.y <- y[label.aa == i.aa]
      common.y <- tmp.y[id.common]
      tmp.y <- tmp.y + common.y     ### Add the old value back.
      tmp.y[id.common] <- -common.y ### Replace back to the old value.
      swap.y <- tmp.y[id.swap]      ### Get a new value.
      tmp.y <- tmp.y - swap.y       ### Substract the new value.
      tmp.y[id.swap] <- -swap.y     ### Replace the negative new value.
      y[label.aa == i.aa] <- tmp.y

      ### Swap y.ci if any.
      if(!is.null(y.ci)){
        tmp.y.ci <- matrix(y.ci[label.aa == i.aa,], ncol = 2)
        common.y.ci <- tmp.y.ci[which(id),]
        tmp.y.ci <- tmp.y.ci + common.y       ### Add the old value back.
        tmp.y.ci[id.common,] <- -common.y.ci  ### Replace back to the old value.
        swap.y.ci <- tmp.y.ci[id.swap,]       ### Get a new value.
        tmp.y.ci <- tmp.y.ci - swap.y.ci      ### Substract the new value.
        tmp.y.ci[id.swap,] <- -swap.y.ci
        y.ci[label.aa == i.aa,] <- tmp.y.ci
      }
    }
  }

  ret <- list(y = y, y.ci = y.ci)
  ret
} # End of adjust.focal.codon().

