insilico.digest.internal <- function(DNAseq, recognition_code, cut_site_5prime, cut_site_3prime)
{
    frag1 <- strsplit(DNAseq, split = recognition_code, fixed = FALSE, perl = FALSE)
    n <- length(frag1)
    for(i in 1:n){
       ni <- length(frag1[[i]])
       if(ni ==1){frag1[[i]] <- frag1[[i]]}
       if(ni > 1){
           for(y in 1:ni){
              if(y == 1){frag1[[i]][1] <- paste(frag1[[i]][1], cut_site_5prime, sep = "")}
              if(y > 1 & y < ni){frag1[[i]][y] <-  paste(cut_site_3prime, frag1[[i]][y], cut_site_5prime, sep = "")}
              if(y == ni){frag1[[i]][y] <- paste(cut_site_3prime, frag1[[i]][y], sep = "")}
           }
       }
    }
    if(n==1 & ni==1){frag2 <- frag1}
    frag2 <- unlist(frag1)
    return(frag2)
}
