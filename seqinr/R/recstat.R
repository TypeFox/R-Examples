##
# This function counts the number of each triplet in a sliding window,
# create a contingency table with triplet composition and then computes
# a correspondence analysis on this table
##
#v.18.08.2011
recstat <- function(seq, sizewin = 90, shift = 30, seqname = "no name")
{
    if (is.character(seq) == FALSE)
    { # class character if only one sequence, class list if more
        print("This file has more than one sequence or class is not character.")
        return()
    }
    if ((shift%%3) != 0)
    { # test if shift give the same reading frame
        print("The windows are not in the same reading frame, please change shift value to a multiple of 3.")
        return()
    }
    if (sizewin%%3 != 0) 
    { 
        print("The length of the window is not a multiple of 3, please change sizewin value.")
        return()
    }
    seqsize <- length(seq) # give the number of 1-mer in the sequence
    v1 <- seq(from = 1, to = seqsize, by = shift) # start vector of window in reading frame 1
    v1 <- v1[1:(which((v1 + sizewin) > seqsize)[1])] # suppression if more than one incomplete window
    v2 <- seq(from = 2, to = seqsize, by = shift)
    v2 <- v2[1:(which((v2 + sizewin)>seqsize)[1])]
    v3 <- seq(from = 3, to = seqsize, by = shift)
    v3 <- v3[1:(which((v3 + sizewin) > seqsize)[1])]
    vdep <- c(v1, v2, v3) # start vector in the 3 reading frames of direct/reverse strand
    vind <- c(rep(1, length(v1)), rep(2, length(v2)), rep(3, length(v3))) # index vector of reading frame for each window
    ##
    ##direct strand##
    ##
    cseq <- c2s(seq)
    vstopd <- c(words.pos("taa", cseq), words.pos("tag", cseq), words.pos("tga", cseq)) # vector of stop codons positions in direct strand
    vinitd <- c(words.pos("atg", cseq)) # vector of start codons positions in direct strand
    resd <- lapply(1:length(vdep), function(x)
    { # calculation on 3 reading frames of direct strand
        seq_tmp <- seq[(vdep[x]):(vdep[x] + sizewin - 1)] # temporary window
        count(seq_tmp, wordsize = 3, start = 0, by = 3) # counting of triplets
    })
    ##
    ##reverse strand##
    ##
    seq_reverse <- rev(comp(seq, ambiguous = TRUE)) # creation of reverse strand
    cseq_reverse <- c2s(seq_reverse)
    vstopr <- c(words.pos("taa", cseq_reverse), words.pos("tag", cseq_reverse), words.pos("tga", cseq_reverse)) # vector of stop codons positions in reverse strand
    vinitr <- c(words.pos("atg", cseq_reverse)) # vector of init codons positions in reverse strand    
    resr <- lapply(1:length(vdep), function(x)
    { # calculation on 3 reading frames of reverse strand
        seq_tmp <- seq_reverse[(vdep[x]):(vdep[x] + sizewin - 1)] # temporary window
        count(seq_tmp, wordsize = 3, start = 0, by = 3) # counting of triplets
    })
    resd <- matrix(unlist(resd), byrow = TRUE, ncol = 64) # conversion vector to contingency table
    resr <- matrix(unlist(resr), byrow = TRUE, ncol = 64)
    ##
    ##CA##
    ##
    resd.coa <- dudi.coa(resd, scannf = FALSE, nf = 4) # CA on direct strand
    resr.coa <- dudi.coa(resr, scannf = FALSE, nf = 4) # CA on reverse strand    
    rec <- list(seq, sizewin, shift, seqsize, seqname, vdep, vind, vstopd, vstopr, vinitd, vinitr, resd, resr, resd.coa, resr.coa)    
    return(rec)
}
