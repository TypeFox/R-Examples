"permutation" <-
function (sequence, modele = "base", frame = 0, replace = FALSE, 
    prot = FALSE, numcode = 1, ucoweight = NULL) 
{
    if (modele == "base") { #modele 1 : permute toutes les bases entre elles
        new <- sample(sequence, replace = replace)
    }
    else if (modele == "position") { #modele 2 : permute les bases en conservant position
        CDSseq <- sequence[(frame + 1):length(sequence)]
        if (length(CDSseq)%%3 != 0) {
            stop("sequence could not be subdivided in codons")
        }
        else {
            seqI <- CDSseq[seq(1, length(CDSseq), by = 3)]
            seqII <- CDSseq[seq(2, length(CDSseq), by = 3)]
            seqIII <- CDSseq[seq(3, length(CDSseq), by = 3)]
            newCDSseq <- as.vector(t(cbind(sample(seqI, replace = replace), 
                sample(seqII, replace = replace), sample(seqIII, 
                  replace = replace))))
        }
        new <- c(sequence[0:frame], newCDSseq)
    }
    else if (modele == "codon") { #modele 3 : permute tous les codons entre eux
        split <- splitseq(sequence, frame) #split in codons (frame=starting position of CDS)
        if (prot == FALSE) {
            new <- c(sequence[0:frame], s2c(c2s(sample(split, 
                replace = replace))), tail(sequence, (length(sequence) - 
                frame)%%3))
        }
        else { #"prot==TRUE" = garder Met et STOP
            l <- length(split)
            new <- c(sequence[0:frame], s2c(c2s(c(split[1], sample(split[2:(l - 
                1)], replace = replace), split[l]))), tail(sequence, 
                (length(sequence) - frame)%%3))
        }
    }
    else if (modele == "syncodon") { #modele 4 : permute/remplace les codons synonymes
        CDSseq <- sequence[(frame + 1):length(sequence)]
        newCDSseq <- synsequence(CDSseq, numcode = numcode, ucoweight = ucoweight)
        new <- c(sequence[0:frame], newCDSseq, tail(sequence, 
            (length(sequence) - frame)%%3))
    }
    return(new)
}
