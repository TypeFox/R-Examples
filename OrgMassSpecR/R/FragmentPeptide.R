FragmentPeptide <- function(sequence, fragments = "by", IAA = TRUE, N15 = FALSE, custom = list()) {

    results_list <- vector("list")
    for(sequence_number in 1:length(sequence)) {

        peptide_vector <- strsplit(sequence[sequence_number], split = "")[[1]]
        peptide_length <- length(peptide_vector)
        if(peptide_length < 2) stop("sequence must contain two or more residues")
        
        C <- 12                                                  # carbon-12
        H <- 1.0078250321                                        # hydrogen-1
        O <- 15.9949146221                                       # oxygen-16
        S <- 31.97207069                                         # sulfer-32    
        N <- ifelse(N15 == TRUE, 15.0001088984, 14.0030740052)   # nitrogen-14 or nitrogen -15
        proton <- 1.007276466
        electron <- 0.00054857990943   # currently unused   

        residueMass <- function(residue) {

            if(residue == "A") mass = C*3 + H*5 + N + O          # alanine
            if(residue == "R") mass = C*6 + H*12 + N*4 + O       # arginine
            if(residue == "N") mass = C*4 + H*6 + N*2 + O*2      # asparagine
            if(residue == "D") mass = C*4 + H*5 + N + O*3        # aspartic acid
            if(residue == "E") mass = C*5 + H*7 + N + O*3        # glutamic acid
            if(residue == "Q") mass = C*5 + H*8 + N*2 + O*2      # glutamine
            if(residue == "G") mass = C*2 + H*3 + N + O          # glycine
            if(residue == "H") mass = C*6 + H*7 + N*3 + O        # histidine
            if(residue == "I") mass = C*6 + H*11 + N + O         # isoleucine
            if(residue == "L") mass = C*6 + H*11 + N + O         # leucine
            if(residue == "K") mass = C*6 + H*12 + N*2 + O       # lysine
            if(residue == "M") mass = C*5 + H*9 + N + O + S      # methionine
            if(residue == "F") mass = C*9 + H*9 + N + O          # phenylalanine
            if(residue == "P") mass = C*5 + H*7 + N + O          # proline
            if(residue == "S") mass = C*3 + H*5 + N + O*2        # serine
            if(residue == "T") mass = C*4 + H*7 + N + O*2        # threonine
            if(residue == "W") mass = C*11 + H*10 + N*2 + O      # tryptophan
            if(residue == "Y") mass = C*9 + H*9 + N + O*2        # tyrosine
            if(residue == "V") mass = C*5 + H*9 + N + O          # valine
             
            if(residue == "C" & IAA == FALSE) mass = C*3 + H*5 + N + O + S    # cysteine, iodoacetylated cysteine
            if(residue == "C" & IAA == TRUE)
                mass <- ifelse(N15 == FALSE, C*5 + H*8 + N*2 + O*2 + S,       # do not include nitrogen-15 in IAA
                               C*5 + H*8 + N + 14.0030740052 + O*2 + S)
            
            # custom residues
            if(length(custom) != 0)
                for(i in 1:length(custom$code)) if(residue == custom$code[i]) mass = custom$mass[i]

            return(mass)

        }
        
        masses <- sapply(peptide_vector, residueMass)

        pm <- sum(masses)
        p1 <- round(pm + H*2 + O + proton, digits = 3)               # precursor ion m/z, z = +1
        p2 <- round((pm + H*2 + O + (2 * proton)) / 2, digits = 3)   # precursor ion m/z, z = +2
        p3 <- round((pm + H*2 + O + (3 * proton)) / 3, digits = 3)   # precursor ion m/z, z = +3

        if(fragments == "by") {
        
            b1 <- vector(mode = "numeric", length = 0)
            b2 <- vector(mode = "numeric", length = 0)
            bs <- vector(mode = "character",length = 0)
            bi <- vector(mode = "integer", length = 0)
            y1 <- vector(mode = "numeric", length = 0)
            y2 <- vector(mode = "numeric", length = 0)
            ys <- vector(mode = "character",length = 0)
            yi <- vector(mode = "integer", length = 0)

            for(i in 1:(peptide_length - 1)) {                       # b-ions
                mass  <- sum(masses[1:i])
                b1[i] <- round(mass + proton, digits = 3)            # +1 charge                         
                b2[i] <- round((b1[i] + proton) / 2, digits = 3)     # +2 charge
                bs[i] <- paste(peptide_vector[1:i], collapse = "")   # b-ion sequence   
                bi[i] <- i                                           # b-ion position
            }
            
            for(j in 2:peptide_length) {                                            # y-ions
                mass    <- sum(masses[j:peptide_length])                  
                y1[j-1] <- round(mass + H*2 + O + proton, digits = 3)               # +1 charge
                y2[j-1] <- round((y1[j-1] + proton) / 2, digits = 3)                # +2 charge
                ys[j-1] <- paste(peptide_vector[j:peptide_length], collapse = "")   # sequence
                yi[j-1] <- peptide_length - j + 1                                   # position
            }

            ms1seq <- rep(sequence[sequence_number], times = ((2*(length(bi))) + (2*(length(yi)))))   # peptide seq column
             
            ms1z1 <- rep(p1, times = ((2*(length(bi))) + (2*(length(yi)))))   # precursor m/z columns
            ms1z2 <- rep(p2, times = ((2*(length(bi))) + (2*(length(yi)))))
            ms1z3 <- rep(p3, times = ((2*(length(bi))) + (2*(length(yi)))))
                       
            ms2seq <- c(rep(bs, times=2), rep(ys, times=2))   # MS2 ion sequence column
                       
            b1.type <- paste("[b", bi, "]1+", sep = "")   # MS2 ion type column
            b2.type <- paste("[b", bi, "]2+", sep = "")
            y1.type <- paste("[y", yi, "]1+", sep = "")
            y2.type <- paste("[y", yi, "]2+", sep = "")
            ms2type <- c(b1.type, b2.type, y1.type, y2.type)
                
            ms2mz <- c(b1, b2, y1, y2)   # MS2 ion m/z column
        
        }

        if(fragments == "cz") {

            c1 <- vector(mode = "numeric", length = 0)
            c2 <- vector(mode = "numeric", length = 0)
            cs <- vector(mode = "character",length = 0)
            ci <- vector(mode = "integer", length = 0)
            z1 <- vector(mode = "numeric", length = 0)
            z2 <- vector(mode = "numeric", length = 0)
            zs <- vector(mode = "character",length = 0)
            zi <- vector(mode = "integer", length = 0)
                
            for(i in 1:(peptide_length - 1)) {                        # c-ions
                mass  <- sum(masses[1:i])
                c1[i] <- round(mass + 3*H + N + proton, digits = 3)   # +1 charge                         
                c2[i] <- round((c1[i] + proton) / 2, digits = 3)      # +2 charge
                cs[i] <- paste(peptide_vector[1:i], collapse = "")    # sequence   
                ci[i] <- i                                            # position
            }
            
            for(j in 2:peptide_length) {                                            # z-ions
                mass    <- sum(masses[j:peptide_length])                  
                z1[j-1] <- round(mass + O - N, digits = 3)                          # +1 charge
                z2[j-1] <- round((z1[j-1] + proton) / 2, digits = 3)                # +2 charge
                zs[j-1] <- paste(peptide_vector[j:peptide_length], collapse = "")   # sequence
                zi[j-1] <- peptide_length - j + 1                                   # position
            }

            ms1seq <- rep(sequence[sequence_number], times = ((2*(length(ci))) + (2*(length(zi)))))   # peptide seq column
             
            ms1z1 <- rep(p1, times = ((2*(length(ci))) + (2*(length(zi)))))   # precursor m/z columns
            ms1z2 <- rep(p2, times = ((2*(length(ci))) + (2*(length(zi)))))
            ms1z3 <- rep(p3, times = ((2*(length(ci))) + (2*(length(zi)))))
                
            ms2seq <- c(rep(cs, times=2), rep(zs, times=2))   # MS2 ion sequence column
               
            c1.type <- paste("[c", ci, "]1+", sep = "")   # MS2 ion type column
            c2.type <- paste("[c", ci, "]2+", sep = "")
            z1.type <- paste("[z", zi, "]1+", sep = "")
            z2.type <- paste("[z", zi, "]2+", sep = "")
            ms2type <- c(c1.type, c2.type, z1.type, z2.type)
                
            ms2mz <- c(c1, c2, z1, z2)   # MS2 ion m/z column

        }
        
        results_list[[sequence_number]] <- data.frame(ms1seq, ms1z1, ms1z2, ms1z3, ms2seq, ms2type, ms2mz)

    }
                  
    return(as.data.frame(do.call("rbind", results_list)))
                
}
