Digest <- function(sequence, enzyme = "trypsin", missed = 0, IAA = TRUE, 
                   N15 = FALSE, custom = list()) {


    ## determine cleavage sites according to enzyme specific rules

    seq_vector <- strsplit(sequence, split = "")[[1]]
    end_position <- length(seq_vector)
    
    if(enzyme == "trypsin") {                                
        if(seq_vector[end_position] == "K" | seq_vector[end_position] == "R") {
           seq_vector[end_position] <- "!"
           seq_string <- paste(seq_vector, collapse = "")
        } else seq_string <- sequence    
        seq_string <- gsub("KP", "!P", seq_string)          # prevent cleavage at K and R if followed by P  
        seq_string <- gsub("RP", "!P", seq_string)     
        seq_vector <- strsplit(seq_string, split = "")[[1]]
        stop <- grep("K|R", seq_vector)
        start <- stop + 1    
    }
    
    if(enzyme == "trypsin.strict") {                                
      if(seq_vector[end_position] == "K" | seq_vector[end_position] == "R") {
        seq_vector[end_position] <- "!"
        seq_string <- paste(seq_vector, collapse = "")
      } else seq_string <- sequence
      seq_vector <- strsplit(seq_string, split = "")[[1]]
      stop <- grep("K|R", seq_vector)
      start <- stop + 1    
    }

    if(enzyme == "pepsin") {
        if(seq_vector[end_position] == "F" | seq_vector[end_position] == "L" |
           seq_vector[end_position] == "W" | seq_vector[end_position] == "Y" | 
           seq_vector[end_position] == "A" | seq_vector[end_position] == "E" | 
           seq_vector[end_position] == "Q") {
               seq_vector[end_position] <- "!"
        } 
        stop <- grep("F|L|W|Y|A|E|Q", seq_vector)       
        start <- stop + 1
    }

   
    ## error checking

    if(enzyme != "trypsin" & enzyme != "trypsin.strict" & enzyme != "pepsin") stop("undefined enzyme, defined enzymes are trypsin, trypsin.strict, and pepsin")                    
    if(length(stop) == 0) warning("sequence does not contain cleavage sites")
    if(missed > length(stop)) stop("number of specified missed cleavages is greater than the maximum possible")


    ## cleave sequence

    cleave <- function(sequence, start, stop, misses) {
        peptide <- substring(sequence, start, stop)
        mc <- rep(misses, times = length(peptide))
        result <- data.frame(peptide, start, stop, mc, stringsAsFactors = FALSE)
        return(result) 
    }

    start <- c(1, start)                           # peptides if 0 missed cleavages
    stop <- c(stop, end_position)
    results <- cleave(sequence, start, stop, 0)

    if(missed > 0) {                               # peptides if missed cleavages > 0
        for(i in 1:missed) {
            start_tmp <- start[1:(length(start) - i)]
            stop_tmp <- stop[(1 + i):length(stop)]
            peptide <- cleave(sequence, start_tmp, stop_tmp, i)
            results <- rbind(results, peptide) 
        } 
    }


    ## calculate m/z values for resulting peptides

    C <- 12.0000000                                         # carbon-12
    H <- 1.0078250321                                       # hydrogen-1
    O <- 15.9949146221                                      # oxygen-16
    S <- 31.97207069                                        # sulfer-32    
    N <- ifelse(N15 == TRUE, 15.0001088984, 14.0030740052)  # nitrogen-14 or -15
    proton <- 1.007276466

    residueMass <- function(residue) {
                                   
        if(residue == "A") mass = C*3 + H*5 + N + O         # alanine
        if(residue == "R") mass = C*6 + H*12 + N*4 + O      # arginine
        if(residue == "N") mass = C*4 + H*6 + N*2 + O*2     # asparagine
        if(residue == "D") mass = C*4 + H*5 + N + O*3       # aspartic acid
        if(residue == "E") mass = C*5 + H*7 + N + O*3       # glutamic acid
        if(residue == "Q") mass = C*5 + H*8 + N*2 + O*2     # glutamine
        if(residue == "G") mass = C*2 + H*3 + N + O         # glycine
        if(residue == "H") mass = C*6 + H*7 + N*3 + O       # histidine
        if(residue == "I") mass = C*6 + H*11 + N + O        # isoleucine
        if(residue == "L") mass = C*6 + H*11 + N + O        # leucine
        if(residue == "K") mass = C*6 + H*12 + N*2 + O      # lysine
        if(residue == "M") mass = C*5 + H*9 + N + O + S     # methionine
        if(residue == "F") mass = C*9 + H*9 + N + O         # phenylalanine
        if(residue == "P") mass = C*5 + H*7 + N + O         # proline
        if(residue == "S") mass = C*3 + H*5 + N + O*2       # serine
        if(residue == "T") mass = C*4 + H*7 + N + O*2       # threonine
        if(residue == "W") mass = C*11 + H*10 + N*2 + O     # tryptophan
        if(residue == "Y") mass = C*9 + H*9 + N + O*2       # tyrosine
        if(residue == "V") mass = C*5 + H*9 + N + O         # valine
         
        if(residue == "C" & IAA == FALSE) mass = C*3 + H*5 + N + O + S   # cysteine, iodoacetylated cysteine
        if(residue == "C" & IAA == TRUE)
            mass <- ifelse(N15 == FALSE, C*5 + H*8 + N*2 + O*2 + S,      # do not include nitrogen-15 in IAA
                           C*5 + H*8 + N + 14.0030740052 + O*2 + S)
       
        if(length(custom) != 0)                             
            for(i in 1:length(custom$code)) if(residue == custom$code[i]) mass = custom$mass[i]
 
        return(mass)
    }


    ## calculate m/z values of peptides

    mz <- vector("list", length = nrow(results))                  
    for(i in 1:nrow(results)) {                               
        peptide_vector <- strsplit(results$peptide[i], split = "")[[1]]
        peptide_mass <- sum(sapply(peptide_vector, residueMass))    
        mz[[i]] <- round((peptide_mass + H*2 + O + (c(1, 2, 3) * proton)) / c(1, 2, 3), digits = 3)
    }
    mz <- as.data.frame(do.call("rbind", mz))
    names(mz) <- c("mz1", "mz2", "mz3")
    results <- cbind(results, mz)

    return(results)

}
