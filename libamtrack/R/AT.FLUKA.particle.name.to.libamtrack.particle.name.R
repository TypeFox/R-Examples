AT.FLUKA.particle.name.to.libamtrack.particle.name <- function( FLUKA.particle.names)
{
    libamtrack.particle.name    <- character(length(FLUKA.particle.names))
    
    # Special cases: P, D, T
    ii                           <- FLUKA.particle.names == "P"
    libamtrack.particle.name[ii] <- rep("1H", sum(ii))
    ii                           <- FLUKA.particle.names == "D"
    libamtrack.particle.name[ii] <- rep("2H", sum(ii))
    ii                           <- FLUKA.particle.names == "T"
    libamtrack.particle.name[ii] <- rep("3H", sum(ii))

    # Special cases II: '*' all isotope of an element (<= O), replace by A = 2Z
    gen.isotope                  <- c("H*", "HE*", "LI*", "BE*", "B*",  "C*",  "N*",  "O*")
    rep.isotope                  <- c("H1", "HE4", "LI6", "BE8", "B10", "C12", "N14", "O16")
    idx.gen.isotope              <- grep('*', FLUKA.particle.names, fixed = TRUE)
    if(length(idx.gen.isotope) > 0){
        jj                           <- idx.gen.isotope %in% 1:length(FLUKA.particle.names)
        for (i in 1:length(gen.isotope)){
            # DEBUG i <- 1 
            ii                           <- FLUKA.particle.names[jj] == gen.isotope[i]
            FLUKA.particle.names[jj][ii] <- rep(rep.isotope[i], sum(ii))
        }
    }

    # Regular cases
    ii                           <- !(FLUKA.particle.names %in% c("P", "D", "T"))

    split.pos                    <- regexpr("[[:digit:]]", FLUKA.particle.names[ii])
    names                        <- substring(FLUKA.particle.names[ii], 1, split.pos-1)
    numbers                      <- substring(FLUKA.particle.names[ii], split.pos, nchar(FLUKA.particle.names[ii]))
    
    jj                           <- nchar(names) == 2
    names[jj]                    <- paste( substring( names[jj], 1, 1), 
                                           tolower(substring(names[jj], 2, 2)), 
                                           sep = "")
    libamtrack.particle.name[ii] <- paste( numbers, names, sep = "")
 
    return(libamtrack.particle.name)
}