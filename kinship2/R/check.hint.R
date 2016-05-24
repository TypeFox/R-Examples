# Automatically generated from all.nw using noweb
check.hint <- function(hints, sex) {
    if (is.null(hints$order)) stop("Missing order component")
    if (!is.numeric(hints$order)) stop("Invalid order component")
    n <- length(sex)
    if (length(hints$order) != n) stop("Wrong length for order component")
    
    spouse <- hints$spouse
    if (is.null(spouse)) hints
    else {
        lspouse <- spouse[,1]
        rspouse <- spouse[,2]
        if (any(lspouse <1 | lspouse >n | rspouse <1 | rspouse > n))
            stop("Invalid spouse value")
        
        temp1 <- (sex[lspouse]== 'female' & sex[rspouse]=='male')
        temp2 <- (sex[rspouse]== 'female' & sex[lspouse]=='male')
        if (!all(temp1 | temp2))
            stop("A marriage is not male/female")
        
        hash <- n*pmax(lspouse, rspouse) + pmin(lspouse, rspouse)
        #Turn off this check for now - is set off if someone is married to two siblings
        #if (any(duplicated(hash))) stop("Duplicate marriage")

        # Break any loops: A left of B, B left of C, C left of A.
        #  Not yet done 
      }
    hints
  }
