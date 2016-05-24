`getZones` <-
function(ints, nofsets, ellipse = FALSE) {
    
    # s - sets; v - version; b - borders; x,y - coordinates
    borders <- read.csv(file.path(system.file("data", package="venn"), "borders.csv.gz"))
    
    # s - sets; v - version; i - intersection; b - border
    ib <- read.csv(file.path(system.file("data", package="venn"), "ib.csv.gz"))
    
    ints <- ints + 1
    
    if (nofsets < 4 | nofsets > 5) {
        ellipse <- FALSE
    }
    
    if (identical(ints, 1)) {
        ints <- seq(2^nofsets)[-1]
    }
    
    if (length(ints) > 1) {
        checkz <- logical(length(ints))
        names(checkz) <- ints
        checkz[1] <- TRUE
        
        result <- list()
        
        while(!all(checkz)) {
            checkz <- checkZone(as.numeric(names(checkz)[1]), ints, checkz, nofsets, ib, ellipse)
            
            result[[length(result) + 1]] <- as.numeric(names(checkz)[checkz])
            ints  <-  ints[!checkz]
            checkz <- checkz[!checkz]
            
            if (length(checkz) > 0) {
                checkz[1] <- TRUE
            }
        }
    }
    else {
        result = list(ints)
    }
    
    
    result <- lapply(result, function(x) {
        
        b <- ib$b[ib$s == nofsets & ib$v == as.numeric(ellipse) & ib$i %in% x]
        
        if (any(duplicated(b))) {
            b <- setdiff(b, b[duplicated(b)])
        }
        
        # print(ib[ib$s == nofsets & ib$v == as.numeric(ellipse) & ib$b %in% b, ])
        
        v2 <- borders[borders$s == nofsets & borders$v == as.numeric(ellipse) & borders$b == b[1], c("x", "y")]
        v2 <- v2[-nrow(v2), ] # get rid of the NAs, we want a complete polygon
        ends <- as.numeric(v2[nrow(v2), ])
        
        checkb <- logical(length(b))
        names(checkb) <- b
        checkb[1] <- TRUE
        
        # print(paste(b, collapse = ","))
        # print(head(v2))
        # print(tail(v2))
        
        counter <- 0
        
        while(!all(checkb)) {
            
            # do.call("rbind", lapply(... ???
            
            for (i in which(!checkb)) {
                
                temp <- borders[borders$s == nofsets & borders$v == as.numeric(ellipse) & borders$b == b[i], c("x", "y")]
                
                flag <- FALSE
                if (all(ends == as.numeric(temp[1, ]))) {
                    v2 <- rbind(v2, temp[-nrow(temp), ])
                    checkb[i] <- TRUE
                }
                else if (all(ends == as.numeric(temp[nrow(temp) - 1, ]))) {
                    temp <- temp[-nrow(temp), ]
                    v2 <- rbind(v2, temp[seq(nrow(temp), 1), ])
                    checkb[i] <- TRUE
                }
                
                if (checkb[i]) {
                    ends <- as.vector(v2[nrow(v2), ])
                }
            }
            
            counter <- counter + 1
            
            if (counter > length(checkb)) {
                # print(checkb)
                cat("\n")
                stop("Unknown error.\n\n", call. = FALSE)
            }
        }
        
        
        return(rbind(v2, rep(NA, 2)))
    })
    
    return(result)
}

    
