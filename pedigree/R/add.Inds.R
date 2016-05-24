`add.Inds` <-
function(ped)
  {
    if(ncol(ped)<3)stop("pedigree should have at least 3 columns")
    head <- names(ped)
    ndams <- match(ped[,2],ped[,1])
    ndams <- as.character(unique(ped[is.na(ndams),2]))
    ndams <- ndams[!is.na(ndams)]
    nsires <- match(ped[,3],ped[,1])
    nsires <- as.character(unique(ped[is.na(nsires),3]))
    nsires <- nsires[!is.na(nsires)]
    nped <- data.frame(matrix(NA,nrow = length(ndams) + length(nsires),
                              ncol = ncol(ped)))
    names(nped) <- names(ped)
    nped[,1] <- c(ndams,nsires)
    ped <- rbind(nped,ped)
    names(ped) <- head
    return(ped)
  }

