makeHaploviewInputFile <-
  function(famid, patid, fid, mid, sex, aff, geno, marker.name,
           marker.position, haploview.pedfile, haploview.infofile) {
    
    aff[is.na(aff)] <- 0
    if (any(is.na(match(aff, c(0, 1, 2))))) {
      stop(paste("Affected status should contain only ", 
                 "0( = unkown), 1(=unaff) or 2(=aff).", sep = ""))
    }
    genot <- alleleRto2(geno)
    genot <- replace(genot, is.na(genot), 0)
    ped.data <- data.frame(famid, patid, fid, mid, sex, aff, genot)
    write.table(ped.data, 
                file = haploview.pedfile,
                quote = FALSE, 
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t")
    rm(ped.data)
    info.data <- data.frame(marker.name, marker.position)
    write.table(info.data, 
                file = haploview.infofile,
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE, 
                sep = "\t")
  }
