write.pop.file <- function(data,
                           file="",
                           digits=2,
                           description="Data from R"
                           )
  {
    which <- sapply(data, is.genotype)
    if(!all(which)) warning("Data contianed non-genotype variables.",
                            " These have been ignored: ",
                            paste(colnames(data)[!which]) )

    data <- data[,which]

    # convert allele names into two or three digit numbers
    for( col in names(data) )
      {
        # first convert to numbers
        a1 <- as.numeric(factor(allele(data[[col]],1)))
        a2 <- as.numeric(factor(allele(data[[col]],2)))

        # convert NA to 0
        a1[is.na(a1)] <- 0
        a2[is.na(a2)] <- 0
        
        # now format to have correct # of digits
        a1 <- formatC( a1, width=digits, flag="0")
        a2 <- formatC( a2, width=digits, flag="0")
        
        # now paste back together
        data[[col]] <- paste(a1,a2,sep="")
      }

    if(file=="")
      f <- stdout()
    else
      f <- file(file,"w")

    # header line
    cat(description, file=f)
    cat("\n", file=f)

    # marker names
    cat(colnames(data),sep=" ", file=f)
    cat("\n", file=f)

    # group token
    cat("POP", file=f)
    cat("\n", file=f)

    # write allele data.  First token is row id, followed by a comma
    # markers are separated by space
    rownames(data) <- paste(rownames(data),",", sep="")
    write.table( data, file=f, sep=" ", quote=FALSE, col.names=F)

    if(file!="")
      close(f)
    
  }

