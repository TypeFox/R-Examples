`read.IPGpedigree` <-
function(pFile = "pedigree.dat",header = TRUE)
  {
    l <- readLines(pFile)
    l <- lapply(l,function(x)unlist(strsplit(x,split = " ")))
    nc <- max(sapply(l,length))
    l <- sapply(l,function(x)
                {
                  line <- character(nc)
                  line[1:length(x)] <- x
                  line
                })
    l <- t(l)
    if(header)
      {
        header <- l[1,]
        l <- l[-1,]
      }
    ped <- data.frame(l)
    if(length(header == nc))names(ped) <- header
    ped$DATE_BIRTH <- as.Date(ped$DATE_BIRTH,"%d-%m-%Y")
    return(ped)
  }

