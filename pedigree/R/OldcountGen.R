OldcountGen <- function(ped)
  {
    findParents <- function(i)
      {
        id <- ped$ID[i]
        s <- match(ped$SIRE[i],ped$ID)
        if(!is.na(s))
          {
            if(is.na(gen[s]))findParents(s)
            genS <- gen[s]+1
          }
        else genS <- 0
                
        d <- match(ped$DAM[i],ped$ID)
        if(!is.na(d))
          {
            if(is.na(gen[d]))findParents(d)
            genD <- gen[d]+1
          }
        else genD <- 0

        gen[i] <<- max(genD,genS)
      }
    if(!is.data.frame(ped))stop("ped should be data.frame")
    if(ncol(ped)!=3)warning("assuming ID SIRE DAM organisation of pedigree")
    ped <- ped[,1:3]
    ped <- ped[orderPed(ped),]
    gen <- rep(NA,nrow(ped))
    for(i in 1:nrow(ped))findParents(i)
    return(gen)
  }
