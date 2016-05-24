msc.level.ind <- function (msLevel, pId, addExtrema = TRUE) 
{


    i <- msLevel$partition == pId
    if(addExtrema){ 
      i[msLevel$min[pId]] <- TRUE
      i[msLevel$max[pId]] <- TRUE
    }
    i
}

