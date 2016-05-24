# Version: 28-06-2013, DF

# Function tested on 28-06-2013, DF

  createGroups <- function(g,desOrder){
    reG <- g
    if(!is.numeric(g))g <- match(g,names(table(g)))
    desOrder <- desOrder
    curClass <- 1
    for(i in 1:length(desOrder))
    {
      tempPos <- g==desOrder[i]
      reG[tempPos] <- curClass
      curClass <- curClass + 1
    }
    as.numeric(reG)
  }  # end of function createGroups

#-----------------------------------------------------------------------------------------------------------------
