
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2011-12-01 14:26:10 +0100 (Thu, 01 Dec 2011) $
# $Rev: 200 $
# ----------------------------------------------------------------

AverageWuxDataFrame <- function(x, INDEX, fun = "mean") {
##
## Args:
##   x: WUX data frame from models2wux dataframe
##   INDEX: factor (data frame column name) over which to aggregate
##   fun: aggregating function (default is arithmetic mean)
##
## History:
##   2011-09-29 | orig code (thm)
  
  ## if any '.run' is passed in INDEX, we have to cut out trailing run
  ## identification from acronyms (eg "cccma_cgcm3_1-r3" will become 
  ## "cccma_cgcm3_1")
  if (length(grep(".run", INDEX)) > 0){
    acronyms <- x[["acronym"]]
    acronym.norun <- sub("-r[0-9]+$", "", acronyms)
    x[["acronym"]] <- acronym.norun
  }
  
  ## factors expanding the dataframe
  id = c("acronym","subreg", "season",  "institute", "gcm", "gcm.run",
    "rcm", "em.scn", "period", "ref.per", "resolution", "corrected")
  ## get all the factors which should prevail the aggrreagation
  id.without.index <- id[!id %in% INDEX]
  ## formula used to aggregate data.frame using reshape::cast
  id.formula <- formula(paste(paste(id.without.index, collapse = " + "),
                              " ~ variable", sep = ""))
  
  ## from wide to long dataframe
  x.melt <- reshape::melt(x, id = id)
  ## aggregate long and get wide again
  x.aggr <- reshape::cast(x.melt, id.formula, fun)

  return(x.aggr)
}

