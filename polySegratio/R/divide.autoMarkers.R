`divide.autoMarkers` <-
function(markers, description=paste("Markers split for",
                                      deparse(substitute(markers))),
                              parent.cols=c(1,2), extra.cols = NULL,
                              cols.drop = c(parent.cols,extra.cols))
{
  
  ## Function: divide.markers
  ## Purpose: read markers (or more correctly dominant 1,0) marker data
  ##          and return list object of containing markers data split
  ##          according to parental alleles, namely 1,0 for each parent and
  ##          1,1 for both parents
  ##
  ## Arguments:
  ## 
  ## markers:  matrix of 1, 0, NA indicating marker alleles
  ##           rownames are markernames, column names are progeny names
  ##           NB: If markers were simulated and stored as an object of
  ##               class simAutoMarkers then  simAutoMarkers$markers
  ##               may need to be split if parental markers misclassified
  ## title:         to be used for printing/plotting
  ## parent.cols:   column(s) for parental markers (default: 1,2)
  ## extra.cols:    extra column(s) to be subsetted (default: NULL)
  ## cols.drop:     columns to be dropped from markers
  ##                before splitting data which
  ##                can be set to NULL if no columns are to be dropped
  ##                (Default: c(parent.cols,extra.cols)
  ##
  ## Values:
  ##     p10, p01, p11  are lists for where the first, second
  ##     are heterozygous for parents 1, 2 and both resp. Each list contains
  ## description:  text containing a description for printing
  ## parent:       label for parent
  ## markers:      markers for specified parental type (including parents etc)
  ## extras:       extra column subsetted if specified
  ## seg.ratios:   segregation ratios as class segRatios
  
  ## index to subset appropriate markers
  index.p10 <- markers[,parent.cols[1]]==1 & markers[,parent.cols[2]]==0
  index.p01 <- markers[,parent.cols[1]]==0 & markers[,parent.cols[2]]==1
  index.p11 <- markers[,parent.cols[1]]==1 & markers[,parent.cols[2]]==1

  ## produce lists
  parent1 <- colnames(markers)[parent.cols[1]]
  parent2 <- colnames(markers)[parent.cols[2]]
  descp <- paste("Parent with 1 is", parent1,"and 0 is",parent2)
  p10 <- list(description=descp, parent.inherited=parent1,
              markers=markers[index.p10,-cols.drop],
              seg.ratios=segregationRatios(markers[index.p10,-cols.drop]))

  descp <- paste("Parent with 0 is", parent1,"and 1 is",parent2)
  p01 <- list(description=descp, parent.inherited=parent2,
              markers=markers[index.p01,-cols.drop],
              seg.ratios=segregationRatios(markers[index.p01,-cols.drop]))

  parent <- paste(parent1,parent2,sep=" & ")
  descp <- paste("Parents both  with 1 -", parent)
  p11 <- list(description=descp, parent.inherited=parent,
                   markers=markers[index.p11,-cols.drop],
              seg.ratios=segregationRatios(markers[index.p11,-cols.drop]))

  if (length(extra.cols) > 0) {
    p10$extras <- markers[index.p10, extra.cols]
    p01$extras <- markers[index.p01, extra.cols]
    p11$extras <- markers[index.p11, extra.cols]
  }
  
  ## return result
  res <- list(description=description, p10=p10, p01=p01,
              p11=p11, time.split=date(), call=match.call())
  oldClass(res) <- "divideAutoMarkers"
  return(res)
}

