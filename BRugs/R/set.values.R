"setValues" <- function(nodeLabel, values)
# set value of node
{
    nodeLabel <- as.character(nodeLabel)
# NA handling, now internal in OpenBUGS?
#    cv <- currentValues(nodeLabel)
#    DoNotSetNA <- is.na(values) & !is.na(cv)
#    if(any(DoNotSetNA))
#        warning("Some NA values formerly had a non-NA value -- left unchanged")
#    values[DoNotSetNA] <- cv[DoNotSetNA]
    nodeSize <- .OpenBUGS(c("BugsRobjects.SetVariable", "BugsRobjects.GetSize"),
                          c("CharArray","Integer"),
                          c(nodeLabel,NA))[[2]]
    if(nodeSize == -1)
        stop(nodeLabel, " is not a node in BUGS model")
    numChains <- getNumChains()
    if(length(values) != nodeSize*numChains)
        stop("length(values) does not correspond to the node size and number of chains")
    .OpenBUGS(c("BugsRobjects.SetVariable", "BugsRobjects.SetValues"),
              c("CharArray","RealArray"),
              list(nodeLabel,as.double(values)))[[2]]
    invisible()
}
