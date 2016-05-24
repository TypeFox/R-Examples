"infoNodeValues" <-  
function(nodeLabel)
# Get current value of node
{
    nodeLabel <- as.character(nodeLabel)
    out <- .OpenBUGS(c("BugsRobjects.SetVariable", "BugsRobjects.GetSize"),
                     c("CharArray","Integer"),
                     list(nodeLabel, NA))
    nodeSize <- out[[2]]
    if(nodeSize == -1) 
        stop(nodeLabel, " is not a node in BUGS model")
    numChains <- getNumChains()
    out <- .OpenBUGS(c("BugsRobjects.SetVariable", "BugsRobjects.GetValues"),
                     c("CharArray","RealArray"),
                     list(nodeLabel, double(nodeSize*numChains)))
    values <- matrix(out[[2]], nrow=nodeSize, ncol=numChains)
    values 
}

infoNodeMethods <- function(nodeLabel)
{
    nodeName <- sQuote(nodeLabel)
    command <- paste("BugsEmbed.SetNode(",nodeName,"); BugsEmbed.Methods");
    .CmdInterpreter(command)
    buffer <- file.path(tempdir(), "buffer.txt")
    result <- read.table(buffer, sep="\t", skip = 1, as.is=TRUE, col.names=c("empty", "Node", "Type", "Size", "Depth"))[,-1]
    for (i in 1:2)
        result[,i] <- gsub(" ", "", result[,i])
    result
}

infoNodeTypes <- function(nodeLabel)
{
    nodeName <- sQuote(nodeLabel)
    command <- paste("BugsEmbed.SetNode(",nodeName,"); BugsEmbed.Types");
    .CmdInterpreter(command)
    buffer <- file.path(tempdir(), "buffer.txt")
    result <- read.table(buffer, sep="\t", skip = 1, as.is=TRUE, col.names=c("empty", "Node", "Type"))[,-1]
    for (i in 1:2)
        result[,i] <- gsub(" ", "", result[,i])
    result
}
