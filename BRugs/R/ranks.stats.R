"ranksStats" <-
function(node)
#   Calculates ranks statistics for vector valued node in OpenBUGS model
{
    if(length(node) > 1 || node == "*")
        stop("node cannot be a vector, nor '*'")
    nodeName <- sQuote(node)
    command <- paste("RanksEmbed.SetVariable(", nodeName, "); RanksEmbed.StatsGuard;",
                             "RanksEmbed.Stats")
    .CmdInterpreter(command)
    buffer <- file.path(tempdir(), "buffer.txt")
    rlb <- readLines(buffer)
    len <- length(rlb)
    if (len > 1) 
        read.table(buffer)
    else
        message(rlb)
}
