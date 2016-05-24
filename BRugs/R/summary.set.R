"summarySet" <-
function(node)
#   Set summary monitor for node in OpenBUGS model
{
    nodeName <- sQuote(node)
    for(i in seq(along=nodeName)){
        command <- paste("SummaryEmbed.SetVariable(", nodeName[i], "); SummaryEmbed.SetGuard;",
                                "SummaryEmbed.Set")
        .CmdInterpreter(command)
        buffer <- file.path(tempdir(), "buffer.txt")
        rlb <- readLines(buffer)
        if(getOption("BRugsVerbose"))
            message("Variable ", nodeName[i], ": ", rlb)
    }
    invisible()    
}
