"summaryStats" <-
function(node)
#   Calculates statistics for summary monitor associated with node in OpenBUGS model
{
    nodeName <- sQuote(node)
    if (is.R())
      result <- data.frame(mean=NULL, sd=NULL, val2.5pc=NULL, 
                           median=NULL, val97.5pc=NULL, sample=NULL)
    else
      result <- data.frame(mean=numeric(), sd=numeric(), val2.5pc=numeric(), 
                           median=numeric(), val97.5pc=numeric(), sample=numeric())
    for(i in seq(along=nodeName)){
        command <- paste("SummaryEmbed.SetVariable(", nodeName[i], "); SummaryEmbed.StatsGuard;",
                         "SummaryEmbed.Stats")
        .CmdInterpreter(command)
        buffer <- file.path(tempdir(), "buffer.txt")
        rlb <- readLines(buffer)
        len <- length(rlb)
        if (len > 1)
            result <- rbind(result, read.table(buffer))
        else{
            if(length(grep("val97.5pc", rlb)))
                message("Variable ", nodeName[i], " has probably not been updated")
            else if(getOption("BRugsVerbose"))
                message("Variable ", nodeName[i], ": ", rlb)
        }
    }
    return(result)
}
