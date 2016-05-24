"dicStats" <-
function()
#   Calculate dic statistics
{
    command <- "DevianceEmbed.SetVariable('*');DevianceEmbed.StatsGuard;DevianceEmbed.Stats"
    .CmdInterpreter(command)
    buffer <- file.path(tempdir(), "buffer.txt")
    rlb <- readLines(buffer)
    len <- length(rlb)
    if (len > 1) {
        writeLines(rlb, buffer)
        read.table(buffer)
    } else {
        message(rlb)
    }
}
