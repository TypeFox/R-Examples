"dicClear" <-
function()
#   Clear monitor for dic
{
    command <- "DevianceEmbed.StatsGuard;DevianceEmbed.Clear"
    invisible(.CmdInterpreter(command))
}
