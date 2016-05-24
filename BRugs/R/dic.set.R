"dicSet" <-
function()
#   Set a monitor for dic
{
    command <- "DevianceEmbed.SetVariable('*');DevianceEmbed.SetGuard;DevianceEmbed.Set"
    .CmdInterpreter(command)
    if(getOption("BRugsVerbose")) 
        buffer()
}
