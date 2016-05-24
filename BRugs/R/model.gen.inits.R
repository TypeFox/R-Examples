"modelGenInits" <-
function()
#   Generate initial values for OpenBUGS model
{
    command <- paste("BugsEmbed.GenerateInitsGuard;",  "BugsEmbed.GenerateInits")
    .CmdInterpreter(command)
    if(getOption("BRugsVerbose")) 
        buffer()
}
