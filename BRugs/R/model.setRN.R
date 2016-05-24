"modelSetRN" <-
function(state)
#   Set the seed of random number generator
{
    if(!state %in% 1:14)
        stop("state must be an integer from 1 to 14")
    state <- as.integer(state)
    command <- paste("BugsEmbed.SetRNGuard; BugsEmbed.SetRNState(", state, ")" )
    invisible(.CmdInterpreter(command))
}
