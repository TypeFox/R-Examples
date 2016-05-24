"modelUpdate" <-
function(numUpdates, thin = 1, overRelax = FALSE)
#   Update the each chain in OpenBUGS model numUpdates * thin time
{
    if(!is.numeric(numUpdates))
        stop("numUpdates ", "must be numeric")
    numUpdates <- as.integer(numUpdates)
    if(!is.numeric(thin))
        stop("thin ", "must be numeric")
    thin <- as.integer(thin)
    if(!is.logical(overRelax))
        stop("overRelax ", "must be logical") 
    command <- paste("BugsEmbed.UpdateGuard",
        ";BugsEmbed.thin := ", thin,
        ";BugsEmbed.overRelax := ", as.integer(overRelax),
        ";BugsEmbed.updates := ", numUpdates,
        ";BugsEmbed.Update")
    .CmdInterpreter(command)
    if(getOption("BRugsVerbose")) 
        buffer()
}
