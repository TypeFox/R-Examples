"samplesSize" <-
function(node)
#   Size of stored sample of single component of OpenBUGS name
{
    sM <- samplesMonitors(node)
    # Doesn't distinguish between nodes not in the model and nodes not monitored
    # so returns 0 for non-existent nodes
    if (any(grep("^no monitor set", sM))) return(0)    
    if (any(grep("^model has probably not yet been updated", sM))) return(0)
    if (any(grep("^inference can not be made", sM))) { warning(sM); return(0) }
    if(length(sM) > 1 || sM != node)
        stop("node must be a scalar variable from the model")
    size <- .OpenBUGS(c(.SamplesGlobalsCmd(node), "SamplesEmbed.SampleSize"),
                      c("CmdInterpreter","Integer"))[[2]]
    size
}
