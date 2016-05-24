"samplesSample" <-
function(node)
#   Get stored sample for single component of OpenBUGS name
{
    if(samplesGetFirstChain() > samplesGetLastChain())
        stop("Number of first chain is larger than last chain!")
    if(length(node) != 1)
        stop("Exactly one scalar node must be given.")
    sM <- samplesMonitors(node)[1]
    if(sM == "model must be initialized before monitors used")
        stop("model must be initialized / updated / monitored before samplesSample is used")
    if(length(grep("^no monitor set for variable", sM)))
        stop(sM)
    nodeSize <- .OpenBUGS(c("BugsRobjects.SetVariable", "BugsRobjects.GetSize"),
                          c("CharArray","Integer"),
                          list(node,NA))[[2]]
    if(nodeSize > 1)
        stop("Only scalar nodes such as ", node, "[1] are allowed.")
    sampleSize <- samplesSize(node)
    sample <- .OpenBUGS(c(.SamplesGlobalsCmd(node), "SamplesEmbed.SampleValues"),
                        c("CmdInterpreter","RealArray"),
                        list(node,double(sampleSize)))[[2]]
    sample
}
