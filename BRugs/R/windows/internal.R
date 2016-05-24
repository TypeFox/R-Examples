### Run a list of OpenBUGS API command strings

.OpenBUGS.platform <- function(cmds, cmdtypes, args)
{
    if (.Platform$r_arch == "x64"){
        out <- .OpenBUGS.helper(cmds, cmdtypes, args)
    }
    else if (.Platform$r_arch == "i386") {
        ncmds <- length(cmds)
        out <- vector(ncmds, mode="list")
        for (i in 1:ncmds) {
            out[[i]] <-  switch(cmdtypes[i],
                                "CmdInterpreter" = {
                                    res <- .C("CmdInterpreter", cmds[i], nchar(cmds[i]), integer(1), PACKAGE="libOpenBUGS")
                                    handleRes(res[[3]])
                                    res
                                },
                                "Integer" = {
                                    values <- .C("Integer", cmds[i], nchar(cmds[i]), integer(1), integer(1), PACKAGE="libOpenBUGS")
                                    handleRes(values[[4]])
                                    as.integer(values[[3]])
                                },
                                "CharArray" = {
                                    values <- .C("CharArray", cmds[i], nchar(cmds[i]), args[[i]], nchar(args[[i]]), integer(1), PACKAGE="libOpenBUGS")
                                    handleRes(values[[5]])
                                    values[[3]]
                                },
                                "RealArray" = {
                                    values <- .C("RealArray", cmds[i], nchar(cmds[i]), args[[i]], length(args[[i]]), integer(1), PACKAGE="libOpenBUGS")
                                    handleRes(values[[5]])
                                    values[[3]]
                                })

        }
    }
    else {
        stop("Unknown architecture ", .Platform$r_arch, " , should be i386 or x64")
    }
    out
}
