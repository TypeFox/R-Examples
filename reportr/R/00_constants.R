#' @import ore
#' @export
OL <- list(Debug=1L, Verbose=2L, Info=3L, Warning=4L, Question=5L, Error=6L, Fatal=7L)

.Defaults <- list(reportrOutputLevel=OL$Info,
                  reportrPrefixFormat="%d%L: ",
                  reportrStderrLevel=OL$Warning,
                  reportrStackTraceLevel=OL$Error,
                  reportrMessageFilterIn=NULL,
                  reportrMessageFilterOut=NULL,
                  reportrStackFilterIn=NULL,
                  reportrStackFilterOut=NULL)

.Workspace <- new.env()
