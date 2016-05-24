justDoItDoE <- function (command)
{
    Message()
    if (!getRcmdr("suppress.X11.warnings")) {
        messages.connection <- file(open = "w+")
        sink(messages.connection, type = "message")
        on.exit({
            sink(type = "message")
            close(messages.connection)
        })
    }
    else messages.connection <- getRcmdr("messages.connection")
    result <- try(eval(parse(text = command),envir = .GlobalEnv))
    ## warnings only, if no error abort
    if (!class(result)[1]=="try-error") checkWarnings(readLines(messages.connection))
    if (getRcmdr("RStudio")) 
        Sys.sleep(0)
     result
}
