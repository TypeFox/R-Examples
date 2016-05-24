#' The reportr message reporting system
#' 
#' Functions for reporting informative messages, warnings and errors. These are
#' provided as alternatives to the \code{\link{message}}, \code{\link{warning}}
#' and \code{\link{stop}} functions in base R.
#' 
#' The \code{reportr} system for reporting messages provides certain useful
#' features over the standard R system, such as the incorporation of output
#' consolidation, message filtering, expression substitution, automatic
#' generation of stack traces for debugging, and conditional reporting based on
#' the current ``output level''. Messages of level at least equal to the value
#' of option \code{reportrStderrLevel} are written to standard error
#' (\code{\link{stderr}}); others are written to standard output
#' (\code{\link{stdout}}).
#' 
#' The output level is set by the \code{setOutputLevel} function, and governs
#' whether a particular call to \code{report} will actually report anything.
#' Output levels are described by the \code{OL} object, a list with components
#' \code{Debug}, \code{Verbose}, \code{Info}, \code{Warning}, \code{Question},
#' \code{Error} and \code{Fatal}. Any call to \code{report} using a level lower
#' than the current output level will produce no output. If \code{report} is
#' called before \code{setOutputLevel}, the output level will default to
#' \code{Info} (with a message).
#' 
#' The \code{flag} function is called like \code{report}, but it stores
#' messages for later reporting, like \code{\link{warning}}, rather than
#' reporting them immediately. Stored messages are reported when \code{report}
#' is next called, at which point multiple instances of the same message are
#' consolidated where possible. The user may also manually force stored
#' messages to be reported by calling \code{reportFlags}, or remove them with
#' \code{clearFlags}. Note that the output level at the time that
#' \code{reportFlags} is called (implicitly or explicitly) will determine
#' whether the flags are printed.
#' 
#' The \code{ask} function requests input from the user, using
#' \code{\link{readline}}, at output level \code{Question}. The text argument
#' then forms the text of the question, and \code{ask} returns the text
#' entered by the user.
#' 
#' The call \code{report(Error,\dots)} is largely similar to \code{stop(\dots)}
#' in most cases, except that a stack trace will be printed if the current
#' output level is \code{Debug}. The "abort" restart is invoked in any case. No
#' other standard conditions are signalled by \code{report}. Stack traces can
#' be generated at lower output levels, if desired, by setting the
#' \code{reportrStackTraceLevel} option.
#' 
#' The \code{withReportrHandlers} function evaluates \code{expr} in a context
#' in which R errors, warnings and messages will be handled by reportr, rather
#' than by the standard R functions.
#' 
#' The \code{prefixFormat} argument to \code{report} and \code{ask} controls
#' how the output message is formatted. It takes the form of a
#' \code{\link{sprintf}}-style format string, but with different expansions for
#' percent-escapes. Specifically, \code{"\%d"} expands to a series of stars
#' indicating the current stack depth; \code{"\%f"} gives the name of the
#' function calling \code{report} or \code{ask}; \code{"\%l"} and \code{"\%L"}
#' give lower and upper case versions of the level of the message,
#' respectively; and \code{"\%p"} expands to the ID of the current R process
#' (see \code{\link{Sys.getpid}}). The default is \code{"\%d\%L: "}, giving a
#' prefix such as \code{"* * INFO: "}, but this default can be overridden by
#' setting the \code{reportrPrefixFormat} option.
#' 
#' A number of other options influence the output produced by reportr.
#' \code{getOutputLevel} and \code{setOutputLevel} get and set the
#' \code{reportrOutputLevel} option, which can be set directly if preferred.
#' The options \code{reportrMessageFilterIn} and \code{reportrMessageFilterOut}
#' can contain a single character string representing a Perl regular
#' expression, in which case only messages which match
#' (\code{reportrMessageFilterIn}) or do not match
#' (\code{reportrMessageFilterOut}) the regular expression will be retained.
#' Likewise, the \code{reportrStackFilterIn} and \code{reportrStackFilterOut}
#' options filter the call stack.
#' 
#' @param level The level of output message to produce, or for
#'   \code{setOutputLevel}, the minimum level to display. See Details.
#' @param \dots Objects which can be coerced to mode \code{character}. These
#'   will be passed through function \code{\link{es}} (from the \code{ore}
#'   package) for expression substitution, and then printed with no space
#'   between them. Options to \code{\link{es}}, such as \code{round}, may also
#'   be given.
#' @param prefixFormat The format of the string prepended to the message. See
#'   Details.
#' @param default A default return value, to be used when the message is
#'   filtered out or the output level is above \code{Question}.
#' @param expr An expression to be evaluated.
#' 
#' @return These functions are mainly called for their side effects, but
#'   \code{getOutputLevel} returns the current output level,
#'   \code{withReportrHandlers} returns the value of the evaluated expression,
#'   and \code{ask} returns a character vector of length one giving the user's
#'   response.
#' 
#' @examples
#' setOutputLevel(OL$Warning)
#' report(Info, "Test message")  # no output
#' setOutputLevel(OL$Info)
#' report(Info, "Test message")  # prints the message
#' 
#' flag(Warning, "Test warning")  # no output
#' flag(Warning, "Test warning")  # repeated warning
#' reportFlags()  # consolidates the warnings and prints the message
#' 
#' \dontrun{name <- ask("What is your name?")
#' report(OL$Info, "Hello, #{name}")}
#' 
#' @seealso \code{\link{es}} (in package \code{ore}) for expression
#'   substitution (which is performed on messages). \code{\link{message}},
#'   \code{\link{warning}}, \code{\link{stop}} and \code{\link{condition}} for
#'   the normal R message and condition signalling framework.
#' @author Jon Clayden
#' 
#' @name reportr
#' @aliases OL
NULL

.resolveOption <- function (name)
{
    value <- getOption(name)
    if (is.null(value))
        value <- .Defaults[[name]]
    return (value)
}

.evaluateLevel <- function (level)
{
    name <- as.character(substitute(level,parent.frame()))
    if (length(name) == 1 && name %in% names(OL))
        return (OL[[name]])
    else
        return (level)
}

#' @rdname reportr
#' @export
setOutputLevel <- function (level)
{
    level <- .evaluateLevel(level)
    if (level %in% OL$Debug:OL$Fatal)
        options(reportrOutputLevel=level)
    invisible(NULL)
}

#' @rdname reportr
#' @export
getOutputLevel <- function ()
{
    if (is.null(getOption("reportrOutputLevel")))
    {
        setOutputLevel(OL$Info)
        report(OL$Info, "Output level is not set; defaulting to \"Info\"", prefixFormat="")
        level <- OL$Info
    }
    else
        level <- getOption("reportrOutputLevel")
    
    names(level) <- names(which(OL == level))
    return (level)
}

.truncate <- function (strings, maxLength)
{
    lengths <- nchar(strings)
    strings <- substr(strings, 1, maxLength)
    lines <- ore.split(ore("\n",syntax="fixed"), strings, simplify=FALSE)
    strings <- sapply(lines, "[", 1)
    strings <- paste(strings, ifelse(lengths>maxLength | sapply(lines,length)>1, " ...", ""), sep="")
    return (strings)
}

#' @rdname reportr
#' @export
withReportrHandlers <- function (expr)
{
    result <- withCallingHandlers(expr, message=function (m) {
        report(OL$Info, ore.subst("\n$","",m$message))
        invokeRestart("muffleMessage")
    }, warning=function (w) {
        flag(OL$Warning, w$message)
        invokeRestart("muffleWarning")
    }, error=function (e) {
        if (is.null(e$call))
            report(OL$Error, e$message)
        else
            report(OL$Error, e$message, " (in \"", as.character(e$call)[1], "(", .truncate(paste(as.character(e$call)[-1],collapse=", "),100), ")\")")
    })
    
    reportFlags()
    return (result)
}

.getCallStack <- function ()
{
    callStrings <- .truncate(as.character(sys.calls()), 100)
    
    handlerFunLoc <- which(callStrings %~% "^withReportrHandlers\\(")
    if (length(handlerFunLoc) > 0)
        callStrings <- callStrings[-seq_len(handlerFunLoc[length(handlerFunLoc)]+1)]
    
    raisingFunLoc <- which(callStrings %~% "^(ask|flag|report|reportFlags|message|warning|stop)\\(")
    if (length(raisingFunLoc) > 0)
        callStrings <- callStrings[-(raisingFunLoc[1]:length(callStrings))]
    
    filterIn <- .resolveOption("reportrStackFilterIn")
    filterOut <- .resolveOption("reportrStackFilterOut")
    if (!is.null(filterIn))
        callStrings <- callStrings[callStrings %~% as.character(filterIn)[1]]
    if (!is.null(filterOut))
        callStrings <- callStrings[!(callStrings %~% as.character(filterOut)[1])]
    
    return (callStrings)
}

.buildPrefix <- function (level, format = NULL)
{
    if (!is.null(format))
        prefix <- as.character(format)[1]
    else
        prefix <- as.character(.resolveOption("reportrPrefixFormat"))[1]
    
    if (prefix == "")
        return (prefix)
    else
    {
        if (prefix %~% "\\%(d|f)")
            stack <- .getCallStack()

        if (prefix %~% "\\%d")
            prefix <- ore.subst(ore("%d",syntax="fixed"), paste(rep("* ",length(stack)),collapse=""), prefix, all=TRUE)
        if (prefix %~% "\\%f")
            prefix <- ore.subst(ore("%f",syntax="fixed"), ore.subst("^([\\w.]+)\\(.+$","\\1",stack[length(stack)]), prefix, all=TRUE)
        if (prefix %~% "\\%l")
            prefix <- ore.subst(ore("%l",syntax="fixed"), tolower(names(OL)[which(OL==level)]), prefix, all=TRUE)
        if (prefix %~% "\\%L")
            prefix <- ore.subst(ore("%L",syntax="fixed"), toupper(names(OL)[which(OL==level)]), prefix, all=TRUE)
        if (prefix %~% "\\%p")
            prefix <- ore.subst(ore("%p",syntax="fixed"), as.character(Sys.getpid()), prefix, all=TRUE)

        return (prefix)
    }
}

.buildMessage <- function (..., round = NULL, signif = NULL)
{
    # This assumes that the environment containing relevant variables is the grandparent of the current one
    message <- es(paste(..., sep=""), round=round, signif=signif, envir=parent.frame(2))
    keep <- TRUE
    
    filterIn <- .resolveOption("reportrMessageFilterIn")
    filterOut <- .resolveOption("reportrMessageFilterOut")
    if (!is.null(filterIn))
        keep <- keep & (message %~% as.character(filterIn)[1])
    if (!is.null(filterOut))
        keep <- keep & (!(message %~% as.character(filterOut)[1]))
    
    if (keep)
        return (message)
    else
        return (NULL)
}

#' @rdname reportr
#' @export
ask <- function (..., default = NULL, prefixFormat = NULL)
{
    outputLevel <- getOutputLevel()
    message <- .buildMessage(...)
    if (!interactive() || outputLevel > OL$Question || is.null(message))
        return (default)
    else
    {
        reportFlags()
        ans <- readline(paste(.buildPrefix(OL$Question,prefixFormat), message, " ", sep=""))
        return (ans)
    }
}

#' @rdname reportr
#' @export
report <- function (level, ..., prefixFormat = NULL)
{
    level <- .evaluateLevel(level)
    outputLevel <- getOutputLevel()
    message <- .buildMessage(...)
    if (outputLevel > level || is.null(message))
        return (invisible(NULL))
    
    reportFlags()
    
    if (level >= .resolveOption("reportrStderrLevel"))
        file <- stderr()
    else
        file <- stdout()
    
    cat(paste(.buildPrefix(level,prefixFormat), message, "\n", sep=""), file=file)
    
    if (outputLevel == OL$Debug)
    {
        if (level >= .resolveOption("reportrStackTraceLevel"))
        {
            stack <- .getCallStack()
            cat("--- Begin stack trace ---\n", file=file)
            for (i in 1:length(stack))
                cat(rep("* ", i), stack[i], "\n", sep="", file=file)
            cat("---  End stack trace  ---\n", file=file)
        }
    }
    
    if (level == OL$Error)
        invokeRestart("abort")
}

#' @rdname reportr
#' @export
flag <- function (level, ...)
{
    level <- .evaluateLevel(level)
    if (getOutputLevel() == OL$Debug)
    {
        if (level >= .resolveOption("reportrStackTraceLevel"))
        {
            report(level, ...)
            return (invisible(NULL))
        }
    }
    
    message <- .buildMessage(...)
    if (is.null(message))
        return (invisible(NULL))
    currentFlag <- list(list(level=level, message=message))
    
    if (!exists("reportrFlags",.Workspace) || is.null(.Workspace$reportrFlags))
        .Workspace$reportrFlags <- currentFlag
    else
        .Workspace$reportrFlags <- c(.Workspace$reportrFlags, currentFlag)
}

#' @rdname reportr
#' @export
reportFlags <- function ()
{
    if (exists("reportrFlags",.Workspace) && !is.null(.Workspace$reportrFlags))
    {
        levels <- unlist(lapply(.Workspace$reportrFlags, "[[", "level"))
        messages <- unlist(lapply(.Workspace$reportrFlags, "[[", "message"))
        
        # This is before the call to report() to avoid infinite recursion
        clearFlags()
        
        for (message in unique(messages))
        {
            locs <- which(messages == message)
            level <- max(levels[locs])
            if (length(locs) == 1)
                report(level, message, prefixFormat="%L: ")
            else
                report(level, paste("[x",length(locs),"] ",message,sep=""), prefixFormat="%L: ")
        }
    }
}

#' @rdname reportr
#' @export
clearFlags <- function ()
{
    .Workspace$reportrFlags <- NULL
}
