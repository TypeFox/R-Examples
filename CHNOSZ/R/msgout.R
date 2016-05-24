# CHNOSZ/msgout.R
# msgout() is a variation on message() from base R
# 20120117 jmd, adapted from...

#  File src/library/base/R/message.R
#  Part of the R package, http://www.R-project.org
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


# a variant of message(), using stdout instead of stderr
# and with a default setting of appendLF = FALSE
# this is now used for informative messages in CHNOSZ 
# (in info() and other functions)
# so that they show up on the console and in Sweave output, but are suppressed
# when running the 'test_that' testing scripts
# (messages produced by cat() would mess up the rows of green dots)
# the appendLF = FALSE is so that the cat("...\n") statements used previously
# can all be replaced with msgout without further modification
msgout <-
function(..., domain = NULL, appendLF = FALSE)
{
    args <- list(...)
    cond <- if (length(args) == 1L && inherits(args[[1L]], "condition")) {
        if(nargs() > 1L)
            warning("additional arguments ignored in message()")
        args[[1L]]
    } else {
        msg <- .makeMessage(..., domain=domain, appendLF = appendLF)
        call <- sys.call()
        simpleMessage(msg, call)
    }
    defaultHandler <- function(c) {
        ## Maybe use special connection here?
        #cat(conditionMessage(c), file=stderr(), sep="")
        cat(conditionMessage(c), file="", sep="")
    }
    withRestarts({
        signalCondition(cond)
        ## We don't get to the default handler if the signal
        ## is handled with a non-local exit, e.g. by
        ## invoking the muffleMessage restart.
        defaultHandler(cond)
    }, muffleMessage = function() NULL)
    invisible()
}
