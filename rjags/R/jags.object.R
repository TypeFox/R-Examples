#  R package rjags file R/jags.object.R
#  Copyright (C) 2007-2011 Martyn Plummer
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License version
#  2 as published by the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

update.jags <- function(object, n.iter = 1, by, progress.bar, ...)
{
    if (!is.numeric(n.iter) || n.iter < 1) {
        stop("Invalid n.iter")
    }

    adapting <- .Call("is_adapting", object$ptr(), PACKAGE="rjags")
    on.exit(object$sync())

    if (missing(progress.bar)) {
        progress.bar <- getOption("jags.pb")
    }
    if (!is.null(progress.bar)) {
        match.arg(progress.bar, c("text","gui","none"))
        if (progress.bar=="none")
            progress.bar <- NULL
    }
    
    do.pb <- interactive() && !is.null(progress.bar) && n.iter >= 100
    if (do.pb) {
        start.iter <- object$iter()
        end.iter <- start.iter + n.iter
        pb <- switch(progress.bar,
                     "text" = txtProgressBar(start.iter, end.iter,
                     initial = start.iter, style=3, width=50,
                     char=ifelse(adapting,"+","*")),
                     "gui" = updatePB(start.iter, end.iter, adapting))
    }
    
    ## Set refresh frequency for progress bar
    if (missing(by) || by <= 0) {
        ##In JAGS 3.x.y there is a memory reallocation bug when
        ##monitoring that slows down updates. Drop refresh
        ##frequency to avoid triggering memory reallocations.
        ##by <- min(ceiling(n.iter/50), 100)
        by <- ceiling(n.iter/50)
    }
    else {
        by <- ceiling(by)
    }
    
    ## Do updates
    n <- n.iter
    while (n > 0) {
        .Call("update", object$ptr(), min(n,by), PACKAGE="rjags")
        if (do.pb) {
            switch(progress.bar,
                   "text" = setTxtProgressBar(pb, object$iter()),
                   "gui" =  setPB(pb, object$iter()))
        }
        n <- n - by
    }
    if (do.pb) {
        close(pb)
    }

    invisible(NULL)
}

adapt <- function(object, n.iter, end.adaptation = FALSE, ...)
{
    if(.Call("is_adapting", object$ptr(), PACKAGE="rjags")) {
        if(n.iter > 0)
			update.jags(object, n.iter, ...)
        ok <- .Call("check_adaptation", object$ptr(), PACKAGE="rjags")
        if (end.adaptation) {
            .Call("adapt_off", object$ptr(), PACKAGE="rjags")
        }
        return(ok)

    }
    else {
        return(TRUE)
    }
}

coef.jags <- function(object, chain = 1, ...) {
    if (!is.numeric(chain) || chain < 1 || chain > object$nchain()) {
        stop("Invalid chain")
    }
    object$state(internal=FALSE)[[chain]]
}

variable.names.jags <- function(object, ...) {
    .Call("get_variable_names", object$ptr(), PACKAGE="rjags")
}
