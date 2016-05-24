
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:               DESCRIPTION:    
#  .fGarchEnv              Create GARCH Environment
#  .setfGarchEnv           Set GARCH Environment       
#  .getfGarchEnv           Get GARCH Environment
################################################################################


.fGarchEnv <- 
new.env(hash = TRUE)


# ------------------------------------------------------------------------------


.setfGarchEnv <-
function(...)
{
    x <- list(...)
    nm <- names(x)
     if (is.null(nm) || "" %in% nm)
        stop("all arguments must be named")
    sapply(nm, function(nm) assign(nm, x[[nm]], envir = .fGarchEnv))
    invisible()
}


# ------------------------------------------------------------------------------


.getfGarchEnv <-
function(x = NULL, unset = "")
{
    if (is.null(x))
        x <- ls(all.names = TRUE, envir = .fGarchEnv)
        ###     unlist(mget(x, envir = .fGarchEnv, mode = "any",
        ###         ifnotfound = as.list(unset)), recursive = FALSE)
    get(x, envir = .fGarchEnv, mode = "any")
}


################################################################################

