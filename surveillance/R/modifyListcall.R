################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Function to modify a list _call_ according to another one similar to
### what utils::modifyList (by Deepayan Sarkar) does for list objects.
###
### Copyright (C) 2012 Sebastian Meyer
### $Revision: 463 $
### $Date: 2012-12-06 17:26:39 +0100 (Don, 06 Dez 2012) $
################################################################################

is.listcall <- function (x)
{
    is.call(x) &&
    as.character(x[[1]]) %in% c("list", "alist")
}

modifyListcall <- function (x, val)
{
    stopifnot(is.listcall(x), is.listcall(val))
    xnames <- names(x)[-1]
    for (v in names(val)[nzchar(names(val))]) {
        xv <- if (v %in% xnames && is.listcall(x[[v]]) && is.listcall(val[[v]]))
            modifyListcall(x[[v]], val[[v]]) else val[[v]]
        x[v] <- list(xv)                # allows for NULL value of val[[v]]
    }
    x
}
