## This function determines the dependencies for the specified
## package, exluding only packages found in "base".
getDependencies <- function (pkgs,
                             dependencies = c("Depends", "Imports", "LinkingTo"),
                             installed=TRUE,
                             available=TRUE,
                             base=FALSE,
                             recommended=FALSE)
{
    pkgs.in = pkgs
    if (is.null(dependencies))
        return(unique(pkgs))
    dep2 <- c("Depends", "Imports", "LinkingTo")


    if(installed && !available)
        all.packages <- installed.packages()
    else if (available && !installed)
        all.packages <- available.packages()
    else
        all.packages <- as.matrix(smartbind(##
                                  installed.packages(),
                                  available.packages()
                                  ))
    rownames(all.packages) <- all.packages[,"Package"]

    p0 <- unique(pkgs)
    miss <- !p0 %in% row.names(all.packages)
    if (sum(miss)) {
        warning(sprintf(ngettext(sum(miss), "package %s is not available (for %s)",
            "packages %s are not available (for %s)"), paste(sQuote(p0[miss]),
            collapse = ", "), sub(" *\\(.*", "", R.version.string)),
            domain = NA, call. = FALSE)
        if (sum(miss) == 1L && !is.na(w <- match(tolower(p0[miss]),
            tolower(row.names(all.packages))))) {
            warning(sprintf("Perhaps you meant %s ?", sQuote(row.names(all.packages)[w])),
                call. = FALSE, domain = NA)
        }
        flush.console()
    }

    ## Whether to exclude base and recommended packages
    if(!base || !recommended)
        {
            priority <- NULL
            if(!base)        priority <- c("base", priority)
            if(!recommended) priority <- c("recommended", priority)
            installed <- installed.packages(priority=priority)
        }
    else
        installed <- installed.packages()[FALSE,]

    p0 <- p0[!miss]

        p1 <- p0
        not_avail <- character()
        repeat {
            deps <- apply(all.packages[p1, dependencies, drop = FALSE],
                1L, function(x) paste(x[!is.na(x)], collapse = ", "))

            res <- .clean_up_dependencies2(
                deps,
                installed=installed,
                all.packages)
            not_avail <- c(not_avail, res[[2L]])
            deps <- unique(res[[1L]])
            deps <- deps[!deps %in% c("R", pkgs)]
            if (!length(deps))
                break
            pkgs <- c(deps, pkgs)
            p1 <- deps
            if (!is.null(dep2)) {
                dependencies <- dep2
                dep2 <- NULL
            }
        }
        if (length(not_avail)) {
            not_avail <- unique(not_avail)
            warning(sprintf(ngettext(length(not_avail),
                                     "dependency %s is not available",
                                     "dependencies %s are not available"),
                            paste(sQuote(not_avail),
                                  collapse = ", ")),
                    domain = NA,
                    call. = FALSE,
                    immediate. = TRUE)
            flush.console()
        }
        pkgs <- unique(pkgs)
        pkgs <- pkgs[pkgs %in% row.names(all.packages)]
        p0 <- pkgs


    p0[ ! p0 %in% pkgs.in  ]
}

