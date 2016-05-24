.onLoad <- function (lib, pkg)
{
	if (.tmpfilesAllowed()) {
		## Clear menus
		.menuClear()
		## ... and create the default one
		.menuFileInit()
		.ctxMenuFileInit()	
	}
}

.onUnload <- function (libpath)
{
	## Clear menus
	if (interactive()) try(.menuClear())
}

.Last.lib <- function (libpath)
{
	## Clear menus
	if (interactive()) try(.menuClear())
}

.packageName <- "svDialogs"

.isJGR <- function () "package:JGR" %in% search()

## Avoid dependency on svMisc for those two functions
.TempEnv <- function ()
{
    pos <-  match("SciViews:TempEnv", search())
    if (is.na(pos)) { # Must create it
        `SciViews:TempEnv` <- list()
        Attach <- function (...) get("attach", mode = "function")(...)
        Attach(`SciViews:TempEnv`, pos = length(search()) - 1)
        rm(`SciViews:TempEnv`)
        pos <- match("SciViews:TempEnv", search())
    }
    return(pos.to.env(pos))
}

.getTemp <- function (x, default = NULL, mode = "any", item = NULL)
{
    if (is.null(item)) Mode <- mode else Mode <- "any"
    if  (exists(x, envir = .TempEnv(), mode = Mode, inherits = FALSE)) {
        dat <- get(x, envir = .TempEnv(), mode = Mode, inherits = FALSE)
        if (is.null(item)) return(dat) else {
            item <- as.character(item)[1]
            if (inherits(dat, "list") && item %in% names(dat)) {
                dat <- dat[[item]]
                if (mode != "any" && mode(dat) != mode) dat <- default
                return(dat)
            } else {
                return(default)
            }
        }
    } else { # Variable not found, return the default value
        return(default)
    }
}

.assignTemp <- function (x, value, replace.existing = TRUE)
    if (replace.existing || !exists(x, envir = .TempEnv(), mode = "any",
		inherits = FALSE))
        assign(x, value, envir = .TempEnv())
