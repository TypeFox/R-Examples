.onLoad <- function (lib, pkg)
{
	## Create .GUI that contains information about the default GUI
	guiAdd(".GUI")
}

#.onUnload <- function (libpath)
#{
#	## We do nothing, because other packages may also use .GUI
#}

.packageName <- "svGUI"

## A copy of TempEnv() from svMisc to avoid a useless dependency (only
## internally used)
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
