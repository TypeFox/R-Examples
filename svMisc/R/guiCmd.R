guiCmd <- function (command, ...)
{
    ## This function sends a command to the GUI client
    ## The actual code is a custom function named .guiCmd in SciViews:TempEnv
	CmdFun <- getTemp(".guiCmd", mode = "function")
    if (!is.null(CmdFun)) return(CmdFun(command, ...)) else return(NULL)
}

guiLoad <- function (...)
{
	## Ask the GUI client to select a .Rdata file to load()
	guiCmd("load", ...)
}

guiSource <- function (...)
{
	## Ask the GUI client to select a .R file to source()
	guiCmd("source", ...)  # TODO: should use sys.source() here
}

guiSave <- function (...)
{
	## Ask the GUI client for a file where to save some data
	guiCmd("save", ...)
}

guiImport <- function (...)
{
	## Ask the client to display a dialog for importing some data
	guiCmd("import", ...)
}

guiExport <- function (...)
{
	## Ask the client to display a dialog for exporting some data
	guiCmd("export", ...)
}

guiReport <- function (...)
{
	## Ask the client to display a dialog for reporting data (send a view...)
	guiCmd("report", ...)
}

guiSetwd <- function (...)
{
	## Ask the GUI client to select a directory to set as active
	guiCmd("setwd", ...)
}
