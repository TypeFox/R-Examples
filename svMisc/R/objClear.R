objClear <- function (id = "default")
{
    ## Clear any reference to a given 'id' object browser
	id <- as.character(id)[1]  # Make sure id is character
	if (id == "") id <- "default"
	Pars <- getTemp(".guiObjParsCache", default = list())
	Pars[[id]] <- NULL
	assignTemp(".guiObjParsCache", Pars)
    List <- getTemp(".guiObjListCache", default = list())
    List[[id]] <- NULL
	assignTemp(".guiObjListCache", List)
	## Also delete corresponding files
	Root <- objDir()
	ParsFile = file.path(Root, paste("Pars_", id, ".txt", sep=""))
	if (file.exists(ParsFile)) unlink(ParsFile)
	ListFile = file.path(Root, paste("List_", id, ".txt", sep=""))
	if (file.exists(ListFile)) unlink(ListFile)
	MenuFile = file.path(Root, paste("Menu_", id, ".txt", sep=""))
	if (file.exists(MenuFile)) unlink(MenuFile)
	invisible(TRUE)
}
