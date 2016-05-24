### tk2reg.R - Functions to interface the reg Tcl package under Windows
### The reg Tcl package is provided with standard install of Tcl in R and in
### ActiveState, but there are no R functions to specifically use it in tcltk
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2007-01-01: fisrt version (for tcltk2_1.0-0)
###
### To do:
### - use a try for tk2reg.get()
### - Change double call for a try() for tk2reg.keys(), tk2reg.type()
###   and tk2reg.values()
### - Add "none to the type of supported formats?

.tk2reg.require <- function ()
{
	## Make sure tcl/tk registry is operational
	if (.Platform$OS.type != "windows")
		stop("This is a Windows-specific function!")
	if (!capabilities("tcltk"))
		stop("This version of R cannot use Tcl/Tk!")
	## This should be installed by default with the tcltk package under Windows
	res <- tclRequire("registry", warn = TRUE)
	if (inherits(res, "tclObj")) res <- tclvalue(res)
	if (res[1] == FALSE)
		stop("Unable to find the 'registry' Tcl/tk package!")
	return(res)  # The package version number
}

tk2reg.broadcast <- function ()
{
    ## Used to warn running apps that something changes in the registry
    ## Use this when you change an environment variable
    .tk2reg.require()
    res <- tclvalue(.Tcl("catch {registry broadcast \"Environment\"}"))
    return(res == "0")  # "0" if OK, "1" otherwise
}

tk2reg.delete <- function (keyname, valuename)
{
    ## Delete a registry value in a key (take care when using this!)
    .tk2reg.require()
    keyname <- as.character(keyname[1])
	valuename <- as.character(valuename[1])
    res <- tclvalue(.Tcl(paste("catch {registry delete {", keyname, "} {",
		valuename, "}}", sep = "")))  # return "0" if OK, "1" otherwise
	return(res == "0")
}

tk2reg.deletekey <- function (keyname)
{
    ## Completely delete a registry key (take care when using this!)
    .tk2reg.require()
    keyname <- as.character(keyname[1])
	res <- tclvalue(.Tcl(paste("catch {registry delete {", keyname, "}}",
		sep = "")))  # Return "0" if OK (even if already deleted) or "1"
	return(res == "0")
}

tk2reg.get <- function (keyname, valuename)
{
    ## Get the content of a key
    .tk2reg.require()
    keyname <- as.character(keyname[1])
	valuename <- as.character(valuename[1])
	## First get the type of this registry key
	Type <- tk2reg.type(keyname, valuename)
	if (is.na(Type)) return(NA)	# The key does not exists
	## The key is found... retrieve its data
	res <- .Tcl(paste("registry get {", keyname, "} {",
		valuename, "}", sep = ""))
	## Convert according to its type...
	res <- switch(Type,
		sz = tclvalue(res),				# A single string
		expand_sz = tclvalue(res),		# This string is NOT expanded!
		multi_sz = as.character(res),	# A vector of strings
		dword = as.numeric(res),		# Numbers,... check very large numbers!
		dword_big_endian = as.numeric(res),	# Is this correct???
		res)  # Other types are probably not handled well!
	return(res)
}

tk2reg.keys <- function (keyname)
{
    ## Get a list of all subkeys in a key
    .tk2reg.require()
    keyname <- as.character(keyname[1])
    ## First check if the command succeeds
    res <- tclvalue(.Tcl(paste("catch {registry keys {", keyname, "}}",
		sep = "")))  # Return "0" if OK, "1" otherwise
	if (res != "0") return(NA)  # Indicate that keyname is probably inexistant
	## Now run the command unprotected
	res <- as.character(.Tcl(paste("registry keys {", keyname, "}", sep = "")))
	return(res)
}

tk2reg.set <- function (keyname, valuename, data,
type = c("sz", "expand_sz", "multi_sz", "dword", "dword_big_endian"))
{
    ## Set a registry key value
    .tk2reg.require()
    keyname <- as.character(keyname[1])
	valuename <- as.character(valuename[1])
    data <- as.character(data)
    if (length(data) > 1)  # Collapse into one string, using {} as separators
    	data <- paste(data, collapse = "\n")
	type <- type[1]
    if (!(type %in% c("sz", "expand_sz", "multi_sz", "dword", "dword_big_endian",
    	"binary", "link", "resource_list", "none")))
		stop("Unrecognized 'type'!")
    res <- tclvalue(.Tcl(paste("catch {registry set {", keyname, "} {",
		valuename, "} {", data, "} {", type, "}}" , sep = "")))
	return(res == "0")  # Because "0" if OK, and "1" otherwise
}

tk2reg.setkey <- function (keyname)
{
    ## Set a registry key
    keyname <- as.character(keyname[1])
	.tk2reg.require()
    res <- tclvalue(.Tcl(paste("catch {registry set {", keyname, "}}",
		sep = "")))  # Return "0" if OK, "1" otherwise
	return(res == "0")
}

tk2reg.type <- function (keyname, valuename)
{
    ## Get the type of a key...
    .tk2reg.require()
    keyname <- as.character(keyname[1])
	valuename <- as.character(valuename[1])
    ## First test it to see if the command succeeds (i.e., if the key exists)
    res <- tclvalue(.Tcl(paste("catch {registry type {", keyname, "} {",
		valuename, "}}", sep = "")))  # return "0" if OK, "1" otherwise
	if (res != "0") return(NA)  # The key is probably missing
	## Run the command unprotected now
	res <- tclvalue(.Tcl(paste("registry type {", keyname, "} {", valuename,
		"}", sep = "")))
	return(res)
}

tk2reg.values <- function (keyname)
{
    ## Get a list of all values in a key
    keyname <- as.character(keyname[1])
	.tk2reg.require()
	## First check if the command succeeds
	res <- tclvalue(.Tcl(paste("catch {registry values {", keyname, "}}",
		sep = "")))  # Returns "0" if OK, "1" otherwise
	if (res != "0") return(NA)  # The key probably does not exist!
	## We issue the command now without protection
    res <- as.character(.Tcl(paste("registry values {", keyname, "}",
		sep = "")))
	return(res)
}
