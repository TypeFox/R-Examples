### tclVarFun.R - A series of additional function to manipulate Tcl variables
### and functions from within R and vice versa
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2007-01-01: fisrt version (for tcltk2_1.0-0)
###
### To do:
### - Add a catch {} in tclFun and handle it
### - A tclFunDispose() function to delete the Tcl equivalent of a function
### - Add a try construct in tclVarExists and tclVarFind
### - better manage catch{} in tclVarName

makeTclNames <- function (names, unique = FALSE)
{
    ## Make valid Tcl variable names (allow_ = TRUE by default in R >= 2.0.0)
    names <- make.names(names, unique = unique)
    ## There is a problem if the variable starts with a dot => prepend it with 'X'
    .names <- grep("^\\.", names)
    names[.names] <- paste("X", names[.names], sep = "")
    ## Although it is accepted, there could be problems with variable names
	## containing dots, so, replace them with '_'
    names <- gsub("\\.", "_", names)
    return(names)
}

### TODO: change this to use closure functions instead!!!
tclFun <- function (f, name = deparse(substitute(f)))
{
    ## Register a simple R function (without arguments) as a callback in Tcl,
	## and give it the same name under Tcl)
    ## Indeed, .Tcl.callback(f) does the job... but it gives criptic names
	## like R_call 0x13c7168

    ## Check that 'f' is a function with no arguments (cannot handle them yet)
    if (!is.function(f)) stop("'f' must be a function!")
    if (!is.null(formals(f))) stop("The function used cannot (yet) have arguments!")
    ## Make sure the name of the function is valid
    if (!is.character(name)) stop("'name' must be a character string!") else
		name <- make.names(name[1])

    res <- .Tcl.callback(f)
    ## Make sure this is correct (R_call XXXXXXXX)
    if (length(grep("R_call ", res) > 0)) {
		## Create a proc with the same name in Tcl
    	.Tcl(paste("proc ", name, " {} {", res, "}", sep = ""))
    }
    ## Return the R_call XXXXXXXX string, as .Tcl.callback() does
    return(res)
    ## Rem: if you delete the R 'f' function, the Tcl 'f' function still works!
    ## You have to explicitly delete the Tcl function
}

tclGetValue <- function (name)
{
    ## Get the value stored in a plain Tcl variable
    if (!is.character(name))
		stop("'name' must be a character!")
	name <- makeTclNames(name[1]) # The usual name conversion

    ## Create a temporary dual variable with tclVar() (name does not mather)
    Temp <- tclVar(init = "")

    ## Copy the content of the var of interest to it
    res <- tclvalue(.Tcl(paste("catch {set ", as.character(Temp), " $", name,
		"}", sep = "")))  # Return "0" if OK, "1" otherwise
	if (res != "0")
		stop(gettextf("Error when getting the value in the '%s' Tcl variable",
			name))

    ## Get the content of the temporary variable
    return(tclvalue(Temp)) # (Temp will be destroyed when the function exits)
}

tclSetValue <- function (name, value)
{
    ## This is the opposite of tclGetValue() and it is a wrapper
	## for 'set name value' Tcl command
    if (!is.character(name))
		stop("'name' must be a character!")
	name <- makeTclNames(name[1]) # The usual name conversion
    
	## Create a temporary dual variable with tclVar() (name does not mather)
    Temp <- tclVar(init = value)
    
    ## Copy the content of this variable to the tcl variable 'name'
	res <- tclvalue(.Tcl(paste("catch {set ", name, " $", as.character(Temp),
		"}", sep = "")))
	if (res != "0")
		stop(gettextf("Error when changing the value of the '%s' Tcl variable",
			name))
	
	## (Temp is destroyed when the function exits)
	return(invisible(name))  # Return the name of the Tcl variable invisibly
}

tclVarExists <- function (name)
    as.integer(tcl("info", "exists", name)) == 1

tclVarFind <- function (pattern)
    as.character(tcl("info", "vars", pattern))

tclVarName <- function (name, init = "", keep.existing = TRUE)
{
    ## tclVar gives names like ::RtclX automatically...
    ## We need to define names ourselve. This is what tclVarName does
    ## If keep existing == TRUE and the variable is already defined, then
    ## we keep its content, instead of initializing it with "init"
    if (!is.character(name)) stop("'name' must be a character!")
    name <- makeTclNames(name[1])  # Make sure the name is correct
    
    ## Temporary save potential content of the Tcl variable elsewhere
	## (catch in case the variable does not exist)
    if (isTRUE(keep.existing))
		.Tcl(paste("catch {set ZZZTempRvariable $", name, "}", sep = ""))

    ## Create the new dual Tcl-R variable
    l <- list(env = new.env())
    assign(name, NULL, envir = l$env)
    reg.finalizer(l$env, function(env) tcl("unset", ls(env)))
    class(l) <- "tclVar"
    tclvalue(l) <- init

    ## Possibly restore the content of the variable, if keep.existing == TRUE
    if (isTRUE(keep.existing)) {
        .Tcl(paste("catch {set", name, "$ZZZTempRvariable}"))
        ## Remove the temporary variable
        .Tcl("unset -nocomplain ZZZTempRvariable")
    }
    return(l)
}
