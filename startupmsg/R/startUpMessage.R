######

# new S3 classes for Startup-message


### constructor of condition "StartupMessage"
StartupMessage <- function (message, call = NULL, pkg = "", type = "version", 
                            endline = FALSE) 
structure(list(message = message, call = call, package = pkg, type = type,
               endline = endline),
               class = c("StartupMessage", "packageStartupMessage", "condition", 
                         "message", "simpleMessage"))

### accessor to slot package
startupPackage <- function(startupmessage) {startupmessage$"package"}
### accessor to slot type
startupType <- function(startupmessage) {startupmessage$"type"}
### accessor to slot endline
startupEndline <- function(startupmessage) {startupmessage$"endline"}


### suppressing Startup messages by a wrapper
suppressStartupMessages<-
function (expr) 
withCallingHandlers(expr, StartupMessage = 
                          function(c) invokeRestart("muffleMessage"))

onlytypeStartupMessages<-
function (expr,atypes="version") 
{withCallingHandlers(expr, StartupMessage = 
                     function(c) {invokeRestart(r = "onlytypeMessage", 
                                                c0 = c, atypes=atypes)}) }

### generating a startupMessage
startupMessage <- function(..., domain=NULL, pkg="", type="version", endline = FALSE) {
    withRestarts( withCallingHandlers(
                       message(..., domain=domain), 
                       message = function(cond)
                              {SM <- StartupMessage(conditionMessage(cond), 
                                                    conditionCall(cond), 
                                                    pkg, type, endline)
                              signalCondition(SM)
                              }      ), 
                  onlytypeMessage = function(c0,atypes){
                               if(startupType(c0) %in% atypes) 
                                  writeLines(conditionMessage(c0),stderr())               
                                                    }, 
                  #as suggested by Seth Falcon:
                  custom = function(c,f) f(c),
                  muffleMessage = function() NULL )
    invisible(NULL)
}

###############################################################
#Utilities for reading the DESCRIPTION file and NEWS file out 
#                     for starting information on the package 
###############################################################

readVersionInformation <- function(pkg, library=NULL){      
# next lines are taken from Valentin Todorov's package "rrcov"
    if(is.null(library)) library <- .Library
    ver <- read.dcf(file.path(library, pkg, "DESCRIPTION"), "Version")
    ver <- as.character(ver)
    title <- read.dcf(file.path(library, pkg, "DESCRIPTION"), "Title")
    title <- as.character(title)
#
    list(ver=ver, title=title)
}

readURLInformation <- function(pkg, library=NULL){      
# next lines are taken from Valentin Todorov's package "rrcov"
    if(is.null(library)) library <- .Library
    URL <- read.dcf(file.path(library, pkg, "DESCRIPTION"), "URL")
    if(is.na(URL)||(is.character(URL)&&length(URL)==0)) return(NULL)
    else return(as.character(URL))
}

pointertoNEWS <- function(pkg, library=NULL){
    if(file.exists(file.path(system.file(package=pkg, lib.loc=library),"NEWS")))
       return(paste("NEWS(\"",pkg,"\")",sep=""))
    else return(NULL)   
}

infoShow <- function(pkg, filename, library=NULL)
   {file.show(file.path(system.file(package = pkg, lib.loc=library), 
    paste(filename,sep=.Platform$file.sep,collapse=.Platform$file.sep)))}
### filename may also be given as a vector of characters with the corresponding
### names of folders i.e. c(folder.1,....,folder.n,filename)
### (to be system-independent) --- for Windows and Linux
### the usual [folder.1/..../folder.n/]filename will do

NEWS<-function(pkg, library=NULL) 
{   ## inspired by Andy Liaw
    infoShow(pkg, filename="NEWS", library=library)
}
#######################################################################

### analogously:
TOBEDONE<-function(pkg, library=NULL)
{   ## inspired by Andy Liaw
    infoShow(pkg, filename="TOBEDONE", library=library)
}
