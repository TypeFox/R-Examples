### R code from vignette source 'iBUGS.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=65)
program="";bugs.directory=""


###################################################
### code chunk number 2: start-iBUGS (eval = FALSE)
###################################################
## ## the GUI will show up once the package is loaded
## library(iBUGS)
## ## or call iBUGS() to generate another GUI
## iBUGS()


###################################################
### code chunk number 3: attach-ex (eval = FALSE)
###################################################
## dat = list(x=1:3,y=rnorm(5))
## ## personally we don't recommend attaching R objects that is a bad habit
## attach(dat)


###################################################
### code chunk number 4: bugs-directory (eval = FALSE)
###################################################
## if (nzchar(prog <- Sys.getenv("ProgramFiles")) && 
##     length(bugs.dir <- list.files(prog, "^(Open|Win)BUGS.*")) && 
##     length(bugs.exe <- dirname(list.files(file.path(prog, bugs.dir), 
##         pattern = "(Open|Win)BUGS.*\\.exe$", full.names = TRUE, 
##         recursive = TRUE)))) {
##     ## if we can find OpenBUGS, use it prior to WinBUGS
##     program = ifelse(length(grep("OpenBUGS", bugs.exe)), "OpenBUGS", 
##         "WinBUGS")
##     ## ignore multiple directories if (several versions of) BUGS installed in multiple places
##     bugs.directory = bugs.exe[grep(program, bugs.exe)][1]
## }


###################################################
### code chunk number 5: gui-demo (eval = FALSE)
###################################################
## library(gWidgets)
## options(guiToolkit = "RGtk2")
## ## create a window and add a button to it
## gw=gwindow('GUI Demo')
## gb=gbutton('Click me!', container=gw, handler=function(h,...){svalue(h$obj)=paste(svalue(h$obj),'haha!')})


