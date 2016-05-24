## CLI torture test, copyright (c) Ph. Grosjean (phgrosjean@sciviews.org)
## GNU GPL => 2 license
## A series of commands to check for R CLI (or console widget)
## Note: test this in English, but also, in internationalized versions of R!

## Copy all text in testCLIcmd.r and paste it in a console window. Copy all
## the content of the console into a document (for instance testCLIres.r)
## Then, run the following code that will do the same, but from a
## processSocket() as if it is issued by a client
## Then, compare both files

## To test in french, $LANG=fr_FR.UTF8 R or
Sys.setenv(LANG = "fr_FR.UTF8")

## Open a clean R session, then, issue these commands and copy/paste the content
## of testCLIcmd.r in it... Copy the results from the console into a *.save file
require(svSocket) || stop("Package 'svSocket' is required")
#source('/Users/phgrosjean/Documents/Pgm/SciViews/svSocket/R/RSocketServer.R')
TempEnv_()
res <- outfile <- out <- i <- cmdfile <- cmd0 <- cmd <- ""
options(width = 80)


## Run this to generate output with processSocket()
setwd('/Users/phgrosjean/Documents/Pgm/SciViews/svSocket/')
cmdfile <- 'testCLIcmd.R'
outfile <- 'testCLIcmd.out'
require(svSocket) || stop("Package 'svSocket' is required")
#source('/Users/phgrosjean/Documents/Pgm/SciViews/svSocket/R/RSocketServer.R')
options(width = 80)

## Read the file with commands
cmd <- readLines(cmdfile)
## Run these commands
out <- file(outfile, "w")
cmd0 <- ""
cat("> ", file = out)
for (i in 1:length(cmd)) {
    cat(cmd[i], "\n", sep = "", file = out)
    if (cmd0 == "") cmd0 <- cmd[i] else
        cmd0 <- paste(cmd0, cmd[i], sep = "<<<n>>>")
    res <- processSocket(cmd0, "", "")
    if (res != "+ ") cmd0 <- ""  # Not a multiline command
    cat(res, file = out)
}
close(out)
## Until here...




## User interrupt
## Here is a long calculation. Hit <esc> or <ctrl-c> in mother R app during
## calculation. See how interrupt is managed
###for (i in 1:10000000) i
## Here is a command that issues a very long output. Let it calculate
## and when output starts, hit <esc> or <ctrl-c> in mother R app
## See how interrupt is managed
###1:100000    # This is like if I have to use flush.console()
## TODO: correct handling of flush.console()

## Interaction at the command line
#scan()  # Read in numbers
## Other tests to be implement
#readline()
#non graphical menu()
#create graph (are there problems between the graph device and Tcl/Tk event loop?)
#test ask = TRUE for graphs and demos
#identify()
#locator()
#lattice graph
#rgl graph
#progress() from the svMisc package
#assignment & ->>
#test tcltk windows and widgets + <Tcl> messages (are they printed correctly?)
#:, ::, :::
#addTaskCallback => do we run it after execution?
#test history and the like
#q("yes"), q("no"), q()
#UTF-8 or other non ASCII characters (input and output)
#break
#debug
#trace
#browser
#test connections
#source
#sink
#capture.output
#try & tryCatch
#on.exit, sys.on.exit
#Cstack_info
#eval, evalq, eval.parent, local
#exists
#test size of console output and allow for automatic adjustment of width (like in Rgui)
#flush + flush.console
#shell & shell.exec
#stderr, stdin, stdout
#interactive()
#utf8ToInt, intToUtf8
#system, system.time
#pause (sm package)
#bringToTop
#dev.interactive
#graphics.off
#stayOnTop under Windows
#choose.dir
#choose.files
#chooseCRANmirror, install & update packages
#Data.entry, dataentry, de, edit, fix, fixInNamespace
#file.edit
#history
#page
#recover
#setWindowTitle under Windows
#winDialog, winDialogString under Windows
#winMenuXXXX under Windows
#gettext, ngettext, bindtextdomain
#readline functionnalities
#test command history with pgup/pgdwn
