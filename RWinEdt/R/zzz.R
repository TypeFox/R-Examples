.onAttach <- function(lib, pkg) {
    ## we have a NAMESPACE now: .First.lib <- function(lib, pkg) {
    if(.Platform$GUI != "Rgui" || .Platform$OS.type != "windows")
        stop("\nR-WinEdt is designed only for *RGui* on Windows!")
    ## we have a NAMESPACE now: library.dynam("RWinEdt", pkg, lib)
    .gW <- getWinEdt()
    if(!.gW$RWinEdtInstalled)
        installWinEdt(.gW$InstallRoot, .gW$ApplData, .gW$WinEdtVersion)
    .gW <<- .gW <- getWinEdt()
    if(!.gW$RWinEdtInstalled)
        warning(paste("Looks like the R-WinEdt installation failed,",
                      "please consider to start R with Administrator privileges the first time you load R-WinEdt.", sep="\n"))
    .rwloc <- file.path(system.file(package="RWinEdt"), "PlugIn")
    if(.gW$RWinEdtVersion != scan(file.path(.rwloc, "R.ver", fsep = "\\"), quiet = TRUE)){
        if("YES" ==
            winDialog(type = "yesno", paste("Looks like you have installed a new Version of R-WinEdt.",
                "\nUser customized settings of R-WinEdt might be lost after upgrade!",
                "\nUpgrade?"))){
                    installWinEdt(.gW$InstallRoot, .gW$ApplData, .gW$WinEdtVersion, force = TRUE)
                    .gW <<- .gW <- getWinEdt()
                    if(.gW$RWinEdtVersion != scan(file.path(.rwloc, "R.ver", fsep = "\\"), quiet = TRUE))
                        warning(paste("Looks like the R-WinEdt upgrade failed,",
                                      "please consider to start R with Administrator privileges",
                                      "the first time you load R-WinEdt after upgrade and try again.", sep="\n"))
        }
    }

    WindowTitle <- getWindowTitle()
    send2RedtLoc <- file.path(.gW$ApplData, "send2R.edt", fsep = "\\")
    if(ismdi()){
        packageStartupMessage("    You are running R in MDI mode which is *not*\n",
            "    supported for non-english translations of RGui.\n",
            "    It is recommended to use R in SDI mode which can be\n",
            "    set in the command line or by clicking in the Menu:\n",
            "    Edit - GUI Preferences: SDI, then Save and restart R.")
        file.copy(file.path(.rwloc, "send2R--mdi.edt"), send2RedtLoc, overwrite = TRUE)
        if(WindowTitle != "RGui"){ # In case we have another Window title such as for REvolution R
            edtFile <- readLines(send2RedtLoc)
            edtFile <- gsub("RGui", WindowTitle, edtFile)
            writeLines(edtFile, con = send2RedtLoc)
        }
    } else{
        file.copy(file.path(.rwloc, "send2R--sdi.edt"), send2RedtLoc, overwrite = TRUE)
        if(WindowTitle != "R Console"){ # In case we have another Window title such as for REvolution R
            edtFile <- readLines(send2RedtLoc)
            edtFile <- gsub("R Console", WindowTitle, edtFile)
            writeLines(edtFile, con = send2RedtLoc)
        }
    }

    winMenuAdd("R-WinEdt")
    winMenuAddItem("R-WinEdt", "Set and start R-WinEdt", "startWinEdt(.gW$InstallRoot, .gW$ApplData, .gW$WinEdtVersion)")
    ## internal pager seems to be better:
    ## winMenuAddItem("R-WinEdt", "Set WinEdt as pager", "options(pager = options('editor')[[1]])")
    winMenuAddItem("R-WinEdt", "Reset R-WinEdt settings", "installWinEdt(.gW$InstallRoot, .gW$ApplData, .gW$WinEdtVersion, force = NULL)")
    options(editor = paste("\"", .gW$InstallRoot, "\\WinEdt.exe\"",
                           if(.gW$WinEdtVersion < 6) " -c=\"R-Editor\" -e=R.ini -V",
                           sep = ""))
    startWinEdt(.gW$InstallRoot, .gW$ApplData, .gW$WinEdtVersion)
}
