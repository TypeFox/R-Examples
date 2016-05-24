startWinEdt <- function(InstallRoot, ApplData, WinEdtVersion, args=NULL){
    if(WinEdtVersion < 6){
        shell(paste('""', InstallRoot, '\\WinEdt.exe" -C="R-WinEdt" -E=', 
                     shQuote(normalizePath(file.path(ApplData, "R.ini"))), '"', sep = ""), 
               wait = FALSE)
    } else {
        shell(paste('""', InstallRoot, '\\WinEdt.exe""', if(!is.null(args)) args, sep = ""), wait = FALSE)
    }
}
