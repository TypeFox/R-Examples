.onAttach <- function(libname, pkgname) {
    Pver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, Pver))
    packageStartupMessage("Insieme di funzioni di supporto al volume")
    packageStartupMessage(sQuote('Laboratorio di Statistica con R'))
    packageStartupMessage("Iacus-Masarotto, MacGraw-Hill Italia, 2006.")
    packageStartupMessage(paste("Si veda",
                   sQuote("library(help=\"labstatR\")"),
                   "per i comandi disponibili."))
    packageStartupMessage("Nota: 'mean.a', 'mean.g' e 'hist.pf' sono state sostituite rispettivamente da 'meana', 'meang' e 'histpf'")
                   
}

