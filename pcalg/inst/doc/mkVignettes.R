## A version that works for "everyone" (not just Markus)
## *after* having installed or already reqquire()d  'pcalg':

manualInst.vignette <- function(fstem, package="pcalg", verbose=FALSE) {
    mkf <- function(ext) paste(fstem, ext, sep=".")
    Rnw <- system.file("doc", mkf("Rnw"), package=package, mustWork = TRUE)
    if(interactive() && verbose)
        file.show(Rnw)
    else message("Using file ", Rnw)
    pd <- packageDescription(package)
    message("Package '", package,"'; built: \"", pd$Built,'";\n',
            "  Loaded from ", dirname(dirname(attr(pd,"file"))))
    ## But we need to work from the *source*package: need *.bib, *.sty extra PDFs, etc
    ## o.wd <- setwd(dirname(Rnw)); on.exit(setwd(o.wd))
    pkgSrcDir <-
        switch(Sys.getenv("USER"),
	       "maechler" = file.path("~/R/Pkgs", package),# or "~/R/Pkgs/pcalg-dev"
               "kalischm" = ".....",    # PATH of pcalg or pcalg-dev
               stop("Must add your (username, pcalg-source) in file {inst/}doc/mkVignettes.R "))

    srcDES <- read.dcf(file.path(pkgSrcDir, "DESCRIPTION"))
    ## now check, if 'Version', 'Date' agree:
    for(key in c("Version", "Date"))
        if(pd[[key]] != unname(srcDES[,key]))
           stop(sprintf("'%s' disagreement: Pkg: '%s',  Source: '%s'",
                       key, pd[[key]], srcDES[,key]))

    o.wd <- setwd(file.path(pkgSrcDir, "vignettes")); on.exit(setwd(o.wd))

    ## Now manually "install" the vignette pdf and R :
    Sweave (Rnw)
    tools::texi2pdf(mkf("tex"))
    Stangle(Rnw)
    ## and test if the code works:
    Rfile <- mkf("R")
    source(Rfile)
    message("\nend of source()ing file ", Rfile,"; wd= ",getwd(),"\n\n")

    (pkgDoc <- file.path(pkgSrcDir, "inst","doc"))
    f3 <- c(Rnw, Rfile, mkf("pdf"))
    message("Now copying 3 files ", paste(f3, collapse=", "),
            "\n --> to ",pkgDoc)
    r <- file.copy(f3, pkgDoc, overwrite=TRUE)
    ## cleanup all but pdf
    file.remove(Rfile, mkf("tex"), # mkf("bbl"),
                mkf("aux"), mkf("log"), mkf("out"), mkf("blg"),
                list.files(patt=paste0(fstem, "-*\\.pdf$")))
    r ## --> TRUE TRUE TRUE  if it works
}

## Now execute it:

## Mlibrary():  for MM, gets "newly installed" version of pkg:
if(exists("Mlibrary", mode="function")) Mlibrary("pcalg") 
stopifnot(r <- manualInst.vignette(fstem = "pcalgDoc"))
