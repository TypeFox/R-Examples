# GIS_LOCK 110814 RSB, suggested by Brian Oney
get.GIS_LOCK <- function() {
    Sys.getenv("GIS_LOCK")
}

set.GIS_LOCK <- function(pid) {
    if (missing(pid)) pid <- round(runif(1, 1, 1000))
    pid <- as.integer(pid)
    stopifnot(!is.na(pid))
    Sys.setenv(GIS_LOCK=pid)
}

unset.GIS_LOCK <- function() {
    Sys.unsetenv("GIS_LOCK")
}

unlink_.gislock <- function() {
    gl <- paste(Sys.getenv("GISDBASE"), Sys.getenv("LOCATION_NAME"),
        Sys.getenv("MAPSET"), ".gislock", sep="/")
    if (file.exists(gl)) unlink(gl)
}

initGRASS <- function(gisBase, home, SG, gisDbase, location, mapset,
    override=FALSE, use_g.dirseps.exe=TRUE, pid) {
    if (nchar(Sys.getenv("GISRC")) > 0 && !override)
      stop("A GRASS location is already in use; to override, set override=TRUE")

    if (nchar(get.GIS_LOCK()) > 0) {
      if(!override)
        stop("A GIS_LOCK environment variable is present; to override, set override=TRUE")
      else unset.GIS_LOCK()
    }

    if (missing(pid)) pid <- round(runif(1, 1, 1000))
    pid <- as.integer(pid)
    stopifnot(!is.na(pid))

    if (!file.exists(gisBase)) stop(paste(gisBase, "not found"))

    SYS <- get("SYS", envir=.GRASS_CACHE) 
    if (SYS == "WinNat") {
# grass63.bat
        Sys.setenv(GRASSDIR=gisBase)
        if (missing(home)) home <- Sys.getenv("USERPROFILE")
        Sys.setenv(HOME=home)
        Sys.setenv(WINGISBASE=Sys.getenv("GRASSDIR"))
        Sys.setenv(GRASS_PROJSHARE=paste(Sys.getenv("GRASSDIR"),
            "\\proj", sep=""))
        Sys.setenv(GRASS_SH=paste(Sys.getenv("GRASSDIR"), 
            "\\msys\\bin\\sh.exe", sep=""))
        Wpath <- Sys.getenv("PATH")
        if (length(grep(basename(Sys.getenv("GRASSDIR")), Wpath)) < 1) {
            Sys.setenv(PATH=paste(Sys.getenv("GRASSDIR"), "\\msys\\bin;", 
                Sys.getenv("PATH"), sep=""))
            Sys.setenv(PATH=paste(Sys.getenv("GRASSDIR"), "\\extrabin;",
                Sys.getenv("GRASSDIR"), "\\extralib;",
                Sys.getenv("PATH"), sep=""))
            Sys.setenv(PATH=paste(Sys.getenv("GRASSDIR"), "\\tcl-tk\\bin;",
                Sys.getenv("GRASSDIR"), "\\sqlite\\bin;",
                Sys.getenv("GRASSDIR"), "\\gpsbabel;", 
                Sys.getenv("PATH"), sep=""))
# etc/Init.bat
            GRASS_addons <- Sys.getenv("GRASS_ADDON_PATH")
            if (GRASS_addons == "")
                Sys.setenv(PATH=paste(Sys.getenv("WINGISBASE"), "\\bin;",
                    Sys.getenv("WINGISBASE"), "\\lib;",
                    Sys.getenv("PATH"), sep=""))
            else 
                Sys.setenv(PATH=paste(Sys.getenv("WINGISBASE"), "\\bin;",
                    Sys.getenv("WINGISBASE"), "\\lib;",
                    GRASS_addons, ";", Sys.getenv("PATH"), sep=""))
            ePyPATH <- Sys.getenv("PYTHONPATH")
            if ((length(grep(basename(Sys.getenv("WINGISBASE")), ePyPATH)) < 1) 
                || nchar(ePyPATH) == 0) {
                GrPyPATH <- paste(Sys.getenv("WINGISBASE"), "/etc/python",
                    sep="")
                if (nchar(ePyPATH) > 0)
                    Sys.setenv(PYTHONPATH=paste(GrPyPATH, ePyPATH, sep=";"))
                else Sys.setenv(PYTHONPATH=GrPyPATH)
            }
            Sys.setenv("PYTHONHOME"=paste(Sys.getenv("WINGISBASE"),
                "Python27", sep="/"))
            Sys.setenv("GRASS_PYTHON"=paste(Sys.getenv("WINGISBASE"),
                "extrabin/python.exe", sep="/"))
            pyScripts <- basename(list.files(paste(Sys.getenv("WINGISBASE"),
                "scripts", sep="/"), pattern="py$"))
            names(pyScripts) <- sub("\\.py", "", pyScripts)
            assign("pyScripts", pyScripts, envir=.GRASS_CACHE)
        }
        Sys.setenv(WINGISRC=paste(Sys.getenv("HOME"), "\\.grassrc6", sep=""))
        if (file.exists(Sys.getenv("WINGISRC")) && !override)
            stop("A GISRC file already exists; to override, set override=TRUE")
        Sys.setenv(GISRC="junk")
        cat("GISDBASE:", getwd(), "\n", file=Sys.getenv("WINGISRC"))
        cat("LOCATION_NAME: <UNKNOWN>", "\n", file=Sys.getenv("WINGISRC"),
            append=TRUE)
        cat("MAPSET: <UNKNOWN>", "\n", file=Sys.getenv("WINGISRC"),
            append=TRUE)
        gisrc <- ifelse (use_g.dirseps.exe, system(paste("g.dirseps.exe -g",
            shQuote(Sys.getenv("WINGISRC"))), intern=TRUE),
            Sys.getenv("WINGISRC"))
        assign("addEXE", .addexe(), envir=.GRASS_CACHE)
        Sys.setenv(GISRC=gisrc)
        if (!missing(gisDbase)) {
            if (!file.exists(gisDbase)) dir.create(gisDbase)
        } else {
            gisDbase <- tempdir()
        }
        gisDbase <- ifelse (use_g.dirseps.exe, system(paste("g.dirseps.exe -g",
            shQuote(gisDbase)), intern=TRUE), gisDbase)
    } else if (SYS == "unix") {
        Sys.setenv(GISBASE=gisBase)
        if (missing(home)) home <- Sys.getenv("HOME")
        ePATH <- Sys.getenv("PATH")
        if (length(grep(basename(Sys.getenv("GISBASE")), ePATH)) < 1) {
            Sys.setenv(PATH=paste(Sys.getenv("GISBASE"), "/bin:",
            Sys.getenv("GISBASE"), "/scripts",
            ifelse(nchar(ePATH) == 0, "", ":"), ePATH, sep=""))
        }
        eLDPATH <- Sys.getenv("LD_LIBRARY_PATH")
        if (length(grep(basename(Sys.getenv("GISBASE")), eLDPATH)) < 1) {
            Sys.setenv(LD_LIBRARY_PATH=paste(Sys.getenv("GISBASE"), "/lib:",
            ifelse(nchar(eLDPATH) == 0, "", ":"), eLDPATH, sep=""))
        }
#FIXME Sys.info()["sysname"] == "Darwin"
        Sys.setenv(GISRC=paste(home, "/.grassrc6", sep=""))
#FIXME
        if (file.exists(Sys.getenv("GISRC")) && !override)
            stop("A GISRC file already exists; to override, set override=TRUE")
        ePyPATH <- Sys.getenv("PYTHONPATH")
        if (length(grep(basename(Sys.getenv("GISBASE")), ePyPATH)) < 1 
            || nchar(ePyPATH) == 0) {
            GrPyPATH <- paste(Sys.getenv("GISBASE"), "etc", "python", sep="/")
            if (nchar(ePyPATH) > 0)
                 Sys.setenv(PYTHONPATH=paste(GrPyPATH, ePyPATH, sep=":"))
            else Sys.setenv(PYTHONPATH=GrPyPATH)
        }
        if (!missing(gisDbase)) {
            if (!file.exists(gisDbase)) dir.create(gisDbase)
        } else {
            gisDbase <- tempdir()
        }
        cat("GISDBASE:", gisDbase, "\n", file=Sys.getenv("GISRC"))
        cat("LOCATION_NAME: <UNKNOWN>", "\n", file=Sys.getenv("GISRC"),
            append=TRUE)
        cat("MAPSET: <UNKNOWN>", "\n", file=Sys.getenv("GISRC"),
            append=TRUE)
    } else stop(paste("Platform variant", SYS, "not supported"))
    set.GIS_LOCK(pid)
    assign("INIT_USED", TRUE, envir=.GRASS_CACHE)
    assign("GIS_LOCK", pid, envir=.GRASS_CACHE)
    system(paste(paste("g.gisenv", get("addEXE", envir=.GRASS_CACHE), sep=""),
        shQuote(paste("set=GISDBASE=", gisDbase))))
    if (missing(location)) location <- basename(tempfile())
    loc_path <- paste(gisDbase, location, sep="/")
    if (!file.exists(loc_path)) dir.create(loc_path)
    if (!file.exists(paste(loc_path, "PERMANENT", sep="/")))
        dir.create(paste(loc_path, "PERMANENT", sep="/"))
    if (missing(mapset)) mapset <- basename(tempfile())
    if (!file.exists(paste(loc_path, mapset, sep="/")))
        dir.create(paste(loc_path, mapset, sep="/"))
    system(paste(paste("g.gisenv", get("addEXE", envir=.GRASS_CACHE), sep=""),
        shQuote(paste("set=GISDBASE", gisDbase, sep="="))))
    system(paste(paste("g.gisenv", get("addEXE", envir=.GRASS_CACHE), sep=""),
        shQuote(paste("set=LOCATION_NAME", location, sep="="))))
    system(paste(paste("g.gisenv", get("addEXE", envir=.GRASS_CACHE), sep=""),
        shQuote(paste("set=MAPSET", mapset, sep="="))))
    system(paste(paste("g.gisenv", get("addEXE", envir=.GRASS_CACHE), sep=""),
        shQuote("set=GRASS_GUI=text")))
    Sys.setenv(GISBASE=gisBase)
    Sys.setenv(GISDBASE=gisDbase)
    Sys.setenv(LOCATION_NAME=location)
    Sys.setenv(MAPSET=mapset)
    gv <- system(paste("g.version", get("addEXE", envir=.GRASS_CACHE),
        sep=""),  intern=TRUE)

    comp <- .compatibleGRASSVersion(gv)
    if ( !comp ){
        stop( attr(comp, "message") )
    }

    assign("GV", gv, envir=.GRASS_CACHE)

    pfile <- paste(loc_path, "PERMANENT", "DEFAULT_WIND", sep="/")
    if (!file.exists(pfile)) {
        mSG <- !missing(SG)
        if (mSG) bb <- bbox(SG)
        if (mSG) gt <- gridparameters(SG)
        cat("proj:       0\n", file=pfile)
        cat("zone:       0\n", file=pfile, append=TRUE)
        cat("north:      ", ifelse(mSG, bb[2, "max"], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("south:      ", ifelse(mSG, bb[2, "min"], 0), "\n",
            sep="", file=pfile, append=TRUE)
        cat("east:       ", ifelse(mSG, bb[1, "max"], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("west:       ", ifelse(mSG, bb[1, "min"], 0), "\n",
            sep="", file=pfile, append=TRUE)
        cat("cols:       ", ifelse(mSG, gt$cells.dim[1], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("rows:       ", ifelse(mSG, gt$cells.dim[2], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("e-w resol:  ", ifelse(mSG, gt$cellsize[1], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("n-s resol:  ", ifelse(mSG, gt$cellsize[2], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("top:        1\n", sep="", file=pfile, append=TRUE)
        cat("bottom:     0\n", sep="", file=pfile, append=TRUE)
        cat("cols3:      ", ifelse(mSG, gt$cells.dim[1], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("rows3:      ", ifelse(mSG, gt$cells.dim[2], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("depths:     1\n", sep="", file=pfile, append=TRUE)
        cat("e-w resol3: ", ifelse(mSG, gt$cellsize[1], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("n-s resol3: ", ifelse(mSG, gt$cellsize[2], 1), "\n",
            sep="", file=pfile, append=TRUE)
        cat("t-b resol:  1\n", sep="", file=pfile, append=TRUE)
    }
    tfile <- paste(loc_path, "PERMANENT", "WIND", sep="/")
    if (!file.exists(tfile)) file.copy(pfile, tfile, overwrite=TRUE)
    tfile <- paste(loc_path, mapset, "WIND", sep="/")
    if (!file.exists(tfile)) file.copy(pfile, tfile, overwrite=TRUE)
    gmeta()
}
