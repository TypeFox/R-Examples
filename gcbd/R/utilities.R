
.dbfile <- "/var/tmp/gcbd.sqlite"

requirements <- function() {

    # are we on Unix ?
    stopifnot(.Platform$OS.type == "unix")

    # are we on a Debian or Ubuntu system ?
    stopifnot(file.exists("/etc/debian_version"))

    # reference blas as minimum standard installed ?
    stopifnot(file.exists("/usr/share/doc/libblas3gf/copyright"))
    stopifnot(file.exists("/usr/share/doc/liblapack3gf/copyright"))

    # wajig frontend to dpkg, apt, ... ?
    stopifnot(file.exists("/usr/bin/wajig"))

    # "littler" scripting frontend to R
    stopifnot(file.exists("/usr/bin/r"))

    # is Goto available -- gotoblas2-help installs this
    stopifnot(file.exists("/etc/default/gotoblas2-helper"))

    # is database available ?
    if ( ! file.exists(.dbfile)) {
        createDatabase(.dbfile)
    }

    invisible(TRUE)
}

installAtlas <- function() {
    res <- system("sudo apt-get -y --force-yes install libatlas3gf-base", intern=TRUE, ignore.stderr=TRUE)
    invisible(res)
}

purgeAtlas <- function() {
    ##-- old local build wajig purge libatlas-base-dev libatlas-dev  libatlas3gf-amd64sse3  libatlas3gf-base
    res <- system("sudo apt-get -y --force-yes purge libatlas3gf-base", intern=TRUE, ignore.stderr=TRUE)
    invisible(res)
}

installMKL <- function() {
    res <- system("sudo apt-get -y --force-yes install revolution-mkl r-revolution-revobase revolution-r", intern=TRUE, ignore.stderr=TRUE)
    invisible(res)
}

purgeMKL <- function() {
    res <- system("sudo apt-get -y --force-yes purge revolution-mkl r-revolution-revobase revolution-r", intern=TRUE, ignore.stderr=TRUE)
    invisible(res)
}

installGoto <- function() {
    res <- system("sudo wajig -y install /var/spool/gotoblas2-helper/archive/gotoblas2_1.13-1_amd64.deb", intern=TRUE, ignore.stderr=TRUE)
    invisible(res)
}

purgeGoto <- function() {
    res <- system("sudo apt-get -y --force-yes purge gotoblas2", intern=TRUE, ignore.stderr=TRUE)
    invisible(res)
}

installAtlas39 <- function() {
    res <- system("sudo wajig -y install /var/spool/atlas39/libatlas39_3.9.25-1_amd64.deb", intern=TRUE, ignore.stderr=TRUE)
    invisible(res)
}

purgeAtlas39 <- function() {
    res <- system("sudo apt-get -y --force-yes purge libatlas39", intern=TRUE, ignore.stderr=TRUE)
    invisible(res)
}

createDatabase <- function(dbfile=.dbfile) {

    dbi <- dbDriver("SQLite")
    dframe <- data.frame(host="", datum="",
                         type="", nobs=NA, nrun=NA,
                         ref=NA, atlas=NA, atl39=NA,
                         mkl=NA, gotob=NA, gpu=NA)
    dftypes <- list(host="text", datum="text", type="text",
                    nobs="integer", nrun="integer",
                    ref="real", atlas="real", atl39="real",
                    mkl="real", gotob="real", gpu="real")
    sql <- dbBuildTableDefinition(dbi, name="benchmark", value=dframe,
                                  field.types=dftypes, row.names=FALSE)
    dbcon <- dbConnect(dbDriver("SQLite"), dbname=dbfile)
    res <- dbGetQuery(dbcon, sql)
    dbDisconnect(dbcon)
}

databaseResult <- function(data,dbfile=.dbfile) {

    newdf <- data.frame(host=Sys.info()["nodename"],
                        datum=format(Sys.Date()),
                        data)
    dbcon <- dbConnect(dbDriver("SQLite"), dbname=.dbfile)
    res <- dbWriteTable(dbcon, "benchmark", newdf, row.names=FALSE, overwrite=FALSE, append=TRUE)
    dbDisconnect(dbcon)

}

isPackageInstalled <- function(package) { 	# Henrik Bengtsson, r-devel, 24 Aug 2010
    path <- system.file(package=package)
    (path != "")
}

#hasMagma <- function() {
#    isPackageInstalled("magma")
#}

hasGputools <- function() {
    isPackageInstalled("gputools")
}

getBenchmarkData <- function(host) {
    dbcon <- dbConnect(dbDriver("SQLite"), dbname=system.file("sql", "gcbd.sqlite", package="gcbd"))
    data <- dbGetQuery(dbcon, paste('select * from benchmark where host="',
                                    host, '" order by nobs', sep=""))
    invisible(dbDisconnect(dbcon))
    invisible(data)
}
