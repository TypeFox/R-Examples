.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is vegdata ",
    utils::packageDescription("vegdata", field="Version"),
    appendLF = TRUE)
#    packageStartupMessage("created:", file.info('/home/jansen/R/x86_64-pc-linux-gnu-library/3.1/vegdata/Meta/package.rds')$ctime
#)
}
