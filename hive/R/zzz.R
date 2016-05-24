## hadoop root directory has to be set

hive <- local({
               .henv <- .hive_default_env()
               function(new){
                 if(missing(new))
                   .henv
                 else
                   .henv <<- new
               }}
              )

.onLoad <- function(libname, pkgname){
    ## initialize hive environment
    hive(.hinit())
    .jpackage(pkgname, lib.loc = libname)
    if( is.environment(hive()) )
    {
        if(!hive_start(hive()))
            warning("Hadoop home exists but no Hadoop cluster was started.")
    }
}

