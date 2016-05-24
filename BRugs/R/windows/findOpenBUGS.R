findOpenBUGS <- function()
{
    dir <- Sys.getenv("OpenBUGS_PATH")
    if(!nchar(dir)){
        deps <- packageDescription("BRugs", fields="SystemRequirements")
        version.req <- gsub(".*OpenBUGS ?\\(>= ?(.+)\\).*", "\\1", deps)

        ob.reg <- try(readRegistry("Software\\OpenBUGS", "HLM", view = "32-bit"), silent = TRUE)
        if (inherits(ob.reg, "try-error")) {
            warning("OpenBUGS ", version.req, " or greater must be installed\n(if so, this indicates missing registry keys of OpenBUGS).\nSetting the environment variable 'OpenBUGS_PATH' in advance of loading 'BRugs' overwrites the path.\nSee ?loadOpenBUGS in order to load OpenBUGS manually.")
            return()
        }
        rnames <- names(ob.reg)
        version.full <- gsub("OpenBUGS ", "", rnames)
        ## remove suffixes from development versions, converts e.g. 3.2.1alpha to 3.2.1
        version.inst <- gsub("(.+[0-9]+)[a-zA-Z]+$","\\1", version.full)

        if(length(version.inst > 1)){
            id <- which(apply(outer(version.inst, version.inst, Vectorize(compareVersion, c("a", "b"))), 1, function(x) all(x >= 0)))
            id <- max(id) # if more than one release with same number, arbitrarily choose last one in registry
            version.inst <- version.inst[id]
            version.full <- version.full[id]
            rnames <- rnames[id]
        }

        if (compareVersion(version.inst, version.req) < 0) {
            warning("Found OpenBUGS version ", version.inst, ".\n Requires ", version.req, " or greater.\nSetting the environment variable 'OpenBUGS_PATH' in advance of loading 'BRugs' overwrites the path.\nSee ?loadOpenBUGS in order to load OpenBUGS manually.")
            return()
        }

        ## OpenBUGS installation location
        dir <- readRegistry(paste("Software","OpenBUGS",rnames,sep="\\"), "HLM", view = "32-bit")[["InstallPath"]]
    } else {
        if(!file.exists(file.path(dir, "libOpenBUGS.dll"))){
            warning("Environment variable OpenBUGS_PATH found but cannot access ", file.path(dir, "libOpenBUGS.dll"))
            return()
        }
        version.inst <- version.full <- NA
    }
    list(dir=dir, version=version.full)
}
