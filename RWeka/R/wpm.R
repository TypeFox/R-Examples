WPM <-
function(cmd, ...)
{
    cmd <- cmd[1L]
    cmds <- c("refresh-cache", "list-packages", "install-package",
              "remove-package", "load-package", "package-info")
    pos <- pmatch(tolower(cmd), cmds)
    if(is.na(pos))
        stop(gettextf("Invalid package manager command '%s'.", cmd),
             domain = NA)
    cmd <- cmds[pos]

    args <- as.character(list(...))

    wpm <- .jnew("weka.core.WekaPackageManager")
    
    if(cmd == "load-package") {
        ## Need to write code ourselves ...
        arg <- args[1L]
        if(is.na(arg))
            stop(gettextf("No package given."),
                 domain = NA)
        ## Avoid repeated loads (for now).
        packages <- Weka_packages_loaded()
        if(!(arg %in% packages)) {
            dir <- .jcall(wpm, "Ljava/io/File;", "getPackageHome")
            dir <- file.path(.jcall(dir, "S", "toString"), arg)
            if(!file.exists(dir)) {
                stop(gettextf("Unavailable package '%s'.", arg),
                     domain = NA)
            }
            .jaddClassPath(Sys.glob(file.path(dir, "*.jar")))
            Weka_packages_loaded(c(packages, arg))
        }
        return(invisible())
    }

    ## * -refresh-cashe
    ## This is
    ##   public static java.lang.Exception
    ##   weka.core.WekaPackageManager.refreshCache(java.io.PrintStream[])
    ## We could try to capture progress into a suitable output stream,
    ## but for now simply allow writing to stdout ... sort of.
    ## Messy ... similar for the others.

    ## Note that as of 8982 the WekaPackageManager class main() method
    ## calls System.exit(0), so using main() is no longer feasible.  We
    ## currently patch the upstream Weka sources for RWekajars, but
    ## should really rewrite the WPM() code to use the appropriate
    ## WekaPackageManager methods directly, instead of main().

    if(cmd == "refresh-cache") {
        .jcall(wpm,
               "Ljava/lang/Exception;",
               "refreshCache",
               .jarray(.jfield("java/lang/System", , "out"),
                       "java/io/PrintStream"))
        return(invisible())
    }
    
    ## Let us use main() for now, but make sure we call it reasonably,
    ## as the Weka 3.7.2 release in fact does System.exit(1) in case of
    ## misuse ...

    switch(EXPR = cmd,
           "list-packages" = {
               arg <- c(args, "all")[1L]
               tab <- c("all", "installed", "available")
               pos <- pmatch(arg, tab)
               if(is.na(pos))
                   stop("Invalid package manager command '%s %s'",
                        cmd, arg)
               args <- c("-list-packages", arg)
           },
           "package-info" = {
               ## This is somewhat silly ...
               ## But we really need 2 arguments here.
               args <- args[c(1L, 2L)]
               tab <- c("repository", "installed", "archive")
               pos <- pmatch(args[1L], tab)
               if(is.na(pos) || is.na(args[2L]))
                   stop("Invalid package manager command '%s %s %s'",
                        cmd, arg[1L], arg[2L])
               args <- c("-package-info", args)
           },
           "install-package" = {
               arg <- args[1L]
               if(is.na(arg))
                   stop(gettextf("No package given."),
                        domain = NA)
               args <- c("-install-package", arg)
           },
           "remove-package" = {
               arg <- args[1L]
               if(is.na(arg))
                   stop(gettextf("No package given."),
                        domain = NA)
               args <- c("-remove-package", arg)
           })

    ## Capture Java output.
    bos <- .jnew("java/io/ByteArrayOutputStream")
    out <- .jfield("java/lang/System", , "out")
    .jcall("java/lang/System", "V", "setOut",
           .jnew("java/io/PrintStream",
                 .jcast(bos,"java/io/OutputStream")))
    err <- .jfield("java/lang/System", , "err")
    .jcall("java/lang/System", "V", "setErr",
           .jnew("java/io/PrintStream",
                 .jcast(bos,"java/io/OutputStream")))

    .jcall(wpm, "V", "main", .jarray(args))

    ## Stop redirecting Java messages.
    .jcall("java/lang/System", "V", "setOut", out)
    .jcall("java/lang/System", "V", "setErr", err)
    ## And display them.
    message(.jcall(bos, "Ljava/lang/String;", "toString"))
}

Weka_packages_loaded <-
local({
    packages <- character()
    function(new) {
        if(!missing(new))
            packages <<- unique(new)
        else
            packages
    }
})

make_Weka_package_loader <-
function(p)
    function() {
        tryCatch(WPM("load-package", p),
                 error =
                 function(e)
                 stop(gettextf("Required Weka package '%s' is not installed.",
                               p),
                      domain = NA))
                 
    }
