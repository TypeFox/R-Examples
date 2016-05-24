packagefile="gdata_2.16.0.tar.gz"
destdir=tempdir()

checkReverseDependencies <- function(packagefile, destdir=tempdir(), cleanup=FALSE )
    {
        if(!file.exists(packagefile))
            stop(packagefile, " does not exist!")

        cat("Using directory '", destdir, "'.  Remember to delete it when done.\n", sep='')

        file.copy(packagefile, destdir)

        package <- gsub("_.*$", "", packagefile)

        rdeps <- tools::package_dependencies(package, db=available.packages(), reverse = TRUE)[[1]]
        cat( length(rdeps), "reverse dependencies:\n")
        print(rdeps)

        tools::check_packages_in_dir(destdir, reverse=list(), Ncpus=6)

        if(cleanup) unlink(destdir, recursive=TRUE, force=TRUE)
    }
