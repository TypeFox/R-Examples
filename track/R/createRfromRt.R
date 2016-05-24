createRfromRt <- function(base, transcript.suffix="Rt", input.suffix="R", output.suffix="Rout.save", location=FALSE) {
    ## If <base>.R doesn't already exist, create it
    ## If <base>.Rout.save exists, use that, otherwise look for <base>.Rt
    ## (and in the latter case, also create <base>.Rout.save)
    if (!file.exists(paste(base, input.suffix, sep=".")) || regexpr("generated automatically", readLines(paste(base, input.suffix, sep="."), 1))) {
        if (file.exists(paste(base, output.suffix, sep="."))) {
            ## Read the the Rout.save file
            ## Note that the Rout.save file must have:
            ## * two extra blank lines at the beginning
            ## * command prompt at the end (empty, with a newline)
            lines <- readLines(paste(base, output.suffix, sep="."), -1)
            if (length(lines)>2 && all(regexpr("^[ \t]*$", lines[1:2])>0))
                lines <- lines[-(1:2)]
            if (length(lines)>1 && all(regexpr("^>[ \t]*$", lines[length(lines)])>0))
                lines <- lines[-length(lines)]
            lines <- paste(lines, "\n", sep="")
        } else if (file.exists(paste(base, transcript.suffix, sep="."))) {
            ## Read the Rt file and create the Rout.save file.
            ## Try to massage the commands so that they will exactly match -
            ## desired output must have:
            ## * two extra blank lines at the beginning
            ## * command prompt at the end (empty, with a newline)
            lines <- paste(readLines(paste(base, transcript.suffix, sep="."), -1), "\n", sep="")
            if (location)
                lines <- c(paste("> # generated automatically in", getwd(), "on", Sys.time(), "\n"), lines)
            else
                lines <- c("> # generated automatically\n", lines)
            cat(c("\n", "\n", lines, "> \n"), sep="", file=paste(base, output.suffix, sep="."))
        } else {
            stop("unable to find ", paste(base, transcript.suffix, sep="."), " or ",
                 paste(base, output.suffix, sep="."))
        }
        lines <- sub("^[+>] ?", "", grep("^[+>] ?", lines, value=TRUE))
        cat(lines, sep="", file=paste(base, input.suffix, sep="."))
    }
    invisible(NULL)
}
