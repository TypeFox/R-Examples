
# Compare sets of student answers from submitted files

# Source() the code from a file and return the values as
# a list;  useful for inspecting student answers
# (without polluting global workspace)
# Only the sought-after objects are returned
sourceFile <- function(filename, modelNames) {
    tempenv <- new.env()
    try(with(tempenv, source(filename, local=TRUE)))
    templist <- as.list(tempenv)
    templist[names(templist) %in% modelNames]
}                       
                       
# The result is a list of "comparison" objects,
# one per model answer
compareFile <- function(filename, modelNames,
                        modelCode=NULL, modelSave=NULL,
                        modelAnswers=NULL,
                        # round may be a list with a different
                        # setting per model answer
                        round=FALSE,
                        ...) {
    # If modelCode and modelSave are both NULL then
    # the assumption is that the
    # model answer objects are either provided or
    # are already in the global workspace
    if (!is.null(modelCode))
        # Allow modelCode to be vector of file names
        for (FILENAME in modelCode)
            source(FILENAME)
    if (!is.null(modelSave))
        for (FILENAME in modelSave)
            load(FILENAME, envir=.GlobalEnv)
    # If modelAnswers not supplied, find them in the global workspace
    if (is.null(modelAnswers))
        modelAnswers <- lapply(modelNames, get,
                               envir=.GlobalEnv, inherits=FALSE)
    # Use local() so as not to despoil global env 
    tempenv <- new.env()
    # Some students' code will not even run
    # NOTE that we keep going though, because maybe they got
    # some answers right before the crash
    result <- try(with(tempenv, source(filename, local=TRUE)))
    # Set up 'round' argument
    if (is.list(round)) {
        # FALSE by default
        roundList <- as.list(rep(FALSE, length(modelNames)))
        names(roundList) <- modelNames
        # Use any named values in 'round'
        namesInCommon <- intersect(names(round), modelNames)
        if (length(namesInCommon) > 0) {
            roundList[namesInCommon] <- round[namesInCommon]
        }
    } else {
        if (is.function(round)) {
            roundList <- vector("list", length(modelNames))
            for (i in seq_along(roundList)) {
                roundList[[i]] <- round
            }
        } else {
            roundList <- as.list(rep(round, length.out=length(modelNames)))
        }
    }
    # Compare with model answers
    results <- mapply(compareName, modelAnswers, as.list(modelNames),
                      round=roundList,
                      MoreArgs=list(compEnv=tempenv, ...),
                      # Force a list result
                      SIMPLIFY=FALSE)
    # FIXME: may need to add some "clean up" code.
    # For example, detach() any attach()es that the source()ed code
    # performed (by detach()ing until the search path is the same
    # length as it was before the source() ?)
    names(results) <- modelNames
    class(results) <- "comparisonList"
    results
}

print.comparisonList <- function(x, ...) {
    print(unclass(x))
}

# The result is a list of lists, one list per file
# Each element of the list is named after the
# corresponding file, unless 'resultNames' is specified
compareFiles <- function(filenames, modelNames,
                         modelCode=NULL, modelSave=NULL,
                         resultNames=filenames, ...) {
    # If modelCode and modelSave are both NULL then
    # the assumption is that the
    # model answer objects are already in the global workspace
    if (!is.null(modelCode))
        # Allow modelCode to be vector of file names
        for (FILENAME in modelCode)
            source(FILENAME)
    if (!is.null(modelSave))
        for (FILENAME in modelSave)
            load(FILENAME)
    # Stick the model answer objects in a list
    # This does two things:
    # (i) keeps individual model answer objects from being
    #     overwritten by assignments to the global workspace
    # (ii) makes it easy to mapply() over all model answers
    modelAnswers <- lapply(modelNames, get,
                           envir=.GlobalEnv, inherits=FALSE)
    results <- lapply(filenames, compareFile, modelNames,
                      modelAnswers=modelAnswers, ...)
    names(results) <- resultNames
    class(results) <- "comparisonListList"
    results
}

print.comparisonListList <- function(x, falseDetails=TRUE, ...) {
    # Collapse to a character matrix combining results
    # plus transforms (if result is TRUE)
    # First step is to find out how many comparisons on each row
    numResults <- max(sapply(x, length))
    report <- do.call("rbind",
                      lapply(x,
                             function(y) {
                                 if (is.null(y)) {
                                     rep("", numResults)
                                 } else {
                                     sapply(y,
                                            function(z) {
                                                paste(z$result,
                                                      if (z$result ||
                                                          falseDetails) {
                                                          paste(z$transform,
                                                                collapse=" ")
                                                      } else {
                                                          ""
                                                      },
                                                      collapse=" ")
                                            })
                                 }
                             }))
    print(report, quote=FALSE)
}

# Utility function to try to clean up a file if the student
# cut-and-pastes an R session, rather than just saving R code.

tidyFile <- function(filename) {
    lines <- readLines(filename)
    # If any lines start with the R prompt, '>'
    if (length(grep("^>", lines)) > 0) {
        # Back up the original file
        file.rename(filename, paste(filename, ".BAK", sep=""))
        # Only keep lines starting with R prompt or
        # continuation prompt
        keptLines <- grep("^(>|[+])", lines)
        # Strip the prompts
        newLines <- gsub("^(>|[+])", "", lines[keptLines])
        writeLines(newLines, filename)
    }
}
