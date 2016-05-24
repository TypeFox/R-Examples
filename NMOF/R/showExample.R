showChapterNames <- function()
    Chapters

showExample <- function(file = "", chapter = NULL,
                        showfile = TRUE, includepaths= FALSE, ...) {

    ChapterDirs <- c("Introduction",
                     "NumAnNutshell",
                     "LinEqsLSP",
                     "FiniteDifferences",
                     "BinomialTrees",
                     "RandomNumberGeneration",
                     "ModelingDependencies",
                     "FinancialSimulations",
                     "CaseStudies",
                     "OptProbFinance",
                     "BasicMethods",
                     "HeuristicsNutshell",
                     "PortfolioOptimization",
                     "EconometricModels",
                     "OptionCalibration")

    path <- system.file(package = "NMOF")
    fpaths <- list.files(paste(path, "/book", sep = ""),
                         recursive = TRUE, full.names = TRUE)

    ## create file names and chapternames
    fnames <- gsub(".*/R/", "", fpaths, ignore.case = TRUE)
    fnames <- gsub(".*ChangeLog.*", "ChangeLog", fnames,
                   ignore.case = TRUE)
    chnames <- gsub(".*/C-([a-zA-Z]+[^/])/.*", "\\1",
                    fpaths, ignore.case = TRUE)
    chnames <- gsub(".*ChangeLog.*", "none", chnames, ignore.case = TRUE)

    filematch <- grepl(file, fnames, ...)

    if (is.null(chapter))
        chapmatch <- rep.int(TRUE, length(fnames))
    else {
        if (is.numeric(chapter))
            tmp <- ChapterDirs[chapter[1]]
        else if (is.character(chapter)) {
            chapmatch <- grepl(chapter, Chapters, ...)
            tmp <- ChapterDirs[chapmatch]
        }
                    chapmatch <- logical(length(fnames))
            for (i in seq_along(tmp))
                chapmatch <- (chapmatch | tmp[i] == chnames)

    }
    results <- chapmatch & filematch
    if (!length(which(results))) {
        message("no matches")
        flist <- data.frame(Chapter = character(0L),
                            File = character(0L))
    } else {
        flist <- data.frame(Chapter = chnames[results],
                            File = fnames[results])
        if (length(which(results)) == 1L)
            file.show(fpaths[results], title = "NMOF",
                      header = fnames[results])
        else
            message("found several files...")
    }
    if (includepaths)
        flist <- cbind(flist, Paths = fpaths[results])
    flist
}
