`print.how` <- function(x, ...) {

    ## only for objects of correct class
    stopifnot(class(x) == "how")

    ## prefix to add to sub-elements
    pfix <- "  "

    cat("\n")
    writeLines(strwrap("Permutation Design:"))
    cat("\n")

    ## Blocks
    writeLines("Blocks:")
    blocks <- getBlocks(x)
    if (is.null(blocks)) {
        writeLines(strwrap("Defined by: none", prefix = pfix))
    } else {
        writeLines(strwrap(paste("Blocks:", x$blocks.name),
                           prefix = pfix))
    }

    cat("\n")

    ## Plots
    writeLines("Plots:")
    plotStr <- getStrata(x, which = "plots")
    plots <- getPlots(x)
    ptype <- getType(x, which = "plots")
    if (is.null(plotStr)) {
        writeLines(strwrap("Defined by: none", prefix = pfix))
    } else {
        writeLines(strwrap(paste("Plots:", plots$plots.name),
                           prefix = pfix))
        writeLines(strwrap(paste("Permutation type:", ptype),
                           prefix = pfix))
        mirrorP <- getMirror(x, which = "plots")
        writeLines(strwrap(paste("Mirrored?:", if(mirrorP) "Yes" else "No"),
                           prefix = pfix))
        if(isTRUE(all.equal(ptype, "grid"))) {
            nr <- getRow(x, which = "plots")
            nr.t <- if(nr > 1) "rows" else "row"
            nc <- getCol(x, which = "plots")
            nc.t <- if(nc > 1) "cols" else "col"
            writeLines(strwrap(paste("Grid dimensions:", nr, nr.t, " ",
                                     nc, nc.t),
                               prefix = pfix))
        }
    }

    cat("\n")

    ## Within plots
    writeLines("Within Plots:")
    wtype <- getType(x, which = "within")
    writeLines(strwrap(paste("Permutation type:", wtype), prefix = pfix))
    mirrorW <- getMirror(x, which = "within")
    constantW <- getConstant(x)
    txt <- "Different permutation within each Plot?:"
    if(isTRUE(wtype %in% c("series", "grid"))) {
        writeLines(strwrap(paste("Mirrored?:", if(mirrorW) "Yes" else "No"),
                           prefix = pfix))
        writeLines(strwrap(paste(txt, if(constantW) "No" else "Yes"),
                           prefix = pfix))
    }
    if(isTRUE(all.equal(wtype, "grid"))) {
        nr <- getRow(x, which = "within")
        nr.t <- if(nr > 1) "rows" else "row"
        nc <- getCol(x, which = "within")
        nc.t <- if(nc > 1) "cols" else "col"
        writeLines(strwrap(paste("Grid dimensions:", nr, nr.t, " ",
                                 nc, nc.t),
                           prefix = pfix))
    }

    cat("\n")

    ## Meta data
    writeLines("Permutation details:")
    writeLines(strwrap(paste("Number of permutations:",
                             getNperm(x)), prefix = pfix))
    writeLines(strwrap(paste("Max. number of permutations allowed:",
                             getMaxperm(x)), prefix = pfix))
    txt <- paste("Evaluate all permutations?:",
                 if(getComplete(x)) "Yes." else "No.",
                 "  Activation limit:", getMinperm(x))
    writeLines(strwrap(txt, prefix = pfix))

}
