outgis <- function (partition)
{
    require(tcltk) || stop("The tcltk support is absent")
    if (is.null(class(partition)) | class(partition) != "nampartition") {
        cat("Argument is not of class 'nampartition' \n")
        return(invisible())
    }
    filename <- tclvalue(tkgetSaveFile(initialfile = "output_gis.txt",
                        defaultextension = ".txt", title = "Save marked point sets into TXT file...",
                        filetypes = "{txt {.txt}} {{All Files} {*.*}}"))
    if (filename == "") return()
    zz <- file(filename, "w")
    cat("LABEL", "X_LONG", "Y_LAT", "NAMCLASS\n", sep = ",", file = zz)
    for (i in 1:nrow(partition$status)) {
        sp <- as.character(partition$status[i,1])
        clas <- partition$status[i,2]
        pts <- partition$occupancy[[sp]]
        write.table(cbind(sp, partition$coords[pts, 1], partition$coords[pts, 2],
                  clas), row.names = FALSE, col.names = FALSE, quote = FALSE,
                  append = TRUE, sep = ",", file = zz)
    }
    close(zz)
}


