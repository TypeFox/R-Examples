outgearth <- function(partition, showlabel = FALSE)
{
    require(tcltk) || stop("The tcltk support is absent")
    if (is.null(class(partition)) | class(partition) != "nampartition") {
        cat("Argument is not of class 'nampartition' \n")
        return(invisible())
    }
    filename <- tclvalue(tkgetSaveFile(initialfile = "output_gearth.kml",
                        defaultextension = ".kml", title = "Save marked point sets into KML file...",
                        filetypes = "{KML {.kml}} {{All Files} {*.*}}"))
    if (filename == "") return()
    zz <- file(filename, "w")
    aux <- split(partition$status[,1], partition$status[,2])
    naux <- names(aux)
    cat("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n", file = zz)
    cat("<kml xmlns = \"http://earth.google.com/kml/2.1\">\n", file = zz)
    cat("<Document>\n", file = zz)
    cat("<name>", "NAM classification", "</name>\n", file = zz)
    a <- "<Placemark>\n<name>"
    b <- "</name>\n<Point>\n<coordinates>"
    c <- ", 0</coordinates>\n</Point>\n<description><![CDATA["
    d <- "]]></description>\n</Placemark>\n"
    for (i in 1:length(aux)) {
      cat("<Folder>\n", file = zz)
      cat("<name>", naux[i],"</name>",sep = "", file = zz)
      if(showlabel) label <- naux[i] else label <- ""
        for(j in aux[[i]]) {
          for(k in partition$occupancy[[j]]) {
            cat(a, label , b, partition$coords[k,1], ", ",
              partition$coords[k,2], c, j, d, sep = "", file = zz)
          }
        }
      cat("</Folder>\n", file = zz)
    }
    cat("</Document>\n</kml>", file = zz)
    close(zz)
}
