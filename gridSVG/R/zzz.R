
.onLoad <- function(libname, pkgname) {
    genAttrTable <- utils::read.table(system.file("svg", "genAttrTable.txt",
                                                  package="gridSVG"),
                                      col.names=c("attribute", "element"),
                                      stringsAsFactors=FALSE)
    assign("genAttrTable", genAttrTable, envir=.gridSVGEnv)
    presAttrTable <- utils::read.table(system.file("svg", "presAttrTable.txt",
                                                   package="gridSVG"),
                                       col.names=c("attribute", "element"),
                                       stringsAsFactors=FALSE)
    assign("presAttrTable", presAttrTable, envir=.gridSVGEnv)
    assign("SVGParList", unique(presAttrTable$attribute), envir=.gridSVGEnv)
}
