## useDynLib() in ../NAMESPACE  loads the C code
.onLoad <- function(lib, pkg) {
    nms <- c("Version", "Date", "Built")
    assign("vlmc.version",
           paste(nms, utils::packageDescription("VLMC")[nms], sep=": "),
	   asNamespace("VLMC"))
}
