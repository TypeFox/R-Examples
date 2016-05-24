.onLoad <-
function(libname, pkgname) {
	library.dynam("maxent",pkgname,libname);
	assign("maximumentropy",Module("maximumentropy",PACKAGE="maxent"),envir=.BaseNamespaceEnv);
	setClass("maxent", representation(model = "character", weights = "data.frame"));
}