.onUnload = function(libpath)
{
	#cat("Unloading library", libpath, "...\n");
	library.dynam.unload("glmmGS", libpath);
}
