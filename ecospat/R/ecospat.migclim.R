ecospat.migclim<-function()
{
message("load the MigClim package")
  packageStartupMessage("initializing ...")
  Sys.sleep(1)
requireNamespace("MigClim")

Sys.sleep(1)
packageStartupMessage("... done")

}