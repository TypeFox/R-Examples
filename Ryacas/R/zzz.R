
.onAttach <- function(lib, pkg) {

   if (.Platform$OS.type != "windows") return()

   yacas.exe <- yacasFile("yacas.exe")
   scripts.dat <- yacasFile("scripts.dat")

   if (!file.exists(yacas.exe) ||  !file.exists(scripts.dat)) {
      packageStartupMessage(yacas.exe, " and/or ", scripts.dat, " not found\n",
         "run yacasInstall() without arguments to installation yacas.")
   }

   invisible()

}
