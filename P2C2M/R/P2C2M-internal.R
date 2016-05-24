.onAttach <-
function(libname = find.package("P2C2M"), pkgname = "P2C2M") {

  # Check if Richard Hudson's `ms` is present in compiled form
  msdir = system.file("msdir", "", package="P2C2M")
  if (!file.exists(paste(msdir, "ms", sep=""))) {
    #packageStartupMessage("-- Richard Hudson's 'ms' (Hudson 2002) is being compiled --")
    system(paste("cd", msdir, ";", "gcc -o ms ms.c streec.c rand1.c -lm"))
    #packageStartupMessage(" DONE -- \n\n")
  }

}
.onLoad <-
function(libname = find.package("P2C2M"), pkgname = "P2C2M") {
    #packageStartupMessage("-- Setting up a package-specific environment --")
    assign("P2C2M_globalVars", new.env(), envir=parent.env(environment()))
}
