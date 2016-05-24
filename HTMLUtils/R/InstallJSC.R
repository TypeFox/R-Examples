InstallJSC <- function#installs the JS components
###prompts the user to install the JS components to the relevant directory, which enables dynamically sortable tables.
(JSCPATH ##<< path to install the jsc directory to. Recommended is the base public html directory.
){
  ## David changed this because wasn't working on beast
  cat(JSCPATH, " does not exist. Boldly copying that directory to the base html directory.");
  system(paste("mkdir ", JSCPATH));
  jscDIR <- paste(system.file(package = "HTMLUtils"), "/jsc/",sep ="")
  system(paste("cp -r ", jscDIR,"* ", JSCPATH,"/",sep=""))

  if (1==0) {
    ## David changed this because wasn't working on beast
    cat(JSCPATH, " does not exist, we recommend that we copy the relevant files, proceed ? (y/n):");
    ans <- readLines(con=file("stdin"),n=1);
    if (ans == "y") {
      system(paste("mkdir ", JSCPATH));
      jscDIR <- paste(system.file(package = "HTMLUtils"), "/jsc/",sep ="")
      system(paste("cp -r ", jscDIR,"* ", JSCPATH,"/",sep=""))
    } else {
      cat("The dynamically sortable tables will likely not work ...\n");
    }  
  }  
}
