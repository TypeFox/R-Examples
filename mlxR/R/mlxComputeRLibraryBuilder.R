#' @importFrom utils install.packages
mlxComputeRLibraryBuilder <- function(lixoftHOME){
  myOS <- Sys.info()['sysname']
  
  # STEP 02: Test if mlxComputeR.so / mlxComputeR.dll already exists ?
  so <- if( myOS == "Windows" ) "dll" else "so"
  mlxComputeFileName <- sprintf("%s/lib/mlxComputeR.%s", lixoftHOME, so)
  
  if (!file.exists(mlxComputeFileName)){
    # STEP 03: If mlxComputeR.so / mlxComputeR.dll does not exist it shall be created
    # STEP 03.1: Create file Makevars (containing compilation directives) in directory mlxLibrary/.../runtime/lib/src
    makeVarsFile <- if (myOS == "Windows"){
      sprintf("%s/lib/mlxComputeR/src/Makevars.win", lixoftHOME)  
    } else {
      sprintf("%s/lib/mlxComputeR/src/Makevars", lixoftHOME)  
    }
    
    fileConn <- file(makeVarsFile, open="w")
    
    writeLines("PKG_CPPFLAGS +=  -DCONNECTOR_R", con = fileConn, sep = "\n", useBytes = FALSE)
    
    writeLines("PKG_LIBS=`Rscript -e \"Rcpp:::LdFlags()\"`", con = fileConn, sep = "\n", useBytes = FALSE)
    
    if (myOS == "Linux"){
      writeLines("PKG_LIBS += -Wl,-rpath,\"$$\"\"ORIGIN\"", con = fileConn, sep = "\n", useBytes = FALSE)
    }
    if (myOS =="Darwin"){
      dirWheremlxComputeRIsInstalled <- sprintf("%s/lib", lixoftHOME)
      pkglibarg <- sprintf("PKG_LIBS += -Wl,-rpath,%s",dirWheremlxComputeRIsInstalled)
      writeLines(pkglibarg, con = fileConn, sep = "\n", useBytes = FALSE)
    }
    
    lineToAppend = sprintf("PKG_LIBS += -L\"%s/lib\" -lmlxCompute_CAPI", lixoftHOME)
    writeLines(lineToAppend, con = fileConn, sep = "\n", useBytes = FALSE)
    close(fileConn)
    
    # STEP 03.2: Define tempCompilationDirectory = mlxLibrary/.../runtime/lib/temp 
    #            the directory where the library will be compiled
    tempCompilationDirectory = sprintf("%s/lib/temp", lixoftHOME);
    dir.create(tempCompilationDirectory)
    
    # STEP 03.3: Launch compilation
    pckName = sprintf("%s/lib/mlxComputeR", lixoftHOME);
    libArgument = sprintf("--library=%s", tempCompilationDirectory);
    install.packages(pckName, tempCompilationDirectory, repos = NULL,INSTALL_opts = c("--no-multiarch", "--no-test-load"), type="source");
    
    # STEP 03.4: Copy mlxComputeR.so / mlxComputeR.dll from mlxLibrary/.../runtime/lib/temp/mlxComputeR/libs to  mlxLibrary/.../runtime/lib
    if (myOS == "Windows"){
      arch <- if( Sys.info()['machine'] == "x86-64" ) "x64" else "i386" 
      fromFile <- sprintf("%s/mlxComputeR/libs/%s/mlxComputeR.dll", tempCompilationDirectory, arch)
      toFile <- sprintf("%s/lib/mlxComputeR.dll", lixoftHOME)
    } else {
      fromFile <- sprintf("%s/mlxComputeR/libs/mlxComputeR.so", tempCompilationDirectory);
      toFile <- sprintf("%s/lib/mlxComputeR.so", lixoftHOME);
    }
    file.rename(fromFile, toFile)
    
    # STEP 03.5: Remove directory tempCompilationDirectory
    unlink(tempCompilationDirectory, recursive = TRUE)
    
    # STEP 03.6: Remove object and lib files created in mlxLibrary/.../runtime/lib/mlxComputeR/src
    dirToClean <- sprintf("%s/lib/mlxComputeR/src", lixoftHOME)
    listFilesToErase <- list.files(path = dirToClean, pattern = "[.](o|so|dll)$", full.names = TRUE);
    unlink(listFilesToErase)
    
  }
  
  if (myOS == "Windows" ){ 
    myOldENVPATH = Sys.getenv('PATH');
    #myNewENVPATH = sprintf("%s;%s/../tools/MinGW/bin;%s/tools/MinGW/bin", myOldENVPATH, lixoftHOME,lixoftHOME);
    myNewENVPATH = sprintf("%s/../tools/MinGW/bin;%s/tools/MinGW/bin;%s", lixoftHOME,lixoftHOME,myOldENVPATH);
    
    Sys.setenv('PATH'=myNewENVPATH)
    dirWheremlxComputeRIsInstalled = sprintf("%s/lib", lixoftHOME)
    owd <- setwd(dirWheremlxComputeRIsInstalled)
    dyn.load(mlxComputeFileName)
    setwd(owd)        
  } else {
    dyn.load(mlxComputeFileName)
  }

}
