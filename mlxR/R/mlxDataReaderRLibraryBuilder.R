#' @importFrom utils install.packages installed.packages
mlxDataReaderRLibraryBuilder <- function (lixoftHOME, myCallingPath)
{

myOS <- Sys.info()['sysname'];
lixoftHOME=normalizePath(lixoftHOME);

#--- install "Rcpp" or check "Rcpp" version
if (!is.element("Rcpp", installed.packages()[,1])){
stop("Please install Rcpp \n",call.="FALSE")
}

# --- Check the Rcpp version installed is the correct one
 #packinfo <- installed.packages ();
  #installedRcppVersion = packinfo[c("Rcpp"), c("Version")]
  #requiredRcppVersion = "0.11.0";

  myOS <- Sys.info()['sysname'];
  #if (myOS == "Windows") {
  #	if (installedRcppVersion < requiredRcppVersion)
   # 		stop("Either install version 0.11.0 or delete file runtime/lib/mlxDataReaderR.dll, install Rtools and recompile") }



#STEP 02: Test if mlxDataReaderR.so / mlxDataReaderR.dll already exists ?
mlxDataReaderFileName = sprintf("%s/lib/mlxDataReaderR.so", lixoftHOME);
if (myOS == "Windows")
{
	mlxDataReaderFileName = sprintf("%s/lib/mlxDataReaderR.dll", lixoftHOME);
}
bFileExists = file.exists(mlxDataReaderFileName)

if (!bFileExists)
{
	# STEP 03: If mlxDataReaderR.so / mlxDataReaderR.dll does not exist it shall be created
	# STEP 03.1: Create file Makevars (containing compilation directives) in directory mlxLibrary/.../runtime/lib/src
	makeVarsFile = sprintf("%s/lib/mlxDataReaderR/src/Makevars", lixoftHOME);
	if (myOS == "Windows")
	{
		makeVarsFile = sprintf("%s/lib/mlxDataReaderR/src/Makevars.win", lixoftHOME);
	}
	fileConn<-file(makeVarsFile, open="w");
	writeLines("PKG_CPPFLAGS += -DCONNECTOR_R", con = fileConn, sep = "\n", useBytes = FALSE);
	writeLines("PKG_LIBS=`Rscript -e \"Rcpp:::LdFlags()\"`", con = fileConn, sep = "\n", useBytes = FALSE);
	if (myOS == "Linux")
	{
		writeLines("PKG_LIBS += -Wl,-rpath,\"$$\"\"ORIGIN\"", con = fileConn, sep = "\n", useBytes = FALSE);
	}
	if (myOS =="Darwin")
	{
	dirWheremlxDataReaderRIsInstalled = sprintf("%s/lib", lixoftHOME);
	pkglibarg = sprintf("PKG_LIBS += -Wl,-rpath,%s",dirWheremlxDataReaderRIsInstalled);
	writeLines(pkglibarg, con = fileConn, sep = "\n", useBytes = FALSE);
	#	writeLines("PKG_LIBS += -Wl,-rpath,\"@\"\"rpath\"", con = fileConn, sep = "\n", useBytes = FALSE);
	}

	lineToAppend = sprintf("PKG_LIBS += -L\"%s/lib\" -lmlxDataReader_CAPI", lixoftHOME);
	writeLines(lineToAppend, con = fileConn, sep = "\n", useBytes = FALSE);
	close(fileConn);

	# STEP 03.2: Define tempCompilationDirectory = mlxLibrary/.../runtime/lib/temp the directory where the library will be compiled
	tempCompilationDirectory = sprintf("%s/lib/temp", lixoftHOME);
	dir.create(tempCompilationDirectory)

	# STEP 03.3: Launch compilation
	pckName = sprintf("%s/lib/mlxDataReaderR", lixoftHOME);
	libArgument = sprintf("--library=%s", tempCompilationDirectory);
	#install.packages(pckName, tempCompilationDirectory, repos = NULL,INSTALL_opts = c("--no-multiarch", "--no-test-load", libArgument), type="source");
	install.packages(pckName, tempCompilationDirectory, repos = NULL,INSTALL_opts = c("--no-multiarch", "--no-test-load"), type="source");

	# STEP 03.4: Copy mlxDataReaderR.so / mlxDataReaderR.dll from mlxLibrary/.../runtime/lib/temp/mlxDataReaderR/libs to  mlxLibrary/.../runtime/lib
	fromFile = sprintf("%s/mlxDataReaderR/libs/mlxDataReaderR.so", tempCompilationDirectory);
	toFile = sprintf("%s/lib/mlxDataReaderR.so", lixoftHOME);
	if (myOS == "Windows")
	{
		if (Sys.info()['machine'] == "x86-64")
		{
			fromFile = sprintf("%s/mlxDataReaderR/libs/x64/mlxDataReaderR.dll", tempCompilationDirectory);
		} else {
			fromFile = sprintf("%s/mlxDataReaderR/libs/i386/mlxDataReaderR.dll", tempCompilationDirectory);
		}
		toFile = sprintf("%s/lib/mlxDataReaderR.dll", lixoftHOME);
	}
	file.rename(fromFile, toFile)

	# STEP 03.5: Remove directory tempCompilationDirectory
	unlink(tempCompilationDirectory, recursive = TRUE)

	# STEP 03.6: Remove object and lib files created in mlxLibrary/.../runtime/lib/mlxDataReaderR/src
	dirToClean = sprintf("%s/lib/mlxDataReaderR/src", lixoftHOME);
	listFilesToErase = list.files(path = dirToClean, pattern = "([.]o){1}", full.names = TRUE);
	for (iFile in 1:length(listFilesToErase))
	{
		file.remove(listFilesToErase[iFile])
	}
	listFilesToErase = list.files(path = dirToClean, pattern = "([.]so){1}", full.names = TRUE);
	if (myOS == "Windows")
	{
		listFilesToErase = list.files(path = dirToClean, pattern = "([.]dll){1}", full.names = TRUE);
	}
	for (iFile in 1:length(listFilesToErase))
	{
		file.remove(listFilesToErase[iFile])
	}
}

# STEP 04: load Rcpp and lxLibrary/.../runtime/lib/mlxDataReaderR.so(.dll)
#library("Rcpp");


if (myOS == "Windows" )
{ 
	myOldENVPATH = Sys.getenv('PATH');
    #    myNewENVPATH = sprintf("%s;%s/../tools/MinGW/bin;%s/tools/MinGW/bin",myOldENVPATH,lixoftHOME,lixoftHOME);
	myNewENVPATH = sprintf("%s/../tools/MinGW/bin;%s/tools/MinGW/bin;%s",lixoftHOME,lixoftHOME,myOldENVPATH);
  Sys.setenv('PATH'=myNewENVPATH);
	dirWheremlxDataReaderRIsInstalled = sprintf("%s/lib", lixoftHOME);
	myCallingPath<-setwd(dirWheremlxDataReaderRIsInstalled);
	dyn.load(mlxDataReaderFileName);
	setwd(myCallingPath);

}

 else {
	dyn.load(mlxDataReaderFileName)
}


}
