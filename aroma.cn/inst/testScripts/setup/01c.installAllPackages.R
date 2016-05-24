library("R.utils");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


verbose && enter(verbose, "Installing additional packages for all tests");

verbose && enter(verbose, "Scanning tests");
path <- system.file("testScripts", package="aroma.cn");
verbose && cat(verbose, "Path to testScripts: ", path);

pattern <- "installPackages.R$";
pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE, recursive=TRUE);

nbrOfTestSets <- length(pathnames);
for (kk in seq(length=nbrOfTestSets)) {
  pathname <- pathnames[kk];
  dropPattern <- sprintf("^%s/", path);
  testSet <- gsub(dropPattern, "", dirname(pathname));
  verbose && enter(verbose, sprintf("Test set #%s ('%s') of %s", kk, testSet, nbrOfTestSets));
  sourceTo(pathname, envir=new.env(), echo=as.logical(verbose));
  verbose && exit(verbose);
} # for (kk ...)

verbose && exit(verbose);


verbose && exit(verbose);
