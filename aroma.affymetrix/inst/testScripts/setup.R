library("aroma.affymetrix");

pathT <- system.file("testScripts", package="aroma.affymetrix");
pathT <- Arguments$getReadablePath(pathT);

path <- "testScripts/";
path <- filePath(path, expandLinks="any");
if (!isDirectory(path)) {
  createLink(target=pathT);
}
path <- Arguments$getReadablePath(path);

pathname <- "menu.Rex";
if (!isFile(pathname)) {
  code <- c(
    'if (interactive()) savehistory();',
    'library("R.menu");',
    'launchMenu("testScripts");'
  );
  cat(file=pathname, code, sep="\n");
}
pathname <- Arguments$getReadablePathname(pathname);
