### This function checks whether a package is installed;
### if not, it installs it. It then loads the package.
safeRequire <- function(packageName, mirrorIndex=NULL) {
  if (!is.element(packageName, installed.packages()[,1])) {
    if (!is.null(mirrorIndex)) {
      chooseCRANmirror(ind=mirrorIndex);
    }
    install.packages(packageName, dependencies=TRUE);
  }
  suppressPackageStartupMessages(require(package = packageName,
                                         character.only=TRUE,
                                         quietly=TRUE));
}