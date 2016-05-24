#### doRUnit.R --- Run RUnit tests
####------------------------------------------------------------------------

### Originally follows Gregor Gojanc's example in CRAN package  'gdata'
### and the corresponding section in the R Wiki:
###  http://wiki.r-project.org/rwiki/doku.php?id=developers:runit

### MM: Vastly changed:  This should also be "runnable" for *installed*
##              package which has no ./tests/
## ----> put the bulk of the code e.g. in  ../inst/unitTests/runTests.R :

### DJS: Added the possibility of different testing levels
### environment variable LEVEL set by make
### default value of 1 is level used in check and install

if(require("RUnit", quietly=TRUE)) {

  ## --- Setup ---

  wd <- getwd()
  pkg <- sub("\\.Rcheck$", '', basename(dirname(wd)))
  level <- Sys.getenv("LEVEL")

  library(package=pkg, character.only=TRUE)

  path <- system.file("unitTests", package = pkg)

  stopifnot(file.exists(path), file.info(path.expand(path))$isdir)

  source(file.path(path, "runTests.R"), echo = TRUE)
}
