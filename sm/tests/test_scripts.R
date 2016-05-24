## Note: R CMD check may run these scripts from an installed package
scripts <- list.files(system.file("scripts", package = "sm"), ".*\\.q$")
## these are interactive
omit2 <- match(c("bissell3.q", "dogs.q"), scripts)
scripts <- scripts[-omit2]
library(sm)
if(.Platform$OS.type == "unix") options(pager="cat") else options(pager="console")
postscript(file="test_scripts.ps")
for(z in scripts) {
    cat("\n============ running script `", z, "' ============\n", sep="")
    set.seed(123)
    source(system.file("scripts", z, package = "sm"), echo=TRUE)
    rm(list = ls(all = TRUE))
}
