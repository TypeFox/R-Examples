## ----setup, include=FALSE------------------------------------------------
# set global chunk options
library(knitr); library(pacman); library(methods)
opts_chunk$set(cache=FALSE, comment=NA)

## Function for embedding high qual text images:
uri_embed <- function(path, add="") {
    uri <- knitr::image_uri(path)
    cat(paste0("<img ", add, " src=\"", uri, "\" />"))
}
opts_knit$set(upload.fun = image_uri, self.contained=TRUE)

## set mirror
options(repos="http://cran.rstudio.com/")

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
uri_embed("r_pacman.png", 
    "width=\"350\", height=\"150\" style=\"display:block; margin-left:auto; margin-right:auto;\"")

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/01_installing_loading_deleting.R")
cat(paste(installing_tab, collapse="\n"))

## ---- eval = FALSE-------------------------------------------------------
#  p_load(..., char, install = TRUE, update = getOption("pac_update"), character.only = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  packs <- c("XML", "devtools", "RCurl", "fakePackage", "SPSSemulate")
#  success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
#  install.packages(names(success)[!success])
#  sapply(names(success)[!success], require, character.only = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  pacman::p_load(XML, devtools, RCurl, fakePackage, SPSSemulate)

## ----eval=FALSE----------------------------------------------------------
#  p_install(dbConnect, qdap, reports)

## ----eval=FALSE----------------------------------------------------------
#  p_install_gh(c("Dasonk/githubSearch", "trinker/regexr", "hadley/httr@v0.4"))

## ----eval=FALSE----------------------------------------------------------
#  p_load_gh("Dasonk/githubSearch", "trinker/regexr", "hadley/httr@v0.4")

## ------------------------------------------------------------------------
p_install_version(
    c("pacman", "testthat"),
    c("0.2.0", "0.9.1")
)

## ----eval=FALSE----------------------------------------------------------
#  p_temp(aprof)
#  
#  p_isinstalled(aprof)
#  p_isloaded(aprof)

## ---- eval = FALSE-------------------------------------------------------
#  p_unload(..., negate = FALSE, char, character.only = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  p_load(lattice)
#  p_isloaded(lattice)
#  p_unload(lattice)
#  p_isloaded(lattice)

## ---- eval = FALSE-------------------------------------------------------
#  p_update()

## ---- eval = FALSE-------------------------------------------------------
#  p_update(FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  p_delete(fakePackage, stats)

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/02_session_information.R")
cat(paste(installing_tab, collapse="\n"))

## ---- eval = FALSE-------------------------------------------------------
#  p_loaded()

## ---- eval = FALSE-------------------------------------------------------
#  p_loaded(all = TRUE)

## ------------------------------------------------------------------------
p_loaded(base, MASS)
p_isloaded(methods, stats)

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/03_local_package_information.R")
cat(paste(installing_tab, collapse="\n"))

## ------------------------------------------------------------------------
p_exists(pacman)  
p_exists(pacman, local = TRUE)
p_exists(I_dont_exist)
## wrapper for `p_exists(local = TRUE)`
p_isinstalled(pacman)

## ------------------------------------------------------------------------
p_depends(lattice)  
p_depends(lattice, local = TRUE)
p_depends(MASS)  
p_depends(MASS, local = TRUE)

## ---- echo=FALSE, results='hide'-----------------------------------------
.pinfo <- names(p_info())
.right_paren <- c(rep("(", length(.pinfo) - 1 ), "and (")
.fields <- paste(.right_paren, letters[1:length(.pinfo)], ") ", .pinfo, sep = "", collapse = ", ")
.fields <- gsub("(c)", "&#40;c)", .fields, fixed = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  p_info(package, ..., fields = NULL)

## ------------------------------------------------------------------------
## Defaults to supply a list of fields with information about R's base package
p_info()
names(p_info())

## ------------------------------------------------------------------------
p_info(pacman, Author)
p_info(pacman, BugReports, URL)
p_info(pacman, fields = "Version")

## ------------------------------------------------------------------------
## without `p_extract`
p_info(MASS, "Depends")  
p_extract(p_info(MASS, "Depends"))
p_extract(p_info(methods, "Imports"))

## ------------------------------------------------------------------------
p_author(pacman)
p_author()

## ------------------------------------------------------------------------
p_cite(pacman)
p_citation()

## ------------------------------------------------------------------------
p_data(lattice)

## ------------------------------------------------------------------------
p_functions(pacman)
p_funs(pacman, all=TRUE)

## ------------------------------------------------------------------------
p_version()
p_ver(pacman)
p_ver(pacman) >= "0.2.0"

## ---- eval=FALSE---------------------------------------------------------
#  p_help(pacman)
#  p_help(pacman, web = FALSE)
#  p_help(pacman, build.pdf = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  p_news()
#  p_news(pacman)
#  ## Grab specfic version subsets
#  subset(p_news(lattice), Version == 0.7)

## ---- eval=FALSE---------------------------------------------------------
#  p_vignette()
#  p_vign(pacman)

## ---- eval=FALSE---------------------------------------------------------
#  p_interactive()

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/04_internet_based_info.R")
cat(paste(installing_tab, collapse="\n"))

## ---- eval=FALSE---------------------------------------------------------
#  p_cran()

## ------------------------------------------------------------------------
length(p_cran())
p_iscran("qdap")

## ---- eval=FALSE---------------------------------------------------------
#  p_sa("color", "package")
#  p_sa("hadley", "author")

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/05_library_information.R")
cat(paste(installing_tab, collapse="\n"))

## ------------------------------------------------------------------------
p_path()

## ---- eval=FALSE---------------------------------------------------------
#  p_lib()
#  p_base()

## ---- eval=FALSE---------------------------------------------------------
#  p_search_library(begins.with = NULL, contains = NULL)

## ---- eval=FALSE---------------------------------------------------------
#  p_sl("pa")
#  p_sl(contains = "man")
#  p_sl(begins.with ="pa", contains = "man")

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/06_pacman_tools.R")
cat(paste(installing_tab, collapse="\n"))

## ----css, echo = FALSE---------------------------------------------------
options(markdown.HTML.stylesheet = "css/style.css")

## ---- echo=FALSE, message=FALSE------------------------------------------
#write.bibtex(file="references.bib")

## ----clean-up, include=FALSE---------------------------------------------
# R compiles all vignettes in the same session, which can be bad
rm(list = ls(all = TRUE))

