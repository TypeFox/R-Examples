## ---- echo = FALSE-------------------------------------------------------
# Grab package version and date from DESCRIPTION file
pkgVersion <- packageDescription("seroincidence")$Version
pkgDate <- packageDescription("seroincidence")$Date
pkgBaseFileName <- paste("seroincidence", pkgVersion, sep = "_")
pkgUrl <- gsub("\n", "", packageDescription("seroincidence")$URL)
pkgPath <- "http://ecdc.europa.eu/en/data-tools/seroincidence-calculator-tool/Documents"
library("knitr")

## ---- include=FALSE, tidy=FALSE------------------------------------------
chunkTemplate <- "```{r, eval = FALSE}
# OPTION A
# Install Windows binary package \"seroincidence\" directly from internet location:
install.packages(
    pkgs = \"{{pkgPath}}/
                {{pkgBaseFileName}}.zip\",
    repos = NULL, type = \"win.binary\")

# OPTION B
# Install source package directly from internet location:
install.packages(
    pkgs = \"{{pkgPath}}/
                {{pkgBaseFileName}}.tar.gz\",
    repos = NULL, type = \"source\")

# OPTION C
# Install Windows binary package from a local file:
#   install.packages(\"[PATH/TO/FILE/]seroincidence_[version].zip\", 
#                       repos = NULL, type = \"win.binary\")
# For instance:
install.packages(pkgs = \"C:/{{pkgBaseFileName}}.zip\", repos = NULL, type = \"win.binary\")

# OPTION D
# Install source package from a local file:
install.packages(pkgs = \"C:/{{pkgBaseFileName}}.tar.gz\", repos = NULL, type = \"source\")
```"

chunkSource <- knit_expand(text = chunkTemplate, pkgBaseFileName = pkgBaseFileName)

## ---- eval = FALSE-------------------------------------------------------
#  # OPTION A
#  # Install Windows binary package "seroincidence" directly from internet location:
#  install.packages(
#      pkgs = "http://ecdc.europa.eu/en/data-tools/seroincidence-calculator-tool/Documents/
#                  seroincidence_1.0.5.zip",
#      repos = NULL, type = "win.binary")
#  
#  # OPTION B
#  # Install source package directly from internet location:
#  install.packages(
#      pkgs = "http://ecdc.europa.eu/en/data-tools/seroincidence-calculator-tool/Documents/
#                  seroincidence_1.0.5.tar.gz",
#      repos = NULL, type = "source")
#  
#  # OPTION C
#  # Install Windows binary package from a local file:
#  #   install.packages("[PATH/TO/FILE/]seroincidence_[version].zip",
#  #                       repos = NULL, type = "win.binary")
#  # For instance:
#  install.packages(pkgs = "C:/seroincidence_1.0.5.zip", repos = NULL, type = "win.binary")
#  
#  # OPTION D
#  # Install source package from a local file:
#  install.packages(pkgs = "C:/seroincidence_1.0.5.tar.gz", repos = NULL, type = "source")

## ---- eval=FALSE---------------------------------------------------------
#  # Load package "seroincidence"
#  library(seroincidence)
#  
#  # Show R help for the package
#  ?seroincidence
#  
#  # Show tutorial for the package
#  vignette(topic = "tutorial", package = "seroincidence")

## ------------------------------------------------------------------------
# Description
packageDescription("seroincidence")

# Citation
citation("seroincidence")

