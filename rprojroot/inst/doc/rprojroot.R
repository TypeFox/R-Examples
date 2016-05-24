## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(error = TRUE)

## ------------------------------------------------------------------------
basename(getwd())

## ------------------------------------------------------------------------
library(rprojroot)

# List all files and directories below the root
dir(find_root("DESCRIPTION"))

## ------------------------------------------------------------------------
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("DESCRIPTION")))

# Find a file relative to the root
file.exists(find_root_file("R", "root.R", criterion = has_file("DESCRIPTION")))

## ------------------------------------------------------------------------
has_file("DESCRIPTION")

## ------------------------------------------------------------------------
as.root_criterion("DESCRIPTION")

## ------------------------------------------------------------------------
criteria

## ------------------------------------------------------------------------
has_license <- has_file("LICENSE")
has_license

is_projecttemplate_project <- has_file("config/global.dcf", "^version: ")
is_projecttemplate_project

## ------------------------------------------------------------------------
# Print first lines of the source for this document
head(readLines(find_package_root_file("vignettes", "rprojroot.Rmd")))

## ------------------------------------------------------------------------
P <- find_package_root_file

# Use a shorter alias
file.exists(P("vignettes", "rprojroot.Rmd"))

## ----error = TRUE--------------------------------------------------------
# Use the has_license criterion to find the root
R <- make_find_root_file(has_license)

# Our package does not have a LICENSE file, trying to find the root results in an error
R()

## ------------------------------------------------------------------------
# Define a function that computes file paths below the current root
F <- make_fix_root_file(is_r_package)

# Show contents of the NAMESPACE file in our project
readLines(F("NAMESPACE"))

## ------------------------------------------------------------------------
# Print the size of the namespace file, working directory outside the project
local({
  oldwd <- setwd("../..")
  on.exit(setwd(oldwd), add = TRUE)
  file.size(F("NAMESPACE"))
})

