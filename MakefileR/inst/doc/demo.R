## ---- echo = FALSE-------------------------------------------------------
library(magrittr)
knitr::opts_chunk$set(
  comment = "#>",
  fig.path = "README-"
)
knit_print.MakefileR <- function(x, options) {
  knitr::knit_print(
    c("``` Makefile",
      format(x) %>% gsub("\t", "⇥      ", .),
      "```") %>%
      paste(">", ., collapse = "\n") %>%
      knitr::asis_output(), options)
}

## ------------------------------------------------------------------------
library(MakefileR)
make_rule("all", c("first_target", "second_target"))

make_rule(".FORCE")

make_rule("first_target", ".FORCE", "echo 'Building first target'")

make_rule("second_target", "first_target",
  c("echo 'Building second target'", "echo 'Done'"))

## ------------------------------------------------------------------------
make_def("R_VERSION", R.version.string)

## ------------------------------------------------------------------------
make_group(make_comment("Definitions")) +
  make_def("R_VERSION", R.version.string) +
  make_def("R_PLATFORM", R.version$platform)

## ------------------------------------------------------------------------
makefile() +
  make_group(
    make_comment("Definitions"),
    make_def("R_VERSION", R.version.string)
  ) +
  make_group(
    make_comment("Universal rule"),
    make_rule("all", c("first_target", "second_target"))
  ) +
  make_group(
    make_comment("Special rule"),
    make_rule(".FORCE")
  ) +
  make_comment(c("============", "Action rules", "============")) +
  make_rule("first_target", ".FORCE", "echo 'Building first target'") +
  make_rule("second_target", "first_target",
    c("echo 'Building second target'", "echo 'Done'"))

