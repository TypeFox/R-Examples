## ---- echo = FALSE, message = FALSE--------------------------------------
library("rprime")
library("knitr")
opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE,
  collapse = TRUE)

## ------------------------------------------------------------------------
library("rprime")
eprime_lists <- FrameList(read_eprime("data/MINP_001L00XS1.txt"))
eprime_lists[[1]]

## ------------------------------------------------------------------------
exp_lines <- read_eprime("data/MINP_001L00XS1.txt")
head(exp_lines)

## ------------------------------------------------------------------------
bad_lines <- read_eprime("data/not_an_eprime_file.txt")
head(bad_lines)

## ------------------------------------------------------------------------
# Experiment aborted on trial 3
aborted <- FrameList(read_eprime("data/MP_Block1_001P00XA1.txt"))

## ------------------------------------------------------------------------
chunks <- extract_chunks(exp_lines)
chunks[[2]]

## ------------------------------------------------------------------------
EprimeFrame(chunks[[2]])

