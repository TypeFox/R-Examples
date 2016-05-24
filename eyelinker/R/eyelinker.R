#' eyelinker: read raw data from Eyelink eyetrackers
#'
#' Eyelink eye trackers output a horrible mess, typically under
#'  the form of an .asc file. The file in question is an assorted collection of
#'  messages, events and raw data. This R package will attempt to make sense of it.
#'
#' The main function in the package is read.asc. 
#' @docType package
#' @name eyelinker
NULL

#' @importFrom stringr str_match str_split str_sub str_trim str_replace str_replace_all str_detect fixed
#' @importFrom readr read_tsv
#' @importFrom plyr llply dlply ldply mutate
#' @importFrom stringi stri_enc_toascii
#' @importFrom magrittr "%>%"
#' @importFrom intervals which_nearest distance_to_nearest Intervals
NULL
