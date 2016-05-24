# =================================================================
# write orthography profile with frequencies for further processing
# =================================================================

write.profile <- function(strings
                         , normalize = NULL
                         , info = TRUE
                         , editing = FALSE
                         , sep = NULL
                         , file.out = NULL
                         , collation.locale = NULL
                         ) {
  
  # set locale
  if (!is.null(collation.locale)) {
    current.locale <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", collation.locale)
    on.exit(Sys.setlocale("LC_COLLATE", current.locale))
  }
  
  # use characters
  if (length(strings) == 1) {
    if (file.exists(strings)) {
    strings <- scan(strings, sep = "\n", what = "character")
    }
  }
  strings <- as.character(strings)
  
  # normalization
  if (is.null(normalize)) {
    transcode <- identity
  } else if (normalize == "NFC") {
    transcode <- stringi::stri_trans_nfc
  } else if (normalize == "NFD") {
    transcode <- stringi::stri_trans_nfd
  } 
  strings <- transcode(strings)
  
  # split using unicode definitions
  # except when 'sep' is specified, then split by sep
  if (is.null(sep)) {
    splitted <- stringi::stri_split_boundaries(strings, type = "character")
  } else {
    splitted <- strsplit(strings, sep)
    # remove empty characters
    splitted <- sapply(splitted, function(x){x[x != ""]}, simplify = FALSE)
  }
    
  # prepare result 
  frequency <- table(unlist(splitted))
  chars <- names(frequency)
  
  # add columns for editing when 'editing = TRUE'
  if (editing) {
    Grapheme <- cbind(  Left = ""
                       , Grapheme = chars
                       , Right = ""
                       , Class = ""
                       , Replacement = chars)
  } else {
    Grapheme <- chars
  }

  # add frequency, codepoints and Unicode names when info = TRUE
  if (info) {    
    codepoints <- sapply(chars, function(x) {
      paste(stringi::stri_trans_general(unlist(strsplit(x,"")), "Any-Hex/Unicode")
            , collapse = ", ")})
    
    names <- sapply(chars, function(x) {
      paste(stringi::stri_trans_general(unlist(strsplit(x,"")), "Any-Name")
            , collapse = ", ")})   
    names <- gsub("\\N{", "", names, fixed= TRUE)
    names <- gsub("}", "", names, fixed = TRUE)
    
    Grapheme <- cbind(Grapheme
                       , Frequency = frequency
                       , Codepoint = codepoints
                       , UnicodeName = names
                      )
  }

  # return result as data frame, or write to file when "file" is specified
  result <- as.data.frame(Grapheme, stringsAsFactors = FALSE)
  rownames(result) <- NULL
  
  if (ncol(result) == 1) {
    colnames(result) <- "Grapheme"
  }
  
  if (is.null(file.out)) {
    
    return(result)
    
  } else {
 
    # check special characters for writing output
    if (sum(grepl("\t", strings)) + 
        sum(grepl("\n", strings)) +
        sum(grepl("\r", strings))) {
      warning("There are tabs and/or newline characters in the input strings. This will lead to problems with the profiles and the tokenization. Consider removing or replacing them in your input strings.")
    }
       
    write.table(result
                , file = file.out
                , quote = FALSE
                , sep = "\t"
                , row.names = FALSE)
    
    return(invisible(result))
  }
}

# ====================================================================
# orthography profiles are just TSV files, so this is just convenience
# ====================================================================

read.profile <- function(profile) {
  
  profile <- read.table(profile
                        , sep = "\t"
                        , quote = ""
                        , header = TRUE
                        , fill = TRUE
                        , colClasses = "character"
                        )
  
  # checking header
  
  if (sum(colnames(profile) == "Grapheme") != 1) {
    stop("There needs to be a column called \'Grapheme\'")
  } else {
    return(profile)
  }

}
