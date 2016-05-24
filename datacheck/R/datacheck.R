#' Check a table against a set of constraints or rules defined in R.
#'
#' The rules can be written in standard R syntax. 
#' A rule must contain the names of 'columns' or variables present in the table 
#' and use R operators or simple functions. If not, the rule will simply be ignored. 
#' Each line must 'test' one rule and return a vector of boolean values as many as the
#' table has rows. Rules must not contain an assignment. The set of rules is simply
#' defined as a set of R statements and can be mixed with empty lines and comments.
#' Comments after a rule will be used for summarizing rule check results in a table and
#' should therefore be short - usually short names. This allows to visually organize
#' rules in a file and also document them. One may put more extensive comments just 
#' before the rule and add a short name or comment on the same line after it. This also
#' allows to use standard R editors for development of the rules.
#' 
#' A simple score is calculated based on the number of rules a datapoint (= table cell)
#' complies with. Like in a school test only the number of correct answers (or rule compliances)
#' are counted. Summaries of scores by row (record) and column (variable) are added to a 
#' score data frame.
#' 
#' The table itself must be a simple dataframe or .csv file.
#' 
#' The package includes a simple graphical user interface as a web page. 
#' This can be started with 
#' \code{run_datacheck()}. This interface shows summaries of the checks by rule and by
#' record. The score table can be 'downloaded'. The user interface is meant as an easy
#' way to get to know the package. All results can be also created using the command line
#' interface of R.
#'
#' The main function and the principal example can be found under \code{datadict_profile}.
#'
#' Several helper functions like \code{is_proper_name} or \code{is_only_lowers} are for convenience
#' and illustration on how to express rules more clearly or succinct.
#'
#' @name datacheck-package
#' @docType package
NA

#library(grDevices)
#library(stringr)
#library(Hmisc)
#library(shiny)

#' Tests for presence of most common punctuation characters
#'
#' These include: .?'!_#$^&*()+=<>,;' - and [tab] and [at]
#'
#' @aliases has.punct
#' @param s a character string
#' @author Reinhard Simon
#' @return boolean TRUE if detects anything
#' @family rule_checks
#' @export
#' @example inst/examples/has_punct.R
has.punct <- function(s) {
  .Deprecated("has_punct")
  stringr::str_detect(s, "[\\.\\?\\'\\!\\_\\@\\#\\$\\^\\&\\*\\(\\)\\+\\=\\<\\>\\,\\:\\;\\'\\ \\\t\\-]")
}

#' Tests for presence of most common punctuation characters
#'
#' These include: .?'!_#$^&*()+=<>,;' - and [tab] and [at]
#'
#' @aliases has_punct
#' @param s a character string
#' @author Reinhard Simon
#' @return boolean TRUE if detects anything
#' @family rule_checks
#' @export
#' @example inst/examples/has_punct.R
has_punct <- function(s) {
  stringr::str_detect(s, "[\\.\\?\\'\\!\\_\\@\\#\\$\\^\\&\\*\\(\\)\\+\\=\\<\\>\\,\\:\\;\\'\\ \\\t\\-]")
}


#' Tests if string is like a proper name with inital letter in upper case
#'
#' @aliases is.properName
#' @param aname a character string
#' @author Reinhard Simon, Jose Francisco Loff
#' @return boolean TRUE if ok
#' @family rule_checks
#' @example inst/examples/is_properName.R
#' @export
is.properName <- function(aname) {
  .Deprecated("is_proper_name")
  no_punct <- !has_punct(aname)
  proper_case <- grepl("^[A-Z]{1}[a-z]+$", aname)
  return(no_punct & proper_case)
}

#' Tests if string is like a proper name with inital letter in upper case
#'
#' @aliases is_proper_name
#' @param aname a character string
#' @author Reinhard Simon, Jose Francisco Loff
#' @return boolean TRUE if ok
#' @family rule_checks
#' @example inst/examples/is_properName.R
#' @export
is_proper_name <- function(aname) {
  no_punct <- !has_punct(aname)
  proper_case <- grepl("^[A-Z]{1}[a-z]+$", aname)
  return(no_punct & proper_case)
}


#' Tests if a string has only lower case letters
#'
#' @aliases is.onlyLowers
#' @param s a character string
#' @author Reinhard Simon
#' @return boolean TRUE if detects anything
#' @family rule_checks
#' @export
#' @example inst/examples/is_onlyLowers.R
is.onlyLowers <- function(s) {
  .Deprecated("is_only_lowers")
  stopifnot(is.character(s))
  grepl("^[a-z]+$", s)
}

#' Tests if a string has only lower case letters
#'
#' @aliases is_only_lowers
#' @param s a character string
#' @author Reinhard Simon
#' @return boolean TRUE if detects anything
#' @family rule_checks
#' @export
#' @example inst/examples/is_onlyLowers.R
is_only_lowers <- function(s) {
  stopifnot(is.character(s))
  grepl("^[a-z]+$", s)
}



#' Tests if a string or 'factor level' is one of a pre-defined set
#'
#' The aset parameter may point to a file with level names.
#' This is useful if there are many levels like in database of world countries.
#' The file path may be an absolute one or relative to the current working directory.
#'
#' The supporting table must have two columns named 'VALUES' and 'LABELS'.
#' The lookup file must be in comma separated format and using the '.csv' extension. It must
#' also be encoded using UTF-8 character set for being able to use foreign characters across
#' operating systems. This is often an issue when using Excel to develop the file.
#'
#' The x parameter may have just one level or multiple levels separated by ';'. Likewise the aset
#' parameter may have just one level or multiple levels separated by ';'. In any case the x parameter
#' must be a subset of aset (or the lookup file): see the example section.
#'
#'
#' @aliases is.oneOf
#' @param x a factor level as character string
#' @param aset a vector of character strings or a path to a custom file (full pathname where necessary)
#' @author Reinhard Simon, Jose Francisco Loff
#' @return boolean TRUE if detects anything
#' @family rule_checks
#' @export
#' @example inst/examples/is_oneOf.R
is.oneOf <- function(x, aset) {
  .Deprecated("is_one_of")
  stopifnot(is.character(x), is.character(aset))
  
  if (stringr::str_detect(aset[1], "(\\.csv)")) {
    if (file.exists(aset)) {
      afile <- aset
      ast <- read.csv(afile, stringsAsFactors = FALSE)
      ast[, c("VALUES")] <- as.character(ast[, c("VALUES")])
      ast <- stringr::str_trim(ast$VALUES)
      aset <- unique(ast)
    } else {
      stop()
    }
  }
  
  res <- x %in% aset
  return(res)
}

#' Tests if a string or 'factor level' is one of a pre-defined set
#'
#' The aset parameter may point to a file with level names.
#' This is useful if there are many levels like in database of world countries.
#' The file path may be an absolute one or relative to the current working directory.
#'
#' The supporting table must have two columns named 'VALUES' and 'LABELS'.
#' The lookup file must be in comma separated format and using the '.csv' extension. It must
#' also be encoded using UTF-8 character set for being able to use foreign characters across
#' operating systems. This is often an issue when using Excel to develop the file.
#'
#' The x parameter may have just one level or multiple levels separated by ';'. Likewise the aset
#' parameter may have just one level or multiple levels separated by ';'. In any case the x parameter
#' must be a subset of aset (or the lookup file): see the example section.
#'
#'
#' @aliases is_one_of
#' @param x a factor level as character string
#' @param aset a vector of character strings or a path to a custom file (full pathname where necessary)
#' @author Reinhard Simon, Jose Francisco Loff
#' @return boolean TRUE if detects anything
#' @family rule_checks
#' @export
#' @example inst/examples/is_oneOf.R
is_one_of <- function(x, aset) {
  stopifnot(is.character(x), is.character(aset))
  
  if (stringr::str_detect(aset[1], "(\\.csv)")) {
    if (file.exists(aset)) {
      afile <- aset
      ast <- read.csv(afile, stringsAsFactors = FALSE)
      ast[, c("VALUES")] <- as.character(ast[, c("VALUES")])
      ast <- stringr::str_trim(ast$VALUES)
      aset <- unique(ast)
    } else {
      stop()
    }
  }
  
  res <- x %in% aset
  return(res)
}

#' Tests if a numeric value is between a minimal and maximum value. Serves as convenience function.
#'
#' @aliases is.withinRange
#' @param val The value to be checked
#' @param min The minimal value (inclusive)
#' @param max The maximum value (inclusive)
#' @author Reinhard Simon
#' @return boolean TRUE if detects anything
#' @family rule_checks
#' @export
#' @example inst/examples/is_withinRange.R
is.withinRange <- function(val, min, max) {
  stopifnot(is.numeric(val), is.numeric(min), is.numeric(max))
  val <- as.numeric(val)
  val >= min & val <= max
}



#' Tests if a numeric value is between a minimal and maximum value. Serves as convenience function.
#'
#' @aliases is_within_range
#' @param val The value to be checked
#' @param min The minimal value (inclusive)
#' @param max The maximum value (inclusive)
#' @author Reinhard Simon
#' @return boolean TRUE if detects anything
#' @family rule_checks
#' @export
#' @example inst/examples/is_withinRange.R
is_within_range <- function(val, min, max) {
  stopifnot(is.numeric(val), is.numeric(min), is.numeric(max))
  val <- as.numeric(val)
  val >= min & val <= max
}


# Applies a rule for a variable in a table (db)
check_data_rule <- function(arule, adf, ruleName = NULL) {
  r <- parse(text = arule)
  res <- rep(NA, nrow(adf))
  dn <- names(adf)[names(adf) %in% ls()]
  rm(list = dn)
  res <- with(adf, {
    eval(r)
  })
  return(res)
}

#' Converts a vector of lines into a rules data frame
#'
#' The rules must be one per line and should evaluate to a vector of TRUE or FALSE.
#'
#' A rule must only refer to one 'column' name.
#' Rule statements may not have assignment operators (= or <-).
#' Rules may be separated by empty or commented lines.
#' A comment after a rule is used to document the specific rule in the summary table.
#'
#' @aliases as.rules
#' @param lines R statements with conditions
#' @author Reinhard Simon
#' @return The results as a datadict object or NA for 'empty' rule lines
#' @family datadict
#' @export
as.rules <- function(lines = "") {
  .Deprecated("as_rules")
  if (!is.character(lines)) 
    return("")
  stopifnot(is.vector(lines))
  
  res <- lines
  res <- res[which(res != "")]
  
  # Filter out comment lines starting with '#'
  res <- stringr::str_trim(res)
  pos <- stringr::str_locate(res, "#")[, 1] > 1
  res <- res[pos | is.na(pos)]
  
  pos <- stringr::str_detect(res, "package =") | stringr::str_detect(res, "package=") | stringr::str_detect(res, " = ") | stringr::str_detect(res, 
    "<-")
  res <- res[!pos]
  
  if (length(res) == 0) 
    return(NA)
  
  # Make a data.frame with two columns: rules, comments
  res <- as.data.frame(res, stringsAsF = F)
  vars <- rep("NA", n = nrow(res))
  coms <- rep("NA", n = nrow(res))
  typs <- rep("character", n = nrow(res))
  res <- cbind(vars, typs, res, coms)
  res <- as.data.frame(res)
  names(res) <- c("Variable", "Type", "Rule", "Comment")
  for (i in 1:ncol(res)) res[, i] <- as.character(res[, i])
  
  for (i in 1:nrow(res)) {
    if (stringr::str_detect(res[i, "Rule"], " #")) {
      cmt <- stringr::str_split(res[i, "Rule"], " #")[[1]][2]
    } else {
      cmt <- "None"
    }
    res[i, "Comment"] <- stringr::str_trim(cmt)
    p <- parse(text = res[i, "Rule"])
    res[i, "Rule"] <- as.character(p)
    res[i, "Variable"] <- all.vars(p)[1]
  }
  
  # Assign a variable type by examination of rules
  vn <- unique(res$Variable)
  for (j in 1:length(vn)) {
    hasVar <- res$Variable == vn[j]
    for (i in 1:nrow(res)) {
      hasTypeI <- stringr::str_detect(res$Rule[i], "is.integer") & res$Variable[i] == vn[j]
      hasTypeD <- stringr::str_detect(res$Rule[i], "is.double") & res$Variable[i] == vn[j]
      hasTypeN <- stringr::str_detect(res$Rule[i], "is.numeric") & res$Variable[i] == vn[j]
      hasTypeL <- stringr::str_detect(res$Rule[i], "is.logical") & res$Variable[i] == vn[j]
      
      if (any(hasTypeI)) 
        res[hasVar, "Type"] <- "integer"
      if (any(hasTypeD)) 
        res[hasVar, "Type"] <- "double"
      if (any(hasTypeL)) 
        res[hasVar, "Type"] <- "logical"
      if (any(hasTypeN)) 
        res[hasVar, "Type"] <- "numeric"
    }
  }
  
  class(res) <- c("data.rules", "data.frame")
  return(res)
}

#' Converts a vector of lines into a rules data frame
#'
#' The rules must be one per line and should evaluate to a vector of TRUE or FALSE.
#'
#' A rule must only refer to one 'column' name.
#' Rule statements may not have assignment operators (= or <-).
#' Rules may be separated by empty or commented lines.
#' A comment after a rule is used to document the specific rule in the summary table.
#'
#' @aliases as_rules
#' @param lines R statements with conditions
#' @author Reinhard Simon
#' @return The results as a datadict object or NA for 'empty' rule lines
#' @family datadict
#' @export
as_rules <- function(lines = "") {
  if (!is.character(lines)) 
    return("")
  stopifnot(is.vector(lines))
  
  res <- lines
  res <- res[which(res != "")]
  
  # Filter out comment lines starting with '#'
  res <- stringr::str_trim(res)
  pos <- stringr::str_locate(res, "#")[, 1] > 1
  res <- res[pos | is.na(pos)]
  
  pos <- stringr::str_detect(res, "package =") | stringr::str_detect(res, "package=") | stringr::str_detect(res, " = ") | stringr::str_detect(res, 
                                                                                                          "<-")
  res <- res[!pos]
  
  if (length(res) == 0) 
    return(NA)
  
  # Make a data.frame with two columns: rules, comments
  res <- as.data.frame(res, stringsAsF = F)
  vars <- rep("NA", n = nrow(res))
  coms <- rep("NA", n = nrow(res))
  typs <- rep("character", n = nrow(res))
  res <- cbind(vars, typs, res, coms)
  res <- as.data.frame(res)
  names(res) <- c("Variable", "Type", "Rule", "Comment")
  for (i in 1:ncol(res)) res[, i] <- as.character(res[, i])
  
  for (i in 1:nrow(res)) {
    if (stringr::str_detect(res[i, "Rule"], " #")) {
      cmt <- stringr::str_split(res[i, "Rule"], " #")[[1]][2]
    } else {
      cmt <- "None"
    }
    res[i, "Comment"] <- stringr::str_trim(cmt)
    p <- parse(text = res[i, "Rule"])
    res[i, "Rule"] <- as.character(p)
    res[i, "Variable"] <- all.vars(p)[1]
  }
  
  # Assign a variable type by examination of rules
  vn <- unique(res$Variable)
  for (j in 1:length(vn)) {
    hasVar <- res$Variable == vn[j]
    for (i in 1:nrow(res)) {
      hasTypeI <- stringr::str_detect(res$Rule[i], "is.integer") & res$Variable[i] == vn[j]
      hasTypeD <- stringr::str_detect(res$Rule[i], "is.double") & res$Variable[i] == vn[j]
      hasTypeN <- stringr::str_detect(res$Rule[i], "is.numeric") & res$Variable[i] == vn[j]
      hasTypeL <- stringr::str_detect(res$Rule[i], "is.logical") & res$Variable[i] == vn[j]
      
      if (any(hasTypeI)) 
        res[hasVar, "Type"] <- "integer"
      if (any(hasTypeD)) 
        res[hasVar, "Type"] <- "double"
      if (any(hasTypeL)) 
        res[hasVar, "Type"] <- "logical"
      if (any(hasTypeN)) 
        res[hasVar, "Type"] <- "numeric"
    }
  }
  
  class(res) <- c("data.rules", "data.frame")
  return(res)
}


#' Reads a file containing rules in data dictionary format.
#'
#' The rules must be one per line and should evaluate to a vector of TRUE or FALSE.
#'
#' A rule must only refer to one 'column' name.
#' Rule statements may not have assignment operators (= or <-).
#' Rules may be separated by empty or commented lines.
#' A comment after a rule is used to document the specific rule in the summary table.
#'
#' @aliases read.rules
#' @param file R file with conditions
#' @author Reinhard Simon
#' @return The results as a datadict object or NA for 'empty' rules file
#' @family datadict
#' @export
read.rules <- function(file = "") {
  .Deprecated("read_rules")
  res <- NA
  if (is.na(file)) 
    return(res)
  if (file == "") 
    return(res)
  
  path <- file.path(file)
  res <- readLines(path)
  
  res <- as_rules(res)
  
  return(res)
}

#' Reads a file containing rules in data dictionary format.
#'
#' The rules must be one per line and should evaluate to a vector of TRUE or FALSE.
#'
#' A rule must only refer to one 'column' name.
#' Rule statements may not have assignment operators (= or <-).
#' Rules may be separated by empty or commented lines.
#' A comment after a rule is used to document the specific rule in the summary table.
#'
#' @aliases read_rules
#' @param file R file with conditions
#' @author Reinhard Simon
#' @return The results as a datadict object or NA for 'empty' rules file
#' @family datadict
#' @export
read_rules <- function(file = "") {
  res <- NA
  if (is.na(file)) 
    return(res)
  if (file == "") 
    return(res)
  
  path <- file.path(file)
  if(!file.exists(path)) stop(paste("This file does not exist:", file))
  res <- NULL
  res <- try(
    readLines(path)
  )
  if(is.null(res)) stop("The file does not seem to contain any rules!")
  
  res <- as_rules(res)
  
  return(res)
}


is_data_rules <- function(x) {
  inherits(x, "data.rules")
}


#' Create a data quality profile (main function)
#'
#' Tests a database against a set of rules (one per line) in a 'data dictionary file'.
#' Rules will be summarized in the returned object: the variable/column, the rule, any comment after the rule,
#' the execution success, the total number of rule violations if any, the record id for any
#' non-compliant records. Rules that can't be executed for any reason will be marked as 'failed'.
#'
#' The rule file must be a simple list of one rule per line. Functions can be used but since
#' they are applied on a 'vector' (the column) they should be used within a sapply statement
#' (see example rule file). Rules may be separated by empty lines or lines with comment character #.
#' Comments after a rule within the same line will be used for display in the summary table and
#' should be short. A rule must only test one variable and one aspect at a time.
#'
#' @aliases datadict.profile
#' @param atable a data.frame
#' @param adictionary a list of rules in rule format
#' @author Reinhard Simon
#' @return a data.profile object or NA
#' @family datadict
#' @export
#' @example inst/examples/datadict_profile.R
datadict.profile <- function(atable, adictionary) {
  .Deprecated("datadict_profile")
  if (is.character(adictionary)) 
    return(NA)
  stopifnot(is.data.frame(atable), is_data_rules(adictionary))
  
  # Progressbar initialization
  pb <- NULL
  steps <- 100
  cat("\n   Checking data dictionary rules:\n")
  pb <- txtProgressBar(0, steps, style = 3, title = "datadict", label = "Checking rules:")
  setTxtProgressBar(pb, 1)
  
  try({
    at <- atable
    dq <- matrix(0, nrow = nrow(at), ncol = ncol(at))
    dq <- as.data.frame(dq)
    names(dq) <- names(at)
    rownames(dq) <- rownames(at)
    ad <- adictionary
    
    # Convert columns in table in expected format
    tn <- names(at)
    n <- length(tn)
    
    for (i in 1:n) {
      try({
        atype <- ad[ad$Variable == tn[i], "Type"][1]
        if (!is.na(atype)) {
          at[, tn[i]] <- as.character(at[, tn[i]])
          if (atype == "integer") 
          at[, tn[i]] <- as.integer(at[, tn[i]])
          if (atype == "numeric") 
          at[, tn[i]] <- as.numeric(at[, tn[i]])
          if (atype == "logical") 
          at[, tn[i]] <- as.logical(at[, tn[i]])
        }
      })
    }
    n <- nrow(ad)
    w <- steps/n
    
    Execution <- rep("ok", n)
    Error.sum <- rep(0, n)
    Error.list <- rep("none", n)
    ad <- cbind(ad, Execution, Error.sum, Error.list)
    ad[, "Error.list"] <- as.character(ad$Error.list)
    ad[, "Execution"] <- as.character(ad$Execution)
    for (i in 1:n) {
      res <- NA
      ers <- which(!res)
      cok <- FALSE
      try({
        res <- check_data_rule(ad[i, "Rule"], at)
        out <- res
        out[out == TRUE] <- 1
        out[out == FALSE] <- 0
        out[is.na(out)] <- 0
        dq[, ad[i, "Variable"]] <- dq[, ad[i, "Variable"]] + out
        ers <- which(!res)
        cok <- TRUE
      }, silent = TRUE)
      if (length(ers) > 0) {
        ad$Error.sum[i] <- length(ers)
        ad$Error.list[i] <- paste(ers, collapse = ",")
        
      }
      if (!cok) {
        ad$Error.sum[i] <- 0
        ad$Error.list[i] <- "NA"
        ad$Execution[i] <- "failed"
      }
      setTxtProgressBar(pb, floor(w * i))
    }
    close(pb)
    ad$Error.sum <- round(ad$Error.sum, 0)
    ad$Error.sum <- as.integer(ad$Error.sum)
    
    dq <- rbind(dq, colSums(dq, na.rm = TRUE))
    rownames(dq)[nrow(dq)] <- "Attribute.score"
    
    rpv <- table(ad$Variable)  # Rules per variable
    rvn <- names(rpv)
    rec <- rep(0, ncol(dq))
    for (i in 1:ncol(dq)) {
      if (names(dq[i]) %in% rvn) {
        rec[i] <- as.integer(rpv[[names(dq[i])]])
      }
    }
    dq <- rbind(dq, rec)
    rownames(dq)[nrow(dq)] <- "Rules.per.variable"
    
    dq <- cbind(dq, rowSums(dq, na.rm = TRUE))
    names(dq)[ncol(dq)] <- "Record.score"
    
    
    res <- list(data = at, checks = ad, scores = dq)
    class(res) <- c("datadict.profile", "list")
    return(res)
  }, silent = TRUE)
  return(NA)
}

#' Create a data quality profile (main function)
#'
#' Tests a database against a set of rules (one per line) in a 'data dictionary file'.
#' Rules will be summarized in the returned object: the variable/column, the rule, any comment after the rule,
#' the execution success, the total number of rule violations if any, the record id for any
#' non-compliant records. Rules that can't be executed for any reason will be marked as 'failed'.
#'
#' The rule file must be a simple list of one rule per line. Functions can be used but since
#' they are applied on a 'vector' (the column) they should be used within a sapply statement
#' (see example rule file). Rules may be separated by empty lines or lines with comment character #.
#' Comments after a rule within the same line will be used for display in the summary table and
#' should be short. A rule must only test one variable and one aspect at a time.
#'
#' @aliases datadict_profile
#' @param atable a data.frame
#' @param adictionary a list of rules in rule format
#' @author Reinhard Simon
#' @return a data.profile object or NA
#' @family datadict
#' @export
#' @example inst/examples/datadict_profile.R
datadict_profile <- function(atable, adictionary) {
  if (is.character(adictionary)) 
    return(NA)
  stopifnot(is.data.frame(atable), is_data_rules(adictionary))
  
  # Progressbar initialization
  pb <- NULL
  steps <- 100
  cat("\n   Checking data dictionary rules:\n")
  pb <- txtProgressBar(0, steps, style = 3, title = "datadict", label = "Checking rules:")
  setTxtProgressBar(pb, 1)
  
  try({
    at <- atable
    dq <- matrix(0, nrow = nrow(at), ncol = ncol(at))
    dq <- as.data.frame(dq)
    names(dq) <- names(at)
    rownames(dq) <- rownames(at)
    ad <- adictionary
    
    # Convert columns in table in expected format
    tn <- names(at)
    n <- length(tn)
    
    for (i in 1:n) {
      try({
        atype <- ad[ad$Variable == tn[i], "Type"][1]
        if (!is.na(atype)) {
          at[, tn[i]] <- as.character(at[, tn[i]])
          if (atype == "integer") 
            at[, tn[i]] <- as.integer(at[, tn[i]])
          if (atype == "numeric") 
            at[, tn[i]] <- as.numeric(at[, tn[i]])
          if (atype == "logical") 
            at[, tn[i]] <- as.logical(at[, tn[i]])
        }
      })
    }
    n <- nrow(ad)
    w <- steps/n
    
    Execution <- rep("ok", n)
    Error.sum <- rep(0, n)
    Error.list <- rep("none", n)
    ad <- cbind(ad, Execution, Error.sum, Error.list)
    ad[, "Error.list"] <- as.character(ad$Error.list)
    ad[, "Execution"] <- as.character(ad$Execution)
    for (i in 1:n) {
      res <- NA
      ers <- which(!res)
      cok <- FALSE
      try({
        res <- check_data_rule(ad[i, "Rule"], at)
        out <- res
        out[out == TRUE] <- 1
        out[out == FALSE] <- 0
        out[is.na(out)] <- 0
        dq[, ad[i, "Variable"]] <- dq[, ad[i, "Variable"]] + out
        ers <- which(!res)
        cok <- TRUE
      }, silent = TRUE)
      if (length(ers) > 0) {
        ad$Error.sum[i] <- length(ers)
        ad$Error.list[i] <- paste(ers, collapse = ",")
        
      }
      if (!cok) {
        ad$Error.sum[i] <- 0
        ad$Error.list[i] <- "NA"
        ad$Execution[i] <- "failed"
      }
      setTxtProgressBar(pb, floor(w * i))
    }
    close(pb)
    ad$Error.sum <- round(ad$Error.sum, 0)
    ad$Error.sum <- as.integer(ad$Error.sum)
    
    dq <- rbind(dq, colSums(dq, na.rm = TRUE))
    rownames(dq)[nrow(dq)] <- "Attribute.score"
    
    rpv <- table(ad$Variable)  # Rules per variable
    rvn <- names(rpv)
    rec <- rep(0, ncol(dq))
    for (i in 1:ncol(dq)) {
      if (names(dq[i]) %in% rvn) {
        rec[i] <- as.integer(rpv[[names(dq[i])]])
      }
    }
    dq <- rbind(dq, rec)
    rownames(dq)[nrow(dq)] <- "Rules.per.variable"
    
    dq <- cbind(dq, rowSums(dq, na.rm = TRUE))
    names(dq)[ncol(dq)] <- "Record.score"
    
    
    res <- list(data = at, checks = ad, scores = dq)
    class(res) <- c("datadict.profile", "list")
    return(res)
  }, silent = TRUE)
  return(NA)
}


#' is.datadict.profile
#'
#' Is this a datadict.profile object
#'
#' @aliases is.datadict.profile
#' @param x The object to be tested.
#' @author Reinhard Simon
#' @return boolean
#' @family datadict
#' @export
#' @example inst/examples/datadict_profile.R
is.datadict.profile <- function(x) {
  .Deprecated("is_datadict_profile")
  inherits(x, "datadict.profile")
}

#' is.datadict.profile
#'
#' Is this a datadict.profile object
#'
#' @aliases is_datadict_profile
#' @param x The object to be tested.
#' @author Reinhard Simon
#' @return boolean
#' @family datadict
#' @export
#' @example inst/examples/datadict_profile.R
is_datadict_profile <- function(x) {
  inherits(x, "datadict.profile")
}


#' Quick check if a rule profile on a table has any errors.
#'
#' @aliases has.ruleErrors
#' @param profile.rules a data.profile object
#' @author Reinhard Simon
#' @return boolean
#' @family datadict
#' @example inst/examples/has_ruleErrors.R
#' @export
has.ruleErrors <- function(profile.rules) {
  .Deprecated("has_rule_errors")
  xx <- sum(profile.rules$Error.sum)
  xx > 0
}

#' Quick check if a rule profile on a table has any errors.
#'
#' @aliases has_rule_errors
#' @param profile_rules a data.profile object
#' @author Reinhard Simon
#' @return boolean
#' @family datadict
#' @example inst/examples/has_ruleErrors.R
#' @export
has_rule_errors <- function(profile_rules) {
  xx <- sum(profile_rules$Error.sum)
  xx > 0
}

#' Prepares a summary table for display in a 'printed' report.
#'
#' Currently reduces the number of displayed record ids to 5
#' and adds a referral.
#'
#' @aliases prep4rep
#' @param rule.checks table in a data.profile object
#' @param txt text to be added after the first 5 record ids
#' @author Reinhard Simon
#' @return the modified rule.checks table
#' @family datadict
#' @export
prep4rep <- function(rule.checks, txt = "... more") {
  rc <- rule.checks
  n <- nrow(rc)
  for (i in 1:n) {
    s <- rc$Error.list[i]
    if (stringr::str_detect(s, ",")) {
      if (stringr::str_count(s, ",") > 5) {
        ss <- stringr::str_split(s, ",")[[1]]
        st <- paste(ss[1:5], collapse = ",")
        st <- paste(st, txt)
        rc$Error.list[i] <- st
      }
    }
  }
  return(rc)
}

#' Get the current version of a package
#'
#' Uses the citation() function.
#'
#' @aliases pkg.version
#' @param pkg the package name
#' @author Reinhard Simon
#' @return a string with the package number
#' @family helper
#' @export
pkg.version <- function(pkg) {
  .Deprecated("pkg_version")
  cit <- citation(pkg)
  stringr::str_extract(cit$note, "[0-9].[0-9].[0-9]")
}

#' Get the current version of a package
#'
#' Uses the citation() function.
#'
#' @aliases pkg_version
#' @param pkg the package name
#' @author Reinhard Simon
#' @return a string with the package number
#' @family helper
#' @export
pkg_version <- function(pkg) {
  cit <- citation(pkg)
  stringr::str_extract(cit$note, "[0-9].[0-9].[0-9]")
}


#' Draws a heatmap based on data quality scores
#'
#' Knows to extract the quality matrix from the profile object and pass it on to the heatmap function.
#' Plots a heatmap.
#' 
#' Currently this function is limited a table siZé of 300 records.
#'
#' @param profile a datadict.profile object
#' @param recLab variable that should be used for labeling the records
#' @param recMax maximum first n records for display
#' @param scoreMax maxim quality score to filter out
#' @param cols color scheme
#' @param ... as in heatmap function
#' @aliases heatmap.quality
#' @author Reinhard Simon
#' @family visuals
#' @export
heatmap.quality <- function(profile, recLab = NULL, recMax = 100, scoreMax = NULL, cols = NULL, ...) {
  .Deprecated("heatmap_quality")
  stopifnot(is_datadict_profile(profile))
  stopifnot(nrow(profile$data) <= 300)
  
  if (is.null(cols)) {
    cols <- colorRampPalette(c("white", "darkgreen"))(max(profile$scores))
  }
  
  
  db <- with(profile, {
    scores[1:(nrow(scores) - 2), ]
  })
  
  if (!is.null(recLab)) {
    rownames(db) <- profile$data[, recLab]
  }
  
  if (!is.null(scoreMax)) {
    db <- db[db[, ncol(db)] <= scoreMax, ]
  }
  
  db <- db[, 1:(ncol(db) - 1)]
  db <- as.matrix(db)
  
  rm <- min(nrow(db), recMax)
  db <- db[1:rm, ]
  
  if (is.matrix(db)) {
    heatmap(db, col = cols, ...)
  }
  
}

#' Draws a heatmap based on data quality scores
#'
#' Knows to extract the quality matrix from the profile object and pass it on to the heatmap function.
#' Plots a heatmap.
#' 
#' Currently this function is limited a table siZé of 300 records.
#'
#' @param profile a datadict.profile object
#' @param recLab variable that should be used for labeling the records
#' @param recMax maximum first n records for display
#' @param scoreMax maxim quality score to filter out
#' @param cols color scheme
#' @param ... as in heatmap function
#' @aliases heatmap_quality
#' @author Reinhard Simon
#' @family visuals
#' @export
heatmap_quality <- function(profile, recLab = NULL, recMax = 100, scoreMax = NULL, cols = NULL, ...) {
  stopifnot(is_datadict_profile(profile))
  stopifnot(nrow(profile$data) <= 300)
  
  if (is.null(cols)) {
    cols <- colorRampPalette(c("white", "darkgreen"))(max(profile$scores))
  }
  
  
  db <- with(profile, {
    scores[1:(nrow(scores) - 2), ]
  })
  
  if (!is.null(recLab)) {
    rownames(db) <- profile$data[, recLab]
  }
  
  if (!is.null(scoreMax)) {
    db <- db[db[, ncol(db)] <= scoreMax, ]
  }
  
  db <- db[, 1:(ncol(db) - 1)]
  db <- as.matrix(db)
  
  rm <- min(nrow(db), recMax)
  db <- db[1:rm, ]
  
  if (is.matrix(db)) {
    heatmap(db, col = cols, ...)
  }
  
}



#' Dotchart of rules per variable
#'
#' Summarizes rule coverage. There should be at least 3x coverage.
#'
#' @aliases ruleCoverage
#' @param profile a datadict profile object
#' @param rLowest lowest acceptable number of rules per variable
#' @param rLow an intermediate low number of rules per variable
#' @param rOk the 'ok' number of rules per variable
#' @param rMax = the maximum number of rules per variable
#' @author Reinhard Simon
#' @family visuals
#' @export
ruleCoverage <- function(profile, rLowest = 1, rLow = 2, rOk = 3, rMax = 10) {
  .Deprecated("rule_coverage")
  stopifnot(is.datadict.profile(profile))
  
  rc <- profile$scores[nrow(profile$scores), 1:(ncol(profile$scores) - 1)]
  rc <- as.data.frame(t(rc))
  ix <- order(rc[, 1])
  x <- rc[ix, ]
  
  avg <- round(mean(rc[, 1]), 1)
  x <- c(x, avg)
  
  rn <- rownames(rc)[ix]
  rn <- c(rn, "Average rules")
  
  dotchart(x, labels = rn, pch = 20, main = "Rules per variable", xlim = c(0, rMax))
  abline(v = rOk, col = "darkgreen")
  abline(v = rLow, col = "orange")
  abline(v = rLowest, col = "red")
}

#' Dotchart of rules per variable
#'
#' Summarizes rule coverage. There should be at least 3x coverage.
#'
#' @aliases rule_coverage
#' @param profile a datadict profile object
#' @param rLowest lowest acceptable number of rules per variable
#' @param rLow an intermediate low number of rules per variable
#' @param rOk the 'ok' number of rules per variable
#' @param rMax = the maximum number of rules per variable
#' @author Reinhard Simon
#' @family visuals
#' @export
rule_coverage <- function(profile, rLowest = 1, rLow = 2, rOk = 3, rMax = 10) {
  stopifnot(is_datadict_profile(profile))
  
  rc <- profile$scores[nrow(profile$scores), 1:(ncol(profile$scores) - 1)]
  rc <- as.data.frame(t(rc))
  ix <- order(rc[, 1])
  x <- rc[ix, ]
  
  avg <- round(mean(rc[, 1]), 1)
  x <- c(x, avg)
  
  rn <- rownames(rc)[ix]
  rn <- c(rn, "Average rules")
  
  dotchart(x, labels = rn, pch = 20, main = "Rules per variable", xlim = c(0, rMax))
  abline(v = rOk, col = "darkgreen")
  abline(v = rLow, col = "orange")
  abline(v = rLowest, col = "red")
}


#' Line chart of cumulative sum of rule scores.
#'
#' A record receives one point per rule which evaluates TRUE.
#' The total number of points is the 'quality score' per record.
#'
#' @aliases scoreSum
#' @param profile a datadict profile object
#' @author Reinhard Simon
#' @family visuals
#' @export
scoreSum <- function(profile) {
  .Deprecated("score_sum")
  stopifnot(is.datadict.profile(profile))
  tt <- table(profile$scores$Record.score)
  xx <- as.data.frame(cumsum(tt))
  dt <- cbind(as.integer(rownames(xx)), xx)
  dt <- dt[-nrow(dt), ]
  names(dt) <- c("Rule score", "Cumulative sum")
  if (is.data.frame(dt) & nrow(dt) > 1) {
    plot(dt, type = "l")
  }
  
}

#' Line chart of cumulative sum of rule scores.
#'
#' A record receives one point per rule which evaluates TRUE.
#' The total number of points is the 'quality score' per record.
#'
#' @aliases score_sum
#' @param profile a datadict profile object
#' @author Reinhard Simon
#' @family visuals
#' @export
score_sum <- function(profile) {
  stopifnot(is_datadict_profile(profile))
  tt <- table(profile$scores$Record.score)
  xx <- as.data.frame(cumsum(tt))
  dt <- cbind(as.integer(rownames(xx)), xx)
  dt <- dt[-nrow(dt), ]
  names(dt) <- c("Rule score", "Cumulative sum")
  if (is.data.frame(dt) & nrow(dt) > 1) {
    plot(dt, type = "l")
  }
  
}


#' Produces a tabular summary of descriptive statistics using the 'Hmisc::describe'
#' function from the Hmisc package.
#'
#' Returns a dataframe of descriptive statistics
#'
#' @aliases shortSummary
#' @param atable a data frame
#' @author Reinhard Simon
#' @family helper
#' @export
shortSummary <- function(atable) {
  .Deprecated("short_summary")
  stopifnot(is.data.frame(atable))
  
  ss <- Hmisc::describe(atable)
  nm <- names(ss)
  
  hd <- c("n", "missing", "unique", "value", "min", "max", "Mean", "sd", ".05", ".10", ".25", ".50", 
    ".75", ".90", ".95")
  dt <- as.data.frame(matrix("", nrow = length(nm), ncol = length(hd)), stringsAsFactors = FALSE)
  names(dt) <- hd
  rownames(dt) <- nm
  for (i in 1:length(nm)) {
    an <- nm[i]
    ar <- ss[[an]]$counts
    ac <- names(ar)
    # dt[an, ac] = ar
    rdx <- which(row.names(dt) == an)
    for (k in 1:length(ar)) {
      if (names(ar)[k] %in% names(dt)) {
        dt[rdx, names(ar[k])] <- ar[k]
      }
      
    }
    
    if (is.numeric(atable[[an]])) {
      amin <- min(atable[[an]], na.rm = T)
      amax <- max(atable[[an]], na.rm = T)
      asd <- round(sd(atable[[an]], na.rm = T), 2)
      dt[an, "min"] <- amin
      dt[an, "max"] <- amax
      dt[an, "sd"] <- asd
    }
    
  }
  dt
}

#' Produces a tabular summary of descriptive statistics using the 'Hmisc::describe'
#' function from the Hmisc package.
#'
#' Returns a dataframe of descriptive statistics
#'
#' @aliases short_summary
#' @param atable a data frame
#' @author Reinhard Simon
#' @family helper
#' @export
short_summary <- function(atable) {
  stopifnot(is.data.frame(atable))
  
  ss <- Hmisc::describe(atable)
  nm <- names(ss)
  
  hd <- c("n", "missing", "unique", "value", "min", "max", "Mean", "sd", ".05", ".10", ".25", ".50", 
          ".75", ".90", ".95")
  dt <- as.data.frame(matrix("", nrow = length(nm), ncol = length(hd)), stringsAsFactors = FALSE)
  names(dt) <- hd
  rownames(dt) <- nm
  for (i in 1:length(nm)) {
    an <- nm[i]
    ar <- ss[[an]]$counts
    ac <- names(ar)
    # dt[an, ac] = ar
    rdx <- which(row.names(dt) == an)
    for (k in 1:length(ar)) {
      if (names(ar)[k] %in% names(dt)) {
        dt[rdx, names(ar[k])] <- ar[k]
      }
      
    }
    
    if (is.numeric(atable[[an]])) {
      amin <- min(atable[[an]], na.rm = T)
      amax <- max(atable[[an]], na.rm = T)
      asd <- round(sd(atable[[an]], na.rm = T), 2)
      dt[an, "min"] <- amin
      dt[an, "max"] <- amax
      dt[an, "sd"] <- asd
    }
    
  }
  dt
}


#' Presents the packages graphical user interface
#'
#' Runs a web server to show the user interface.
#'
#' @aliases runDatacheck
#' @param port the port where to listen; 1971 by default.
#' @author Reinhard Simon
#' @family interface
#' @export
runDatacheck <- function(port = 1971L) {
  .Deprecated("run_datacheck")
  shiny::runApp(system.file("www", package = "datacheck"), port = port)
}

#' Presents the packages graphical user interface
#'
#' Runs a web server to show the user interface.
#'
#' @aliases run_datacheck
#' @param port the port where to listen; 1971 by default.
#' @author Reinhard Simon
#' @family interface
#' @export
run_datacheck <- function(port = 1971L) {
  shiny::runApp(system.file("www", package = "datacheck"), port = port)
}
