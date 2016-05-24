# ================================================================
# visual help function: expand values for combination of attributes
# =================================================================

expandValues <- function(attributes, data) {
  combination <- expand.grid(
    sapply( attributes, function(x){ c(levels(data[,x]),NA) }, simplify = FALSE )
  )
  combination <- apply(combination,1,function(x){paste(x, collapse = " + ")})
  names(combination) <- 1:length(combination)
  return(as.list(combination))
}

# ===================
# write YAML-template
# ===================

write.recoding <- function(attributes, data, file, yaml = TRUE) {
  
  # prepare the template for one attribute
  makeTemplate <- function(attribute, data) {
    if (length(attribute) > 1) {
      originalValues <- expandValues(attribute, data)
    } else {
      originalValues <- levels(data[,attribute])
    }
    return(list(
      recodingOf = attribute,
      attribute = NULL,
      values = list(NULL,NULL),
      link = NULL,
      originalValues = originalValues,
      comments = NULL
    ))
  }
  
  # combine all templates
  attributes <- as.list(sapply(attributes,function(x){colnames(data)[x]}))
  result <- list(
    title = NULL,
    author = NULL,
    date = format(Sys.time(),"%Y-%m-%d"),
    original_data = NULL,
    recoding = sapply(attributes, function(x) { makeTemplate(x, data) }, simplify = FALSE)
  )
  
  # return the result, defaults to a yaml-file
  if (yaml) {
    if (is.null(file)) {
      stop("please specify file")
    }
    yaml <- yaml::as.yaml(result)
    yaml <- gsub("\n- recodingOf:","\n# ==========\n- recodingOf:",yaml)
    cat(yaml, file = file)
  } else {
    return(result)
  }
}

# ========================================
# Read YAML files, and normalize shortcuts
# ========================================

read.recoding <- function(recoding, file = NULL, data = NULL) {
  
  # recodings can be a file as input
  # remember any metadata already included
  if (is.character(recoding)) {
    infile <- yaml::yaml.load_file(recoding)
    meta <- infile[-which(names(infile)=="recoding")]
    recoding <- infile$recoding  
  } else {
    if (!is.null(recoding$recoding)) {
      meta <- recoding[-which(names(recoding)=="recoding")]
      recoding <- recoding$recoding      
    } else {
      meta <- NULL
    }
  }
  
  # Allow for various shortcuts in the writing of recodings
  # The following lines normalise the input to the cannonical form
  
  reallabels <- c(
    "recodingOf", "attribute", "values", "link",
    "originalValues", "doNotRecode", "comments")
  remove <- c()
  
  for (i in 1:length(recoding)) {
    
    # write labels in full
    names(recoding[[i]]) <- reallabels[pmatch(names(recoding[[i]]),reallabels)]    
    
    # when doNotRecode is specified, you're ready to go
    # if not: then 
    if (is.null(recoding[[i]]$doNotRecode)) {
      # recodingOf is necessary, otherwise break
      if (is.null(recoding[[i]]$recodingOf)) {
        stop(paste("Specify **recodingOf** for recoding number", i, sep = " "))
      }
      # with no link, add doNotRecode
      if (is.null(recoding[[i]]$link)) { 
        recoding[[i]] <- list(doNotRecode = recoding[[i]]$recodingOf)
      } else {
        recoding[[i]]$link <- as.integer(recoding[[i]]$link)
        # make attribute and value names if necessary
        if (is.null(recoding[[i]]$attribute)) {
          recoding[[i]]$attribute <- paste("Att", i, sep = "")
        }
        if (is.null(unlist(recoding[[i]]$values))) {
          recoding[[i]]$values <- paste("val", 1:length(recoding[[i]]$link), sep = "")
        }
      }
    } else {
      if (!is.null(recoding[[i]]$link)) {
        stop(paste("Both doNotRecode and link specified in recoding number", i, sep = " "))
      }
    }
    
    # when data is specified, add names of original attributes and original values
    # this leads to nicer documentation of the recoding
    if (!is.null(data)) {
      if (is.numeric(recoding[[i]]$recodingOf)) {
        recoding[[i]]$recodingOf <- colnames(data)[recoding[[i]]$recodingOf]
      }
      if (length(recoding[[i]]$recodingOf) == 1) {  
        recoding[[i]]$originalValues <- levels(data[,recoding[[i]]$recodingOf])
      }
      if (length(recoding[[i]]$recodingOf) > 1) {
        recoding[[i]]$originalValues <- expandValues(recoding[[i]]$recodingOf, data)
      }
      if (is.numeric(recoding[[i]]$doNotRecode)) {
        recoding[[i]]$doNotRecode <- colnames(data)[recoding[[i]]$doNotRecode]
      }    
    }
    # put everything in the same order
    recoding[[i]] <- recoding[[i]][reallabels]
    recoding[[i]] <- recoding[[i]][na.omit(names(recoding[[i]]))]
    
    # merge sequences of doNotRecode
    if (i > 1) {
      if (!is.null(recoding[[i]]$doNotRecode) & !is.null(recoding[[i-1]]$doNotRecode)) {
        recoding[[i]]$doNotRecode <- c(recoding[[i-1]]$doNotRecode, recoding[[i]]$doNotRecode)
        remove <- c(remove,(i-1))
      }
    }
  }
  
  # remove superflous recodings because of contractions of doNotRecode
  if (!is.null(remove)) {
    recoding <- recoding[-remove]
  }
  
  # return result
  if (is.null(file)) {
    return(recoding)
  } else {
    # add metadata and write out as yaml
    if (!("date" %in% names(meta))) {
      meta <- c(list(date = format(Sys.time(),"%Y-%m-%d")), meta)
    }
    if (!("author" %in% names(meta))) {
      meta <- c(list(author = NULL), meta)
    }
    if (!("title" %in% names(meta))) {
      meta <- c(list(title = NULL), meta)
    }
    outfile <- c(meta, list(recoding = recoding))
    yaml <- yaml::as.yaml(outfile)
    yaml <- gsub("\n- recodingOf:","\n# ==========\n- recodingOf:",yaml)
    yaml <- gsub("\n- doNotRecode:","\n# ==========\n- doNotRecode:",yaml)
    cat(yaml, file = file)
  }
}