# ================
# tokenize strings
# ================

tokenize <- function(strings
                      , profile = NULL
                      , transliterate = NULL
                      , method = "global"
                      , ordering = c("size","context","reverse")
                      , sep = " "
                      , sep.replace = NULL
                      , missing = "\u2047"
                      , normalize = "NFC"
                      , regex = FALSE
                      , silent = FALSE
                      , file.out = NULL
                    ) {
  
  # ---------------
  # preprocess data
  # ---------------

  if (length(strings) == 1) {
    if (file.exists(strings)) {
    strings <- scan(strings, sep = "\n", what = "character")
    }
  }
  strings <- as.character(strings)
  
  # option gives errors, so removed for now
  case.insensitive = FALSE
  
	# separators
	internal_sep <-  intToUtf8(1110000)
	user_sep <- sep

	# normalization
	if (normalize == "NFC") {
	  transcode <- stri_trans_nfc
	} else if (normalize == "NFD") {
	  transcode <- stri_trans_nfd
	} else {
    warning("Only the normalization-options NFC and NFD are implemented. No normalization will be performed.")
	  transcode <- identity
	}
	
	# keep original strings, and normalize NFC everything by default
	originals <- as.vector(strings)
	strings <- transcode(originals)

  # collapse strings for doing everything at once
	NAs <- which(is.na(strings))
	strings[NAs] <- ""
  all <- paste(strings, collapse = internal_sep)
  all <- paste0(internal_sep, all, internal_sep)
  
  # --------------------
  # read or make profile
  # --------------------
  
  # read orthography profile (or make new one)
  if (is.null(profile)) {
    
    # make new orthography profile
    if (normalize == "NFC") {
      profile  <- write.profile(strings
                                , normalize = normalize
                                , sep = NULL
                                , info = FALSE
                                )  
    } else {
      profile  <- write.profile(strings
                                , normalize = normalize
                                , sep = ""
                                , info = FALSE
                                ) 
    }
    
  } else if (is.null(dim(profile))) {
    
    if (length(profile) > 1) {
      # assume that the profile are graphemes
      profile <- data.frame(Grapheme = profile, stringsAsFactors = FALSE)
    } else {      
      if (file.exists(profile)) {
        # read the linked profile
        profile <- read.profile(profile)
      } else {
        # there is a one-string profile :-0
        profile <- data.frame(Grapheme = profile, stringsAsFactors = FALSE)
        warning("Profile provided consists of a single string. Are you sure?")
      }
    }
  } else {
    # assume the profile is a suitable R object
    profile <- profile
  }

  # first-pass reordering, only getting larger graphemes on top
  # ordering by grapheme size, if specified
  # necessary to get regexes in right order
  if (sum(!is.na(pmatch(ordering,"size"))) > 0) {
    size <- nchar(stri_trans_nfd(profile[,"Grapheme"]))
    profile <- profile[order(-size), ,drop = FALSE]
  }
  
  # normalise characters in profile, just to be sure
  graphs <- transcode(profile[,"Grapheme"])
  if (!is.null(transliterate)) {
    trans <- transcode(profile[,transliterate])
  }
  
  # is there contextual information?
  l_exists <- sum(colnames(profile) == "Left") == 1
  r_exists <- sum(colnames(profile) == "Right") == 1
  c_exists <- sum(colnames(profile) == "Class") == 1

  # then normalise them too
  if (l_exists) { 
    left <- transcode(profile[,"Left"]) 
  } else {
    left  <- ""
  }
  if (r_exists) {
    right <- transcode(profile[,"Right"])
  } else {
    right  <- ""
  }
  
  # -----------------------------------------
  # prepare regexes with context from profile
  # -----------------------------------------
  
  if (!regex) {
    
    contexts <- graphs
    
  } else {
    
    # replace regex boundaries with internal separator
    tmp <- intToUtf8(1110001)
    
    right <- gsub("\\$", tmp, right, fixed = TRUE)
    right <- gsub("\\$$", internal_sep, right)
    right <- gsub(tmp, "\\$", right, fixed = TRUE)
    
    left <- gsub("^\\^", internal_sep, left)
    left <- gsub("([^\\[])\\^", paste0("\\1",internal_sep), left)
    
    graphs <- gsub("\\$", tmp, graphs, fixed = TRUE)
    graphs <- gsub("\\$$", internal_sep, graphs)
    graphs <- gsub(tmp, "\\$", graphs, fixed =  TRUE)
    
    graphs <- gsub("^\\^", internal_sep, graphs)
    graphs <- gsub("^\\.", paste0("[^", internal_sep, "]"), graphs)
    
    # make classes if there is anything there
    if (c_exists && sum(profile[,"Class"] != "") > 0) {
      
      classes <- unique(profile[,"Class"])
      classes <- classes[classes != ""]
      groups <- sapply(classes,function(x){
        graphs[profile[,"Class"] == x]
      }, simplify = FALSE)
      classes.regex <- sapply(groups,function(x){
        paste( "((", paste( x, collapse = ")|(" ), "))", sep = "")
      })
      
      for (i in classes) {
        left <- gsub(i, classes.regex[i], left, fixed = TRUE)
        right <- gsub(i, classes.regex[i], right, fixed = TRUE)
        graphs <- gsub(i, classes.regex[i], graphs, fixed = TRUE)
      }
    }
    
    # add lookahead/lookbehind syntax and combine everything together
    left[left != ""] <- paste("(?<=", left[left != ""], ")", sep = "")
    right[right != ""] <- paste("(?=", right[right != ""], ")", sep = "")
    
    # replace dot in context with internal separator
    left <- gsub("(?<=."
                 , paste0("(?<!", internal_sep, ")(?<=" )
                 , left
                 , fixed =  TRUE
    )
    right <- gsub("(?=."
                 , paste0("(?!", internal_sep, ")(?=" )
                 , right
                 , fixed =  TRUE
    )
    
    contexts <- paste0(left, graphs, right)
    
  }
  
  # -----------------
  # reorder graphemes
  # -----------------
  
  if (is.null(ordering)) {
    
    graph_order <- 1:length(graphs)
  
  } else {
    
    # ordering by grapheme size
    if (sum(!is.na(pmatch(ordering,"size"))) > 0) {
      size <- nchar(stri_trans_nfd(graphs))
    } else {
      size <- rep(T, times = length(graphs))
    }
 
    # ordering by existing of context
    if (regex && (l_exists || r_exists)) {   
      context <- (left != "" | right != "") 
    } else {
      context <- rep(T, times = length(graphs))        
    }
    
    # reverse ordering
    if (sum(!is.na(pmatch(ordering,"reverse"))) > 0) {
      reverse <- length(graphs):1
    } else {
      reverse  <- 1:length(graphs)
    }
    
    # ordering by frequency of occurrence
    if (sum(!is.na(pmatch(ordering,"frequency"))) > 0) {
      frequency <- stri_count_regex(all
                         , pattern = contexts
                         , literal = !regex
                         , case_insensitive = case.insensitive
                         )
    } else {
      frequency <- rep(T, times = length(graphs)) 
    }
      
    # order according to dimensions chosen by user in "ordering"    
    dimensions <- list(  size = - size # largest size first
                       , context = - context # with context first
                       , reverse = reverse # reverse first
                       , frequency = frequency # lowest frequency first
                       )
    graph_order <- do.call(order, dimensions[ordering])

  }
  
  # change order
  graphs <- graphs[graph_order]
  contexts <- contexts[graph_order] 
  if (!is.null(transliterate)) {
    trans <- trans[graph_order]
  }
  
  # --------------
  # regex matching
  # --------------
  
  if (!regex) {
    
    matches <- stri_locate_all_fixed(
                      all
                      , pattern = contexts
                      , overlap = TRUE
                      , case_insensitive = case.insensitive
                      )
    matched_parts <- stri_extract_all_fixed(
                      all
                      , pattern = contexts
                      , overlap = TRUE
                      , case_insensitive = case.insensitive
                      )
    
  } else {
  
    matches <- stri_locate_all_regex(
                      all
                      , pattern = contexts
                      , case_insensitive = case.insensitive
                      )
    matched_parts <- stri_extract_all_regex(
                      all
                      , pattern = contexts
                      , case_insensitive = case.insensitive
                      )
    
  }
    
  # --------------------------------------
  # tokenize data, either global or linear
  # --------------------------------------

	if (!is.na(pmatch(method,"global"))) {

    # =================
		# function to check whether the match is still free
    # and insert graph into "taken" when free
    test_match <- function(context_nr) {
      
      m <- matches[[context_nr]]
      
      # check whether match is not yet taken
      not.already.taken <- apply(m, 1, function(x) {
                            if (is.na(x[1])) { NA } else {
                            prod(is.na(taken[x[1]:x[2]])) == 1
                            }})
      free <- which(not.already.taken)
      if (length(free) > 0) {
        no.self.overlap <- c(TRUE
                         , head(m[free,,drop = FALSE][,2],-1) < 
                           tail(m[free,,drop = FALSE][,1],-1)
                         )
        free <- free[no.self.overlap]
      }
      
      # check whether graph is regex with multiple matches
      different_graphs <- unique(matched_parts[[context_nr]])
      is.regex <- length(unique(different_graphs)) > 1
      
      # take possible matches
      for (x in free) {
        r <- m[x,]
        if (!is.regex) {
          taken[r[1]:r[2]] <<- different_graphs
        } else {
          taken[r[1]:r[2]] <<- matched_parts[[context_nr]][x]
        }
      }
      
      return(m[free, , drop = FALSE])
    }
		# =================
    
    # preparation
    taken <- rep(NA, times = nchar(all))
    
    # select matches
    selected <- sapply(1:length(matches), test_match, simplify = FALSE)
    
    # count number of matches per rule
    matched_rules <- sapply(selected, dim)[1,]
      
    # insert internal separator
    where_sep <- stri_locate_all_fixed(all, internal_sep)[[1]][,1]
    taken[where_sep] <- internal_sep
    
    # remaining NAs are missing parts
		missing_chars <- sapply(which(is.na(taken))
		                        , function(x) { stri_sub(all, x, x) }
		)
		taken[is.na(taken)] <- missing
    
    # transliteration
    if (!is.null(transliterate)) {
      transliterated <- taken
      sapply(1:length(selected), function(x) {
        apply(selected[[x]], 1, function(y) {
          transliterated[y[1]:y[2]] <<- trans[x]
        })
      })
    }
    
		# =================
    # functions to turn matches into tokenized strings
    reduce <- function(taken) {
 
      # replace longer graphs with NA, then na.omit
      sapply(selected, function(x) {
        apply(x, 1, function(y) {
          if (y[1] < y[2]) {
            taken[(y[1]+1) : y[2]] <<- NA
          }
        })
      })
      result <- na.omit(taken)
      
      return(result)
    }
      
    postprocess <- function(taken) {
      
      # replace separator
      if (!is.null(sep.replace)) {
        taken[taken == user_sep] <- sep.replace
      }
      
      # bind together tokenized parts with user separator
      taken <- paste(taken, collapse = user_sep)
      
      # remove multiple internal user separators
      taken <- gsub(paste0(user_sep,"{2,10}"), user_sep, taken)
      
      # Split string by internal separator
      result <- strsplit(taken, split = internal_sep)[[1]][-1]
      
      # remove user_sep at start and end
      result <- substr(x = result
                       , start = nchar(user_sep)+1
                       , stop = nchar(result)-nchar(user_sep)
                       )
      
      return(result)
    }
		# =================
    
    # make one string of the parts selected
    tokenized <- postprocess(reduce(taken))

    # make one string of transliterations
		if (!is.null(transliterate)) {
      transliterated <- postprocess(reduce(transliterated))
		}
    
  # ---------------------------------------------------------
  # finite-state transducer behaviour when parsing = "linear"
  # ---------------------------------------------------------	
	
	} else if (!is.na(pmatch(method,"linear"))) {
		
    # preparations
		all.matches <- do.call(rbind,matches)[,1]
		position <- 1
		tokenized <- c()
    transliterated <- c()
    missing_chars <- c()
    matched_rules <- rep.int(x = 0, times = length(contexts))
		where_sep <- stri_locate_all_fixed(all, internal_sep)[[1]][,1]
		
		graphs_match_list <- unlist(matched_parts)
    contexts_match_list <- rep(1:length(matches)
                               , times = sapply(matches, dim)[1,]
                               )   
		if (!is.null(transliterate)) {
			trans_match_list <- rep(trans
                              , times = sapply(matches, dim)[1,]
                              )
		}
		
		# loop through all positions and take first match
		while(position <= nchar(all)) {
			
      if (position %in% where_sep) {
        tokenized <- c(tokenized, internal_sep)
        if (!is.null(transliterate)) {
          transliterated <- c(transliterated, internal_sep)
        }
        position <- position +1
      } else {
      
  			hit <- which(all.matches == position)[1]
  			if (is.na(hit)) {
  				tokenized <- c(tokenized, missing)
          missing_chars <- c(missing_chars
                             , substr(all, position, position)
                             )
  				if (!is.null(transliterate)) {
            transliterated <- c(transliterated, missing)
  				}
  				position <- position + 1
  			} else {
  				tokenized <- c(tokenized, graphs_match_list[hit])
  				if (!is.null(transliterate)) {
  				  transliterated <- c(transliterated, trans_match_list[hit])
  				}
  				position <- position + nchar(graphs_match_list[hit])
          rule <- contexts_match_list[hit]
          matched_rules[rule] <- matched_rules[rule] + 1
  			}
      }
		}

    # =============
		postprocess <- function(taken) {
		  
		  # bind together tokenized parts with user separator
		  taken <- paste(taken, collapse = user_sep)
		  
		  # Split string by internal separator
		  result <- strsplit(taken, split = internal_sep)[[1]]
		  
		  # remove user_sep at start and end
		  result <- substr(result, 2, nchar(result)-1)
		  result <- result[-1]
		  
		  return(result)
		}
		# =============
 
    # postprocessing 
    tokenized <- postprocess(tokenized)
		if (!is.null(transliterate)) {
      transliterated <- postprocess(transliterated)
		}
    
	} else {
    stop(paste0("The tokenization method \"",method,"\" is not defined"))
	}
	
  # ----------------------
  # preparation of results
  # ----------------------
  
  tokenized[NAs] <- NA
  
  if (is.null(transliterate)) {
    
    strings.out <- data.frame(
      cbind(originals = originals
            , tokenized = tokenized
            )
      , stringsAsFactors = FALSE
      )
    
  } else {
    
    transliterated[NAs] <- NA
    
    strings.out <- data.frame(
      cbind(originals = originals
            , tokenized = tokenized
            , transliterated = transliterated
            )
      , stringsAsFactors = FALSE
      )
  }
    
  # Make a list of missing and throw warning
  whichProblems <- grep(pattern = missing, x = tokenized)
  problems <- strings.out[whichProblems, c(1,2)]
  colnames(problems) <- c("originals", "errors")

  if ( nrow(problems) > 0) {
    
    # make a profile for missing characters
    problemChars <- write.profile(missing_chars)
    
    if ( !silent ) {
      warning("\nThere were unknown characters found in the input data.\nCheck output$errors for a table with all problematic strings.")
      
    }
  } else {
    problems <- NULL
    problemChars <- NULL
  }
  
  # Reorder profile according to order and add frequency of rule-use
  # frequency <- head(frequency, -1)
  
  profile.out <- data.frame(profile[graph_order,]
                        , stringsAsFactors = FALSE
                        )
  if (ncol(profile.out) == 1) {colnames(profile.out) <- "Grapheme"}
  profile.out <- cbind(matched_rules, profile.out)
  
  # --------------
  # output as list
  # --------------
  
  result <- list(strings = strings.out
                 , profile = profile.out
                 , errors = problems
                 , missing = problemChars
                  )
  
  if (is.null(file.out)) {
    
    return(result)
    
  } else {   
    
  # ---------------
  # output to files
  # ---------------
    
    # file with tokenization is always returned
    write.table(  strings.out
                  , file = paste(file.out, "_strings.tsv", sep = "")
                  , quote = FALSE, sep = "\t", row.names = FALSE)
    
    # file with orthography profile
    write.table(  profile.out
                  , file = paste(file.out, "_profile.tsv", sep="")
                  , quote = FALSE, sep = "\t", row.names =  FALSE)
    
    # additionally write tables with errors when they exist
    if ( !is.null(problems) ) {      
      write.table(  problems
                    , file = paste(file.out, "_errors.tsv", sep = "")
                    , quote = FALSE, sep = "\t", row.names = TRUE)   
      write.table(  problemChars
                    , file = paste(file.out, "_missing.tsv", sep = "")
                    , quote = FALSE, sep = "\t", row.names = FALSE)
    }
    
    return(invisible(result))
    
  } 
}
