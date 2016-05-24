#################################################################
# 
# File:         read.px.R
# Purpose:      reads a PC-Axis file into R
#
# Created:      20110618
# Authors:      fvf, cjgb, opl
#
# Modifications: 
#		20111210, cjgb: in the data string, "-" may represent the value 0
#		20111210, cjgb: fixing the strsplit when the split character is contained in the data part
#		20120329, cjgb: number strings in the DATA part can contain ";" as separators.
#				Although deprecated, cases still lurk.
#               20130228, cjgb: There can be ; inside quoted strings that create chaos
#               20130608   fvf: Ability to read files with keys in data area.
#                               ":"  added to defaut na.string (EuroStat files)
#               20130624:  use str_split (line 91) to read DATA area
#               20130917, cjgb: changes to prevent errors with EOL characteres
#               20131115, cjgb: some files do not have heading (or stub): only one of 
#                               them is really required
#               20131118, cjgb: fixed a bug happening when missing (i.e. "..") was the last value in DATA
#                               fixing it required that the last quote was not eliminated (same for first quote)
#################################################################

read.px <- function(filename, encoding = "latin1", 
                    na.strings = c('"."', '".."', '"..."', '"...."', '":"')) {

    ## auxiliary functions ##

    clean.spaces <- function(x){
        gsub("^[[:space:]]+|[[:space:]]+$", "", x) # discards heading|trailing whitespace
    }

    get.attributes <- function(x){
        x <- gsub( "([A-Z-]*)\\((.*)\\).*", "\\1;\\2", x ) ## separates label-attribute with ";"
        x <- ldply(strsplit(x, ";"), 
                   function(y) c(y, "value")[1:2])
    }

    break.clean <- function(x) {
        x <- clean.spaces( strsplit(x, split = '\\"')[[1]] )    ## breaks by '"'
        x[! x %in% c("," , "")]                                 ## and drops spurious seps
    }


    ## end: auxiliary functions ##

    a <- scan(filename, what = "character", sep = "\n", quiet = TRUE, fileEncoding = encoding)

    # modification by  fvf: 130608 
    a <- paste(a, collapse = "\n")       # Se mantienen "CR/LF luego se quitaran selectivamente

    tmp <- strsplit( a, "DATA=" )[[1]]
    tmp[1] <- gsub("\n", " ", tmp[1])    # fvf[130608]: elimina CR de la cabecera
    
    tmp[2] <- gsub(";.*", "", tmp[2])    # removing ";" (and everything following it) within DATA number strings
    a <- paste(tmp[1], "DATA=", tmp[2], sep = "")

    ## modification by cjgb, 20130228 concerning line separators within quoted strings
    ## ; is the logical line end in px files
    ## so we should do:
    ## a <- unlist(strsplit(a, ";"))	
    ## but there might be ; inside quoted strings
    ## so we need the following workaround:

    punto.coma <- str_locate_all(a, ";")[[1]][,1]	# where the ";" are
    comillas   <- str_locate_all(a, '"')[[1]][,1]	# where the '"' are

    ## ";" not after an odd number of '"'
    ## these are the proper "cuts"
    cortes     <- Filter( function(x) sum(comillas < x) %% 2 == 0, punto.coma )		
    
    a <- str_sub(a, c(1, cortes + 1), c(cortes - 1, str_length(a)))
    a <- a[!is.na(a)]
    a <- a[a != ""]
    
    ## end of modification by cjgb, 20130228 concerning line separators within quoted strings
    
   
    # change strsplit by str-split. In big px-files:
    #  "Error: C stack usage is too close to the limit"
    a <- do.call(rbind, str_split(a, "=", n = 2))   
                                                        
    a <- data.frame(cbind(get.attributes(a[, 1]), a[, 2]))
    colnames(a) <- c("label", "attribute", "value")
  
    ## build a px object: list with px class attribute ##  
    
    a$label     <- make.names(clean.spaces(a$label))
    a$attribute <- make.names(clean.spaces(gsub('\\"', "", a$attribute)))
    
    # need to avoid that quotes are removed in DATA part because of a bug:
    # a case was reported where the data part ended in ".." and the last quote was erased
    # and this affected the scan function below
    a.data                     <- as.character(a[a$label == "DATA", "value"])
    a.value                    <- gsub('^\\"|\\"$', "", a$value)   # removes " at beginning / end
    a.value[a$label == "DATA"] <- a.data
    names(a.value)             <- a$attribute
    
    px <- tapply(a.value, a$label, as.list)    

    ## these metadata keys contain vectors (comma separated)
    ## we need to split them (and clean the mess: extra spaces, etc.)
    px$STUB$value    <- if(!is.null(px$STUB))    make.names(break.clean(px$STUB$value))
    px$HEADING$value <- if(!is.null(px$HEADING)) make.names(break.clean(px$HEADING$value))

    px$VALUES <- lapply(px$VALUES, break.clean)
    px$CODES  <- lapply(px$CODES,  break.clean)

    
    #### read the data part into a 'melted' dataframe ###
    
    ## there are two cases: files with/without KEYS keyword
    ## which need to be processed independently

    # fvf[130608]: add to to read files with keys in data area 
    
    if ("KEYS" %in% a$label ){
      
      ## read the whole block
      tc <- textConnection(px$DATA$value); on.exit( close(tc) )
      raw <- read.table(tc, sep = ",", colClasses = "factor")
      
      ## extract and process the data part (the numbers)
      data.part <- as.character( raw[, ncol(raw)] )          # numbers (last column of the data.frame)
      data.part <- gsub('"-"', 0, data.part)                 # 0's might be encoded as "-"
      data.part <- scan(text = data.part, na.strings = na.strings, quiet = T)
      
      ## extract and process the keys part (it needs to be staked a number of times, 
      ##  as many as there are entries in the data vector in each row in the block)
      keys.part <- raw[, -ncol(raw), drop = FALSE]    
      keys.part <- keys.part[ rep(1:nrow(keys.part), each = length(data.part) / nrow(keys.part) ), ]
      colnames(keys.part) <- names(px$KEYS)
      
      ## change CODES (if any) in keys part to VALUES (consistency issue)
      for (col.name in colnames(keys.part)[unlist(px$KEYS) == "CODES"])
        keys.part[[col.name]] <- mapvalues(keys.part[[col.name]], 
                                           from = px$CODES[[col.name]], 
                                           to   = px$VALUES[[col.name]])
      
      
      ## extract and process the variables that are not keys
      no.keys.part <- px$VALUES[ setdiff( names(px$VALUES), names(px$KEYS) ) ]                       
      no.keys.part <- expand.grid(rev(no.keys.part))
      
      ## put everything together & cleanup
      px$DATA$value <- data.frame( keys.part, 
                                   no.keys.part, 
                                   value = data.part,
                                   row.names = NULL)
    }  
    else
    {
      tmp <- gsub('"-"', 0, px$DATA$value)        # 0 can be encoded as "-"
      tmp <- gsub("\n", " ", tmp)                 # delete CR/LF of DATA area fvf[130608]

      tc  <- textConnection(tmp); on.exit( close(tc) )
      raw <- scan(tc, na.strings = na.strings, quote = NULL, quiet = TRUE)
      #raw <- scan(tc, na.strings = na.strings, quiet = TRUE)
      
      names.vals <- c( rev(px$HEADING$value), rev( px$STUB$value ) )
      output.grid <- data.frame(do.call(expand.grid, px$VALUES[names.vals]))
      
      # sanity check: avoids the problem of "reclycling" of values if
      # the ratio of lenghts of variables and values is an exact integer      
      if (nrow(output.grid) != length(raw))
        stop( "The input file is malformed: data and varnames length differ" )
      
      px$DATA$value           <- data.frame(output.grid, raw)
      colnames(px$DATA$value) <- c(names.vals, "value")    

    }
    
    class(px) <- "px"
    px
}

