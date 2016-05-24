
# Split up <expression> into symbols.
# Returns a vector of symbols.
scan <- function(expression, lowerCase=FALSE)
{
  if (lowerCase)
    expression <- tolower(expression)
  
  # add extra whitespace to brackets and operators
  
  expression <- gsub("(\\(|\\)|\\[|\\]|\\||\\&|\\!|\\=|\\.\\.|[a-zA-Z0-9_]*)"," \\1 ", expression)

  # strip multiple whitespace characters
  expression <- gsub("[ ]+", " ", expression)
  expression <- gsub("^[ ]+", "", expression)
  expression <- gsub("[ ]+$", "", expression)
  
  # split up at whitespace positions
  res <- strsplit(expression, " ", fixed=TRUE)[[1]]
  return(res)
}

# Parse a Boolean function in <expression>,
# and build a corresponding parse tree.
parseBooleanFunction <- function(expression, varNames=c(), lowerCase=FALSE)
{
  minArgCounts <- c("&"=2,"|"=2, "all" = 1, "any" = 1, "maj" = 1, 
                    "sumis" = 2, "sumgt" = 2, "sumlt" = 2,
                    "timeis" = 1, "timelt" = 1, "timegt" = 1)
  maxArgCounts <- c("&"=Inf,"|"=Inf,"all" = Inf, "any" = Inf, "maj" = Inf, 
                    "sumis" = Inf, "sumgt" = Inf, "sumlt" = Inf,
                    "timeis" = 1, "timelt" = 1, "timegt" = 1)
  
  isBooleanArgument <- function(arg)
  {
    return (arg$type %in% c("atom","operator") ||
            arg$type == "constant" && arg$value %in% c(0,1))
  }
  
  checkPredicate <- function(pred)
  {
    if (pred$operator %in% c("&","|","maj"))
    {
      for (el in pred$operands)
      {
        if (!isBooleanArgument(el))
        {
          stop(paste("Invalid argument to \"",pred$operator,"\": ", stringFromParseTree(el),sep=""))
        }
      }
    }
    else
    if (pred$operator %in% c("sumis","sumgt","sumlt"))
    {
      for (el in pred$operands[-length(pred$operands)])
      {
        if (!isBooleanArgument(el))
        {
          stop(paste("Invalid argument to \"",pred$operator,"\": ", stringFromParseTree(el),sep=""))
        }
      }
      if (pred$operands[[length(pred$operands)]]$type != "constant")
        stop(paste("Last argument of \"",pred$operator,"\" must be an integer value!"))
    }
    else
    if (pred$operator %in% c("timeis"))
    {
      if (pred$operands[[1]]$type != "constant" || pred$operands[[1]]$value < 0)
        stop(paste("Argument of \"",pred$operator,"\" must be a non-negative integer value!"))
    }
  }
  
  # internal function to step forward one symbol
  advanceSymbol <- function(errorIfEnd = TRUE)
  {
    pos <<- pos + 1
    #print(symbols[-(1:pos-1)])    
    if (pos > length(symbols))
      if (errorIfEnd)
        stop("Unexpected end of expression!")
      else
        return(NA)
    else
      return(symbols[pos])
  }
  
  getCurrentSymbol <- function()
  {
    #print(symbols[-(1:pos-1)])
    if (pos > length(symbols))
      return(NA)
    else
      return(symbols[pos])
  }
  
  # internal function to preview the next symbol
  previewSymbol <- function()
  {
    if (pos + 1 > length(symbols))
      return(NA)
    else
      return(symbols[pos + 1])
  }
  
  getSymbolPos <- function()
  {
    return(pos)
  }
  
  setSymbolPos <- function(newPos)
  {
    pos <<- newPos
  }
  
  # internal function to parse a sub-expression in brackets
  parseExpression <- function(temporalVariables = c())
  {
    operators <- c()
    children <- c()
    symb <- advanceSymbol()
    
    while (TRUE)
    {
      if (symb == "(")
        # a new sub-expression
        children[[length(children)+1]] <- parseExpression(temporalVariables=temporalVariables)
      else
      if (symb == "!")
        # a negation
        children[[length(children)+1]] <- parseNegation(temporalVariables=temporalVariables)
      else
      if (regexpr("^[a-zA-Z0-9_]*$",symb) == -1)
        # something unexpected happened
        stop(paste("Unexpected symbol:",symb))
      else
        children[[length(children)+1]] <- parseAtomOrPredicate(symb, temporalVariables=temporalVariables)
      
      symb <- advanceSymbol(FALSE)  
      
      if (is.na(symb) || symb %in%  c(")",","))
        # end of expression was reached      
        break
        
      if (!symb %in% c("&","|"))
        stop("Operator expected!")
      
      operators <- c(operators,symb)
      
      symb <- advanceSymbol()
    }
    
    if (length(children) == 1)
    {
      # an operator with only one child
      return(children[[1]])
    }
    else
    if (length(unique(operators)) == 1)
    {
      # an n-ary operator
      res <- list(type="operator",operator=operators[1],
                  negated=FALSE, operands=children)
                  
      checkPredicate(res)                  
      return(res)
    }
    else
    # a mixture of multiple operators => ensure correct precedence
    {
      i <- 1
      startPos <- NA
      operators <- c(operators, "|")
      operands <- list()

      while (i <= length(operators))
      # identify AND operators and move them to a subtree
      {
        if (operators[i] == "&" && is.na(startPos))
          # start of a conjunction
          startPos <- i
        if (operators[i] == "|" && !is.na(startPos))
        # end of a conjunction => create subtree
        {
          subOp <- list(type="operator", operator="&",
                        negated=FALSE, operands=children[startPos:i])
          operands[[length(operands)+1]] <- subOp
          startPos <- NA
        }
        else
        if (is.na(startPos))
          operands[[length(operands)+1]] <- children[[i]]
        i <- i + 1  
      }
      res <- list(type="operator",operator="|",
                  negated=FALSE, operands=operands)
      checkPredicate(res)                
      return(res)
    }
  }
  
  # internal function to parse atoms or temporal predicates
  parseAtomOrPredicate <- function(name, temporalVariables=c())
  {
    if (regexpr("^[0-9]*$",name) != -1 || name %in% c("true","false"))
    # a (Boolean or integer) constant
    {
      if (name == "true")
        value <- 1L
      else
      if (name == "false")
        value <- 0L
      else
        value <- as.integer(name)
         
      return(list(type="constant",
                  value=value,
                  negated=FALSE))
    } 
    if (tolower(name) %in% names(minArgCounts))
    # a temporal predicate
    {
      name <- tolower(name)
      # translate all into & and any into |
      operator <- switch(name,
                         all = "&",
                         any = "|",
                         name)

      nextSym <- previewSymbol()
      if (!is.na(nextSym) && nextSym == "[")
      # a temporal iterator definition follows
      {
        # discard bracket
        advanceSymbol()
      
        varName <- advanceSymbol()
        
        if (regexpr("^[a-zA-Z0-9_]*$",varName) == -1)
          # something unexpected happened
          stop(paste("Unexpected symbol:",varName))
        
        symb <- advanceSymbol()   
        if (symb != "=")
        {
          stop(paste("Expected \"=\" instead of \"",symb,"\" in temporal predicate!"))
        }
        
        # parse start time point
        startTime <- parseTimePoint(temporalVariables=temporalVariables, TRUE)
        symb <- getCurrentSymbol()
        
        if (symb != "..")
        {
          stop(paste("Expected \"..\" instead of \"",symb,"\" in temporal predicate!"))
        }

        # parse end time point
        endTime <- parseTimePoint(temporalVariables=temporalVariables, TRUE)        
        symb <- getCurrentSymbol()
                    
        if (symb != "]")
            stop(paste("Unexpected symbol:",symb))
        
      }
      else
      # no temporal iterator defined
      {
        startTime <- 0
        endTime <- 0
        varName <- NULL
      }            

      symb <- advanceSymbol()      
      if (symb != "(")
          stop(paste("Unexpected symbol:",symb))

      subElements <- c()
            
      while(TRUE)
      # parse parameters of predicate
      {
        oldPos <- getSymbolPos()
        
        for (time in startTime:endTime)
        # expand arguments according to temporal range
        {
          tempVars <- temporalVariables
          
          if (!is.null(varName))
          {
            if (length(tempVars) > 0 && !is.na(tempVars[varName]))
            {
              stop(paste("Duplicate variable ",varName,"! Please use only unique variable names!",sep=""))
            }
            tempVars[varName] <- time
          }
          # re-parse the argument with a new time point => reset parser position
          setSymbolPos(oldPos)          
          subElement <- parseExpression(temporalVariables=tempVars)
          
          subElements <- c(subElements, list(subElement))  
          if (subElement$type == "constant")
            break
        }
        
        subElements <- unique(subElements)
        
        symb <- getCurrentSymbol()
        if (is.na(symb))
          stop("Unexpected end of expression!")
        
        if (symb == ")")
          break
          
        if (symb != ",")
          stop ("Expected \",\" or \")\"!")
        
      }
      
      # check if number of argument matches
      if (length(subElements) < minArgCounts[name] || length(subElements) > maxArgCounts[name])
        stop(paste("Invalid number of arguments to ",name,": ",length(subElements),sep=""))
      
      res <- list(type="operator",operator=operator,
                  negated=FALSE, operands=subElements)

      # check if argument types match
      checkPredicate(res)
                  
      return(res)
    }
    else
    # an atom with or without temporal information
    {
      nextSym <- previewSymbol()
      if (!is.na(nextSym) && nextSym == "[")
      # the atom has a time point specification
      {
        symb <- advanceSymbol()
        time <- parseTimePoint(temporalVariables=temporalVariables)
        symb <- getCurrentSymbol()        
                  
        if (symb != "]")
          stop(paste("Unexpected symbol:",symb))
      }
      else
      # no time point specified => default to previous time point
        time <- -1L
        
      if (length(varNames) > 0)
      # check if variable is known
      {
        varIndex <- which(varNames == name)
        if (length(varIndex) == 0)
          stop(paste("Unknown variable:",name))
        negated <- FALSE
        
        return(list(type="atom",
                    name=name,
                    index=as.integer(varIndex),
                    time=as.integer(time), 
                    negated=negated))
          
      }
      else
      {
          return(list(type="atom",
                      name=name,
                      time=as.integer(time), 
                      negated=FALSE))
      }

    }
  }
  
  # internal function to parse a time specification
  parseTimePoint <- function(temporalVariables=c(), allowPositive=FALSE)
  {
    symb <- advanceSymbol()
    time <- 0
    while (TRUE)
    {
      sign <- 1
      
      if (symb == "+")
      # ignore positive sign
        symb <- advanceSymbol()
      else
      if (symb == "-")
      # set multiplier to negative sign
      {
        sign <- -1
        symb <- advanceSymbol()
      }
      
        
      if (symb %in% names(temporalVariables))
      # a variable was found => add its value
      {
        time <- time + sign * temporalVariables[symb]
        symb <- advanceSymbol()
      }      
      else
      # this has to be an integer time value
      {
        suppressWarnings(addTime <- as.integer(symb))
        
        if (is.na(addTime))
        {
          stop(paste("Invalid time specification: ",symb, 
                     ". Time specifications must be integers <= -1)", sep=""))
        }
        time <- time + sign * addTime
        symb <- advanceSymbol()
      }      
      if (symb %in% c("]",".."))
        break
    }
    
    if (time > -1 && !allowPositive)
    {
      stop(paste("Invalid time specification: ",time, 
                 ". Time specifications must be integers <= -1)", sep=""))
    }
    return(as.integer(time))
  }
  
  # Internal function to parse a negation
  parseNegation <- function(temporalVariables=c())
  {
    symb <- advanceSymbol()
    if (symb == "(")
    {
      res <- parseExpression(temporalVariables)
    }
    else
    if (symb %in% c(")","&","|","]","[","="))
      stop(paste("Unexpected symbol:",symb))
    else
      res <- parseAtomOrPredicate(symb, temporalVariables=temporalVariables)
    res$negated <- !res$negated
    return(res)
  }

  tryCatch(
  {
    symbols <- scan(expression, lowerCase=lowerCase)
    
    pos <- 0
    res <- parseExpression()
  },
  error = function(e)
  {
    stop(sprintf("Error parsing expression \"%s\": %s", expression, e$message))
  })
  if (!isBooleanArgument(res))
    stop(paste("Non-Boolean formula ",stringFromParseTree(res),"!",sep=""))
  
  return(res)
}

# Regenerate the expression string 
# from the parse tree <tree>
stringFromParseTree <- function(tree)
{
  res <- switch(tree$type,
    operator = 
    {
      if (tree$operator %in% c("|","&"))
      {
        paste({if (tree$negated) "!" else ""},
              "(",
              paste(sapply(tree$operands,stringFromParseTree), 
                    collapse=paste(" ",tree$operator," ",sep="")),
              ")", sep="")
      }
      else
      {
        paste({if (tree$negated) "!" else ""},
              tree$operator,"(",paste(sapply(tree$operands,stringFromParseTree), collapse=", "),")",sep="")
      }
    },    
    atom = paste({if (tree$negated) "!" else ""},
                 tree$name,
                 {if (tree$time != -1) paste("[",tree$time,"]",sep="") else ""},
                 sep=""),    
    constant = paste({if (tree$negated) "!" else ""},
                     tree$value,
                     sep="")
  )                     
  return(res)  
}

