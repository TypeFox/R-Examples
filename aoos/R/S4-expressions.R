ExpressionTree <- function(.mc, .where) {
  # Function to construct the central language object. Will be processed to be
  # used to construct argument lists for S4s set-functions (setClass,
  # setGeneric, setMethod)
  # The systax is as follows:
  # [name1 : ... : nameN-1 : ] nameN([<argList>]) %<fun>% expr
  # name1 to nameN will be the names of the generic / type and inheritance
  #   structure
  # <argList> the arguments. Can be usual argument expressions (argName =
  # defaultValue) or 'type expressions' (argName ~ typeName).
  # expr is the function body. 
  
  .seperate <- function(x, delim) {
    x <- sub(delim, "#", x) # this is saver because of args with defaults
    lArgs <- lapply(x, . %>% splitTrim("#"))
    args <- sapply(lArgs, . %>% .[2])
    names(args) <- argNames
    args
  }
  
  .processClassUnions <- function(nameTypeExpr) {
    sapply(nameTypeExpr, USE.NAMES = FALSE, function(nte) {
      if (grepl("#", nte)) {
        classes <- splitTrim(nte, "#")
        nameClassUnion <- paste(classes, collapse = "OR")
        try(silent = TRUE,
            setClassUnion(nameClassUnion, classes, topenv(.where))
        )
        nameClassUnion
      } else {
        nte
      }
    })
  }
  
  .collapseMatchingParen <- function(x) {
    # collapses elements of x, if they have matching parenthesis. There is only
    # nice code comming in because the R parser won't allow for anything else.
    parenComplete <- function(x) {
      sapply(strsplit(x, ""), function(x) sum(x == "(") == sum(x == ")"))
    }
    
    collapseTwoElements <- function(x) {
      pos <- Position(Negate(parenComplete), x)
      x[pos] <- paste(x[pos:(pos + 1)], collapse = ", ")
      x[-(pos+1)]
    }
    
    if (all(parenComplete(x))) {
      x
    } else {
      Recall(collapseTwoElements(x))
    }
  }
  
  .lhs <- deparse(.mc$lhs) %>% paste(collapse = "") %>% sub("\\n", "", .)
  body <- deparse(.mc$rhs)
  names <- deleteInParan(.lhs) %>% gsub("\\|", "#", .) %>% splitTrim(":")
  names <- deleteQuotes(names) %>% rev %>% .processClassUnions
  args <- deleteBeforeParan(.lhs) %>% 
    deleteEnclosingParan %>% splitTrim(",") %>% .collapseMatchingParen
  argNames <- sapply(args, . %>% splitTrim("=|~") %>% .[1], USE.NAMES = FALSE)
  argDefaults <- args %>% .seperate("=")
  argClasses <- args %>% .seperate("~") %>% deleteQuotes
  argClasses <- gsub("\\|", "#", argClasses) %>% .processClassUnions
  
  retList("ExpressionTree")
  
}

ClassExpressionTree <- function(.exprTree, where) {
  # Constructs the argument list for setClass
  
  .makeExpression <- function(funName, args) {
    parse(text = paste0(
      funName, "(", paste(args, collapse = ", "), ")"
    ))
  }
  
  .protoIsGood <- function(proto) 
    proto@dataPart || !identical(proto@slots, character())
  
  .localEval <- function(expr) eval(expr, envir = where)
  
  .getExplicitClass <- function(.exprTree) {
    slots <- .exprTree$argClasses
    slots[is.na(slots)] <- "ANY"
    slots <- slots[!(names(slots) == ".Data")]
    slots
  }
  
  .getPrototype <- function(.exprTree) {
    args <- .exprTree$argDefaults[!is.na(.exprTree$argDefaults)]
    args <- paste(names(args), args, sep = "=")
    args %<>% sub(".Data( )?\\=", "", .)
    .localEval(.makeExpression("prototype", args))
  }
  
  .mergeProtoClasses <- function(slots, proto) {
    protoClasses <- sapply(attributes(proto@object), . %>% class %>% `[`(1))
    slots[names(protoClasses)] <- protoClasses[names(protoClasses)]
    slots[slots == "name"] <- "ANY"
    slots
  }
  
  slots <- .getExplicitClass(.exprTree)
  proto <- .getPrototype(.exprTree)
  slots <- .mergeProtoClasses(slots, proto)
  if (!.protoIsGood(proto)) rm("proto")
  
  Class <- .exprTree$names[1]
  contains <- if (is.na(.exprTree$names[2])) 
    character() else 
      .exprTree$names[-1]
  
  retList("ClassExpressionTree")
  
}

InitMethodExpressionTree <- function(.exprTree, where) {
  # Constructs the argument for setMethod for the generic 'initialize'
  
  .wrapInCurlyBraces <- function(x) {
    if (length(x) == 1) c("{", gsub("^\\{|\\}$", "", x), "}")
    else x
  }
  
  .body <- .wrapInCurlyBraces(.exprTree$body)
  .body[1] <- .body[1] %p0% "\n.Object <- callNextMethod()"
  
  f <- getGeneric("initialize", where = where)
  signature <- c(.Object = .exprTree$names[1])
  definition <- makeFunDef(c(".Object", "..."), .body, where)
  
  retList("MethodExpressionTree")
  
}

ConstExpressionTree <- function(.exprTree, envir) {
  # Constructs the argument list to construct the constructor function for the
  # type
  .getConstArgs <- function(.exprTree) {
    args <- .exprTree$argDefaults[names(.exprTree$argDefaults) %without% ".Data"]
    args <- ifelse(
      is.na(args), 
      names(args), 
      paste(names(args), .exprTree$argDefaults, sep = "=")
    )
    c(args, "...")
  }
  
  .getConstArgsNew <- function(.exprTree) {
    args <- names(.exprTree$argDefaults) %without% ".Data"
    c("'" %p0% .exprTree$names[1] %p0% "'", paste(args, args, sep = "="), "...")
  } 
  
  args <- .getConstArgs(.exprTree)
  body <- c("new(", paste(.getConstArgsNew(.exprTree), collapse = ", "), ")")
  
  retList("ConstExpressionTree")
  
}

GenericExpressionTree <- function(.mc, where) {
  # .mc is a match.call() from %g% and where the parent frame
  # This function will construct a list of arguments for a call to setGeneric
  
  .exprTree <- ExpressionTree(.mc, where)
  name <- .exprTree$names[1]
  valueClass <- if (is.na(.exprTree$names[2])) character() else .exprTree$names[2]
  def <- makeFunDef(.exprTree$args, .exprTree$body, where)
  
  retList("GenericExpressionTree")
  
}

MethodExpressionTree <- function(.mc, where) {
  # Constructs the argument list for method::setMethod
  .exprTree <- ExpressionTree(.mc, where)
  f <- eval(parse(text = .exprTree$names[1]), envir = where)
  .genericArgNames <- names(formals(f)) %without% "..."
  
  signature <- .exprTree$argClasses[!is.na(.exprTree$argClasses)]
  
  .args <- ifelse(
    .exprTree$argNames %in% .genericArgNames, 
    .exprTree$argNames,
    .exprTree$args
  )
  
  definition <- makeFunDef(.args, .exprTree$body, where)
  
  retList("MethodExpressionTree")
  
}

makeFunDef <- function(args, body, envir) {
  # args and body are expected to be character which will then be parsed to
  # R-Code
  args <- if (is.character(args)) 
    "(" %p0% paste(args, collapse = ", ") %p0% ")" else 
      stop(args, "is a ", class(args))
  
  body <- if (is.character(body)) 
    paste(body, collapse = "\n") else 
      stop(body, "is a", class(body))
  
  defCall <- "function" %p0% args %p0% body
  eval(parse(text = defCall), envir = envir)
  
}

# Helpers:
deleteQuotes <- . %>% sub("^[\"\']", "", .) %>% sub("[\"\']$", "", .)

deleteBeforeParan <- . %>% splitTrim("\\(") %>% { .[1] <- ""; . } %>% 
  paste0(collapse = "(")

deleteEnclosingParan <- . %>% sub("\\)$", "", .) %>% sub("^\\(", "", .)

deleteInParan <- . %>% gsub("\\(.*\\)", "", .)

splitTrim <- function(x, pattern) {
  strsplit(x, pattern) %>% unlist %>% trimws
}

"%p0%" <- function(lhs, rhs) paste0(lhs, rhs)

"%p%" <- function(lhs, rhs) paste(lhs, rhs)

"%without%" <- function(lhs, rhs) lhs[!(lhs %in% rhs)]