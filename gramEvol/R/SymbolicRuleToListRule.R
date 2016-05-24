 is.call.primitive <- function(x) {
 # returns true if the given input (a "call" object) is a .primitive
   if (length(x) != 1) {
     return(FALSE)
   }
   
   str = deparse(x)
   if (length(str) != 1) {
     return(FALSE)
   }
   
   # it is not a string if it contains paranthesis
   if (length(grep("(", str, fixed=TRUE)) != 0) {
     return(FALSE)
   }
   
   return(is.primitive(try(eval(x), silent=TRUE)))
 }

SymbolicRuleToListRule <- function(ruleDef) {
  # ruleDef: Receives the symbolically defined rule
  # returns: a grammar rule in string list format with gt and lt escaped
  
  # extract rule names
  ruleDefIndex = names(ruleDef)
  
  # create rule name escape list
  escape_rule <- function(x) as.call(parse(text = paste0(".GE_lt.", as.character(x), ".GE_gt.")))
  ruleEscapeList = lapply(ruleDefIndex, escape_rule)
  names(ruleEscapeList) = ruleDefIndex
  
  
  # for every rule in expression
  for (i in 1:length(ruleDef)) {
    
    if (!"GERule" %in% class(ruleDef[[i]])) {
      stop("Invalid Grammar Ruleset")
    }
    
    # if the input was not a string to start with
    if (!"GEStringRule" %in% class(ruleDef[[i]])) {
      # replace the rule with string version
      for (j in 1:length(ruleDef[[i]])) {
        ruleDef[[i]][[j]] = SymbolicRuleToString(ruleDef[[i]][[j]], ruleDefIndex, ruleEscapeList)
      }
    } 
    
    # use list(name, list(rule1, rule2, ...)) structure
    ruleDef[[i]] = list(ruleDefIndex[i] , unclass(ruleDef[[i]]))
  }
  
  return (ruleDef)
}

SymbolicRuleToString <- function(current_rule, ruleDefIndex, ruleEscapeList) {

  # the escaping function that replaces rules in situ
  escape_expression_rules <- function(expr) eval(substitute(substitute(e, ruleEscapeList), list(e = expr)))
  
  # do not process primitives
  if (is.call.primitive(current_rule)) {
    current_rule = paste0('`', escape.gt.lt(current_rule), '`')
  } else {
    # add ".GElt." and ".GEgt." around parsed blocks to mark them as future substitutables
    current_rule = escape_expression_rules(current_rule)
    
    # convert current rule to string properly
    current_rule = escape.gt.lt(GERule.EscapeDot(current_rule))
    
    # unescape expression, by replacing ".GE_lt." and ".GE_gt." with < and >
    for (r in ruleDefIndex) {
      rx = paste0(".GE_lt.", r, ".GE_gt.()")
      ry = paste0("<", r, ">")
      current_rule = gsub(rx, ry, current_rule, fixed=TRUE)
    }

    # concatenate if multiline
    current_rule = paste(current_rule, sep="", collapse="\n")
  }
  
  return (current_rule)
}
