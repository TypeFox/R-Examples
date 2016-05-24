
.replaceVarsForContent <- function(text, vars){
  # Replace the variables in vars (in text are preceded by a $ symbol) by its content
  # in the text AT ONCE (it works independently of the variable
  # names and the content of each variable for the replacement)
  #
  # text The original text of the template
  # vars A list in which each element have the format "name_of_variable=text_to_replace"
  
  # Just in case some variable in vars is not present in text, we remove it previously from vars.
  nameVars <- paste0("\\{\\{",names(vars),"\\}\\}")
  valid <- sapply(nameVars,FUN=function(v) grepl(v,text))
  vars <- vars[valid]
  nameVars <- nameVars[valid]
  
  # FOR DEBUG PURPOSE: If there are variables in text that are not present in vars, we warn it
  varsInTextMatching <- gregexpr("\\{\\{[^\\{\\}]+\\}\\}",text)[[1]]
  varsInText <- mapply(FUN = function(ind,len){
    substr(text,ind+2,ind+len-3)
  }, varsInTextMatching, attr(varsInTextMatching,"match.length"))
  if(any(!(varsInText %in% names(vars))))
    warning(.callWarningMessage("varsInTemplateNotExist", toString(varsInText[!(varsInText %in% names(vars))])))
  
  # For vars whose field is NULL, we assign its var name surrounded by '{{' '}}' as field
  isNull <- sapply(vars, is.null)
  vars[isNull] <- paste0("{{",names(vars)[isNull],"}}")
  
  # We extract the occurrences of each variable in the text and where they are
  indicesByVars <- sapply(names(vars), simplify = FALSE, FUN =function(v){ 
    res <- gregexpr(sprintf("\\{\\{%s\\}\\}",v), text)[[1]]  # Find indices of matches
    attr(res,"name") <- v # We assign the name of the variable to the element (we will need it later)
    res
  })
  # We extract the length of each occurrence (we know it by the length of names(vars), but we want them also
  # in the same order as indicesByVars), so we obtain them with this sencence
  longitudByVars <- sapply(indicesByVars,FUN=function(v){attr(v,"match.length")})
  # We extract the var name of each occurrence (in the same order as indicesByVars)
  namesByVars <- sapply(indicesByVars,FUN=function(v){rep(attr(v,"name"),length(v))})
  
  # Now we transform the previous data structure un plain arrays, and ordered by index occurrence.
  indices <- unlist(indicesByVars)
  newOrder <- order(indices)
  indices <- indices[newOrder]
  longitudes <- unlist(longitudByVars)
  longitudes <- longitudes[newOrder]
  names <- unlist(namesByVars)
  names <- names[newOrder]
  
  # We will split the text in one substrings per variable occurrence in the text.
  # Each one starts just afrer and ends just before a variable name.
  
  # We obtain the first and last index for each substring
  inicio <- c(1,indices+longitudes)
  fin <- c(indices,1+nchar(text))-1
  
  # We build a data.frame to deal with it in the next apply function.
  intervals <- data.frame(inicio=inicio, fin=fin)
  
  # Now split the text in the desired substrings
  substrings <- apply(intervals,1,FUN=function(v){substr(text,v[1],v[2])})
  
  paste0(paste0(substrings[-length(substrings)],unlist(vars[names]),collapse=""),substrings[length(substrings)])
  # The way it works, we have tu put a %s just after each substring unless the last one
  #formatedString <- paste0(paste0(substrings[-length(substrings)],rep("%s",length(substrings)-1),collapse=''),substrings[length(substrings)])
  
  # Finally, we use the sprintf function to generate the output of this function.
  #do.call(function(...) sprintf(formatedString,...), vars[names])
}

.sanitizeLatexCode <- function(text){
  newText <- text
  
  # Single and double quotes in latex are represented by `' and ``". Because
  # we represent that internally as '' and "", we need to change them when they
  # are use for latex code
  newText <- gsub("'([^']*)'","`\\1'",newText)
  newText <- gsub("\\\"(([^\\][^\"])*)\\\"","``\\1\\\"",newText)
  
  # > and < symbols are used as arrows. For that, we transform -> into 
  # $\rightarrow$ and <- into $\leftarrow$
  newText <- gsub("<-","$\\\\leftarrow$",newText)
  newText <- gsub("->","$\\\\rightarrow$",newText)
  
  newText
}

