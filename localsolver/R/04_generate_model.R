#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------

generate.model <- function(lsp, inp.append) {
  result <- list()
  
  for(funcName in names(lsp$functions)) {
    funcDesc <- lsp$functions[[funcName]]

    funcStartLN <- inp.append(sprintf("%s {", funcDesc$decl)); 
    inp.append(funcDesc$body)
    
    if (funcName == "model") {
      inp.append("\n\t__auto_model();\n")
    } else if (funcName == "param") {
      inp.append("\n\t__auto_param();\n")
    }
        
    funcEndLN <- inp.append("}\n\n")
    result[[funcName]] <- list(start.line = funcStartLN, end.line = funcEndLN)
  }
  
  return(result)
}
