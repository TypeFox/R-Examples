#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------


generate.param <- function(lsp, inp.append){
  if (!is.null(lsp$functions$param)) {
    inp.append("\n\nfunction __auto_param() {\n")
    inp.append("\t// this function is called at the very end of provided param function\n")
  } else {
    inp.append("\n\nfunction param() {\n")
  }
  
  for(paramName in names(lsp$params)){
    if(paramName == "indexFromZero") { next }
    
    paramValue <- lsp$params[[paramName]]
    if(!is.null(paramValue)){      
      if(length(paramValue) == 1){
        inp.append(sprintf("\t%s = %d;\n", paramName, paramValue))          
      }else{
        paramValueVector <- paste(paramValue, collapse = ", ")
        inp.append(sprintf("\t%s = {%s};\n", paramName, paramValueVector))
      }
    }
  }

  inp.append("}\n\n")
}
