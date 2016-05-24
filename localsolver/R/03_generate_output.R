#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------

generate.output <- function(lsp, inp.append, indexFromZero) {
  range_str <- function(length) {
    return(sprintf("%d..%d", ifelse(indexFromZero, 0, 1), ifelse(indexFromZero, length - 1, length)))
  }
  
  
  inp.append("function __auto_model() {\n")
  inp.append("\t// this function is called at the very end of model function\n")
  for(varName in names(lsp$output.exprs)) {
    desc <- lsp$output.exprs[[varName]]
    
    if(length(desc$dimensions) == 3) {
      # 3-dim array
      inp.append(sprintf('\t%s[i in %s][j in %s][k in %s] <- %s[i][j][k];\n',
                         varName, 
                         range_str(desc$dimensions[1]), range_str(desc$dimensions[2]), range_str(desc$dimensions[3]), 
                         desc$expr))        
    } else if(length(desc$dimensions) == 2){
      # matrix
      inp.append(sprintf('\t%s[i in %s][j in %s] <- %s[i][j];\n', 
                         varName, range_str(desc$dimensions[1]), range_str(desc$dimensions[2]), desc$expr))
      
    } else if (length(desc$dimensions) == 1) {
      if(desc$dimensions > 1){
        # vector
        inp.append(sprintf('\t%s[i in %s] <- %s[i];\n', varName, range_str(desc$dimensions), desc$expr))
      }else{
        # single number
        inp.append(sprintf("\t%s <- %s;\n", varName, desc$expr))
      }
    } else {
      stop(sprintf("Invalid output expression(%s) dimensions", desc$expr))
    }
  }
  inp.append("}\n\n")
  
  inp.append("function output() {\n")  
  for(varName in names(lsp$output.exprs)) {
    desc <- lsp$output.exprs[[varName]]
    
    if(length(desc$dimensions) == 3) {
      # 3-dim array
      inp.append(sprintf('\tprint("$:"); print[i in %s][j in %s][k in %s](getValue(%s[i][j][k]), " "); println();\n',
                         range_str(desc$dimensions[1]), range_str(desc$dimensions[2]), range_str(desc$dimensions[3]), varName))        
    } else if(length(desc$dimensions) == 2){
      # matrix
      inp.append(sprintf('\tprint("$:"); print[i in %s][j in %s](getValue(%s[i][j]), " "); println();\n',
                         range_str(desc$dimensions[1]), range_str(desc$dimensions[2]), varName))        
    } else if (length(desc$dimensions) == 1) {
      if(desc$dimensions > 1){
        # vector
        inp.append(sprintf('\tprint("$:"); print[i in %s](getValue(%s[i]), " "); println();\n', 
                           range_str(desc$dimensions), varName))
      }else{
        # single number
        inp.append(sprintf('\tprintln("$:" + getValue(%s));\n', varName))
      }
    } else {
      stop(sprintf("Invalid output expression(%s) dimensions", desc$expr))
    }
  }  

  inp.append("}\n\n")
}
