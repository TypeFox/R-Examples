generateModelFromPkModel <-  function(parameter,output)
{
  if (is.data.frame(parameter)){
    p_name <- names(parameter)
    i0 <- which(p_name %in% c("id", "group"))
    p_name <- p_name[-i0]
  } else if (is.vector(parameter))
    p_name <- names(parameter)
    else
    p_name   = parameter$name
  str1 <- paste(p_name,collapse=",")
  model_txt="
          [LONGITUDINAL]
          input = {param.list}
          EQUATION:
          output = pkmodel(param.list)
          "
  model_txt = gsub("param.list",str1,model_txt)
  if(length(output$name)==1){
    str2 = output$name[[1]]
  }else{
    str2 = paste("{",output$name[[1]],",",output$name[[2]],"}")
  }    
  model_txt = gsub("output",str2,model_txt)
  model     = inlineModel(model_txt)
  
  return(model)
}