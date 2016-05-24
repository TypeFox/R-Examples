pirls.table <- 
  function(variable, by, data, export=FALSE, name= "output", folder=getwd()) {
    output <- intsvy.table(variable, by, data, config=pirls_conf)
    
    if (export)  {
      write.csv(output, file=file.path(folder, paste(name, ".csv", sep="")))
    }
    
    return(output)
}
