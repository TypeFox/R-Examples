ElliOut <-
function(my.results, setup.file, out.file, use.ratio = FALSE){
  #Read section orientatiion file
  my.data <- read.table(file = setup.file, header = TRUE, sep = '\t', comment.char = "#")
  
  #Loop through each analysis to compile .elli file
  for(i in 1:length(my.results)){
    #Match row from orientation data to section name
    row.n <- which(my.data$section == my.results[[i]]@sectionName)
    #Error message for match fail
    if(as.numeric(row.n) == 0){stop("Cannot find section title.\nCheck file nomenclature")}
if(use.ratio){
  #Write header on first iteration
      if(i == 1){
        sink(file = out.file)
    cat("#", "strike","dip","rake","shape ratio ","","","","","","","","\n",sep="\t")
      }
  #Write data line
      cat(i, "\t", my.data$strike[row.n], "\t", my.data$dip[row.n], "\t", my.results[[i]]@vectorMean, "\t", my.results[[i]]@strainRatio, "\t\t", 1, "\n", sep="")

} else{
  #Write header on first iteration
      if(i == 1){
        sink(file = out.file)
    cat("#", "strike","dip","rake","long axis","short axis","","","","","","","\n",sep="\t")
      }
  #Write data line
      cat(i, "\t", my.data$strike[row.n], "\t", my.data$dip[row.n], "\t", my.results[[i]]@vectorMean, "\t", my.results[[i]]@rsAxes[1], "\t", my.results[[i]]@rsAxes[2], "\t", 1, "\n", sep="")
  
}
  }
  
  #Close connection
  sink()
}
