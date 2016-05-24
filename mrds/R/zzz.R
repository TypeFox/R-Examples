.onAttach<-function(...){
  # uses packageStartupMessage which can then be
  # surpressed

  version <- utils::packageVersion("mrds")
  built <- utils::packageDescription("mrds",fields="Built")

  hello <- paste("This is mrds ",version,"\nBuilt: ",built,sep="")
  packageStartupMessage(hello)
}
