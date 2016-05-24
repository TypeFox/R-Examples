.onAttach<-function(...){
  # uses packageStartupMessage which can then be
  # surpressed
  version <- utils::packageVersion("dsm")
  built <- utils::packageDescription("dsm",fields="Built")

  hello <- paste0("This is dsm ",version,"\nBuilt: ",built)
  packageStartupMessage(hello)
}
