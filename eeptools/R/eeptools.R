.onAttach<-function(...){
  if (!interactive() ) return()
packageStartupMessage("Welcome to eeptools for R version 0.9.1!", appendLF=TRUE)
packageStartupMessage("Developed by Jared E. Knowles 2012-2015", appendLF=TRUE)
packageStartupMessage("for the Wisconsin Department of Public Instruction", appendLF=TRUE)
packageStartupMessage("Distributed without warranty.", appendLF=TRUE)

}