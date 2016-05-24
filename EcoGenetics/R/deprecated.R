# Functions deprecated in EcoGenetics version 1.2.0-2

aue.idig <- function(...) {
  stop("aue.idig is deprecated for Ecogenetics 1.2.0-2. Use eco.format")
}

aue.char2num <- function(...) {
  stop("aue.char2num is deprecated for Ecogenetics 1.2.0-2. Use eco.format")
}

eco.2columns <- function(...) {
  stop("eco.2columns is deprecated for Ecogenetics 1.2.0-2. Use eco.convert")
}

eco.append <- function(...) {
  stop("eco.append is deprecated for Ecogenetics 1.2.0-2. 
        The function has been improved with accessors. Use the 
        accessor <ecoslot.OUT>.
       See help(\"EcoGenetics accessors\"), Details and Examples")
}