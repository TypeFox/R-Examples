### This file contains function to show all possible options in phyclust.

.show.option <- function(){
  my.cat("Options available in phyclust:", "\n")
  my.cat("boundary method: ", paste(.boundary.method, collapse = ", "), "\n")
  my.cat("code type: ", paste(.code.type, collapse = ", "), "\n")
  my.cat("edist model: ", paste(.edist.model, collapse = ", "), "\n")
  my.cat("em method: ", paste(.em.method, collapse = ", "), "\n")
  my.cat("identifier: ", paste(.identifier, collapse = ", "), "\n")
  my.cat("init method: ", paste(.init.method, collapse = ", "), "\n")
  my.cat("init procedure: ", paste(.init.procedure, collapse = ", "), "\n")
  my.cat("standard code: \n")
  my.print(as.matrix(.nucleotide))
  my.print(as.matrix(.snp))
  my.cat("substitution model: \n")
  my.print(as.matrix(.substitution.model))
  invisible()
} # End of .show.option()
