# is.validsyntax
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

is.validsyntax <- function(reaction){
  if (grepl("([[:digit:]][[:blank:]][[:digit:]][[:blank:]])+",reaction) ||
      grepl("[[:blank:]](<?)-(>?)[[:blank:]]",reaction) ||
      grepl("[[:alnum:]]\\+[[:alnum:]]", reaction) ||
      (!grepl("[[:blank:]]",reaction)) ||
      grepl("[[:blank:]]\\-[[:alnum:]]",reaction)){
    return(FALSE)
  } else{
    return(TRUE)
  }
}



