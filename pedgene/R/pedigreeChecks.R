## Purpose: pedigree checks on data.frame with pedigree data
## Author: Dan Schaid 5/2013

pedigreeChecks <- function(pedigree, male.code=1, female.code=2)
{
  valid = FALSE
  
  pedSize <- nrow(pedigree)
  if(length(unique(pedigree$person)) != pedSize)
    {
      stop("Person ID's not unique")
    }
  
  if(any(pedigree$father == 0 & pedigree$mothher !=0)){
    stop("Father = 0, but mother != 0")
  }
  
  if(any(pedigree$father != 0 & pedigree$mothher ==0)){
    stop("Father != 0, but mother = 0")
  }
 
  parent <- ifelse(pedigree$father != 0 & pedigree$mother != 0, TRUE, FALSE)
  father.id <- unique(pedigree$father[parent])
  mother.id <- unique(pedigree$mother[parent])

  sex <- ifelse(pedigree$sex == male.code, 1, NA)
  sex <- ifelse(pedigree$sex == female.code, 2, sex)

  if(any(is.na(sex))){
    stop("Unknown sex codes")
  }
  
  father.index <- match(father.id, pedigree$person)
  if(any(is.na(father.index))){
    stop("Person info missing for a father")
  }
  if(any(sex[father.index]!=1)){
    stop("Sex of father not correct")
  }

  mother.index <- match(mother.id, pedigree$person)
  if(any(is.na(mother.index))){
    stop("Person info missing for a mother")
  }
  if(any(sex[mother.index]!=2)){
    stop("Sex of mother not correct")
  }

  valid = TRUE
  return(valid)
}
