FamiliasPedigree <- function(id, dadid, momid, sex) {
   if (missing(id) | missing(dadid) | missing(momid) | missing(sex))
      stop("All four parameters must be specified.")
   npers <- length(id)
   if (length(unique(c(id, NA))) < npers + 1)
      stop("The id parameter should be a vector of unique person identifiers.")
   if (length(dadid)!=npers | length(momid)!=npers | length(sex)!=npers)
      stop("All four parameters must be vectors of the same lengths.")
   if (any(!(dadid %in% c(id, NA))) | any(!(momid %in% c(id, NA))))
      stop("The parameters dadid and momid must consist of identifiers present in the id parameter, or NA.")
   if (any(sex!="male" & sex!="female"))
      stop("The sex parameter can have values 'male' and 'female' only.")
   male <- sex=="male"
   dad <- match(dadid, id, nomatch = 0)
   mom <- match(momid, id, nomatch = 0)
   if (any(!(male[dad[dad>0]])))
      stop("All persons indicated as fathers in the dadid parameter must be male.")
   if (any(male[mom[mom>0]]))
      stop("All persons indicated as mothers in the momid parameter must be female.")
   #Check acyclicity: 
   #The following function takes as input a person x and a vector v indicating a path ending at x. 
   #It returns true if the path can be extended so that it intersects itself. 
   hasIntersection <- function(x, v) {
          (mom[x]>0 && (mom[x] %in% v || hasIntersection(mom[x], c(v, mom[x])))) || 
          (dad[x]>0 && (dad[x] %in% v || hasIntersection(dad[x], c(v, dad[x])))) }
   for (i in 1:npers)
      if (hasIntersection(i,i))
         stop("The dadid and momid parameters result in a cyclic pedigree.")
   result <- list(id=id, findex=dad, mindex=mom, sex=sex)
   class(result) <- "FamiliasPedigree"
   result
}