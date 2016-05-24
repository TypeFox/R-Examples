FamiliasPrior <- function (pedigrees, generationsParameter = 1, inbreedingParameter = 1, partnerParameter = 1, maxGenerations) 
{
   if (missing(pedigrees) || length(pedigrees)<1)
      stop("The pedigrees parameter must be an object of type 'pedigree' or 'FamiliasPedigree', or a list of such.") 
   if (class(pedigrees)=="pedigree" | class(pedigrees)=="FamiliasPedigree") pedigrees <- list(pedigrees)
   if (class(pedigrees)!="list")
      stop("The pedigrees parameter must be an object of type 'pedigree' or 'FamiliasPedigree', or a list of such.") 
   for (i in pedigrees) {
       if (class(i)!="pedigree" && class(i)!="FamiliasPedigree")
          stop("The pedigrees parameter must be an object of type 'pedigree' or 'FamiliasPedigree', or a list of such.") 
   }
   npeds <- length(pedigrees)
   firstped <- pedigrees[[1]]
   persons <- firstped$id
   for (i in pedigrees[-1]) persons <- persons[persons %in% i$id]
   npers <- length(persons)
   if (npers<2) 
     stop("The function is meaningless unless there are at least two persons common to all pedigrees.")

   indextable <- matrix(0, npers, npeds)
   
   for (i in 1:npeds)
      for (j in 1:npers) {
         indextable[j,i] <- match(persons[j], pedigrees[[i]]$id)
         if (i>1) {
            if (pedigrees[[i]]$sex[indextable[j,i]] != firstped$sex[indextable[j,1]])
               stop("Persons common to all pedigrees must have the same sex in all pedigrees!")
         }
      } 

   #NewFamilias()
   .C("NewFamilias")

   #Input all persons with data
   for (j in 1:npers) {
      #AddPerson(!(firstped$sex[indextable[j,1]]=="female"))
      result <- .C("AddPerson", as.integer(!(firstped$sex[indextable[j,1]]=="female")), 
		as.integer(-1), as.integer(FALSE), index = integer(1), error = integer(1))
      if (result$error>0) 
         stop("ERROR: Problems with list of persons common to all pedigrees.")	
   }


   #Input all pedigrees
   for (i in pedigrees) {
      nPersons <- length(i$sex)
      neworder <- rep(0, nPersons)
      nExMales <- nExFemales <- 0
      for (j in 1:nPersons) {
         mm <- match(i$id[j], persons, nomatch=0)
         if (mm>0) 
            neworder[j] <- mm
         else if (i$sex[j]=="female") {
            nExFemales <- nExFemales + 1
            neworder[j] <- nExFemales
         } else { #People of unknown sex become males. 
            nExMales <- nExMales + 1
            neworder[j] <- nExMales
         }
      }
      for (j in 1:nPersons) {
         if (!(i$id[j]%in%persons)) { 
            if (i$sex[j]=="female") 
               neworder[j] <- neworder[j] + npers
            else
               neworder[j] <- neworder[j] + npers + nExFemales
         }
      }

      #index <- AddPedigree(nExFemales, nExMales) 
      result <- .C("AddPedigree", as.integer(nExFemales), 
		   as.integer(nExMales), 
		   index = integer(1), 
		   error = integer(1))
      if (result$error>0) 
         stop("ERROR: Problem adding pedigree.")			  
     
         index <- result$index + 1

      for (j in 1:nPersons) {
         if (i$findex[j]>0) {
            #AddRelation(neworder[i$findex[j]], neworder[j], index)
            result <- .C("AddRelation", as.integer(neworder[i$findex[j]]-1), 
		         as.integer(neworder[j]-1), 
		         as.integer(index-1), 
		         error = integer(1))
            if (result$error==1) 
               stop(paste("ERROR: Problem in pedigree", index))	
            else if (result$error==2)
               stop(paste("ERROR: Problem in pedigree", index,": Illegal relation based on Year-of-birth or is-Child data."))
            else if (result$error==3)
               stop(paste("ERROR: Problem in pedigree", index,": Cycle in the pedigree or duplicate parent."))
         }
         if (i$mindex[j]>0) { 
            #AddRelation(neworder[i$mindex[j]], neworder[j], index)
            result <- .C("AddRelation", as.integer(neworder[i$mindex[j]]-1), 
		         as.integer(neworder[j]-1), 
		         as.integer(index-1), 
		         error = integer(1))
            if (result$error==1) 
               stop(paste("ERROR: Problem in pedigree", index))	
            else if (result$error==2)
               stop(paste("ERROR: Problem in pedigree", index,": Illegal relation based on Year-of-birth or is-Child data."))
            else if (result$error==3)
               stop(paste("ERROR: Problem in pedigree", index,": Cycle in the pedigree or duplicate parent."))
         }
      }
   }

   if (missing(maxGenerations)) 
      maxGenerations <- -1 
   if (generationsParameter < 0 | inbreedingParameter < 0 | partnerParameter < 0) 
      stop("ERROR: The parameters cannot be negative.")
   result <- .C("GetProbabilities", 
       	  as.double(generationsParameter), 
          as.integer(maxGenerations), 
    	  as.double(inbreedingParameter), 
    	  as.double(partnerParameter), 
    	  as.integer(FALSE),
	  as.double(0), 
	  redundant = integer(npeds), 
	  probabilities = double(npeds), 
	  likelihoods = double(0), 
	  error = integer(1))
   if (result$error==1) 
      stop("ERROR: Wrong input.")
   if (result$error==2)
      stop("ERROR: All pedigrees have probability zero.")
   if (any(result$redundant))
      stop("Error: some pedigrees are duplicate, remove duplicates.")

   #NewFamilias()
   .C("NewFamilias")

   result$probabilities
}
