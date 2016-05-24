test.data.objects <-
function(admix.gen=NULL, loci.data=NULL,
                            parental1=NULL, parental2=NULL){
  if (dim(parental1)[1] != dim(parental2)[1])
    stop ("the number of loci is not the same for both parental populations")
  else if (dim(parental1)[1] != dim(admix.gen)[1])
    stop (paste("the number of loci is not the same for the admixed and parental populations:", dim(parental1)[1], dim(admix.gen)[1]))
  else if (dim(parental1)[1] != dim(loci.data)[1])
    stop ("loci.data does not contain the appropriate number of loci")
  return (NULL)
}

