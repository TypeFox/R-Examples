famMix <- function (x, R, id.U, id.V = NULL,partialmarker=NULL,theta=0,mutationRateFemale = 0, mutationRateMale = 0, mutationModelFemale = "stable", mutationModelMale = "stable", mutationRangeFemale = 0.1, mutationRangeMale = 0.1, silentFrequency = 0,check=TRUE) {
   if(silentFrequency==0) silentFrequency=NULL
   stopifnot(x$nMark==0)
   if (is.null(partialmarker)) m = stop("Partialmarker must be supplied") else m = partialmarker
   allgenos = allGenotypes(attr(m, 'nalleles'))
   allgenos_in_R = which(allgenos[, 1] %in% R & allgenos[, 2] %in% R)

   all_typed = which(m[,1] != 0 | m[,2] != 0)
   contrib_typed = setdiff(all_typed, id.V)
   contrib_untyped = id.U
   K = unique(c(m[contrib_typed, ]))
   R_not_masked = setdiff(R, K)

   #if(check) stopp=.checkInput(x,R,id.U,id.V,partialmarker=m,all_typed,K,R_not_masked)
   if(check) stopp=checkInput(x,R,id.U,id.V,all_typed,K,R_not_masked)
   set=generate(R,K,length(contrib_untyped))
   
   if(is.null(set)) x=addMarker(x,m)
   else {
      dd=dim(set)[2]/2
      for (i in 1:dd) x=addMarker(x,m)
      x=modifyMarkerMatrix(x, contrib_untyped, set)
  	}
    m2=x$markerdata[[1]]
    x$available=which(m2[,1] != 0 | m2[,2] != 0)
    x=relabel(x,1:x$nInd,x$orig.ids)
    .preFamilias2(x,new=TRUE,mutationRateFemale = mutationRateFemale, mutationRateMale = mutationRateMale , mutationModelFemale = mutationModelFemale, mutationModelMale = mutationModelMale, mutationRangeFemale =  mutationRangeFemale, mutationRangeMale = mutationRangeMale, silentFrequency = silentFrequency )
    PE=GetProbabilities(kinship=theta)$likelihoodsPerSystem
	list(x=x,likelihood = sum(PE), allLikelihoods = PE[PE>0])
}


.preFamilias2=function (x, new = TRUE, mutationRateFemale = 0, mutationRateMale = 0, 
    mutationModelFemale = "stable", mutationModelMale = "stable", 
    mutationRangeFemale = 0.1, mutationRangeMale = 0.1, silentFrequency = NULL) 
{
    if(is.null(silentFrequency)) hack=0 else hack=silentFrequency
    if (new) {
        NewFamilias()
        for (i in x$orig.ids) AddPerson(x$pedigree[i, 4] == 1)
        for (i in 1:x$nMark) {
		    afreq=attr(x$markerdata[i][[1]], which = "afreq")
			afreq=(1-hack)*afreq
		    AddAlleleSystem(afreq, mutationRateFemale = mutationRateFemale, 
            mutationRateMale = mutationRateMale, mutationModelFemale = mutationModelFemale, 
            mutationModelMale = mutationModelMale, mutationRangeFemale = mutationRangeFemale, 
            mutationRangeMale = mutationRangeMale, silentFrequency = silentFrequency)
		}	
        for (i in 1:x$nMark) for (j in x$available) AddDNAObservation(j, 
            i, x$markerdata[[i]][j, 1], x$markerdata[[i]][j, 
                2])
    }
    ped = AddPedigree()
    if (max(x$ped[, 2:3]) > 0) {
        for (i in x$orig.ids) {
            Father = x$pedigree[i, 2]
            Mother = x$pedigree[i, 3]
            Child = x$pedigree[i, 1]
            if (Father > 0) 
                AddRelation(Father, Child, ped)
            if (Mother > 0) 
                AddRelation(Mother, Child, ped)
        }
    }
    ped
}
