
ecipex<-function(formulas, isoinfo=nistiso, limit=1e-12, id=FALSE, sortby="abundance"){

sortby<-match.arg(sortby, c("abundance", "mass"))

# determine which elements are present in each molecule and the corresponding copy number
elementCount <- lapply(formulas, CHNOSZ::count.elements)

# these lists will be used as indices later
elementList <- lapply(elementCount, names)
copyList <- lapply(elementCount, format, scientific=FALSE, trim=TRUE)

# treat them as vectors for now
elements <- unlist(lapply(elementCount, names))
copies <- as.numeric(unlist(elementCount))

uniqueElements <- unique(unlist(elements))

# define a list with an entry for each unique element, stating the unique set of copy numbers that will be required for the formulas provided
compositionList<-lapply(uniqueElements, function(x){
	sort(unique(copies[elements==x]))
})
names(compositionList)<-uniqueElements

# determine the maximum copy number for each element in all of the formulas provided. Needed to determine dimensions of isofft array.
maxElements<-sapply(compositionList, max)

# define a list with an entry for each unique element containing the isotopic expansions for all copy numbers required by the formulas
pureIsos<-lapply(seq(uniqueElements), function(x){
	
	abundance<-isoinfo$abundance[isoinfo$element==uniqueElements[x]]
	mass<-isoinfo$mass[isoinfo$element==uniqueElements[x]]
	nucleons<-isoinfo$nucleons[isoinfo$element==uniqueElements[x]]
	
	#in case abundances of zero were included
	mass<-mass[abundance>0]
	abundance<-abundance[abundance>0]

	# dimensions of isofft array will be isoLength^isoDimension
	isoLength<-nextn(maxElements[x]+1) # may increase length of isofft but makes it more "composite" so that fft can more efficiently be taken
	isoDimension<-length(abundance)-1

	# for elements with only one isotopic variant
	if(isoDimension==0){ 
		abundanceOut<-lapply(seq(length(compositionList[[x]])), function(y){1})
		massOut<- lapply(seq(length(compositionList[[x]])), function(y){compositionList[[x]][y]*mass})
		
		names(abundanceOut)<-compositionList[[x]]
		names(massOut)<-compositionList[[x]]
		
		if(!id){
			return(list("mass"=massOut, "abundance"=abundanceOut))	
		}
		if(id){
			speciesOut<-lapply(seq(length(compositionList[[x]])), function(y){s<-matrix(compositionList[[x]][y]); colnames(s)<-paste(nucleons, uniqueElements[x], sep=""); return(s)})
			names(speciesOut)<-compositionList[[x]]
			return(list("mass"=massOut, "abundance"=abundanceOut, "isoSpecies"= speciesOut))
			}
	}
	
	# define array containing abundances for a single atom of the current element. Its dimensions must be able to accomodate the appropriate maxElements-fold convolution without confounding overlap
	abundanceIn<-array(rep.int(0, isoLength^isoDimension), dim=rep(isoLength, isoDimension))
	abundanceIn[1]<-abundance[1]
	abundanceIn[1+isoLength^seq(0,isoDimension-1)]<-abundance[seq(2, isoDimension+1)]

	# The output of the fft is stored in a hypercube with isoLength^isoDimension entries. We require only the entries in a series of the simplexes it contains. We therefore define an array of the same dimensions as abundanceIn (isoLength^isoDimension) which is enumerated in such a way that it is easy to extract the simplexes and save memory
	simplexIndex<-seq(0,isoLength-1, by=1)
	i=1
	while(i < isoDimension){
		simplexIndex<-outer(simplexIndex, seq(0, isoLength-1, by=1), FUN="+")	
		i=i+1
	}
		
	# perform the actual multidimensional fft convolution to obtain the desired abundances. Also remove abundance-elements lower than than the specified limit to save memory.
	fftIn<-fft(abundanceIn)
	convolOut<-lapply(compositionList[[x]], function(n){	
		arrayOut<-Re(fft(fftIn^n, inverse=TRUE))/(isoLength)^isoDimension
		species<-which((arrayOut>=limit)&(simplexIndex<=n), arr.ind=TRUE) # these indices enable us to determine which isotopic species we are dealing with
		return(list("abundance"=arrayOut[species], "species"=species))
	})

	# obtain the corresponding masses and remove abundance-elements lower than than the specified limit to save memory.
	massDiffs<-outer(mass[-1]-mass[1], seq(0, isoLength-1))
	output<-massDiffs[1,]
	if(length(mass)>2){
			for(i in 2:nrow(massDiffs)){
				output<-outer(output, massDiffs[i,], FUN="+")
			}
		}
	massOut <- lapply(seq_along(compositionList[[x]]), function(n){
		mass[1]*compositionList[[x]][n]+output[convolOut[[n]]$species]
		})
	
	# sort entries in order of descending abundance 
	descendingOrder<-lapply(seq_along(compositionList[[x]]), function(y){order(convolOut[[y]]$abundance, decreasing=TRUE)})
	abundanceOut<-lapply(seq_along(compositionList[[x]]), function(y){convolOut[[y]]$abundance[descendingOrder[[y]]]})
	massOut<-lapply(seq_along(compositionList[[x]]), function(y){massOut[[y]][descendingOrder[[y]]]})	

	names(abundanceOut) <- format(compositionList[[x]], scientific=FALSE, trim=TRUE)
	names(massOut) <- format(compositionList[[x]], scientific=FALSE, trim=TRUE)
		
	if(!id){
		return(list("mass"=massOut, "abundance"=abundanceOut))
	}
	if(id){
		speciesOut<-lapply(seq_along(compositionList[[x]]), function(y){ # determine the isotopic species of the values listed in massOut and abundanceOut
			specMat<-cbind(compositionList[[x]][y]-apply(as.matrix(convolOut[[y]]$species[descendingOrder[[y]],])-1, 1, sum), convolOut[[y]]$species[descendingOrder[[y]],]-1); 
			colnames(specMat)<-paste(nucleons, uniqueElements[x], sep="")
			return(specMat)
			})
		names(speciesOut) <- format(compositionList[[x]], scientific=FALSE, trim=TRUE)
	
		return(list("mass"=massOut, "abundance"=abundanceOut, "isoSpecies"= speciesOut))
	}
})
names(pureIsos)<-uniqueElements


# list that will contain the full isotopic expansion for all formulas
fullIsos<-lapply(seq(formulas), function(x){
	
	# extract highest abundance for each atomic species in current formula
	maxProbs<-sapply(seq(elementList[[x]]), function(i){pureIsos[[elementList[[x]][i]]]$abundance[[copyList[[x]][i]]][1]})	
	
	# for each atomic species we immediately discard abundances which are less than the limit after multiplication by highest abundances in all other species. Similarly for masses, and isotopic species
	ma<-lapply(seq(elementList[[x]]), function(i){
		z<-pureIsos[[elementList[[x]][i]]]$mass[[copyList[[x]][i]]][pureIsos[[elementList[[x]][i]]]$abundance[[copyList[[x]][i]]]*prod(maxProbs[-i])>limit]
		if(length(z)==0){return(pureIsos[[elementList[[x]][i]]]$mass[[copyList[[x]][i]]][which.max(pureIsos[[elementList[[x]][i]]]$abundance[[copyList[[x]][i]]])])}
		return(z)
		})	
	ab<-lapply(seq(elementList[[x]]), function(i){
		z<-pureIsos[[elementList[[x]][i]]]$abundance[[copyList[[x]][i]]][pureIsos[[elementList[[x]][i]]]$abundance[[copyList[[x]][i]]]*prod(maxProbs[-i])>limit]
		if(length(z)==0){return(max(pureIsos[[elementList[[x]][i]]]$abundance[[copyList[[x]][i]]]))}
		return(z)
		})
	if(id){
	is<-lapply(seq(elementList[[x]]), function(i){
		if(ncol(pureIsos[[elementList[[x]][i]]]$isoSpecies[[copyList[[x]][i]]])==1){return(pureIsos[[elementList[[x]][i]]]$isoSpecies[[copyList[[x]][i]]])} # if there is only 1 isotope
		z<-as.matrix(pureIsos[[elementList[[x]][i]]]$isoSpecies[[copyList[[x]][i]]][pureIsos[[elementList[[x]][i]]]$abundance[[copyList[[x]][i]]]*prod(maxProbs[-i])>limit,])
		if(ncol(z)<=1){return(t(as.matrix(pureIsos[[elementList[[x]][i]]]$isoSpecies[[copyList[[x]][i]]][which.max(pureIsos[[elementList[[x]][i]]]$abundance[[copyList[[x]][i]]]),])))}
		return(z)
		})
	}
	
	# these vectors will contain full isotopic expansions, but initially they contain the first elemental expansion
	moleculeMass<-ma[[1]]
	moleculeAbundance<-ab[[1]]
	if(id){moleculeSpecies<-is[[1]]}
	
	# fold them an appropriate number of times
	# each time we fold, a number of abundances will fall below limit and these are excluded. Furthermore, we can exclude abundances that are less than the limit after they have been multiplied by the product of the highest abundances in the remaining atomic species
	remainingMaxAbundance<-c(cumprod(maxProbs[c(length(maxProbs):1)])[c(length(maxProbs):1)], 1) # compare with e.g. prod(maxProbs[c(3:5)])
	i<-2
	while(i <= length(elementList[[x]])){
		moleculeMass<-outer(moleculeMass, pureIsos[[elementList[[x]][i]]]$mass[[copyList[[x]][i]]], FUN="+")
		moleculeAbundance<-outer(moleculeAbundance, pureIsos[[elementList[[x]][i]]]$abundance[[copyList[[x]][i]]])
		keepers<-which(moleculeAbundance>limit/remainingMaxAbundance[i+1], arr.ind=TRUE)
				
		moleculeMass<-moleculeMass[keepers]
		moleculeAbundance<-moleculeAbundance[keepers]		

		if(id){

			moleculeSpecies<-cbind(matrix(moleculeSpecies[keepers[,1],], ncol=ncol(moleculeSpecies), dimnames=list(NULL, colnames(moleculeSpecies))), matrix(pureIsos[[elementList[[x]][i]]]$isoSpecies[[copyList[[x]][i]]][keepers[,2],], ncol=ncol(pureIsos[[elementList[[x]][i]]]$isoSpecies[[copyList[[x]][i]]]), dimnames=list(NULL, colnames(pureIsos[[elementList[[x]][i]]]$isoSpecies[[copyList[[x]][i]]]))))
				
			}		
		i<-i+1
	}
	
	if(sortby=="abundance"){
		dex<-order(moleculeAbundance, decreasing=TRUE)
	}
	if(sortby=="mass"){
		dex<-order(moleculeMass, decreasing=FALSE)
	}
	
	if(id){
		return(cbind(data.frame("mass"= moleculeMass[dex], "abundance"=moleculeAbundance[dex]), matrix(moleculeSpecies[dex,], ncol=ncol(moleculeSpecies), dimnames=list(NULL, colnames(moleculeSpecies)))))
		
	}
	if(!id){
		return(data.frame("mass"= moleculeMass[dex], "abundance"=moleculeAbundance[dex]))		
	}
	
})
names(fullIsos)<-formulas

return(fullIsos)	
}
