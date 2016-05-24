CombineCaseControlC <-
function(cases, controls) {
	## Argument: Cases, Controls point processes
	## Return: combX, combZ, combL
	## Combines the cases and controls point processes
	## Produce a process realization combL containing distinct read positions
	## and X_i be total number of reads at the ith distinct read position
	## and Z_i be number of reads at the ith distinct read position from case
	## and L_i be the label (genome location) at the ith distinct read position
	cases = sort(cases)
	controls = sort(controls)
	combCC = .Call("CombineSortedVectorC", as.numeric(cases), as.numeric(controls), PACKAGE="seqCBS")
	combL = .Call("FindUniqueInSortedArrayC", as.numeric(combCC), PACKAGE="seqCBS")
	combZX = .Call("CombineToUniqueValueC", as.numeric(cases), as.numeric(controls), as.numeric(combL), PACKAGE="seqCBS")
	combZ = as.numeric(combZX[,1])
	combX = as.numeric(combZX[,2])
	return(list(combX=combX, combZ=combZ, combL=combL))
}

