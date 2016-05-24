"realign.pims" <-function(dm){
#  Value:
#
#  new.indices - a vector of new indices for the old PIM values.  The old
#                PIM values are 1:length(new.indices) and the new index is
#                the corresponding value.  For example, new.indices=c(1,1,2,2)
#                renumbers the structure 1,2,3,4 such that 1,2 are now 1
#                and 3,4 are now 2.
#
#  Get all the rows in the design matrix and paste all the values
#  in each row together.
#
	allvals=apply(dm,1,paste,collapse="")
#
#  Get all the unique rows in the design matrix and paste all the values
#  in each row together.
#
	uniquevals=unique(allvals)
#
#  Find the corresponding sets of indices by matching allvals into uniquevals
#
	new.indices=match(allvals, uniquevals)
	

	return(new.indices)
}
