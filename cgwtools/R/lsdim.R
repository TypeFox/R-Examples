# a toy to return dim() when it's sensible and length() elsewise
# items must be collection or list of character strings, e.g. output of ls()
lsdim<- function(items) {
	dims<-vector('list',length(items))
	names(dims)<-items
	for(thing in seq(1,length(items))) {
		if (is.null(dim(get(items[thing])))) {
	# I wanted to do dims$(thenameitself), but I fear only a paste
	# eval(parse(text=paste("dims$",names,'<-length(get(names))',sep="")))
	# would work. Simpler just to use numerical loop index
			dims[[thing]]<-length(get(items[thing]))
			} else{
				#load with dim()
				dims[[thing]]<-dim(get(items[thing]))
				}
		}
	return(dims)
	}