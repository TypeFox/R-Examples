lstype <-
function(type='closure'){
#simple command to get only one type of object in current environment
# Note: if you foolishly create variables named 'c' ,'q' ,'t' or the like,
# this will fail because typeof finds the builtin function first
	inlist<-ls(.GlobalEnv)
	if (type=='function') type <-'closure'
	typelist<-sapply(sapply(inlist,get),typeof)
	return(names(typelist[typelist==type]))
}
