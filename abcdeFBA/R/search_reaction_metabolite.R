SEARCH_reaction<-function(react_name,fba_object)
{	
	if(is.character(react_name)){
	print(fba_object$reaction_list[grep(react_name,fba_object$reaction_list,ignore.case=TRUE)])
	print(grep(react_name,fba_object$reaction_list,ignore.case=TRUE))
	}
	
	if(is.numeric(react_name)){print(fba_object$reaction_list[react_name])}
}

SEARCH_metabolite<-function(metabolite_name,fba_object)
{
	if(is.character(metabolite_name)){
	print(cbind(fba_object$metabolite_name[grep(metabolite_name,fba_object$metabolite_name,
	ignore.case=TRUE)],grep(metabolite_name,fba_object$metabolite_name,ignore.case=TRUE),
	fba_object$comp_name[fba_object$compartment[grep(metabolite_name,fba_object$metabolite_name,
	ignore.case=TRUE)]]))}

	if(is.numeric(metabolite_name)){print(fba_object$metabolite_name[metabolite_name])}
}
