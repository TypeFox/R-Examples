Sybil_2_FBA_obj<-function(Sybil_S4_object)
{
	FBA_obj=list()
	FBA_obj$mat<-Sybil_S4_object@S
	FBA_obj$obj<-Sybil_S4_object@obj_coef
	FBA_obj$rhs<-Sybil_S4_object@rhs
	FBA_obj$types<-rep("C",dim(Sybil_S4_object@S)[2])
	FBA_obj$dir<-rep("==",dim(Sybil_S4_object@S)[1])
	FBA_obj$bounds<-list(lower=list(ind=1:length(Sybil_S4_object@lowbnd),val=Sybil_S4_object@lowbnd),upper=list(ind=1:length(Sybil_S4_object@uppbnd),val=Sybil_S4_object@uppbnd))
	FBA_obj$max<-TRUE
	FBA_obj$reaction_list<-Sybil_S4_object@react_name
	FBA_obj$metabolite_name<-Sybil_S4_object@met_name
	FBA_obj$sub_system<-Sybil_S4_object@subSys
	FBA_obj$compartment<-Sybil_S4_object@met_comp
	FBA_obj$comp_name<-Sybil_S4_object@mod_compart
	FBA_obj$all_genes<-Sybil_S4_object@allGenes
	FBA_obj$gpr<-Sybil_S4_object@gpr
	FBA_obj$rxnGeneMat<-Sybil_S4_object@rxnGeneMat

	return(FBA_obj)
}


