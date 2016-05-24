CHANGE_RXN_BOUNDS<-function(reaction_number=NULL,fba_object,lb=0,ub=0)
{
	if(length(reaction_number)>0)
	{
	fba_object$bounds$lower$val[reaction_number]=lb
	fba_object$bounds$upper$val[reaction_number]=ub
	}
	return(fba_object)
}

CHANGE_OBJ_FUNCTION<-function(obj_reaction=NULL,fba_object,new_wt=1,old_wt=0)
{
	if(length(obj_reaction)>0)
	{
		if(old_wt>0)
		{fba_object$obj[which(fba_object$obj==1)]=old_wt}
		if(old_wt==0)
		{fba_object$obj[which(fba_object$obj==1)]=old_wt}

		fba_object$obj[obj_reaction]=new_wt
	}
	return(fba_object)
}
