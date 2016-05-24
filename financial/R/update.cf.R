"update.cf" <-
function (object,flows=NULL,i=NULL,safe=NULL,rein=NULL,...)
{

	if (is.null(flows)) flows=object$cf;
	if (is.null(i)) i=object$tab[,1];
	if (is.null(safe)) safe=object$mirr[,1];
	if (is.null(rein)) rein=object$mirr[,2];

	return(cf(flows,i,safe,rein));

}

