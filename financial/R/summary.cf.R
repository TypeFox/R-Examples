"summary.cf" <-
function (object,flows=2:length(object$cf),...) 
{
	for (k in flows) {

	cat("\n#",k,"\n");

	print(cf(object$cf[1:k],object$tab[,1]));

	}

}

