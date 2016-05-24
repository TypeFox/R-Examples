"makeSample" <-
function(file, nbPatient,disXi,disP)
{
#disXi is a 3 elements list : rangen <=>  a random positive variable generator
#       ; nbparam <=> number of parameter of this distribution and param <=> a list of parameters
#disP is a 3 elements list : disfun <=> a distribution function
#       ; nbparam <=> number of parameter of this distribution and param <=> a list of parameters

       	#simulation of variables Xi and Zi
       	if(disXi$nbparam==1) x<-disXi$rangen(nbPatient,disXi$param[[1]]);
       	if(disXi$nbparam==2) x<-disXi$rangen(nbPatient,disXi$param[[1]],disXi$param[[2]]);
       	if(disXi$nbparam==3) x<-disXi$rangen(nbPatient,disXi$param[[1]],disXi$param[[2]],disXi$param[[3]]);
	   	xCumsum<-cumsum(x);
	
	   	p_x<-0;
	
	   	if(disP$nbparam==1) p_x<-sapply(x,disP$param[[1]],FUN=disP$disfun);
       	if(disP$nbparam==2) p_x<-sapply(x,disP$param[[1]],disP$param[[2]],FUN=disP$disfun);
       	if(disP$nbparam==3) p_x<-sapply(x,disP$param[[1]],disP$param[[2]],disP$param[[3]],FUN=disP$disfun);
       	if(disP$nbparam==4) p_x<-sapply(x,disP$param[[1]],disP$param[[2]],disP$param[[3]],disP$param[[4]],FUN=disP$disfun);

       	z<-rbinom(nbPatient,1,p_x);
	 	write(c("no","sejour","evt"), file=file,ncolumns=3);
        write(t(matrix(c(1:nbPatient,x,z),ncol=3)), file=file,ncolumns=3,append=TRUE);
        
		#    Normally 'writeLines' is used with a text connection, and the
		#    default separator is converted to the normal separator for that
		#    platform (LF on Unix/Linux, CRLF on Windows, CR on Classic MacOS).
		#    For more control, open a binary connection and specify the
		#    precise value you want written to the file in 'sep'.  For even
		#    more control, use 'writeChar' on a binary connection.

        return(NULL)
}

