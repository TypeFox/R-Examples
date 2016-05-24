"makeSample2" <-
function(file, nbPatient,disXi,disYi)
{
#disXi is a 3 elements list : rangen <=>  a random positive variable generator
#       ; nbparam <=> number of parameter of this distribution and param <=> a list of parameters
#disYi is a 3 elements list : rangen <=> a random positive variable generator
#       ; nbparam <=> number of parameter of this distribution and param <=> a list of parameters

       	#simulation of variables Xi, Yi and Zi
       	if(disXi$nbparam==1) x<-disXi$rangen(nbPatient,disXi$param[[1]]);
       	if(disXi$nbparam==2) x<-disXi$rangen(nbPatient,disXi$param[[1]],disXi$param[[2]]);
       	if(disXi$nbparam==3) x<-disXi$rangen(nbPatient,disXi$param[[1]],disXi$param[[2]],disXi$param[[3]]);
       
       	if(disYi$nbparam==1) y<-disYi$rangen(nbPatient,disYi$param[[1]]);
       	if(disYi$nbparam==2) y<-disYi$rangen(nbPatient,disYi$param[[1]],disYi$param[[2]]);
       	if(disYi$nbparam==3) y<-disYi$rangen(nbPatient,disYi$param[[1]],disYi$param[[2]],disYi$param[[3]]);
       
	   	xCumsum<-cumsum(x);
	   	yCumsum<-cumsum(y);
			
	 
	   	z<-array(0,nbPatient);
		
		for (i in 1:nbPatient)
		{
			if(i>1)
			{
				if(y[i]<x[i])
				{
					z[i]<-1;
				}
				else
				{
					z[i]<-0;
				}
			}
		}
        
        
        
	 	write(c("no","sejour","expo","evt"), file=file,ncolumns=4);
        write(t(matrix(c(1:nbPatient,x,y,z),ncol=4)), file=file,ncolumns=4,append=TRUE);
	        
		#    Normally 'writeLines' is used with a text connection, and the
		#    default separator is converted to the normal separator for that
		#    platform (LF on Unix/Linux, CRLF on Windows, CR on Classic MacOS).
		#    For more control, open a binary connection and specify the
		#    precise value you want written to the file in 'sep'.  For even
		#    more control, use 'writeChar' on a binary connection.

            return(NULL)
}

