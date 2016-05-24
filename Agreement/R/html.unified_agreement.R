html.unified_agreement <-
function(result,dir=getwd(),file="report",CSS="R2HTML",...){
	k <- result$Inputs[1];
	m <- result$Inputs[2];
	n <- result$Inputs[3];
	str <- paste("<center>","For k=",k,", n=", n, ", and m=", m,"</center>", sep="");
	target  <-  HTMLInitFile(outdir=dir,filename=file, Title="Unified Agreement Statistics",...);
  HTMLCSS(file=target, CSSfile=CSS);
	HTML("<center>Unified Agreement Statistics</center>",file=target);
	HTML(summary(result),file=target, digits=as.numeric(result$Inputs[6]));
	HTML(str, file=target);
	HTML("<center>*: The relative bias squared (RBS) must be less than 1 or 8 for CP_a of 0.9 or 0.8, respectively, in order for the approximated TDI and CP values to be valid.</center>", file=target);
	HTMLEndFile();
}

