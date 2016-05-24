bestcod<-function(quantile,continuous,boxcox,cutpoint,pval.naive)
{
posminp<-which.min(pval.naive)
if (missing(cutpoint))
{
	if (class(quantile)=="NULL")
	{
					if (missing(boxcox))
					{
						if (continuous=="TRUE")
					
						{
						transf<-"original continuous variable"

										}
						else{
						transf<-"Coding definition is missing"
	
						}
					}
					else{
						if (continuous=="TRUE")
						{ 
						if(posminp==1){transf<-"original continuous variable"}
						if(posminp>1){transf<-boxcox[posminp-1]}

						}
						else{
						transf<-boxcox[posminp]
						}
					}
	}

	else{
					if (missing(boxcox))
					{
						if (continuous=="TRUE")
						{
						if(posminp<nrow(quantile)+1){transf<-quantile[posminp,]}
						if(posminp>nrow(quantile)){transf<-"original continuous variable"}
					
										}
						else{
						transf<-quantile[posminp,]
										}
					}
					else{
						if (continuous=="TRUE")
						{
						if(posminp<nrow(quantile)+1){transf<-quantile[posminp,]}
						if(posminp==nrow(quantile)+1){transf<-"original continuous variable"}
						if(posminp>nrow(quantile)+1){transf<-boxcox[posminp-nrow(quantile)+1]}


						}
						else{
						if(posminp<nrow(quantile)+1){transf<-quantile[posminp,]}
						if(posminp>nrow(quantile)){transf<-boxcox[posminp-nrow(quantile)]}

						}
				
					}

	}
}
else{
	if (class(quantile)=="NULL")
	{
					if (missing(boxcox))
					{
						if (continuous=="TRUE")
					
						{
						if(posminp<nrow(cutpoint)+1){transf<-cutpoint[posminp,]}
						if(posminp>nrow(cutpoint)){transf<-"original continuous variable"}

										}
						else{
						transf<-cutpoint[posminp,]
	
						}
					}
					else{
						if (continuous=="TRUE")
						{ 
						
						if(posminp<nrow(cutpoint)+1){transf<-cutpoint[posminp,]}
						if(posminp==nrow(cutpoint)+1){transf<-"original continuous variable"}
						if(posminp>nrow(cutpoint)+1){transf<-boxcox[posminp-nrow(cutpoint)+1]}

						}
						else{
						if(posminp<nrow(cutpoint)+1){transf<-cutpoint[posminp,]}
						if(posminp>nrow(cutpoint)){transf<-boxcox[posminp-nrow(quantile)]}

						}
					}
	}

	else{
					if (missing(boxcox))
					{
						if (continuous=="TRUE")
						{
						if(posminp<nrow(quantile)+1){transf<-quantile[posminp,]}
						if( (posminp>nrow(quantile)) & (posminp<nrow(quantile)+nrow(cutpoint+1)) ){transf<-cutpoint[posminp-nrow(quantile),]}
						if(posminp==(nrow(quantile)+nrow(cutpoint)+1)){transf<-"original continuous variable"}
					
										}
						else{
						if(posminp<nrow(quantile)+1){transf<-quantile[posminp,]}
						if( (posminp>nrow(quantile)) & (posminp<nrow(quantile)+nrow(cutpoint+1)) ){transf<-cutpoint[posminp-nrow(quantile),]}

										}
					}
					else{
						if (continuous=="TRUE")
						{
						if(posminp<nrow(quantile)+1){transf<-quantile[posminp,]}
						if( (posminp>nrow(quantile)) & (posminp<nrow(quantile)+nrow(cutpoint+1)) ){transf<-cutpoint[posminp-nrow(quantile),]}
						if(posminp==(nrow(quantile)+nrow(cutpoint)+1)){transf<-"original continuous variable"}
						if(posminp>(nrow(quantile)+nrow(cutpoint)+1)){transf<-boxcox[posminp-(nrow(quantile)+nrow(cutpoint)+1)]}



						}
						else{
						if(posminp<nrow(quantile)+1){transf<-quantile[posminp,]}
						if( (posminp>nrow(quantile)) & (posminp<nrow(quantile)+nrow(cutpoint+1)) ){transf<-cutpoint[posminp-nrow(quantile),]}
						if(posminp>(nrow(quantile)+nrow(cutpoint))){transf<-boxcox[posminp-(nrow(quantile)+nrow(cutpoint))]}

						}
				
					}

	}
}
return(transf)
}
