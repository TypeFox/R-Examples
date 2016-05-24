reportEquivalentVariables <-
function (object,pvalue=0.05,data,variableList,Outcome="Class", type = c("LOGIT", "LM","COX"),eqFrac=0.9,description=".") 
{
    type <- match.arg(type);
	cthr = abs(qnorm(pvalue));
  
	varsList <- as.list(attr(terms(object),"variables"))
	vnames <- as.vector(variableList[,1]);
    
	orgsize <- length(varsList);
    

	outCome = paste(varsList[2]," ~ ");
	startlist = 3;
	auxmodel <- object;
	indexzVAL = startlist - 1;
	
	outvarnames <- vector();
	outrelated <- vector();
	outzVAL <- vector();

	if (type != "LM")
	{
		orgzVAL <- as.vector(getVar.Bin(object,data,Outcome,type)$testData.z.IDIs);
	}
	else
	{
		orgzVAL <- -qnorm(as.vector(getVar.Res(object,data,Outcome,type)$FP.value),0,1);
	}

	print (orgzVAL,digits=3);
	orgzVAL <- eqFrac * orgzVAL;
	print (orgzVAL,digits=3);

	 for ( i in startlist:length(varsList))
	{
		outvarnames <- append(outvarnames,as.character(varsList[i]));
		cat(as.character(varsList[i]),"\n");	
		namelist ="{"; 
		zVALlist ="{";
		for (j in 1:length(vnames))
		{
			frm1 = outCome;
			for ( n in startlist:length(varsList))
			{
				if (n == i)
				{
					frm1 <- paste(frm1,paste(" + ",vnames[j]));
				}
				else
				{
					frm1 <- paste(frm1,paste(" + ",varsList[n]));
				}
			}
			ftmp <- formula(frm1);
			auxmodel <- modelFitting(ftmp,data,type)
      
			if (orgsize == length(as.list(attr(terms(auxmodel),"variables"))))
			{
				if (type != "LM")
				{
					zVAL <- as.vector(getVar.Bin(auxmodel,data,Outcome,type)$testData.z.IDIs);
				}
				else
				{
					zVAL <- -qnorm(as.vector(getVar.Res(auxmodel,data,Outcome,type)$FP.value),0,1);
				}
#		        cat(frm1," : ");
#		        print(zVAL,digits = 2 );
				if ((zVAL[i-indexzVAL] > cthr) && ( zVAL[i-indexzVAL] > orgzVAL[i-indexzVAL] ))
				{
				  namelist <- paste(namelist,vnames[j]);
				  if (description != ".") namelist <- paste(namelist,"[",variableList[j,description],"]");
				  zVALlist <- paste(zVALlist,sprintf("%5.3f",zVAL[i-indexzVAL]));
				}        
			}
		}
		namelist <- paste(namelist,"}");
		zVALlist <- paste(zVALlist,"}");
		cat(namelist,"\n\n");
		outrelated <- append(outrelated,namelist);
		outzVAL <- append(outzVAL,zVALlist);
	}
	result <- cbind(outvarnames);
	result <- cbind(result,outrelated);
	result <- cbind(result,outzVAL);

    return (result)
}
