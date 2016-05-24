power.plot <-
function(powers, #named list
		num.noise = seq(from=0.1,to=3,by=0.1) #or a matrix of MI values
		,mains = c("Linear","Quadratic","Cubic",expression("Sine: period 4" * pi),expression("Sine: period 16" * pi),"X^(1/4)","Circle","Step function","Torus")
		,col = c("black","red","blue","green","cyan","brown","pink")
		,labels = TRUE #what to label the tickmarks
		,which = 1:nrow(powers[[1]])
		,show.legend = "bottomright"
){
	pchs = (1:length(powers)) %% 25  #pch > 25 are not implemented 
	for(i in which){
		if(length(dim(num.noise)) > 0){
			if(is.null(rownames(num.noise))){
				plot(0,type="n",ylim=c(0,1),xlim=range(num.noise[,i]),main=mains[i],xlab="Amount of Gaussian noise",ylab="power",xaxt="n") 
				axis(1,at=num.noise[,i],labels=labels,las=2)
			}else{
				plot(0,type="n",ylim=c(0,1),xlim=range(as.numeric(rownames(num.noise))),main=mains[i],xlab="Mutual Information",ylab="power",xaxt="n")#,log="x")
				axis(1,at=as.numeric(rownames(num.noise)),labels=labels,las=2)
			}
		}else{
			plot(0,type="n",ylim=c(0,1),xlim=c(0,max(num.noise)),main=mains[i],xlab="noise",ylab="power")
		}
		sapply(1:length(powers),function(ind){
					power = powers[[ind]]
					if(length(dim(num.noise)) > 0){
						if(is.null(rownames(num.noise))){
							points(num.noise[,i],power[i,1:length(num.noise[,i])],pch=pchs[ind],col=col[ind],type="b")
						}else{
							points(as.numeric(rownames(num.noise)),power[i,1:length(num.noise[,i])],pch=pchs[ind],col=col[ind],type="b")
						}
					}else{
						points(num.noise,power[i,1:length(num.noise)],pch=pchs[ind],col=col[ind],type="b")
					}
					
				})
		if(!is.null(show.legend)){
			legend(show.legend,names(powers),col=col,pch=pchs,bg="white")
		}
	}
}
