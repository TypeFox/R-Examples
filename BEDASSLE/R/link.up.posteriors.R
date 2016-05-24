link.up.posteriors <-
function(MCMC.output1,MCMC.output2,linked.up.output.file.name){ 
# recover()
	load(MCMC.output1)
	parameters <- objects()	
		for(i in 1:length(parameters)){
			assign(sprintf("tmp.%s",parameters[i]),get(parameters[i]))	
		}	
	load(MCMC.output2)	
		for(i in 1:length(parameters)){
			if(length(get(parameters[i])) > 1 && grepl(pattern="last.params",parameters[i])==FALSE){		
				if(class(get(parameters[i])) == "numeric"){
					if(grepl("moves",parameters[i]) || grepl("accept",parameters[i])){
						assign(parameters[i],
							c(get(sprintf("tmp.%s",parameters[i])),
								get(parameters[i]) +
								get(sprintf("tmp.%s",parameters[i]))[length(get(sprintf("tmp.%s",parameters[i])))]))
					} else{
						assign(parameters[i],c(get(sprintf("tmp.%s",parameters[i])),get(parameters[i])))				
					}
				}
				if(class(get(parameters[i])) == "matrix"){
					assign(parameters[i],cbind(get(sprintf("tmp.%s",parameters[i])),get(parameters[i])))
				}			
			}
		}
	rm(list=unique(c(objects(pattern="tmp."),objects(pattern="MCMC.output"))))
	save(list=setdiff(ls(all.names=TRUE),"linked.up.output.file.name"),file=paste(linked.up.output.file.name,".Robj",sep=""))
}
