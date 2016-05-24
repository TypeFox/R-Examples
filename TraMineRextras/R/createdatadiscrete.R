# TODO: Add comment
# 
# Author: nmuller
###############################################################################


createdatadiscrete <- function(ids, data, vars, agemin, agemax,supvar=NULL) {
	
	nvars <- length(vars)
	datadiscrete <- list("data"=list(), "vars"=list(), "ids"=ids)

		for(i in 1:nvars) {
			message("Creating data set for variable ", vars[i])
			vardep <- vars[i]
			varindep <- vars[-i]
  #    print(vardep)
  #    print(varindep)
			datadiscrete[["data"]][[vars[i]]] <- as.data.frame(do.call("cdd", list(data, ids, varindep, vardep, agemin, agemax,supvar)))
		}

 return(datadiscrete)
}

cdd <- function(data, ids, varindep, vardep, agemin, agemax, supvar) {
	res <- matrix(nrow=0, ncol=length(varindep)+3+length(supvar))
 # print(data$IDPERS)
  idv <- ids

  data <- as.matrix(data)

	if(is.null(supvar)) { dimnames(res) <- list(c(), c("IDPERS", "AGE", vardep, varindep)) }
	else { dimnames(res) <- list(c(), c("IDPERS", "AGE", vardep, varindep, supvar)) }

	for(i in 1:length(idv)) {
		
		if(i%%1000==0) { message(round((i/length(idv))*100,digits=2), "% completed") }
		
		tmp <- data[data[,"IDPERS"]==idv[i],]
		
		
		if(!is.na(agemin[agemin$IDPERS==idv[i],"AGE"]) && !is.na(agemax[agemax$IDPERS==idv[i],"AGE"])) {
	
			agemini <- agemin[agemin$IDPERS==idv[i],"AGE"]
			agemaxi <- agemax[agemax$IDPERS==idv[i],"AGE"]
			
			if(tmp[paste(vardep,"ST",sep="")]==1) { 
				agemaxi <- as.numeric(tmp[paste(vardep)])
			}

			nl <- as.integer((agemaxi-agemini)+1)
			
			if(is.null(supvar)) { 
				rest <- matrix(data=0,nrow=nl,ncol=length(varindep)+3)
				dimnames(rest) <- list(c(), c("IDPERS", "AGE", vardep, varindep)) 
			}
			else { 
			
				rest <- matrix(data=0,nrow=nl,ncol=(length(varindep)+4))			
		
				dimnames(rest) <- list(c(), c("IDPERS", "AGE", vardep, varindep, supvar)) 
				
				sv <- as.integer(tmp[supvar]) 
				rest[1:nl,supvar] <- sv

			}
			rest[1:nl,"IDPERS"] <- idv[i]
			rest[1:nl,"AGE"] <- agemini:agemaxi
			rest[nl,vardep] <- ifelse(tmp[paste(vardep,"ST",sep="")]==1,1,0)	
			
			for(k in 1:length(varindep)) {
				if(tmp[paste(varindep[k],"ST", sep="")]==1 & tmp[varindep[k]] <= agemaxi) { 
					nlmi <- as.integer((tmp[varindep[k]]-agemini) + 1)
					rest[nlmi:nl,varindep[k]] <- 1				
				}
			}

        if(is.na(agemini) | is.na(agemaxi)) {
				stop("agemini/agemaxi NA id:", ids[i])
			}
		}
		res <- rbind(res, rest)
	}
	
	return(res)
}

