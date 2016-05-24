`fact2dummy` <-
function (data, all=TRUE, lab="x") 
{
	dum.fcn <- function(x, all=TRUE){
        fine <- model.matrix(~x-1)
		colnames(fine) <- levels(x)
		if(!all) fine <- fine[,-1, drop=FALSE]
		fine
    }
	
	if(is.null(dim(data))){
		if(class(data)[1]=="numeric" || class(data)[1]=="integer" || class(data)[1]=="logical") oo <- cbind(1*data)
		else{
			oo <- dum.fcn(data, all=all)
			colnames(oo) <- paste(lab, colnames(oo), sep="")
		}
	}
	else{
        if(is.matrix(data)) oo <- data
        else{
            p <- ncol(data)
            n <- nrow(data)
		    vtype <- lapply(data, class)
		    out <- as.list(rep(NA,p))
		
		    for(i in 1:p){
                if(vtype[[i]][1]=="numeric" || vtype[[i]][1]=="integer" || vtype[[i]][1]=="logical") {
				    aa <- matrix(1*data[,i])
				    colnames(aa) <- names(data)[i]
				    out[[i]] <- aa
                }
	            else{
                    aa <- dum.fcn(data[,i], all=all)
				    colnames(aa) <- paste(names(data)[i], colnames(aa), sep="")
				    out[[i]] <- aa
                }
		    }
            oo <- unlist(out)
            oo <- matrix(oo, nrow=n)
            dimnames(oo) <- list(row.names(data), unlist(lapply(out, colnames)))
		}
	}	
	oo
}

