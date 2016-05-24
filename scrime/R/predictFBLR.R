`predictFBLR` <-
function(file, bin, kmax=10, int.level=2){
    	compare <- function(x,y) (all(x == y))*1

    	eval.logic.single <- function(logic,bin){
    		size <- logic[1]
    		if(size == 0) return(NULL)
    		lvars <- logic[2:(size+1)]
    		norm <- (lvars > 0)*1       # norm means that variables are not negated
    		if(size == 1)  
			return( (bin[,abs(lvars)] == norm)*1)
    		else 
			apply(bin[,abs(lvars)],1,compare,y = norm)
    	}	

    	eval.logic <- function(model,bin){
    		k <- model[1]
    		if(k == 0) 
			return(as.matrix(rep(1,n)))
    		else{
    			logics <- matrix(model[2:((int.level+1)*k+1)],nrow=int.level+1)
    			cbind(rep(1,n),apply(logics,2,eval.logic.single,bin=bin)) 
		}
    	}

    	predictBLR <- function(model,bin){
    		pnorm(eval.logic(model[1:modlen],bin)%*%model[(modlen+3):(modlen+3+model[1])])
	}

    	modlen <- ((int.level+1)*kmax+1)
    	n <- dim(bin)[1]
    	erg <- read.table(file=file)
    	test <- apply(erg,1,predictBLR,bin=bin)
    	pdachbay <- apply(test,1,sum)/dim(erg)[1]
	pdachbay
}

