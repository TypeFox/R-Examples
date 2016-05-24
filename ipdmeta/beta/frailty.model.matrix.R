frailty.model.matrix <- 

function(random,data){

	process.formula <- function(f){
   	   str <- as.character(f)[2]
   	   str <- strsplit(strsplit(str,"\\|")[[1]][1],"\\(")[[1]][2]
   	   as.formula(paste(c("~",str),collapse=""))
	}

	   vars <- all.vars(random)
	   n <- length(vars)

	   group <- factor(data[,vars[n]])
	   
	   f.random <- process.formula(random)

	   if(attr(terms(f.random),"intercept")){
	      data$int <- 1
	      vars <- c("int",vars) 
	   }

	   if(length(levels(group))==1){
	      Z.matrix <- data[,vars[-length(vars)]]
	      Z.matrix <- unlist(Z.matrix)
	   }
	   else{

	      Z <- function(x){model.matrix(~-1+x:group)}

	      Z.matrix <- sapply(vars[-length(vars)],function(var){
	   	  Z(data[,var])
	      })
	   }

   Z.matrix <- matrix(as.vector(Z.matrix),nrow=nrow(data))

   return(Z.matrix)
}
