################################################################################
# check consistency of input to din-method (data, q.matrix, ...)               #
################################################################################

check.input <- function( data , q.matrix , conv.crit = 0.001 , maxit = 100 ,
                    constraint.guess = NULL , constraint.slip = NULL ,
                    guess.init = rep(.2 , ncol(data) ) , slip.init = guess.init ,
                    weights = rep( 1 , nrow( data ) ) ,  rule = "DINA" ,
                    progress = TRUE ){

# Call: from din()
# Input: cf. din()
# Output: if possible cleaned arguments 
#  	else a warning message which leads to the termination of the procedure.
                    
################################################################################
# check consistency of data object                                             #
################################################################################

	# check for data classes matrix and data.frame
	if ((data.class(data) != "matrix") && (data.class(data) != "data.frame"))
    	   return(warning("data must be matrix or data frame"))
	data <- as.matrix(data)

	# check for data entries being dichotomous or missing
	gt <- data[ is.na( data ) == F ]
	# gt <- data[ ! is.na( data) ]

#	if(any(gt == 9||gt == 99||gt == .99)){
	if( sum(gt == 9 ) + sum(gt == 99 ) + sum(gt == .99) > 0 ){
  	return(warning("Recode your data! Only responses with values 0 or 1 (or NA) are valid.",
  	"\nMaybe missing values coded as 9, 99, .99.\n"))
  	}
	
  # return all response pattern not containing NA
	gt <- unique( gt[ gt %in% c(0,1) == F ] )
	if(length(gt) > 0){ 
		return(warning("Recode your data! Only responses with values 0 or 1 (or NA) are valid.\n")) 
					}
	
	# check for provision of row- and colnames
	if(is.null(rownames(data))) rownames(data) <- 1:nrow(data)
	if(is.null(colnames(data))) colnames(data) <- paste("Item",1:ncol(data),sep="")
  
################################################################################
# check consistency of q.matrix object                                         #
################################################################################

	# check for data classes matrix and data.frame
	if ((data.class(q.matrix) != "matrix") && (data.class(q.matrix) != "data.frame"))
       	return(warning("data must be matrix or data frame"))
	att.lbl <- attributes(q.matrix)$skill.labels
  q.matrix <- as.matrix(q.matrix)

  
	# return all response pattern not containing NA
	# gt_q <- data[ is.na( q.matrix ) == F ]
	gt_q <- q.matrix[ ! is.na( q.matrix ) ]
	
	# gt_q <- unique( gt_q[ gt_q %in% c(0,1) == F ] )
	gt_q <- setdiff( unique( gt_q ) , c(0,1) )
	if(length(gt) > 0){ return(warning("Check your Q-matrix! Only values 0 or 1 are valid.\n")) }
	
 	# check if q.matrix obtains same number of items as the data set
	if(nrow(q.matrix)!=ncol(data)){ return(warning("Check your Q-matrix! Number of assigned items (rows)
  		must fit the number of items in the data (columns).\n")) }

	# return warning message if there is a Zero-Row in the q.matrix
	rq <- rowSums(q.matrix)
	if (min(rq) == 0){
    	return(warning("Check your Q-matrix! The following items are not related to attributes:"
      		, "\n" , "Items " , paste( (1:(nrow(q.matrix)))[ rq == 0] , collapse=" , " ) , "\n" ))}

	# check for provision of row- and colnames
	if(is.null(rownames(q.matrix))) rownames(q.matrix) <- paste("Item",1:nrow(q.matrix),sep="")
	if(is.null(colnames(q.matrix))) colnames(q.matrix) <- paste("Skill",1:ncol(q.matrix),sep="")
	
  # check for provision of skill labels
  if(is.null(att.lbl)) attr(q.matrix, "skill.labels") <- colnames(q.matrix)
  if(length(att.lbl) != ncol(q.matrix) & length(att.lbl) != 0){
    attr(q.matrix, "skill.labels") <- colnames(q.matrix)
    warning("Unreasonable number of skill labels; skill labels replaced by colnames of q.matrix")
  }else{
    attr(q.matrix, "skill.labels") <- att.lbl
  }
  
################################################################################
# check consistency of arguments for parameter estimation routine              #
################################################################################

if(!is.numeric(conv.crit)||!is.numeric(maxit)) return(warning("Check your routine criteria"))
if(conv.crit<=0||maxit<1) return(warning("Check your routine criteria"))



################################################################################
# check consistency of constraint arguments for parameter boundaries           #
################################################################################

	# slip constraints see help files
	if(!is.null(constraint.slip)){                                                  #NULL permitted
		
	  if (any(is.na(constraint.slip))||any(!is.numeric(constraint.slip))||         #numeric values only
	    (!is.vector(constraint.slip) && (data.class(constraint.slip) != "matrix") && 
	    (data.class(constraint.slip) != "data.frame"))||                            #object typ
	    (length(constraint.slip %% 2 != 0) && ncol(constraint.slip)!=2)){				#two columns!
	       return(warning("check your error parameter constraints. See Help-files."))
	    }
	    if(is.vector(constraint.slip)) 
#	      try(constraint.slip <- matrix(constraint.slip, ncol=2, byrow=T))                 
	    if(data.class(constraint.slip) == "data.frame"){ 
			onstraint.slip <- as.matrix(constraint.slip) }
	    
	    if(any(duplicated(constraint.slip[,1]))||                                   #no duplicates
	    any(!constraint.slip[,1]%in%1:ncol(data))||                            		#first column may only be indicees
	    all(!(constraint.slip[,2]>=0 && constraint.slip[,2]<=1))){                  #all entries between 0 and 1
	       return(warning("check your error parameter constraints. See Help-files."))
	    }   
	}
	
	# guessing constraints see help files
	if(!is.null(constraint.guess)){                                                 #NULL permitted
	  if(any(is.na(constraint.guess))||any(!is.numeric(constraint.guess))||         #numeric values only
	    (!is.vector(constraint.guess) && (data.class(constraint.guess) != "matrix") 
	    && (data.class(constraint.guess) != "data.frame"))||          				#object typ
	    (length(constraint.guess)%%2!=0 && ncol(constraint.guess)!=2)){            	#two columns!                            
	       return(warning("check your error parameter constraints. See Help-files."))
	    }
#	    if(is.vector(constraint.guess)) 
#	      try(constraint.guess <- matrix(constraint.guess, ncol=2, byrow=T))                 
	    if(data.class(constraint.guess) == "data.frame"){ 
			constraint.guess <- as.matrix(constraint.guess)
									}
	    
	    if(any(duplicated(constraint.guess[,1])) ||                                 #no duplicates
	    any(!constraint.guess[,1] %in% 1:ncol(data))||                           		#first column may only be indicees
	    all(!(constraint.guess[,2]>=0&&constraint.guess[,2]<=1))){                  #all entries between 0 and 1
	       return(warning("check your error parameter constraints. See Help-files."))
	    }   
	}

################################################################################
# check consistency of init arguments                                           #
################################################################################

	# slipping initialization see help files
#	try({slip.init <- as.vector(slip.init)
#	     guess.init <- as.vector(guess.init)}, silent=T)
	if(!is.null(slip.init)){     
	if(any(is.na(slip.init))||
	  !all(is.numeric(slip.init))||
	  !all(slip.init>=0&&slip.init<=1)|| 
	  (length(slip.init)!=ncol(data)))
	 return(warning("Check your initial error parameter values. See Help-files."))
	}

	# guessing initialization see help files
	if(!is.null(guess.init)){     
	if(any(is.na(guess.init))||
	  !all(is.numeric(guess.init))||
	  !all(guess.init>=0&&guess.init<=1)|| 
	  (length(guess.init)!=ncol(data)))
	 return(warning("Check your initial error parameter values. See Help-files."))
	}
	
	
################################################################################
# check consistency of weight argument                                         #
################################################################################

	# weight see help files	
#	try(weights <- as.vector(weights), silent=T)     
	if(any(is.na(weights)) || is.null(weights) || !all(is.numeric(weights)) ||
	  !all(weights>0)|| (length(weights)!=nrow(data)))
	 return(warning("Check your specificated weights of the response patterns. See Help-files."))

################################################################################
# check consistency of rule argument                                           #
################################################################################

	# rule specification see help files	
	if(length(rule)!=1 && length(rule)!=ncol(data)){
		return(warning("Check the condensation rule for parameter estimation. The character string has
		to be of length 1 or of length ncol(data)."))
	}
#	try(if(!all(unique(rule)%in%c("DINA", "DINO")))  
	if(!all(unique(rule)%in%c("DINA", "DINO"))){  
		return(warning("Check the condensation rule for parameter estimation. Only \"DINA\" and \"DINO\" possible.")) }

################################################################################
# check consistency of progress argument                                       #
################################################################################

	# progress see help files	
	if(!(is.logical(progress))){
		return(warning("Check specification whether or not the estimation progress should be printed."))
					}
 
	return(list("data"=data, "q.matrix"=q.matrix,
	  "conv.crit"=conv.crit, "maxit"=maxit, "constraint.guess"=constraint.guess,
	  "constraint.slip"=constraint.slip, "guess.init"=guess.init,
	  "slip.init"=slip.init, "weights"=weights, "rule"=rule, "progress"=progress))
}