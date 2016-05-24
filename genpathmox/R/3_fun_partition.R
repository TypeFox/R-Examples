#' @title Defining binary partitions given a segmentation variable (factor). 
#' @details
#' Internal function. \code{partition} is called by \code{all.part.pls}.
#' @param x  single factor or a data.frame of segmenation variables.
#' @param \dots Further arguments passed on to \code{\link{partition}}.
#' @return list of matrices containing the all possibiles binary partions given a segmenation variable.
#' @keywords internal
#' @export

partition <- function(x,...) 
{
	cat		=	sapply(x, is.factor)
	x[cat] 	=	lapply(x[cat],factor)
	
	nc 		= ncol(x)                                         		
	split	= list()
	
	length(split) = nc
	
	for (j in 1:nc) 
	{                                       		
		
		j.pval.opt	= NULL                            		
		modnum 		= as.numeric(x[,j])               
		nmod		= length(levels(x[,j]))            
		
		if (nmod == 1) next                              
		
		mod = 1:nmod                                	
		if (nmod == 2) {splits = matrix(c(1,2),1,2)}          
		
		if (nmod > 2) 
		{                                      			
			if (class(x[,j])[1] == "ordered") 
			{       	
				splits = matrix(NA,(nmod-1),nmod)         
				np = 0                                     					
				for (m in (2:nmod)) 
				{                      
					np 	= np+1 
					two = rep(1,nmod)                    			
					two[mod>=m] = 2                      		    
					splits[np,] = two                                  
				}
			}
			
			else
			{
				npart		= 2^(nmod-1)-1                     			
				np 			= 0                                   					
				breakflag 	= FALSE                       			
				splits 		= matrix(NA,npart,nmod)          		
				for (k in (1:floor(nmod/2))) 
				{           			
					comb = comb(nmod,k)         	
					for (i in 1:nrow(comb)) 
					{
						np = np + 1
						if (np > npart) breakflag = TRUE    
						if (breakflag) break
						two = rep(2,nmod)              				
						for (l in 1:ncol(comb)) {two[mod==comb[i,l]] = 1}       
						splits[np,] = two             					
					} 
				}
			}
		}
		
		split[[j]] = splits
	} 
	list(split=split)  
}

