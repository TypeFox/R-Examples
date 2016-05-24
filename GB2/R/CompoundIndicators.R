# At-risk-of-poverty threshold
arpt.cgb2 <- function(prop, shape1, scale, shape2, shape3, pl0, pl, decomp="r"){
     median <- qcgb2(0.5, shape1, scale, shape2, shape3, pl0, pl, decomp)   # scaled median
     return(prop*median)
}  
         
# At-risk-of-poverty rate
arpr.cgb2 <- function(prop, shape1, shape2, shape3, pl0, pl, decomp="r"){
     return(pcgb2(arpt.cgb2(prop, shape1, 1, shape2, shape3, pl0, pl, decomp),shape1, 1, shape2, shape3, pl0, pl, decomp))
}  

# Relative median poverty gap
rmpg.cgb2 <- function(arpr, shape1, shape2, shape3, pl0, pl, decomp="r"){
 return(1-qcgb2(arpr/2, shape1, 1, shape2, shape3, pl0, pl, decomp)/qcgb2(arpr, shape1, 1, shape2, shape3, pl0, pl, decomp))   
}

# Quintile share ratio
qsr.cgb2 <- function(shape1, shape2, shape3, pl0, pl, decomp="r") {
	q20 <- qcgb2(0.2, shape1, 1, shape2, shape3, pl0, pl, decomp)
	q80 <- qcgb2(0.8, shape1, 1, shape2, shape3, pl0, pl, decomp)
	return((1-incompl.cgb2(q80, 1, shape1, 1, shape2, shape3, pl0, pl, decomp))/incompl.cgb2(q20, 1, shape1, 1, shape2, shape3, pl0, pl, decomp))
	}
	
# The four indicators and the median
main.cgb2 <- function(prop, shape1, scale, shape2, shape3,pl0,pl,decomp="r"){
	arpr <-arpr.cgb2(prop, shape1,shape2, shape3, pl0, pl,decomp=decomp)
	main <- c(qcgb2(0.5, shape1, scale, shape2, shape3,pl0, pl,decomp=decomp),
	          moment.cgb2(1,shape1, scale, shape2, shape3, pl0, pl,decomp=decomp),
		  arpr ,
		  rmpg.cgb2(arpr,shape1,shape2, shape3, pl0, pl,decomp=decomp),
		  qsr.cgb2(shape1,shape2, shape3, pl0, pl,decomp=decomp))
	names(main) <- c("median","mean","arpr","rmpg","qsr")
	return(main)
}