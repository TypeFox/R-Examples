

################################################
# create interactions
create_interactions <- function (y_,xobs_,xall_,index_int_,min_int_cor_,maxcols_){ 
	.Call("create_interactions_cpp", 
		y_,xobs_,xall_,index_int_,min_int_cor_,maxcols_, 
		PACKAGE = "miceadds")
					}


# res <- create_interactions( y_=y[ry] , xobs_=as.matrix( x[ry,] ) , 
#			    xall_=as.matrix(x) ,index_int_ = as.matrix(dfr), 
#			    min_int_cor_= min.int.cor , maxcols_= min(nrow(dfr),5000) )
					# total number of interactions