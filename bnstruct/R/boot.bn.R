# boot.bn <- function(x, node.sizes, B = 100, cont.nodes = c(),
#   ess = 1, verbose = FALSE, seed = 0, k.impute = 0, layering = NULL, 
#   method = c("sm","mmhc"), max.fanin = NULL, max.fanin.layers = NULL,
#   chi.th = 0.05, layer.struct = NULL )
# {
# 	set.seed(seed) 
# 	n.nodes <- ncol(x)
# 	n.cases <- nrow(x)
#   method <- match.arg(method)
#   
# 	boot.sample <- matrix(sample.int(n.cases,size = B*n.cases,replace=TRUE),n.cases,B)
# 	finalPDAG <- matrix(0,n.nodes,n.nodes)
# 	
# 	for( i in seq_len(B) )
# 	{
# 		if( verbose ) print(i)
#     
#     data <- x[boot.sample[,i],]
# 		if( verbose ) print("impute")
#     if( k.impute > 0 )
#       data <- knn.impute( data, k.impute, setdiff(1:n.nodes,cont.nodes) )
#     
#     # write.table( data, paste(i,".txt",sep=""), quote=F, row.names=F, col.names=F )
#     
#     if( method == "mmhc" )
#     {
#       if( verbose ) print("mmpc")
#       cpc.mat <- mmpc( data, node.sizes, cont.nodes, chi.th, layering, 
#                        layer.struct )
#       if( verbose ) print("hc")
#       dag <- hc( data, node.sizes, cpc.mat, cont.nodes, ess )
#     }
#     else if( method == "sm")
#     {
#       if( verbose ) print("sm")
#       dag <- sm(data, node.sizes, cont.nodes, max.fanin, layering,
#                 max.fanin.layers, ess)
#     }
# 
#     # print(dag)
# 		if( verbose ) print("cpdag")
#     finalPDAG <- finalPDAG + dag.to.cpdag( dag, layering )
#     # print( dag.to.cpdag( dag, layering ) )
# 	}
# 	
# 	return(finalPDAG)	
# }