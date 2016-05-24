sim.reg<-function(n.nodes,n.edges)
{
z<-.C("Rregsim",
	as.integer(n.nodes),
	as.integer(n.edges),
	mat=integer(n.nodes*n.nodes), PACKAGE="brainwaver")

mat<-z$mat
mat<-matrix(mat,n.nodes,n.nodes)

return(mat)


}








