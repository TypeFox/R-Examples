print.ds <-
function(x,type='', ...)
{
cat("\nCommand:\n ")
 print(x$Call)
#
	if(x$tipo=="O"){Type<-c("Ordinay Dual Scaling")}
	if(x$tipo=="A"){Type<-c("Forced Classification")}
	#if(type=="B"){Type<-c("Ignoring Criterion Item")}
	
cat("\nType of Analysis:\n", Type,"\n")
#
cat("\nResults:\n")
#
switch(x$tipo,
	O = {
		if(is.null(x$Out_O)){stop("\n Your ds object needs 'A', 'B' or 'C' printing modes!")}
		print(x$Out_O,row.names=FALSE,digits=3)
		cat("\nDistribution of Information Over",x$N.Comp, "Components:\n")
		print(x$Inf_O,row.names=TRUE)
		for(i in 1:x$N.Comp){
			cat("\nInter Item Correlation for Component",i,":\n")
			print(data.frame(q=x$Rij_O[,,i]),row.names=x$cn,digits=3)
			}
			},
    A = {print(x$Out_A,row.names=FALSE,digits=3) 
       	cat("\nDistribution of Information Over",x$N.Comp, "Components:\n")
		print(x$Inf_A,row.names=TRUE)
		for(i in 1:x$SolFC){
			cat("\nInter Item Correlation for Component",i,":\n")
			print(data.frame(q=x$Rij_A[,,i]),row.names=x$cn,digits=3)
			}

			})
#
if(type=='B'){print(x$Out_B,row.names=FALSE, digits=3)}
#if(type=='C'){print(x$Out_C,row.names=FALSE, digits=3)}
}
