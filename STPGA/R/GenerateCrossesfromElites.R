GenerateCrossesfromElites <-
function(Elites, Candidates, npop, mutprob){
newcrosses<-lapply(1:npop, FUN=function(x){
	x1<-Elites[[sample(1:length(Elites),1)]]
	x2<-Elites[[sample(1:length(Elites),1)]]
	return(makeonecross(x1=x1,x2=x2,Candidates=Candidates,mutprob=mutprob))
})
	return(newcrosses)
}
