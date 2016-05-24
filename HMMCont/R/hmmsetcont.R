hmmsetcont <-
function (Observations, Pi1=0.5, Pi2=0.5, A11=0.7, A12=0.3, A21=0.3, A22=0.7, Mu1=5, Mu2=(-5), Var1=10, Var2=10) 
{
	Parameters<-matrix(c(Pi1, Pi2, A11, A12, A21, A22, Mu1, Mu2, Var1, Var2), byrow=T, ncol=10, nrow=1)
	colnames(Parameters) <- c("Pi1", "Pi2", "A11", "A12", "A21", "A22", "Mu1", "Mu2", "Var1", "Var2")
	Results<-matrix(NA, ncol=6, nrow=1)
	colnames(Results) <- c("Pal", "Pbe", "Pxi", "AIC", "SBIC", "HQIC")
	Viterbi<-matrix(0, nrow=(length(Observations)), ncol=1)
	B<-matrix(NA, nrow=(length(Observations)), ncol=2)
	hmm<-list(Observations=(Observations), Parameters=(Parameters), Results=(Results), Viterbi=(Viterbi), B=(B))
	class(hmm)<-"ContObservHMM"
	return(hmm)
}
