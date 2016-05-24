#' @title Compute the most probable path of states
#' @description This function is used by the MMCNS interface to perfom the \code{\link[HMM]{viterbi}} algorithm,
#'  which computes the most probable path of states for a sequence of symbols and a given Hidden Markov Model.
#' @param nbStates number of states (not including the NA state).
#' @param nbSymbols number of symbols (not including the NA symbol).
#' @param transProbs matrix containing the transition probabilities between the states.
#' @param emissionProbs matrix containing the emission probabilities of the states.
#' @param startProb vector of initial probability distribution of the states.
#' @param symbolSeq vector of symbol sequence.
#' @return MarkovViterbiEstimation returns a list containing:
#' \item{estimatedHMM}{A HMM, which is a list of 5 elements (see \code{\link[HMM]{initHMM}}).}
#' \item{stateSeq}{The estimated state sequence.}
#' @importFrom HMM viterbi initHMM
#' @seealso \code{\link[HMM]{viterbi}} \code{\link[HMM]{initHMM}}

.MarkovViterbiEstimation<-function(nbStates,nbSymbols,transProbs,emissionProbs,startProb,symbolSeq){

	#Initialisation des Variables et constantes
	matEmi=emissionProbs[-1,-1];
	matTrans=transProbs[-1,-1];
	stateNames=c(1:nbStates);

	#Construction des symboles au format caractere
	symbolNames=c()
	for(i in 1:nbSymbols){
		symbolNames=c(symbolNames,paste("G",i,sep=""))
	}

	#Calcul du modele
	estimatedHMM=initHMM(States=stateNames, Symbols=symbolNames, startProbs=startProb, transProbs=matTrans, emissionProbs=matEmi)

	#Prediction de la classification sur de nouvelles donnees
	symbolNamesSeq=symbolNames[symbolSeq-1]
	calculViterbi=viterbi(estimatedHMM,symbolNamesSeq)

	#Ajout des donnees manquantes dans le vecteur d'etats predit
	toRemove<-symbolSeq==1 
	stateSeq=rep(0,length(symbolSeq))
	stateSeq[!toRemove]=calculViterbi
	stateSeq=stateSeq+1
	
	return(list(estimatedHMM=estimatedHMM,stateSeq=stateSeq))
	
}
