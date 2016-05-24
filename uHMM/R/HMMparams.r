#' @title Hidden Markov Model parameter estimation
#' @description This function is used by the \code{\link{uHMMinterface}} to estimate parameters of a Hidden Markov Model.
#' @param stateSeq a numeric vector of state sequencing.
#' @param symbolSeq a numeric vector of symbol sequencing.
#' @return HMMparams returns a list containing :
#' \item{trans}{The transition matrix.}
#' \item{emis}{The emission matrix.}
#' \item{startProb}{The vector of initial probability distribution (initial states are supposed equiprobable).}
#' @seealso \code{\link{transitionMatrix}} \code{\link{emissionMatrix}}
#' @export

HMMparams<-function(stateSeq,symbolSeq){

  matTransition<-transitionMatrix(stateSeq)
  matEmission<-emissionMatrix(stateSeq,symbolSeq)

  # Proba initiales
  if (sum(which(stateSeq==1))==0){  # si aucune donnee dans la classe 1 (on suppose que l'etat initial ne peut pas etre l'etat des NA's)
    Nstates<-length(unique(stateSeq))
  }else{
    Nstates<-length(unique(stateSeq))-1
  }
  startProb=rep((1/Nstates),Nstates);
  
  return(list(trans=matTransition,emis=matEmission,startProb=startProb))
  
}
  
#' @title Transition matrix estimation
#' @description This function estimates the transition matrix of a (Hidden) Markov Model from a vector of state sequencing.
#' @param states a numeric vector of state sequencing.
#' @return Estimated transition matrix.
#' @seealso \code{\link{HMMparams}}
#' @examples 
#' states<-c(1,1,3,2,1,2,1,3)
#' A<-transitionMatrix(states)
#' A
#' @export


transitionMatrix<-function(states){

  if (sum(which(states==1))==0){  # si aucune donnee dans la classe 1
    Nstates<-length(unique(states))+1
  }else{
   Nstates<-length(unique(states))
  }
  
  M<-matrix(0,nrow=Nstates,ncol=Nstates)
  
  #Construction de la matrice en comptant le passage d'un etat a l'autre
  for(i in 1:(length(states)-1)){
    M[states[i],states[i+1]]=M[states[i],states[i+1]]+1;
  }
  
  #Calcul de la somme de chaque ligne de la matrice
  sumByRow<-rowSums(M)
  
  #Reduction de la matrice pour que tout soit compris entre 0 et 1
  for(i in 1:Nstates){
    if(sumByRow[i]!=0){
      M[i,]=M[i,]/sumByRow[i];
    }  
  } 
  
  return(M)
  
}

#' @title Emission matrix estimation
#' @description This function estimates the emission matrix of a Hidden Markov Model from vectors of state and symbol sequencing.
#' @param states a numeric vector of state sequencing.
#' @param symbols a numeric vector of symbol sequencing.
#' @return Estimated emission matrix.
#' @seealso \code{\link{HMMparams}}
#' @examples
#' states<-c(1,1,3,2,1,2,1,3)
#' symbols<-c(4,1,3,1,4,4,4,2)
#' B<-emissionMatrix(states,symbols)
#' B
#' @export

emissionMatrix<-function(states,symbols){
  
  if (sum(which(states==1))==0){  # si aucune donnee dans l'etat 1
    Nstates<-length(unique(states))+1
  }else{
    Nstates<-length(unique(states))
  }
  
  if (sum(which(symbols==1))==0){  # si aucune donnee dans le symbole 1
    Nsymbols<-length(unique(symbols))+1
  }else{
    Nsymbols<-length(unique(symbols))
  }
  M<-matrix(0,nrow=Nstates,ncol=Nsymbols)
  
  #Construction de la matrice en comptant le nombre de fois que l'on est dans un etat et dans un groupe
  for(i in 1:(length(states))){
    M[states[i],symbols[i]]=M[states[i],symbols[i]]+1;
  }
  #Reduction de la matrice pour que la somme des termes vaille 1
  M=M/sum(M);
  
  return(M)
  
}
	

	




