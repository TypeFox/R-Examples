"KM" <-
function(fileName,toplot=TRUE,header=TRUE)
{
    echantillon<-read.table(fileName,header)
#     A data frame is a list of variables of the same length with unique
#     row names, given class '"data.frame"'.
    
    nb_patient<-length(echantillon["no"][[1]]);
    durees<-echantillon["sejour"][[1]]; #X_i
    evt<-echantillon["evt"][[1]]; #Z_i
      
    
    # Initializations
	U <- array(42,1);
	D <- array(42,1);
	N <- array(42,1);
	S <- array(42,1);
	H <- array(42,1);
	TH <- array(42,1);
	evt_trie <- array(42,1);
	
    durees_triees <- (sort(durees, index.return=TRUE));
    idurees <- durees_triees$ix;
    durees_triees <- durees_triees$x;
    
    #durees$ix is the vector of sorted indexes
    #construction of evt_trie
    for (i in 1:length(durees)){
	    evt_trie[i] <- evt[idurees[i]];    
	}
    
	#construction of D : number of side effects
    tmp <- 0;
	for (i in 1:length(durees)){
	    if (evt_trie[i]==1){
	        D[i] <- tmp+1; 	           
        }
        else {
	        D[i] <- tmp;  
	    }
	    tmp <- D[i]
	}
	print("D vaut");
	print(D);
	
	
	#construction of N : number of patients that have not had side effets arisen  	
	tmp <- length(durees)+1;
	for (i in 1:length(durees)){
	    N[i] <- tmp-1;
	    tmp <- N[i];
	}
	print("N vaut");
	print(N);
	
	
	#construction of H : instantaneous risk
	for (i in 1:length(durees)){
	    H[i] <- D[i]/N[i];
	}
	print("H vaut");
	print(H);
	
	lambda_p <- (length(durees_triees)-1)/sum(durees_triees);
	print("estimateur de lambda si exp exp")
	print(lambda_p)
	
	tmp <- 1
	#construction of S : survival function
	for (i in 1:length(durees)){
		if (evt_trie[i]==1){  #if a side effect has been reported
	        S[i] <- tmp * (1-H[i]);
	        TH [i]<- exp((-1)*lambda_p*durees_triees[i]);
	        tmp <- S[i];
        }
        else {   #if no side effect has been reported
	        S[i] <- tmp; 
	        # tmp <- S[i]  #inutile
	        TH [i]<- exp((-1)*lambda_p*durees_triees[i]);   
	    }
	}
	
	
	
	
	# plots the survival function
	print("les durees triees valent")
	print(durees_triees);
	print("la fonction de survie vaut")
	print(S);
	plot(durees_triees,S);
	lines (durees_triees,S,col="green");
	lines (durees_triees,TH,col="red");
	
	    
}

