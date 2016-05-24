# Generate initial population of DE when searching for lower and upper bounds of each alpha-cut,
# subject to constraints specified in an interval matrix
.generateInitialPopulation<-function(numvars, lowerBoundsMatrix, upperBoundsMatrix)
{    
    n = length(lowerBoundsMatrix[1,]);
    NP = 10*numvars;
    resultado = matrix(0,NP,numvars);    
    # Accumulate lower bound of every row
    acumExtrIzquierdos = rowSums(lowerBoundsMatrix);
    occupied = matrix(0,n,n); # First position of each row: number K of occupied elements (indicated in columns 2 to K+1)    
                              # Rest of the row: indices of columns of lowerBoundsMatrix or upperBoundsMatrix containing non-zero elements
    for(i in 1:n){
      index = 2;
      for(j in 1:n){
        if(lowerBoundsMatrix[i,j] > 0 || upperBoundsMatrix[i,j] > 0){
          occupied[i,index] = j;
          index = index + 1;
        }
      }
      occupied[i,1] = index-2;
    }
        
    for(fila in 1:NP){
        index = 1;        
        for(i in 1:n){  # for each row we generate a distribution
          if(occupied[i,1] > 1){            
            probabilities = rep(0,occupied[i,1]);            # example: row i of occupied: occupied[i,] = [3, 8, 10, 15] -> 3 occupied in columns 8, 10 and 15
            selected = sample(occupied[i,1],occupied[i,1]); # example: there are 3 occupied, so "selected" is a permutation of them, e.g. [3,1,2]
            terminar = FALSE;
            while(!terminar){ # now assign the free probability mass after every variable has taken at least its lower bound
              sobrante = 1 - acumExtrIzquierdos[i];
              if(sobrante <= 0.0001){
                terminar = TRUE;
                for(k in 1:length(probabilities)){
                  fil = i; col = occupied[i, 1 + k];
                  probabilities[k] = lowerBoundsMatrix[i, col];
                }
              }
              else{
                for(k in 1:length(selected)){
                  fil = i; col = occupied[i,1 + selected[k]]; # columns of occupied start in column 2 of the row of occupied matrix
                  amplitud = upperBoundsMatrix[fil,col] - lowerBoundsMatrix[fil,col];
                  if(k < length(selected)){
                    aleat = runif(1, 0,min(amplitud, sobrante));
                    sobrante = sobrante - aleat;
                    probabilities[selected[k]] = lowerBoundsMatrix[fil,col] + aleat;
                  }else{  # last one: add all the probability mass that is has not been assigned yet
                    probabilities[selected[k]] = lowerBoundsMatrix[fil,col] + sobrante;
                    if(probabilities[selected[k]] < upperBoundsMatrix[fil,col]) terminar = TRUE;
                    # if this instruction fails, we cannot accumulate all the remaining mass at once -> retry this complete row
                  }
                }
              }
            }
            for(k in 1:(occupied[i,1]-1)){
           
              # The last one is not stored as a variable but we are sure it falls inside the feasible interval              
              # It is necessary to sort them from left to right, from smallest to the greatest column they are in lowerBoundsMatrix
              # (i.e. we have to keep the original ordering of occupied[i,]
              resultado[fila,index] = probabilities[k];
              index = index + 1;
            }
          }
        }
        stopifnot(index-1 == numvars);
    }
    return(resultado);
}

.generateInitialPopulationQuadraticRegression<-function(NP,centerpoint, side){
  
  numvars = 2;
  resultado = matrix(0,NP,numvars);
  b0 = centerpoint;
  if(side == "left"){
    for(fila in 1:NP){
      a = runif(1,-500, 500);
      if(a < 0){ 
        lder=(1-a)/b0;
        lizq = lder - 100;
      }else{
        lder = -2*sqrt(a/b0^2) - (2*a)/b0;
        lizq = lder - 100; 
      }
      b = runif(1,lizq, lder);
      resultado[fila,1] = a;
      resultado[fila,2] = b;
    }
  }
  else{ # side == "right"
    for(fila in 1:NP){
      decidir = runif(1,0,1);
      if(decidir < 0.5){
        lizqa = -1/(-1+b0);
        ldera = 1/(1-2*b0+b0^2);
        a = runif(1,lizqa,ldera);
        lizqb = (-1+a-a*b0^2)/(-b0+b0^2);
        lderb = (1-a)/b0;
        b = runif(1,lizqb,lderb);
      }
      else{
        lizqa = 1/(1-2*b0+b0^2);
        ldera = lizqa+100;
        a = runif(1,lizqa,ldera);
        lizqb = 2*sqrt(a/(b0^2)) - (2*a)/b0;
        lderb = (1-a)/b0;
        b = runif(1,lizqb, lderb);
      }
      resultado[fila,1] = a;
      resultado[fila,2] = b;
    }
  }
  return(resultado);
}