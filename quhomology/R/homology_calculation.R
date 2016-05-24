#first define the subfunctions

#this is used to calculate the left elementary matrix X
findX <- function(A){
  X <- diag(nrow(A)) #create Identity matrix
  G_X <- cbind(A,X) #combine identity matrix with original matrix for Gaussian elimination
  GX_res <- GaussianElimination(G_X) #Gaussian elimination
  X <- GX_res[, (ncol(GX_res)-nrow(A)+1):ncol(GX_res)] #extract the left matrix, X
  return(X)
}


push_down <- function(D){
  for(i in 1:(length(D) - 1)){
    a <- D[i]
    b <- D[i + 1]
    if(a!=1||b!=1){
      d <- numbers::GCD(a,b)
      if(ifelse(T==all.equal(a,0),T,F) || ifelse(T==all.equal(b,0),T,F)){
        d <- ifelse(ifelse(T==all.equal(a,0),T,F),abs(b),abs(a))
      }
      if(d!=0){
        alpha <- a/d
        D[i] <- d;
        D[i + 1] <- -(b * alpha)
      } else{
        D[i] <- 0
        D[i+1] <- 0
      }
    }
  }
  return(D)
}



check_more_push <- function(D){
  for(i in 2:length(D)){
    if(D[i] < D[i - 1]){
      if(D[i]!=0){
        return(TRUE)
      } else{
        return(FALSE)
      }
    }
  }
  return(FALSE)
}



#this calculates the smith normal form based on the calculation of the Hermite Normal Form.
#In particular, HNF((HNF(A))^T)
smith <- function(S){
  S <- numbers::hermiteNF(S)$H
  S <- t(numbers::hermiteNF(t(S))$H)
  D <- diag(S)
  D <- push_down(D)
  for(i in 1:length(D)){
    D[i] <- ifelse(D[i]<0, -D[i], D[i])
  }
  j <- 1
  more_push <- check_more_push(D)
  while(more_push){
    D <- push_down(D)
    for(i in 1:length(D)){
      D[i] <- ifelse(D[i]<0, -D[i], D[i])
    }
    print(j)
    j <- j+1
    more_push <- check_more_push(D)
  }
  for(i in 1:length(D)){
    D[i] <- ifelse(D[i]<0, -D[i], D[i])
  }
  diag(S) <- D
  return(S)
}



matrix_rank <- function(A){
  A <- GaussianElimination(A)
  A <- unique(A)
  return(nrow(A) - 1)
}



row_space <- function(B){
  B <- t(numbers::hermiteNF(t(B))$H)
  return(B)
}




#here is the main function to calculate the homology
homology <- function(degree, k, quandle=TRUE,return_values = FALSE){
  if(degree < 2){
    print(paste0("we cannot calculate the degenerate homology groups H_",degree,". Please choose a higher group."))
    return(NA)
  }
  boundary_F <- boundary_matrix(degree + 1, k, quandle)
  boundary_G <- boundary_matrix(degree, k, quandle)
  rho <- matrix_rank(boundary_G) #first, this calculates the rank of the matrix G. This removes the need to calculate D and Y later.
  q <- nrow(boundary_G)
  q_rho <- q - rho
  boundary_G <- findX(boundary_G)
  boundary_G <- boundary_G[(rho + 1):q, ] #only take the rows that map to zero (i.e. only the boundaries)
  boundary_F <- row_space(boundary_F) #identify the cycles, i.e. remove the boundaries via Gaussian elimination
  boundary_G <- round(boundary_F %*% MASS::ginv(boundary_G)) #calculate N. Details in documentation
  boundary_G <- smith(boundary_G)
  Delta <- diag(boundary_G) #Extract the values necessary for the output.
  
  if(return_values){
    return(Delta)
  } else{
    #the following is the output depending on the number of zeroes.
    output_results(ifelse(quandle,"quandle","rack"),Delta,degree,k)
  }
  return(NULL)
}




degenerate_homology <- function(degree, k, return_values = FALSE){
  if(degree < 3){
    print(paste0("we cannot calculate the degenerate homology group H_",degree,". Please choose a higher group."))
    return(NA)
  }
  boundary_F <- boundary_matrix_degenerate(degree + 1, k)
  boundary_G <- boundary_matrix_degenerate(degree, k)
  rho <- matrix_rank(boundary_G) #first, this calculates the rank of the matrix G. This removes the need to calculate D and Y later.
  q <- nrow(boundary_G)
  r <- ncol(boundary_G)
  q_rho <- q - rho
  X <- findX(boundary_G)
  Z <- X[(rho + 1):q, ] #only take the rows that map to zero (i.e. only the boundaries)
  B <- row_space(boundary_F) #identify the cycles, i.e. remove the boundaries via Gaussian elimination
  N <- round(B %*% MASS::ginv(Z)) #calculate N. Details in documentation
  S <- smith(N)
  Delta <- diag(S) #Extract the values necessary for the output.
  if(return_values){
    return(Delta)  
  } else{
    #the following is the output depending on the number of zeroes.
    output_results("degenerate",Delta,degree,k)
  }
  return(NULL)
}

output_results <- function(hom_type,Delta,degree,k){
  s <- length(Delta)
  l <- 0
  ones <- 0
  output <- c()
  for(i in 1:s){
    if(ifelse(T==all.equal(Delta[i],1),T,F)){
      ones <- ones + 1 #count number of ones
    } else if(Delta[i]!=0){
      l <- l + 1
      output <- append(output,Delta[i]) #count and extract nonzero and non-one values in the diagonal
    } else{
      break
    }
  }
  if(s>l+ones){#check if there are any values not equal to one or zero
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))), " ",hom_type, " homology group of R_",k," is isomorphic to Z^", s-(l+ones)," plus the following:"))
  } else{
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))), " ",hom_type, " homology group of R_",k," is isomorphic to the following:"))
  }
  if(l>0){ #if so, print out the resulting Z_n groups
    for(i in 1:l){
      print(paste0("Z_",Delta[ones+i],ifelse(i!=l," plus","")))
    }
  } else{
    print("0")
  }
  return(0)
}