
zscore <- function (sequence, simulations = NULL, modele, exact = FALSE, alphabet = s2c("acgt"), ... ){
  if (is.null(simulations)){
    if (modele=="base"){
      uni <- count(sequence, 1, freq = TRUE, alphabet = alphabet)
      di <- count(sequence, 2, freq = TRUE, alphabet = alphabet)
      rep1 <- rep(uni, 4)
      rep2 <- rep(uni, each = 4)
      rho <- di/(rep1*rep2)
      if (exact==FALSE){
        zscore <- ((rho-1)/sqrt(((1-rep1)*(1-rep2))/((length(sequence))*rep1*rep2)))
      }
      else if (exact==TRUE){
        rm(uni,di,rep1,rep2)
        uni <- count(sequence, 1, freq = FALSE, alphabet = alphabet)
        di <- count(sequence, 2, freq = FALSE, alphabet = alphabet)
        rep1 <- rep(uni, 4)
        rep2 <- rep(uni, each = 4)
        n <- length(sequence)
        e <- c(rep(c(0,rep(n/(n-1),4)),3),0)
        for (i in seq(1,4)){
          e[1+5*(i-1)] <- (n*(uni[i]-1))/(uni[i]*(n-1))
        }
        v=NULL
        for (i in seq(1,16)){
          if (i==1 | i==6 | i==11 | i==16){
            v[i] <- ((n^2*(rep1[i]-1))/(rep1[i]^2*(n-1)^2))*(1-rep1[i]+(n/rep1[i])*(1+2*(rep1[i]-2)/(n-1)+((n-4)*(rep1[i]-2)*(rep1[i]-3))/((n-2)*(n-3))))
          }
          else{
            v[i] <- ((n^3)/((n-1)^2*rep1[i]*rep2[i]-n^2/(n-1)^2+(n^3*(n-4)*(rep1[i]-1)*(rep2[i]-1))/(rep1[i]*rep2[i]*(n-2)*(n-3)*(n-1)^2)))
          }
        }
        zscore <- ((rho-e)/sqrt(v))
      }
    }
    else if (modele=="codon"){
      split <- splitseq(sequence)
      n <- length(split)
      pos <- sort(c(seq(3,length(sequence)-1,by=3),seq(4,length(sequence),by=3))) #position3-1 codons
      xy <- table(factor(splitseq(sequence[pos],word=2),levels=words(2)))
      n1 <- sapply(s2c('acgt'), function(x){length(which(substring(split,3,3)==x))})
      n1 <- rep(n1,each=4)
      n2 <- sapply(s2c('acgt'), function(x){length(which(substring(split,1,1)==x))})
      n2 <- rep(n2,4)
      n3 <- sapply(s2c('acgt'), function (x){sapply(s2c('acgt'), function(y){length(which(substring(split,1,1)==y & substring(split,3,3)==x))})})
      n3 <- as.vector(n3)
      e <- (n1*n2-n3)/n
      v <- (e+((1/(n*(n-1)))*(2*n3*(n1+n2-n1*n2-1)+n1*n2*(n1-1)*(n2-1))))-e^2
      zscore <- ((xy-e)/sqrt(v))
    }
    else{
      stop("analytical solution not implemented for this model")
    }    
  }
  else { #simulations = nb d'iterations
    rhopermut <- sapply(seq(simulations),function(x){rho(permutation(sequence =
  sequence, modele = modele, ...))})
    mean <- sapply(seq(dim(rhopermut)[1]),function(x){mean(rhopermut[x,])})
    var <- sapply(seq(dim(rhopermut)[1]),function(x){var(rhopermut[x,])})
    zscore <- ((rho(sequence)-mean)/sqrt(var))
  }
  return(zscore)
}
