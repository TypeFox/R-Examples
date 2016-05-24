# internal function from package crossdes, slightly modified
"myGF" <-
function(p,n){

  primen100 <- c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,
      61,67,71,73,79,83,89,97)
  if (!(p %in% primen100)){stop("p is not a prime number less than 100.")}

  # Now obsolete: require(conf.design); if( (p%%1) || (p < 2) ){stop("p is not a prime number.")}; if( primes(p)[length(primes(p))]!=p ){stop("p is not a prime number.")}

  if( (n%%1) || (n < 1) ){stop("n is not a positive integer.")}

  ord <- p^n                                     # order of the field 
  
  if( n==1 )                                     # trivial case
   { elemente <- c(0:(p-1))
     primpol <- NULL
     out<-list(elemente,primpol)}
  else
   {  
     a <- factor.comb(p,n)                       # The monic polynomials are constructed using a full factorial plan 
     g <- cbind( rep(1,ord), a )  
                                                 # evaluate polynomials for 0,1,...,p-1
                                                 # rows where h is all non-zero correspond to irreducible polynomials 
                                                 # get them into irr      
     h <- matrix(0, ord, p)
    
     for (i in 0:(p-1)){                                 
       a <- numeric(ord)
       for (j in 1:n){ 
         a <- a + g[,j]*(i^(n+1-j))  
       }
       h[,i+1] <- (a + g[,n+1])%%p
     }
     
     irr <- g[!apply(h==0, 1, any), , drop = FALSE]

     # Get the primitive roots. For all polynomials compute : p^1,p^2,...,p^(ord-1).
     # Cycle is completed after at most ord-1 powers. If ord-1 powers are needed, a primitive root is found.
     # Look for smallest power of p that has coefficients equal to (0,...,0,0,1).
     # First power is always equal to (0,...,0,1,0), don't need to check that.
    
     nirr <- nrow(irr)
     z1 <- numeric(nirr)
    
     for (i in 1:nirr){
       dummy <- c(numeric(ord-3),1,0)            # 1st power of (1,0); ord-3>0 since p,n>2. 
       j <- 0
       while( z1[i]==0 ){                        # ord-1st to 2nd power of (1,0)
         j <- j+1
         dummy <- myredu.modp(c(dummy,0), irr[i,], p)        
                                                 # division modulo p
         if( all(dummy==c(numeric(n-1),1)) ) {z1[i] <- j+1}       
           # dummy is a polynomial of degree n-1.
           # all(...)==TRUE if the result by division modulo p is (0,1), i.e. the number one.
           # Here rest(x^n : polynom) = rest(x*rest(x^n-1 : polynom))
           # get smallest power that equals (0,1)
       }
     }
    
     primpol <- irr[which(z1==(ord-1)),,drop=FALSE]            

     elemente <- matrix(0,ord,n,byrow=TRUE)
     elemente[2,n] <- 1                          # 1st element is =0,...,0; 2nd element is 0,...0,1;
     elemente[3,n-1] <- 1                        # 3rd element ist =0,..,0,1,0.
      
     dummy <- c(numeric(ord-3),1,0)              # 1st power of (1,0); ord-3>0 since p,n>2. 
     for( i in 4:ord ){                          # 2nd ord-2nd power of (1,0)
        dummy <- myredu.modp(c(dummy,0), primpol[1,], p)        
        elemente[i,] <- dummy
     }    
    
     out<-list(elemente,primpol)
    
   }  
  out
}
