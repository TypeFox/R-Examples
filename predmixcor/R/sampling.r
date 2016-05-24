
##############################################################################

training <- function (
     data,alpha,p,
     iter,iters.phi,ratio.phi,iters.lblphi,
     phi.range,r.shape,r.rate,no.r.bank,
     correction,no.phi.cor,commonr,approxim,
     no.r1=50)
              
{   #preliminary
    #the possible value of ``r'' can take
    rspace <- gen.alpha.IGamma(no.r.bank,r.shape,r.rate) 
    rspace1 <- gen.alpha.IGamma(no.r1,r.shape,r.rate)
    oiter <- 1
    I1 <- I2 <- phi <-  matrix(0, iter + oiter, ncol(data))
    label <- matrix (0, nrow( data ), iter + oiter)
    N1 <- N2 <- rep(0, iter + oiter)

    #r is a matrix used to store the value of r and r1(for only response) 
    #in the paper
    r <- rep (0, iter + oiter)
    r1 <- r

   
    r1[1] <- r[1] <- 10
    phi[1,] <- runif( ncol(data),0.1,0.9 )
    
    label [, 1] <- sample( x = c(1, 2), size = nrow( data ), replace = TRUE )
    n.obs <- nrow( data )
    N1[1] <- sum ( label[, 1] == 1 )
    N2[1] <- n.obs - N1[1]
    I1[1,] <- apply ( data[label[, 1 ] == 1,,drop=FALSE ], 2, sum )
    I2[1,] <- apply ( data[label[, 1 ] == 2,,drop=FALSE ], 2, sum )
    
    record <- matrix(0, 1, 6) 
    if(correction){
        
         #calculating B_alpha, which records the bound of 
	 #n11(col 2 and col3) 
	 #for each
         #n01(col1)
         indexalpha <- bound(alpha,sum(data[, 1] == 0),
	                     sum(data[,1] == 1))
         
         prob.cor.sum.sch <- function(indexalpha,N,phi.range,r,n)
         {             
             #Searching from the past values of r   
             pp <- prob.sch(N = N,r = r,record = record)
             if ( pp == FALSE )
             {
                 pp <- prob.cor.sum(indexalpha = indexalpha,N = N,
                                    phi.range = phi.range, r = r, 
                                    n = n)            
                 record <<- rbind( record, c(r, pp) ) 
             }
             pp
         }
         #the function used to calculate the correction factor for 
	 #updating ``r'',
         #the output is like that of fr, a vector containing the 
	 #likelihood value at
         #each possible ``r'' in ``rspace''
         lf.r.md <- function(rspace,p,indexalpha,phi.range,N,n)
         {
             logp <- rep (0, length( rspace ) )
             for ( i in 1:length( rspace ) )
                logp [i] <- p * log( prob.cor.sum.sch( 
		                     indexalpha = indexalpha, 
                                     N = N, r = rspace[i], 
				     phi.range = phi.range,
                                     n = n))
             logp
         }
    } 
    if(commonr) ix.r <- -1 else ix.r <- 1:ncol(data)  
    ######################################################################
    #starting markov chain
    for (i in oiter+(1:iter))
    {
     	label[, i] <- label[, i - 1]
     	N1[i] <- N1[i - 1]
     	N2[i] <- N2[i - 1]
     	I1[i,] <- I1[i - 1, ]
     	I2[i,] <- I2[i - 1, ]
     	r[i] <- r[i-1]
     	r1[i] <- r1[i-1]
     	
	for(i_lbl in 1:iters.lblphi){
      		#updating phi 
      		for ( j in 1:ncol(data) )
      		{    if(j > 1) thisr <- r[i] else thisr <- r1[i]
        	       phi[i,j] <- mh1.unif( f = lf.phi, x0 = phi[i - 1, j], 
	                         iter = iters.phi,
                                 r = thisr, I1 = I1[i,j],I2 = I2[i,j], 
                                 N1 = N1[i], N2 = N2[i],x.range=phi.range,
                                 ratio = ratio.phi)
           
      	        }
      		#updating labels
      		for (nr in 1:nrow(data))
      		{
        		if( label[nr,i] == 1) {
             		N1[i] <- N1[i] - 1 
	     		I1[i,] <- I1[i,] - data[nr,] 
	 		}
         		else{ 
	     		N2[i] <- N2[i] - 1 
	     		I2[i,] <- I2[i,] - data[nr,] 
	 		}
         		theta2 <- theta1 <- rep(0,ncol(data))
      
         		theta1[-1] <- ( r[i] * phi[i,-1] + I1[i,-1] ) / ( r[i] + N1[i] )
         		theta2[-1] <- ( r[i] * phi[i,-1] + I2[i,-1] ) / ( r[i] + N2[i] )
         		theta1[1] <- (r1[i] * phi[i,1] + I1[i,1]) / (r1[i] + N1[i])
         		theta2[1] <- (r1[i] * phi[i,1] + I2[i,1]) / (r1[i] + N2[i])

         		#counting n0.c1,n1.c1,n0.c2,n1.c2
         		label[nr, i] <- 1
         		N <- countN( data[, 1], label[, i] )	 
         		if(!approxim & correction)
           		prob.z1 <- prob.cor.sum.sch( indexalpha = indexalpha, N = N, 
                        	            r = r[i], phi.range = phi.range,
                                	    n = no.phi.cor)
				    
         		else prob.z1 <- 1			    
                                 
         		log.prob1 <- sum(dbinom( data[nr, ],1,theta1,log = TRUE)) +
                        				     p * log( prob.z1 ) + 
			     				     log((N1[i] + 1)/(N1[i] + N2[i] + 2))

         		label[nr, i] <- 2
         		N <- countN( data[, 1], label[, i] )
	 		if(!approxim & correction)
           		prob.z2 <- prob.cor.sum.sch( indexalpha = indexalpha, N = N, 
                                            r = r[i], phi.range = phi.range,
                                            n = no.phi.cor)
         		else prob.z2 <- 1                             
         		log.prob2 <- sum(dbinom( data[nr, ], 1, theta2, log = TRUE )) +
                        			p * log( prob.z2 ) + 
		        			log( ( N2[i] + 1 ) / ( N1[i] + N2[i] + 2 ) )
				
         		probratio <- exp( log.prob2 - log.prob1 )
         		label[nr, i] <- 1 + as.numeric( runif(1) > 1/(1 + probratio) )

         		if( label[nr,i] == 1){
            	 		 N1[i] <- N1[i] + 1
	    	  		I1[i,] <- I1[i,] + data[nr,] 
	 		}
         		else { 
	   	  		N2[i] <- N2[i] + 1
	    	  		I2[i,] <- I2[i,] + data[nr,] 
			}
	 
      		}
        }
     	#updating r
     	N <- countN(data[, 1], label[, i])

     	log.prob.r <- lf.r(rspace,phi[i,ix.r],I1[i,ix.r],I2[i,ix.r], 
                        N1[i], N2[i]) 
     	if(correction) {
        	log.prob.r <- log.prob.r + 
        	              lf.r.md(rspace = rspace, p = p, 
	                      indexalpha = indexalpha,
                              N = N,phi.range = phi.range,
		              n = no.phi.cor) 
         
    	}
      
     	r[i] <- sample(rspace, 1, 
                    prob = exp(log.prob.r - LOG.sum.exp(log.prob.r)))
     	#updating r1
     	if(commonr) r1[i] <- r[i]
     	else{
             log.prob.r1 <- lf.r(rspace1,phi[i,1],I1[i,1],I2[i,1],N1[i],N2[i]) 
             r1[i] <- sample(rspace1,1,
	                 prob = exp(log.prob.r1 - LOG.sum.exp(log.prob.r1)))
        }
   }
   list(label=label,I1=I1,I2=I2,N1=N1,N2=N2,phi=phi,r=r,r1=r1,alpha_set= 
   rspace,alpha0_set=rspace1)
}

#############################################################################
#this function generates no.alpha alpha's from Inverse Gamma distribution 
#with shape and rate parameters, by taking the quantiles corresponding 
#to the the probabilities spaced equally from w to 1-w, 
#where w = 1/no.alpha/2
gen.alpha.IGamma <- function(no.alpha,shape,rate)
{    w <- 1/no.alpha/2
     1/qgamma(seq(1 - w, w, length = no.alpha),shape,rate)

}

#########################################################################
#this function is used to search whether the probability we want has been
#calculated previously to save time
prob.sch = function(N,r,record)
{
    N1 = N
    N2 = c(N[3:4],N[1:2])
    k = nrow(record)
    found = FALSE
    while ( (!found) & k >= 1 )
    {
        if ( all( record[k,1:4] == N1, r == record[k,5] )  ||
             all( record[k,1:4] == N2, r == record[k,5] )
           )
        {
            found <- record[k,6]
            break
        }
        else k <- k - 1
    }   
    found
}

##########################################################################
#The function for calculating the bound of n11 given n01
#n01 is stored in the first column of the output matrix, the corresponding
#starting and ending points of n11 are stored in the second and the 
#third column.
bound <- function(alpha,n0,n1)
{

   n <- n0+n1
   x.bar <- n1 / n
   x.std <- sqrt(x.bar^2*n0 + (1- x.bar)^2 * n1)
   #the function used to calculate the correlation for 
   #given n01 and n11, in
   #case that the bound determined by equation is possible to 
   #have an error of 1
   calcorl=function(n01,n11)
   {
      if(n01+n11 == 0 || n01+n11 == n)
      0
      else
      abs((-x.bar*n01+(1-x.bar)*n11) / 
          ( x.std*sqrt((n01+n11)-(n01+n11)^2/n) ) )

   }

   indexalpha=matrix(0,n0+1,3)
   indexalpha[,1]=0:n0
   #specially for n01=0
   indexalpha[1,2]=0
   indexalpha[1,3] <- min(
        floor( 1 / ( 1/n + ( (1-x.bar) / ( alpha*x.std ) )^2 ) ) , n1)

   for(n01 in 1:n0)
   {
    tmp.m=(1-x.bar)/n01
    tmp.v=(alpha*x.std/n01)^2
    upper=min(floor ( 1/( tmp.m+0.5*tmp.v -
            sqrt( (1/4*tmp.v+tmp.m-1/n)*tmp.v ) )-n01 + 1e-4 ),
              n1)

    indexalpha[n01+1,3] <- upper - ( calcorl(n01,upper) > alpha )

    lower=max(
                ceiling ( 1/( tmp.m+0.5*tmp.v +
                sqrt( (1/4*tmp.v+tmp.m-1/n)*tmp.v ) )-n01-1e-4)
            ,0)
    indexalpha[n01+1,2] <- lower + ( calcorl(n01,lower) > alpha )
   }

    indexalpha
}
#########################################################################
#This function counts n0.c1,n1.c1,n0.c2,n1.c2
countN <- function(data,label)
{
    c ( sum( (label == 1)*(data == 0) ),
        sum( (label == 1)*(data == 1) ),
        sum( (label == 2)*(data == 0) ),
        sum( (label == 2)*(data == 1) )
      )
}

#########################################################################
#the function used to calculate the log-likelidhood on each ``r'' 
#in ``rspace''
lf.r <- function(rspace,phi,I1,I2,N1,N2)
{
    logp <- rep(0,length(rspace))

    for( i in 1:length(rspace) )
       logp[i] <- sum(
            lgamma( rspace[i] * phi + I1 ) +
            lgamma( rspace[i] * phi + I2 ) +
            lgamma( rspace[i] * ( 1 - phi ) + N1 - I1 ) +
            lgamma( rspace[i] * ( 1 - phi ) + N2 - I2 ) -
            2 * lgamma( rspace[i] * phi ) -
            2 * lgamma( rspace[i] * (1-phi) )
            ) + length(phi) * ( 2 * lgamma( rspace[i] ) -
                lgamma( N1 + rspace[i] ) - lgamma( N2 + rspace[i] ) )
    logp
}
########################################################################
#defining the function of ``p'' for M-H sampling
lf.phi <- function( phi,r,I1,I2,N1,N2 )
{
    
    lgamma( r * phi + I1 ) +
    lgamma( r * phi + I2 ) +
    lgamma( r * (1-phi) + N1 - I1 ) +
    lgamma( r * (1-phi) + N2 - I2 ) -
    2 * lgamma( r * phi ) - 
    2 * lgamma( r * ( 1-phi ) )
}
