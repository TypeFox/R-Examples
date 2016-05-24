d011ch <- function(z, d, K, konst, 
          identical = rep(0, length(z)), maxiter = 49, error = 0.00001) 
{
# Written by Mai Zhou (mai@ms.uky.edu), and Kun Chen(chen@ms.uky.edu). 
# Last revision July 18,1999. Gnu GPL2 copyright reserved. 
#
# Chang, M. N. (1990). Weak convergence in doubly censored data.
# Ann. Statist. 18, 390-405
# Chen, K and Zhou, M (2000). Nonparametric confidence interval and testing
#                             hypothesis for doubly censored data.
#
# maxiter and error are used to control the EM iterations. The iteration 
# will stop whenever the difference between 2 consecutive iteration is 
# less then error, or number of iteration exceeds maxiter.
#
# Before calling the function, make sure z is a vector contain only 
# real numbers and d is a vector contain only 2's or 1's or 0's.
# Finally, z and d must be of same length. I did not write a lot of built-in
# checkings, to keep things simple and save some cpu time.
#
# Output: time (of z according only to d=0,1, -1; where a -1 means those
#               added z(i) for (0,2) pattern, in that case, I took the 
#               mid point as the time)
#         status (the indicator, = 0, 1, or -1.)
#         surv, jump, (these are self explaining)
#         conv (the # of iteration and diff between last 2 iterations)

       # check the input of z(observed times) and d(censoring indicator)
       N <- n <-length(z)
       if(n <= 2) stop ("Need more observation")
       if(length(d) != n) stop ("length of z and d must agree")
       if(any((d!=0)&(d!=1)&(d!=2)))
        stop ("d must be 0(right-censored), 1(uncensored) or 2(left-censored)")
       if(!is.numeric(z))
         stop("z must be numeric values --- observed times") 
       

       niceorder <- order(z,  - d)   # order them as z increase. When z 
       sortz <- z[niceorder]         # has ties, use d=2, 1, 0 to order 
       sortd <- d[niceorder]         # the tied z values. 
       dupsortz <- duplicated(sortz)  #see if there is dup in z. But even
       argdiff <- c(1, diff(sortd))  #z's tie, if d is diff, it's not a tie
       dupsortz[argdiff != 0] <- FALSE   #seems I need not dupsortz ==T &
       dupsortz[identical != 0] <- FALSE #also, do not collaps if identical !=0

       sortz <- sortz[dupsortz != TRUE] # get the unique values of sortz and 
       sortd <- sortd[dupsortz != TRUE] # sortd. the weight or duplicatility
                         # w[i] will need to be computed later in C function

       count <- (1:length(dupsortz))[dupsortz != TRUE]
       weight <- diff( c(count, (1+length(dupsortz)) ) )

# change the last and first obs to uncensored as nesessory. 
# This is tedious but must be done in order to get an estimator that is a 
# proper distribution. I did not move this part into C-code since I feel
# there may be other ways of redefining the last right censored
# observations and first left censored observation. (flexibility!)

        m<-length(sortd) 
        d01 <- sortd[sortd < 1.5]
        last <- length(d01)
        if(d01[last] == 0) {
            i <- m 
            z01 <- sortz[sortd < 1.5]
            while(sortd[i] != 1 && i > 0) {    
#
# one more condition in while: sortz[i] >= z01[last] ?
#
                if(sortd[i] == 0 && sortz[i] == z01[last]) 
                               sortd[i] <- 1
                i <- i - 1 
            }
        }
        d12 <- sortd[sortd > 0.5]        
        if(d12[1] == 2) {
             i <- 1 
             z12<-sortz[sortd > 0.5] 
             while(sortd[i] != 1 && i <= m ) { 
                  if(sortd[i] == 2 && sortz[i] == z12[1]) 
                  sortd[i] <- 1 
                  i <- i + 1 
             } 
        } 

        if( all(sortd == 1) ) {   ## i.e. no censoring at all ##
               nn <- length(sortz)
               time <- sortz
               status <- rep(1, nn)
               surv <- ((nn-1):0)/nn
               dist <- 1 - surv
               jump <- rep(1/nn, nn)
               exttime <- time
               extweight <- weight
               extstatus <- status
               extsurv <- surv
               extjump <- jump
               if(K < min(sortz))
               {
                  konstjump <- rep((1-konst)/nn,nn)
                  konstdist <- cumsum(c(konst,konstjump))  #changed 4/2008  
                  konstdist <- konstdist[2:nn]
               }
               else if(K > max(sortz))
               {
                  konstjump <- rep(konst/nn,nn)
                  konstdist <- cumsum(c(0,konstjump))
                  konstdist <- konstdist[2:nn]
               }
               n1 <- length(sortz[sortz<=K])
               n2 <- nn - n1 
               konstjump1 <- rep(konst/n1,n1)
               konstjump2 <- rep((1-konst)/n2,n2)
               konstdist1 <- cumsum(c(0,konstjump1))
               konstdist1 <- konstdist1[2:n1]
               konstdist2 <- cumsum(c(konst, konstjump2)) 
               konstjump <- c(konstjump1, konstjump2)
               konstdist <- c(konstdist1, konstdist2)
	       loglik1 <- sum(log(jump)) 

               loglik2 <- sum(log(konstjump)) 
 
   # Compute the -2loglikehood ratio using { loglikelihood 
   # maximized under (Ho+H1) } - { loglikelihood maximized 
   # under the constraint of Ho }. 
               llratio <- 2*(loglik1 - loglik2)
               maxit <- 1 
        }  
        else {
               sur <- rep(9, length(sortz))
               jum <- sur
               wext <- dext <- rep(9, length(sortd)+length(sortd[sortd>1.5]))
               konstdist <- konstjump <- zext <- wext 
               llratio <- 0.1 
               tes <- .C("d011ch",
                          as.double(sortz),
                          as.integer(sortd),
                          as.character(dupsortz),
                          as.double(sur),
                          as.double(jum),
                          as.integer(maxiter),
                          as.double(error),
                          as.integer(length(dupsortz)),
                          as.integer(length(sortd)),
                          as.integer(length(sortd[sortd>1.5])),
                          as.double(zext),
                          as.integer(dext),
                          as.double(wext),
                          as.double(K),
                          as.double(konst),
                          as.double(konstdist),
                          as.double(konstjump),
                          as.double(llratio),
                          PACKAGE = "dblcens")

               time <- tes[[1]][1:tes[[9]]]
               status <- tes[[2]][1:tes[[9]]]
               surv <- tes[[4]][1:tes[[9]]]
               jump <- tes[[5]][1:tes[[9]]]
               extstatus <- tes[[12]][1:tes[[8]]]
               exttime   <- tes[[11]][1:tes[[8]]]
               extweight <- tes[[13]][1:tes[[8]]]

               extjump <- rep(0, tes[[8]])
               extjump[extstatus != 2] <- jump
               extsurv <- 1- cumsum(extjump)

               konstdist <- tes[[16]][1:tes[[8]]]
               konstjump <- tes[[17]][1:tes[[8]]]

               llratio <- tes[[18]]
               maxit <- tes[[6]]
       }
  
       list(time=time,
            status=status,
            surv=surv,
            jump=jump,
            exttime=exttime,
            extstatus=extstatus,
            extjump=extjump,
            extsurv.Sx=extsurv,
            konstdist=konstdist,
            konstjump=konstjump,
            konsttime=K,
            theta=konst,
            "-2loglikR"=llratio,
            maxiter=maxit)
}

