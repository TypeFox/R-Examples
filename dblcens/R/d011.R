####################################################################
# d011      compute the NPMLE of distribution function from doubly #
#           censored data. Now returns censoring distribution too. #
#           There is an option to return influence functions       #
#           (this may take long time and large memory).            #
# copyright The software is GNU GPL2 copyrighted                   #
#           Please send comments, bug report etc. to               # 
#           mai@ms.uky.edu  or                                     #
#                                                                  #
#           Mai Zhou                                               #
#           Department of Statistics                               #
#           University of Kentucky                                 #
#           Lexington, KY 40506-0027                               # 
####################################################################
d011 <- function(z, d, identical = rep(0, length(z)), 
      maxiter = 49, error = 0.00001, influence.fun = FALSE) {
# Written by Mai Zhou (mai@ms.uky.edu) and Li Lee (most of the C part). 
# Last revision Aug.18,1999. 
#
# New version returns the NPMLE of the 2 censoring distributions too.
# As an option, it will also try to compute influence functions for Sx, 
# 3 of them and the estimated var of Sx(t). 
#
# Turnbull (1976) The empirical distribution function with
# arbitrarily gruped, censored and truncated data. JRSS B 290-295.
#
# Chang, M. N. (1990). Weak convergence in doubly censored data.
# Ann. Statist. 18, 390-405
#
# Inputs: 
# z -- a vector of length n denoting observed times, (ties permitted)
# d -- a vector of length n that contains censoring indicator: 
#      d= 2 or 1 or 0, (according to z being left, no, right censoring)
# identical, maxiter and error, influence.fun are optional input. 
#
# maxiter and error are used to control the EM iterations. The iteration 
# will stop whenever the difference between 2 consecutive iteration is 
# less then error, or number of iteration exceeds maxiter.
#
# identical is a vector of length n that has values either 0 or 1.
# identical[i]=1 means even if (z[i],d[i]) is identical with
# (z[j],d[j]), some j!=i, they still stay as 2 observations, 
# (not 1 obs. with weight 2, which only happen if identical[i]=0 
# and identical[j] =0). One reason for this may be because they have 
# different covariates not shown here so that it has more flexibility 
# for regression applications.
#
# influence.fun is a logical flag, indicating if you want estimate of influence
# fun. of the NPMLE of lifetime dist. estimator. Default is FALSE.
#
# Before calling the function, make sure z is a vector contain only 
# real numbers and d is a vector contain only 2's or 1's or 0's.
# Finally, z and d must be of same length. 
#
# Output: time (of z according only to d=0,1, -1; where a -1 means those
#               added z(i) for (0,2) pattern, in that case, I took the 
#               mid point as the time)
#         status (the indicator, = 0, 1, or -1.)
#         surv, jump, (these are self explaining)
#         conv (the # of iteration and diff between last 2 iterations)
#
# Bug: When compute influence function, the number of nodes cannot 
#      be too large (due to difficulty of finding inverse matrix by solve(), 
#      below 300 is more reasonable.) We use censored observations as nodes.
#
# Example: if z=(1,2,3,4,5); d=(1,0,2,2,1) then d011(z,d) gets:
#            $time: 
#            [1] 1.0 2.0 2.5 5.0    (notice the times, (3,4), coresponding 
#                                   to d=2 are removed, and time 2.5 added
#            $status:               since there is a (0,2) pattern at
#            [1]  1  0 -1  1        times 2, 3. The status indicator of -1
#                                   show that it is an added time )
#            $surv:
#            [1] 0.5000351 0.5000351 0.3333177 0.0000000
#       
#            $jump:
#            [1] 0.4999649 0.0000000 0.1667174 0.3333177
#
#            $conv:
#            [1] 3.300000e+01 8.788214e-06           (32 iterations done)
#       The true MLE of surv is (1/2, 1/2, 1/3, 0) at times (1,2,2.5,5). 

       N <- n <-length(z) 
       if(n < 2) stop("need more observations")
       niceorder <- order(z,  - d)   # order them as z increase. When z 
       sortz <- z[niceorder]         # has ties, use d=2, 1, 0 to order 
       sortd <- d[niceorder]         # the tied z values. 
       dupsortz <- duplicated(sortz)  #see if there is dup in z. But even if
       argdiff <- c(1, diff(sortd))  #z's tie, if d is diff, it's not a tie.
       dupsortz[argdiff != 0] <- FALSE  #seems I need not dupsortz ==T &
       dupsortz[identical != 0] <- FALSE #also, do not collaps if identical !=0

       sortz <- sortz[dupsortz != TRUE] # get the unique values of sortz and 
       sortd <- sortd[dupsortz != TRUE] # sortd. the weight or duplicatility
                        # w[i] will need to be computed later in C function?
                        # may be I should try to compute the weight here
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


        if( all(sortd == 1) ) { 
               nn <- length(sortz)
               list( time = sortz, 
                     status = rep(1, nn), 
                     surv = ((nn-1):0)/nn, 
                     jump = rep(1/nn, nn), 
                     exttime = sortz,
                     extstatus = rep(1, nn),
                     extweight = weight,
                     extjump = rep(1/nn, nn),
                     extsurv.Sx = ((nn-1):0)/nn, 
                     surv0.Sy = 1,
                     jump0 = 0, 
                     surv2.Sz = 0,
                     jump2 = 0,
                     conv = c("no censoring" , 0 ) , 
                     VarFt = (((nn-1):0)/nn)*((1:nn)/nn)/nn )   }
        else {

        sur <- rep(9, length(sortz))
        jum <- sur 
        zext <- wext <- dext <- rep(9, length(sortd)+length(sortd[sortd>1.5]))
        tes <- .C("urnew010",
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
                PACKAGE = "dblcens")

       status <- tes[[2]][1:tes[[9]]]
       surv <- tes[[4]][1:tes[[9]]]
       jump <- tes[[5]][1:tes[[9]]]

       extstatus <- tes[[12]][1:tes[[8]]]
       exttime   <- tes[[11]][1:tes[[8]]]
       extweight <- tes[[13]][1:tes[[8]]]

       extjump <- rep(0, tes[[8]])
       extjump[extstatus != 2] <- jump
       extsurv <- 1- cumsum(extjump)


       dzero <- extstatus + 1          # Now compute the surv function of
       dzero[dzero > 1.5] <- 0         # censoring distributions 
       dzero[dzero < 0.5] <- 0

       jumpzero <- (extweight * dzero)/(n * extsurv)
####    jumpzero[which.na(jumpzero)] <- 0  # let 0/0=0. 
       jumpzero[is.na(jumpzero)] <- 0      # (changed for R compatibility)
       survzero <- 1 - cumsum(jumpzero)    # do I need to do 
                                           # if(length(which.na(jumpzero)))...?

       dtwo <- extstatus - 1
       dtwo[dtwo < 0.5] <- 0

       jumptwo <- (extweight * dtwo)/(n*(1 - extsurv))  
####       jumptwo[which.na(jumptwo)] <- 0    # let 0/0=0 
       jumptwo[is.na(jumptwo)] <- 0           # (use is.na() for R compatible)  
       survtwo <- rev(cumsum(rev(jumptwo)))

       vec1 <- extjump/(survtwo - survzero)   # since jump>0, I changed sign.
####       vec1[which.na(vec1)] <- 0        # vec1 is for use in computing IC
       vec1[is.na(vec1)] <- 0               # use is.na for R compatible
       vec1 <- cumsum(vec1) 
       IC1tu2<-statusofnode<-nodes<-var<-IC1tu<-IC2tu<-IC3tu<-NA
#####################################################    
# the following computation needs (when influence.fun=TRUE): 
#    nodes,  nodesurv,   nodesurvzero,  nodesurvtwo,  nodejump, 
#    nodejumpzero, nodejumptwo, vec1, n (now is length of nodes) 
#####################################################
  if(influence.fun) { 
########          if( no left censor) { use KM formular }
########          else { 

  nodestatus <- extstatus != 1                  #exclude d= 1 or -1
  nodestatus[extstatus == -1] <- FALSE 
  nodes <- exttime[nodestatus]                   
  nodejump <- extjump[nodestatus]
  nodesurv <- extsurv[nodestatus] 
  nodesurvzero <- survzero[nodestatus]
  nodesurvtwo <- survtwo[nodestatus]
  nodejumpzero <- jumpzero[nodestatus]
  nodejumptwo <- jumptwo[nodestatus]
  statusofnode <- extstatus[nodestatus] 
  vec1 <- vec1[nodestatus] 
  n <- length(nodejump)

         # set up the matrix of Fi functions, those are
         #                (Fi(t1,u1)...Fi(t1,un))
         #      Fi(t,u) = (   .      .     .    ) 
         #                (Fi(tn,u1)...Fi(tn,un))

  mmm <- 1/(nodesurvzero - nodesurvtwo)          # what to do with inf?

  F1tu <- matrix(-mmm, n, n, byrow = TRUE)
  F1tu[row(F1tu) < col(F1tu)] <- 0       #make it lower triangle  

      # creat a triangle matrix of   int_ti^tj  dSx/(Sy-Sz)

   integ <- matrix(vec1,n,n)           # vec1 comes from extjump
   integ <- integ - t(integ) 
   integ[row(integ) <= col(integ)] <- 0    #make it lower triangle
       
   fff1 <- t(matrix(1/nodesurv,n,n))    #what if inf happen?
   F2tu <- fff1 * integ 

  tjjj <- matrix(vec1,n,n)
  jjj <- t(tjjj)
  jjj[row(jjj) < col(jjj)] <- tjjj[row(tjjj) < col(tjjj)] 
  fff <- t(matrix(1/(1-nodesurv),n,n))                      #what if inf?
  F3tu <- fff * jjj
 
   # after the F funtion (matrix) been set, the following is the 
   # solution,    relatively simple. Question: how many points ? 

  gts <- F3tu %*% diag(nodejumptwo) - F2tu %*% diag(nodejumpzero) #notice sign
  IK <- diag(n) - gts 
 
  IKinv <- solve(IK)              # what if inverse do not exist?

  IC1tu <- IKinv %*% F1tu  
  IC2tu <- IKinv %*% F2tu
  IC3tu <- IKinv %*% F3tu 

#####################################################
# to compute the IC1(t, u) at u= those delta =1 points
# we need only to create F1tu for those u
#####################################################
  nodes1 <- exttime[extstatus == 1] 
  n1 <- length(nodes1)
  two <- survtwo[extstatus == 1] 
  zero <- survzero[extstatus == 1]
  F1tu2 <- matrix(1/(two - zero), n, n1, byrow=TRUE)
  F1tu2[matrix(nodes,n,n1) < matrix(nodes1,n,n1,byrow=TRUE)] <- 0
 
  IC1tu2 <- IKinv %*% F1tu2 
#
# Now compute the var of NPMLE at t= the nodes.
#
    sum1 <- rowSums(IC1tu2)/N       ##  sum1 <- apply(IC1tu2,1,sum)/N
    sum2 <- rowSums((IC1tu2)^2)/N   ## apply((IC1tu2)^2,1,sum)/N   # do I need to add weight?

    sum3 <- rowSums(as.matrix(IC2tu[,statusofnode == 0]))  ## apply(as.matrix(IC2tu[,statusofnode == 0]),1,sum)/N
    sum4 <- rowSums(as.matrix(IC2tu[,statusofnode == 0])^2)/N

    sum5 <- rowSums(as.matrix(IC3tu[,statusofnode == 2]))/N 
    sum6 <- rowSums(as.matrix(IC3tu[,statusofnode == 2])^2)/N 

   var <- sum2 + sum4 + sum6 - (sum1 + sum3 + sum5)^2

  } 

         if(influence.fun) 
         list(time = tes[[1]][1:tes[[9]]],
              status = status,
              surv = surv, 
              jump = jump, 
              exttime = exttime,
              extstatus = extstatus, 
              extweight = extweight, 
              extjump = extjump,
              extsurv.Sx = extsurv, 
              surv0.Sy = survzero, 
              jump0 = jumpzero, 
              surv2.Sz = survtwo,
              jump2 = jumptwo, 
              conv = c(tes[[6]], tes[[7]]),
              Nodes = nodes,
              NodeStatus = statusofnode, 
              IC1tu = IC1tu,
              IC1tu2 = IC1tu2, 
              IC2tu = IC2tu,
              IC3tu = IC3tu,
              VarFt = var )
          else 
          list(time = tes[[1]][1:tes[[9]]],
              status = status,
              surv = surv,
              jump = jump,
              exttime = exttime,
              extstatus = extstatus,
              extweight = extweight,
              extjump = extjump,
              extsurv.Sx = extsurv,
              surv0.Sy = survzero,
              jump0 = jumpzero,
              surv2.Sz = survtwo,
              jump2 = jumptwo,
              conv = c(tes[[6]], tes[[7]]))
    }
 } 
######################################################################
#The new version d011(). It adds two main features to 
#d009() (at Statlib): 1. It now computes the NPMLE of the two censored
#distributions too (in addition to the NPMLE of life distribution.)
#2. It now has an option: influence.fun. When you set it to TRUE,
#it will TRY to compute the Influence Function (three of them) 
#for the NPMLE of life distrbution, and the estimated var of it.
#(default is FALSE). When there is no censoring, or only right, or only
#left censoring this may not work.
#
#Due to the difficulty of inverting big matrix, the second feature
#works very slow for data with more then 500 censored (both right and left)
#observations. For data with 100 censored points, it is pretty fast.
#For example if you have a sample of size 1000, but among the 1000
#there are only 100 censored observationa, then it is OK.
#This should improve once S has a better solve() function.     
######################################################################
