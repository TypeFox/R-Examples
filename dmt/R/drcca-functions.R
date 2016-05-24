drCCAcombine <- function(datasets, reg = 0, nfold = 3, nrand = 50)
{

  mat <- datasets # list of data matrices
  #reg: regularization value

  n <- nfold # number of cross validations

  iter <- nrand # number of iterations for random matrix generation

  if(reg < 0){reg <- 0}

  if(nfold <= 0){n <- 3} #use default value

  if(nrand <= 0){iter <- 50} #use default value

  m <- length(mat)  #number of  data matrices


##########################################################
#########################################################



## internal subroutine 1 , part
part <- function(data, nfold)
{
    mat <- as.matrix(data) # read the data file

    r <- nrow(mat)

    n <- nfold #the number for partition

   parts <- array(list(),n) #divide the data in n parts

   for(i in 1:n)
   {
    parts[[i]] <- as.matrix(mat[((((i-1)*r)%/%n) + 1): ((i*r)%/%n), ])

     #print(dim(parts[[i]]))
   }

   test <- array(list(),n)
   train <- array(list(),n)

   for(i in 1:n)
   {
      test[[i]] <- parts[[i]]

      for(j in 1:n)
      {
      if(j != i)
      train[[i]] <- rbind(train[[i]],parts[[j]])
      }
   }

     #partition <- new.env()
     #partition$train <- train
     #partition$test <- test
     #as.list(partition)

     partition <- list(test = test, train = train)
     return(partition);
}


###subroutine 2
## supporting script for subroutine 2

 drop <- function(mat,parameter)
{
      mat <- as.matrix(mat)

      p <- parameter


      fac <- svd(mat)

      len <- length(fac$d)
      ord <- order(fac$d) #sort in increasing order,ie smallest eig val first

      ord_eig <- fac$d[ord] #order eig vals in increasing order

     # cumulative percentage of Frobenius norm

     portion <- cumsum(ord_eig^2)/sum(ord_eig^2)

     ind <- which(portion < p) #ind r those who contribute less than p%

     #print('original number of dimensions')
     #print(length(portion))
     #print('number of dimensions dropped')
     #print(length(ind))
     #print(ind)

      #drop eig vecs corrsponding to small eig vals dropped
       num <- len -length(ind)

      fac$v <- as.matrix(fac$v[,1:num])

      fac$u <- as.matrix(fac$u[,1:num])

      mat_res <- (fac$v %*% sqrt(diag(1/fac$d[1:num],nrow = length(fac$d[1:num])))%*% t(fac$u))

      return(mat_res);

}


##Subroutine 2.
# To calculate gCCA for training sets, Regularized or unregularized gCCA
# input 1. concatenated training matrices
# input 2. Number of features in each training sets
# input 3. Regularization parameter, 0 will indicate unregularized gCCA

   reg_cca <- function(list,reg= 0)
{

   mat <- list # list of mats

   reg <- reg

   m <- length(mat)

   covm <- array(list(),m) #covariance matrices

   for(i in 1:m)
   {

    mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})

    covm[[i]] <- cov(mat[[i]])
   }

   #whitening of data

    white <- array(list(),m) #whitening matrices

    whiten_mat <- array(list(), m) #whitened data
    
    if (reg > 0) 
       {
            for (i in 1:m) 
                {
                whiten_mat[[i]] <- mat[[i]] %*% drop(covm[[i]],reg)
                white[[i]] <- drop(covm[[i]],reg)
                }
        }else 
           {
              for (i in 1:m) 
               {
                whiten_mat[[i]] <- mat[[i]] %*% SqrtInvMat(covm[[i]])
                white[[i]] <- SqrtInvMat(covm[[i]])
               }
           }
   
        tr_a <- concatenate(whiten_mat)
        tr_z <- cov(tr_a)
        eig <- svd(tr_z)
        x <- list(eigval = eig$d, eigvecs = eig$v, white = white)
        return(x)

   

}


### subroutine 6 

separate <- function(data,fea)
{

  data <- data #concatenated data

  fea <- fea #vector of column numbers

  m <- length(fea)

  mat <- array(list(), m)

       j <- 0
       for (i in 1:m) 
           {
            mat[[i]] <- data[, (j + 1):(j + fea[i])]
            j <- j + fea[i]
           }

    return(mat)

}


##suboutine 3.
#SUBROUTINE FOR TEST DATA

eigtest <- function(test,train,fea,reg=0)
{

   test <- test #concatenated test set

   train <- train #concatenated train set

    fea <- fea

    reg <- reg

   ##gcca of training data

     #create list, subroutine 6
     tr_mat <- separate(train,fea)

     cca <- reg_cca(tr_mat,reg) #subroutine 5.

     vecs <- cca$eigvecs #tranining eig vecs

     white <- cca$white #whitening matrices

     #print('check')
     #print(dim(white[[1]]))
    
    ## test matrices 

      mat <- separate(test,fea) #list of test data

    #whitening of data with training whitening matrices

    k <- length(fea)

    whiten_mat <- array(list(),k)  # array of matrices after whitening
 
    for(i in 1:k)
     {

       mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})

       whiten_mat[[i]] <- mat[[i]] %*% white[[i]]
     }    

     ##concatenating the matrices

     con <- concatenate(whiten_mat)

     #full covariance matrix
      z  <- cov(con)

      # it is now zy = lambda y
      # we use traning eig vecs in place of y & calculate eig vals for test

      v <- ncol(vecs)

    lambda <- matrix(,v) # create a matrix for lambdas

      for(i in 1:v)
      {
       lambda[i,] <- (t(vecs[,i])) %*% z %*% (vecs[,i])
      
      }

    return(lambda);

}





##SUBROUTINE 4.
#To calculate random matrices from training data and their eigen values
#Inputs :concatenated training matrix; feature vector; number of iterations
random_eig <- function(train,fea,iter,row)
{
     mat <- train #concatenated training matrices

     iter <- iter #number of iterations

     fea <- fea #  feature vector column

     row <- row #number of samples in original data

     #print(fea)

     dims <- sum(fea) #total number of features

     k <- length(fea) # number of matrices

      #first separate all matrices


      mats <- separate(mat,fea)

       #now generate random matrices

      ################################################################

      #create random matrix from multivariate normal distribution
      # using training data

      #################################################################


   #will create n samples of random matrices from training data
   #will calculate eigen values for each sample
   #store eigen values in a matrix, one row per sample

    mat_eigen <- matrix( ,iter,dims) # eig vals for each sample as each row
                                     # of this matrix


    for(n in 1:iter)  # n is running for number of samples

    {
    ran1 <- array(list(),k) #first set of random matrices

       for(i in 1:k)
          {
           ran1[[i]] <- random(mats[[i]],row)      
          
          # print(dim(ran1[[i]]))
          }

   #eig <- reg_cca(ran1)
   #mat_eigen[n,] <- eig$eigval #eig vals getting stored in a matrix

   #cross validation on each random set

           ranfull <- concatenate(ran1)
           
   #create training and test

          r_train <- ranfull[1:nrow(mat),]

          r_test <- ranfull[(nrow(mat)+1):row,] 
           

    # test eigen value

           mat_eigen[n,] <- eigtest(r_test,r_train,fea,0)

    }
  
  mat_eigen; #returning eig vals in a matrix

  return(mat_eigen)

}


#subroutine

# create random matrices

random <- function(mat,row)

  {
       mat <- mat

       row <- row #samples to generate

       v <- apply(mat,2,var)

       sigma <- diag(v,ncol = length(v))

        mu <- c(rep(0,length(v)))

       rand <- mvrnorm(row,mu,sigma)

        return(rand);

  }

## SUBROUTINE 5.
# will compare each test eigvals with whole random eigvals
# return the number of test eigen values which are greater than all random
#eigen values

#input: each row of test eigvals and corresponding matrix of random eigvals

 comp <- function(test_eigen_value,ran_eigen_value)
  {

     test_eig <- test_eigen_value  #one row

     ran_eig <- ran_eigen_value    #random rig val matrix

     p <- length(test_eig)

     r <- nrow(ran_eig)

      #print(dim(ran_eig))

      count <- 0
      for(i in 1:p)
      {
        if( (score <- length(which(test_eig[i] >= ran_eig))) >  (98*r*p)/100)

          {count <- count +1}else{break}
       }
      return(count);
   }

#SUBROUTINES END HERE#
################################################################
################################################################

  
  ##PRINT THE NAME OF FILES IN THE INPUT FILE
   print("Number of input matrices")
   print(m)

   fea <- c(rep(0,m)) # Number of features in each matrix

  # make the data matrices zero mean

    for(i in 1:m)
    {

     mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})


     fea[i] <- ncol(mat[[i]])
    }

  ## concatenate all matrices columnwise

     a <- concatenate(mat)    #use internal subroutine


  # CREATING TEST AND TRAINING SETS for concatenated data for NFOLD Cross
  #Validation
   
     res <- part(a,n)  #calling subroutine 1.

     test <- res$test  # all test sets

     train <- res$train # all training sets

     print("Number of test and training matrices created")

     print( len <-length(test))
 
  #----------------------------------------------#
 
  # EIGEN VALUES FOR TEST DATA AND EIG VECS OF TRAINING DATA

       p <- sum(fea)

   ##!!! have to check if m/n >> p, else there will 
   ##    be degeneracy in solution

       test_eigvals <- matrix(,n,p)

       for(i in 1:n)
        {
           test_eigvals[i,] <- eigtest(test[[i]],train[[i]],fea,reg)
        }

      #print(test_eigvals[1])
      #print("Dimension of the matrix of test eigen values")
      #print(dim(test_eigvals))


    #----------------------------------------------------#


 

  # CREATE RANDOM MATRICES for each Training set and calculate Canonical
   # Correlation

        ran_eig <- array(list(),n) #each matrix here will contain eigvals
                                   # of random matrix generated from
                                   # each train[[i]]

                                   #no of rows = number of iterations

       for(i in 1:n)     # n runs for nfold validation, ie number of train[[i]]
        {

        ran_eig[[i]] <- random_eig(train[[i]],fea,iter,nrow(a)) #subroutine 4.

        }

   # concatenate all ran_eig rowwise and use this matrix for 
   # comparison

       random_all <- ran_eig[[1]]
       if(n > 1)
            {
               for( i in 2:n)
               {
                 random_all <- rbind(random_all,ran_eig[[i]])
               }
            }


#-----------------------------------------------------#

  #NOW COMPARE EIGENVALS OF TEST SETS WITH THE EIGVALS OF RANDOM MATRICES
  # take the average of eigen values of all test data sets 

        avg_eigs <- apply(test_eigvals,2,mean)


        counts <- comp(avg_eigs,random_all) 

        ################################################################
        #Calculating drCCA for given data and returning the projected###
        #data with recommended features(counts)

         cca_data <- regCCA(mat,reg) ##calculating using the function regCCA
         projection <- cca_data$proj    
         print(dim(projection))
         proj <- projection[,1:counts]               
         return(list(proj=proj, n=counts));




}


########################


##########################################################################
# Regularized CCA
# We will drop those eigen values which contributes less than a threshold 
# in the frobenius norm of the matrix
#########################################################################

###INPUTS###

#input 1. list of  data matrices
#input 2. regularization parameter (0 to 1), default is 0

# regularization <- a value p  will drop smaller eigen values 
# contributing less than p% of total frobenius norm or 1-norm
# "A value of zero mean no regularization chosen"

###OUTPUT VALUES###

# A list containing three components
# 1. x$eigval <- canonical correlations 
# 2. x$eigvecs <- canonical variates
# 3. x$proj <- Projected data on canonical variates
# 4. x$meanvec <- array of columnwise means for each data matrix
# 5. x$white <- array of whitening matrices for each data matrix

 "regCCA" <- 
 function(datasets,reg = 0)
{
 
    mat <- datasets #list of data matrices
  
    m <- length(mat)
  
    #print('total number of files\n')
    #print(m);

    para <- reg # between 0 and 1

    if(para < 0 || para > 1)
    stop(" Value of 'reg' should be greater than zero and less than one")


##################################################
####SUBROUTINES####


drop <- function(mat,parameter)
{
      mat <- as.matrix(mat)

      p <- parameter

      fac <- svd(mat)

      len <- length(fac$d)
 
      ord <- order(fac$d) #sort in increasing order,ie smallest eig val first
   
      ord_eig <- fac$d[ord] #order eig vals in increasing order


     # cumulative percentage of Frobenius norm
     
     portion <- cumsum(ord_eig^2)/sum(ord_eig^2)

     #portion <- cumsum(ord_eig)/sum(ord_eig)


     ind <- which(portion < p) #ind r those who contribute less than p%


     #drop eig vecs corrsponding to small eig vals dropped
       num <- len -length(ind)
  
      fac$v <- as.matrix(fac$v[,1:num])

      fac$u <- as.matrix(fac$u[,1:num])

   mat_res <- (fac$v %*% sqrt(diag(1/fac$d[1:num], nrow = length(fac$d[1:num])))%*% t(fac$u))

     w <- list(mat_res = mat_res, num = num);

     return(w);
}



########SUBROUTINES END HERE
##################################################
    # mat is pointer to array of matrices

       covm <- array(list(),m) #array of covariance matrices

       fea <- c(rep(0,m)) # numebr of dimensions in each data matrix

      mean_m <- array(list(),m) # array of columnwise mean vector for each mat

        for(i in 1:m)
         {

          mean_m[[i]] <- apply(mat[[i]],2,mean)

           #subtracting columnwise mean
           mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})

           
           covm[[i]]<- cov(mat[[i]]); #covariance matrix

           fea[i] <- ncol(mat[[i]])
         }


# whitening of data...Regularized or Normal

     whiten_mat <- array(list(),m) # array of matrices after whitening

     cov_whiten <- array(list(),m) # covariance of whitened matrices

      white <- array(list(),m) #array of matrices,which does whitening

     d <- c(rep(0,m)) # to store the number of dimensions we kept in each
                      # data matrix after regularization

    # approximation of covariance matrices for regularization


      if(para > 0)
      {
             
               for(i in 1:m)
             {

               dummy <- drop(covm[[i]], para)

              whiten_mat[[i]] <- mat[[i]] %*% dummy$mat_res 
              #  print(det(cov(whiten_mat[[i]])))
              d[i] <- dummy$num #dimensions kept for whitening

              white[[i]] <- dummy$mat_res #storing whitening matrix
              }
      }else
      {
 
           for(i in 1:m)
           {
            whiten_mat[[i]] <- mat[[i]] %*% SqrtInvMat(covm[[i]])

            d[i] <- ncol(mat[[i]]) #all dimensions, no whitening

            white[[i]] <- SqrtInvMat(covm[[i]])

            #print('dim of sqrt_mat inverse mat')
            # print(det(cov(whiten_mat[[i]])))
            #print(dim(sqrt_mat_inv(covm[[i]])))
           }
      }


            if(sum(d) ==  sum(fea))
            {
                 print('Normal gCCA')
            }else{print('Regularized gCCA')}

            #print(fea)
            #print(d)
           # print(sum(d))

 
    #concactenating the whitened matrices

           a <- concatenate(whiten_mat) #internal function
   
           #print('concatenate')
           #print(dim(a))


   # z is covariance matrix of concatenated data

           z <- cov(a)
       
           eig <- svd(z) #use svd, eigen is not so good here

        # projected data

           proj_data <- a %*% eig$v


##-------------------------------------------##

     # whitening matrix %*% eig$v

       
     eig_wh <- array(list(),m) # whitened eig vecs for eaxh matrix

      if(para > 0)
      {
             
               j <- 0

               for(i in 1:m)
             {

               dummy <- drop(covm[[i]], para)

              eig_wh[[i]] <- dummy$mat_res %*% eig$v[(j+1) : (j+fea[i]),]

              j <- j + fea[i]
              #print(dim(eig_wh[[i]]))

              }
      }else
      {
 
           j <- 0
           for(i in 1:m)
           {
            
            eig_wh[[i]] <- SqrtInvMat(covm[[i]]) %*% eig$v[(j+1) : (j+fea[i]),]

              j <- j + fea[i]
              
            }
      }


       
        #fullwh <- eig_wh[[1]] #rbind all eig_wh
        #if(m >1)
        #{
           #for(i in 2:m)
           #{
            # fullwh <- rbind(fullwh,eig_wh[[i]])           
 
           #}
        #}

        #fullmat <- concatenate(mat)

        #wproj <- fullmat %*% fullwh
          
            


##-------------------------------------------##





     # make list of eigen values, eigen vectors and projected data

      gencorr <- list(eigval = eig$d,eigvecs = eig_wh, proj = proj_data, meanvec = mean_m, white = white)
      
      return(gencorr)

   

}


#######################


# this program is like sharedVar.R but it will calculate the variances
# with in each data sets. This finds the within data variation captured by
# cca and pca projections.

#input 1. datasets : list of data sets
#input 2. regcca : output of regCCA
#input 3. dim : number of dimensions to be used for back projection


#input 4. pca = FALSE : default parameter

# output : a list of 2 vectors containing within data variation
#          for cca projected data and for pca projected data

specificVar <- function(datasets, regcca, dim, pca = FALSE)
{


#########################################################
#####subroutines start heremi/doc/highlights/drafts2006/


      reverse <- function(proj,dir,dim)
      {

        proj <- proj #projected data
        dir <- dir # direction matrix

        dim <- dim # number of dimension to use

              if(dim > ncol(proj)) stop("Number of input dimension is greater than the number of dimensions in the projected data")



        p <- proj[,1:dim]

        d <- dir[,1:dim]

        inv_d <- solve(t(d)%*% d)%*% t(d)
             
        x <- p  %*% inv_d # same as dimension of proj
      
       return(x)
       }

#subroutine for pca

    withinpca <- function(mat,dim)
    {

     mat <- mat # list of data matrices

     dim <- dim # number of dimensions of projected data to be used to
              # regenerate original data sets

     m <- length(mat)


     fea <- c(rep(0,m)) # numebr of dimensions in each data matrix
     for(i in 1:m)
         {

           fea[i] <- ncol(mat[[i]])
         }

        #concatenate all matrices

          x <- concatenate(mat)


      #performing pca#

         
          pcax <- prcomp(x) #pca of x


          z <- pcax$x  # z is pca rotated data
          eig_z <- pcax$rotation # pca directions
 
      # reverse the pca projection, z[,1:dim]

          x_z <- reverse(z,eig_z,dim)

      # separating all matrices


             arr_xz <- array(list(),m) #storing matrices form x_z

             c <- 0
             for(i in 1:m)
             {

               t <- x_z[,(c+1):(c+fea[i])] #ith data set

               #storing matrices in a array
                arr_xz[[i]] <- t

                c <- c + fea[i]
              }
 
         
       # calculating within data  variation

            pc_var <- c(rep(0,m))

           for(i in 1:m)
           {
       
            var <- cov(arr_xz[[i]])

            dummy <- svd(var)

            pc_var[[i]] <- sum(dummy$d)

          }

           return(pc_var)

        }

#########################################################
# MAIN PROGRAM STARTS HERE



     mat <- datasets  # list of matrices

     cc <- regcca # output of regCCA
     eigV <- cc$eigvecs #, whitened projections directions
     proj <- cc$proj #, projected data from regCCA
     white <- cc$white #, whitening matrices from regCCA
     
    dim <- dim # number of dimensions to be used for back projection


     m <- length(mat)


           covm <- array(list(),m) #covariance matrices
           fea <- c(rep(0,m)) # number of dimensions in each data matrix

           for(i in 1:m)
           {
             
            mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})
           
            covm[[i]]<- cov(mat[[i]])  #covariance matrix

            fea[i] <- ncol(mat[[i]])
           }


        ##inverse whitening for projection directions
       
          for(i in 1:m)
          {
            eigV[[i]] <- solve(white[[i]]) %*% eigV[[i]]
          }

       ## and concatenating all eigV's rowwise

           eigfull <- eigV[[1]]
           for(i in 2:m)
           {
             eigfull <- rbind(eigfull,eigV[[i]])
           }

      # reverse process for regcca proj data for a given 'dim'

         fullrev_cca <- reverse(proj,eigfull,dim)

      # separate all matrices

          rev_cca <- array(list(),m) # list of regenerated data sets from
                                     # proj for given 'dim
          
             c <- 0
             for(i in 1:m)
             {

               t <- fullrev_cca[,(c+1):(c+fea[i])] #ith data set

               t <- t %*% solve(white[[i]])

               #storing matrices in a array
                rev_cca[[i]] <- t

                c <- c + fea[i]
              }

         # within data variation for regCCA projected data

           cc_var <- c(rep(0,m)) #variance for each data set

                 for(i in 1:m)
                 {

                   var <- cov(rev_cca[[i]])

                   dummy <- svd(var)

                   cc_var[[i]] <- sum(dummy$d)
                 }

         # within data variance for original data sets

           orig <- c(rep(0,m)) #for each data set

                 for(i in 1:m)
                 {
                    y <- svd(covm[[i]])

                    orig[i] <- sum(y$d)
                 }

         # mean of within data variance for  drCCA

           mcc <- mean(cc_var/orig)


    if(pca) {
      pc_var <- withinpca(mat,dim)

      mpc <- mean(pc_var/orig)

      result <- list(cc = sum(cc_var)/sum(orig), pc = sum(pc_var)/sum(orig), mcca = mcc, mpca=mpc)  

     return(result)
   } else{ 
     return(cc = sum(cc_var)/sum(orig), mcca = mcc)  
   }

  

}

####################################

   
#It will be used to calculate cca and pca and then to reverse the 
# the projected data back to original one using first few projected
# components; then calculating pairwise mutual information and comparing
# the result for drcca and PCA

#input 1. file: a file containing nanmes of the matrices
#inout 2. regcca_output : output of regCCA
#input 3. dimension : number of dimensions in projected data to be used to
#                     regenarate original data sets
#input 4. pca : if pca = TRUE , only then the pca variances will be calculated
#output : returns a list of three matrices, one for original data sets, one for cca and other for pca
#         entries are pairwise mutual information



"sharedVar" <- 
function(datasets,regcca,dimension,pca=FALSE)
{

        mat <- datasets # list of data matrices

        cc <- regcca  # projection direction from regCCA

        #print(length(cc))
        #print(dim(cc[[3]]))
        eigV <- cc$eigvecs #, whitened projections directions
        proj <- cc$proj #, projected data from regCCA
        white <- cc$white #, whitening matrices from regCCA
        xmean <- cc$meanvec #, array of columnwise mean for each data matrix


        dim <- dimension #dimensions to be used for back projection


        m <- length(mat)

#############################################################
#####subroutines start here


#subroutine 1
reverse <- function(proj,dir,dim)
{

     proj <- proj #projected data
     dir <- dir # direction matrix

     dim <- dim # number of dimension to use

   if(dim > ncol(proj)) stop("Number of input dimension is greater than the number of dimensions in the projected data")



     p <- proj[,1:dim]

     d <- dir[,1:dim]

     inv_d <- solve(t(d)%*% d)%*% t(d)
             
     x <- p  %*% inv_d # same as dimension of proj
      
     return(x)

}

### subroutine 2  ###

 xvar <- function(index,array)
 {

       ind <- index

       arr <- array

       len <- length(arr)

        var <- c(rep(0,(len))) # to store pairwise mi with each index

        for(i in (ind+1):len)
        {

           dummy <- t(arr[[ind]]) %*% arr[[i]]

           scov <- svd(dummy)

           var[i] <- sum(scov$d)

           #den1 <- svd(cov(arr[[ind]]))
           #d1 <- sum(den1$d)
           #den2 <- svd(cov(arr[[i]]))
           #d2 <- sum(den2$d)
           #var[i] <- sum(scov$d)/sqrt(d1*d2)
           
        }

       return(var)

 }

##pca subroutine 3

pwisepca <- function(mat,dim)
{

    mat <- mat # array of matrices
 
    dim <- dim 

    m <- length(mat)

    fea <- c(rep(0,m)) # numebr of dimensions in each data matrix
    for(i in 1:m)
         {

            mat[[i]] <- mat[[i]]
           #mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})

           fea[i] <- ncol(mat[[i]])
         }

        #concatenate all matrices

          x <- concatenate(mat)
        

      #performing pca#

         
         pcax <- prcomp(x) #pca of x


         z <- pcax$x  # z is pca rotated data
         eig_z <- pcax$rotation # pca directions
 
      # reverse the pca projection, z[,1:dim]

         x_z <- reverse(z,eig_z,dim)

      # separating all matrices


             arr_xz <- array(list(),m) #storing matrices form x_z

             c <- 0
             for(i in 1:m)
             {

               t <- x_z[,(c+1):(c+fea[i])] #ith data set

               

               #storing matrices in a array
                arr_xz[[i]] <- t

                c <- c + fea[i]
              }
 
         
       # calculating pairwise variation

        cross_pca <- matrix(0,(m-1),m)

       for(i in 1:(m-1))
       {
       
       cross_pca[i,] <- xvar(i,arr_xz)

       #print(offdiag_pca[[i]])
       }

       return(cross_pca)

}


################################################################




          
           fea <- c(rep(0,m)) # number of dimensions in each data matrix
           for(i in 1:m)
           {
             
            mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})
           

            fea[i] <- ncol(mat[[i]])

           }

                 
       ##inverse whitening for projection directions
       
          for(i in 1:m)
          {

            eigV[[i]] <- solve(white[[i]]) %*% eigV[[i]]
          }


       ## and concatenating all eigV's rowwise

           eigfull <- eigV[[1]]
           for(i in 2:m)
           {
             eigfull <- rbind(eigfull,eigV[[i]])
           }


       # reverse process for each data set for a given 'dim'

         fullrev_cca <- reverse(proj,eigfull,dim)

       # separate all matrices


          rev_cca <- array(list(),m) # list of regenerated adat sets from
                                     # y for given 'dim
          
             c <- 0
             for(i in 1:m)
             {

               t <- fullrev_cca[,(c+1):(c+fea[i])] #ith data set

               t <- t %*% solve(white[[i]])

               #storing matrices in a array
                rev_cca[[i]] <- t

                #print(rev_cca[[i]][1,1:4])

                c <- c + fea[i]
              }

      
    ########## CALCULATING PAIRWISE CROSS VARIATION ##########

    # pairwise cross variation for drcca, rev_cca

      cross_cca <- matrix(0,(m-1),m) 

      for(i in 1:(m-1))
      {

       cross_cca[i,] <- xvar(i,rev_cca) # var of x1%*%t(xj)

       #print(offdiag_cca[[i]])
      }

     #pairwise cross variation for original data sets

      orig <- matrix(0,(m-1),m)
      for(i in 1:(m-1))
      {

       orig[i,] <- xvar(i,mat)

      }

     # calculating the mean of pairwise variation

      mean_cc <- mean(cross_cca[cross_cca !=0]) # for drCCA
        


      mean_oo <- mean(orig[orig !=0]) # for original

           mcc <- mean_cc/mean_oo  #normalized mean for drCCA
           


   if(pca)
   {
      cross_pca <- pwisepca(mat,dim)

      mean_pc <- mean(cross_pca[cross_pca !=0]) # for PCA
      mpc <- mean_pc/mean_oo  #normalized mean for PCA

      result <- list(oo = orig, cc= cross_cca, pc = cross_pca, mcca =mcc, mpca=mpc)            

      return(result)
   }else{ return(list(oo = orig, cc =cross_cca, mcca =mcc))}

}



  




# this program will run sharedVar.R and specificVar.R for a given set of
# dimensions which will be used to project the data back for regCCA and PCA.

# it will also have an option to plot the variations for two programs

#input 1. datasets <- a list of data matrices
#input 2. regcca <- output of regCCA
#input 2. dimVector <- a vectors of integers, dimensions which will be used for
                       # for back projection
#input 3. plot = FALSE , if it is TRUE, a plot will be generated

plotVar <- function(datasets,regcca, dimVector, plot = FALSE)
{

        mat <- datasets        

        regcca <- regcca

        dims <- dimVector # list of values for dimensions

        m <- length(dims)

   ### pwiseVar.R

         pwcca <- c(rep(0,m)) # sum of pwise variation for each dims[i]

         pwpca <- c(rep(0,m)) # sum of pwise variation for each dims[i]

         for(i in 1:m)
         {

           dummy <- sharedVar(mat,regcca,dims[i],pca = TRUE)

           pwcca[i] <- sum(dummy$cc)/sum(dummy$oo) 

           pwpca[i] <- sum(dummy$pc)/sum(dummy$oo)
         }



   ### specificVar.R 

         withcca <- c(rep(0,m)) # sum of within data variation for each dims[i]

         withpca <- c(rep(0,m)) # sum of within data variation for each dims[i]

        for(i in 1:m)
         {

           dummy <- specificVar(mat,regcca,dims[i], pca = TRUE)

            withcca[i] <- sum(dummy$cc)


            withpca[i] <- sum(dummy$pc)
         }

   ###plotting the graphs 

  if(plot)
  {
   pw <- rbind(pwcca,pwpca)
   wn <- rbind(withcca,withpca)

   #split.screen(c(2,1))
   par(mfrow=c(2,1))
   plot(pw[1,],t='l',xlab='Dimensionality of the projection', ylab='Shared variation', main='Shared variation, drCCA vs PCA(red)', ylim=c(min(pw),max(pw)))
   lines(pw[2,], col ='red')
   #screen(2)
   plot(wn[1,],t='l',xlab='Dimensionality of the projection', ylab='Data-specific variation', main='Data-specific variation, drCCA vs PCA(red)', ylim=c(min(wn),max(wn)))
   lines(wn[2,], col='red')
   #close.screen(all = TRUE)
  }

 return(list(pw_cca = pwcca, pw_pca = pwpca, within_cca = withcca, within_pca = withpca))

   

}