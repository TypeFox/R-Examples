# Main function for calling the Bayesian Model Selection

bayesfac <- function(indnr, paramnr, SEx, SEy, xv, yv, chx, chy, zv, chz, SEz, vv, chv, SEv)
{   
  if (indnr == 2)
  {
    nmodelterms = paramnr 
    
    # calling selectterms function
    tmp = selectterms(indnr, paramnr, SEx, SEy) 
    Mx <- tmp[[1]]
    My <- tmp[[2]]
    
    # creating empty vectors to be filled in the next step with selected models 
    # based on te SEtestx/y results, out of all possible model combinations
    Mx_allvars = c()
    My_allvars = c()
    
    count = 1       
    for (ii in 1:nmodelterms)
    {
      M <- combs(1:17, ii) 
      Mx_allvars <- c(Mx_allvars, M[Mx[ii],])
      My_allvars <- c(My_allvars, M[My[ii],])
      count = count + ii
    }    
    print(c("Selected model terms (dx):"), quote=F)
    print(Mx_allvars)
    print(c("Selected model terms (dy):"), quote=F)
    print(My_allvars)
    
    # creating empty vectors to be filled in the next step with Bayes factors 
    # and best parameters calling polyfitbayes function
    Bf1 <- c()
    parambest1 <- c()
    Bf2 <- c()
    parambest2 <- c()
    
    M1 <- Mx_allvars
    count = 1   
    for (j in 1:nmodelterms)
    {
      tmp = polyfitbayes(indnr, xv, yv, chx, M1[count:(count+j-1)]) 
      bestm <- tmp[[1]]
      indexbestm <- tmp[[2]]
      
      Bf1[j] <- bestm     
      parambest1[j] <- indexbestm
      count = count + j;
    }  
    
    M2 <- My_allvars
    count = 1    
    for (j in 1:nmodelterms)
    {
      tmp = polyfitbayes(indnr, xv, yv, chy, M2[count:(count+j-1)])
      bestm <- tmp[[1]]
      indexbestm <- tmp[[2]]
          
      Bf2[j] <- bestm 
      parambest2[j] <- indexbestm
      count = count + j;
    }
    
    # printing the Bayes Factor for the best models for each nmodelterms
    print(c("Bayes factors of the best models per number of modelterms for dx:", Bf1), quote=F)
    print(c("Bayes factors of the best models per number of modelterms for dy:", Bf2), quote=F)   
  }
  
  ############################# indnr == 2 ends, indnr == 3 begins #################################
  
  if (indnr == 3)
  {
    nmodelterms = paramnr
    # calling selectterms function
    tmp = selectterms(indnr, paramnr, SEx, SEy, SEz) 
    Mx <- tmp[[1]]
    My <- tmp[[2]]
    Mz <- tmp[[3]]
#     Mxy <- tmp[[5]]  ## the user might want to comment lines 83-88 to obtain
#     Mxz <- tmp[[7]]  ## Bayes Factors of models with only two indicators in a three-variabe
#     Myx <- tmp[[6]]  ## model setting to see if e.g. variabel x and y are sufficient to model 
#     Myz <- tmp[[8]]  ## changes in dx or if the third variable z is significantly 
#     Mzx <- tmp[[9]]  ## improving the bayes factor. That is the user might want to 
#     Mzy <- tmp[[10]] ## comment this and other related parts out if he/she wants to 
                       ## compare Bayes Factors of models with two and three indicators 
                       ## commenting out these parts required uncommenting related parts
                       ## in selectterms (follow comment instruction in selectterms.R)
    
    # creating empty vectors to be filled in the next step with selected models 
    # based on te SEtestx/y results, out of all possible model combinations
    Mx_allvars = c()
    My_allvars = c()
    Mz_allvars = c()
    ## comment out lines 100-105 for comparing Bayes Factor of 
    ## two-variable and three-variable-models (see descriptions above)
#     Mxy_vars = c()   
#     Mxz_vars = c()   
#     Myx_vars = c()
#     Myz_vars = c()
#     Mzx_vars = c()
#     Mzy_vars = c()
    
    count = 1
    for (ii in 1:nmodelterms)
    {
      M <- combs(1:39, ii) 
      Mx_allvars <- c(Mx_allvars, M[Mx[ii],])
      My_allvars <- c(My_allvars, M[My[ii],])
      Mz_allvars <- c(Mz_allvars, M[Mz[ii],])
      ## comment out lines 116-121 for comparing Bayes Factor of 
      ## two-variable and three-variable-models (see descriptions above)
#       Mxy_vars <- c(Mxy_vars, M[Mxy[ii],])
#       Mxz_vars <- c(Mxz_vars, M[Mxz[ii],])
#       Myx_vars <- c(Myx_vars, M[Myx[ii],])
#       Myz_vars <- c(Myz_vars, M[Myz[ii],])
#       Mzx_vars <- c(Mzx_vars, M[Mzx[ii],])
#       Mzy_vars <- c(Mzy_vars, M[Mzy[ii],])
      count = count + ii
    }
    
    print(c("Selected model terms (dx):"), quote=F)
    print(Mx_allvars)
    print(c("Selected model terms (dy):"), quote=F)
    print(My_allvars)
    print(c("Selected model terms (dz):"), quote=F)
    print(Mz_allvars)
    ## comment out lines 130-135 for comparing Bayes Factor of 
    ## two-variable and three-variable-models (see descriptions above)
#     print(Mxy_vars)
#     print(Mxz_vars)
#     print(Myx_vars)
#     print(Myz_vars)
#     print(Mzx_vars)
#     print(Mzy_vars)
    
    # creating empty vectors to be filled in the next step with Bayes factors 
    # and best parameters calling polyfitbayes function
    Bf1 <- c()
    parambest1 <- c()
    Bf2 <- c()
    parambest2 <- c()
    Bf3 <- c()
    parambest3 <- c()
    
    M1 <- Mx_allvars
    count = 1   
    for (j in 1:nmodelterms)
    {
      tmp = polyfitbayes(indnr, xv, yv, chx, M1[count:(count+j-1)], zv)
      bestm <- tmp[[1]]
      indexbestm <- tmp[[2]]
      
      Bf1[j] <- bestm
      parambest1[j] <- indexbestm
      count = count + j;
    }  

    ## comment out lines 161-187 for comparing Bayes Factor of 
    ## two-variable and three-variable-models (see descriptions above)     
#     Bf1a <- c()
#     parambest1a <- c()
#     M1a <- Mxy_vars
#     count = 1    
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chx, M1a[count: (count+j-1)], zv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf1a[j] <- bestm
#       parambest1a[j] <- indexbestm
#       count = count + j;
#     } 
#
#     Bf1b <- c()
#     parambest1b <- c()    
#     M1b <- Mxz_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chx, M1b[count: (count+j-1)], zv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf1b[j] <- bestm
#       parambest1b[j] <- indexbestm
#       count = count + j;
#     } 
    
    M2 <- My_allvars
    count = 1    
    for (j in 1:nmodelterms)
    {
      tmp = polyfitbayes(indnr, xv, yv, chy, M2[count:(count+j-1)], zv)
      bestm <- tmp[[1]]
      indexbestm <- tmp[[2]] 
      
      Bf2[j] <- bestm
      parambest2[j] <- indexbestm
      count = count + j;
    } 
    
    ## comment out lines 204-230 for comparing Bayes Factor of 
    ## two-variable and three-variable-models (see descriptions above)
#     Bf2a <- c()
#     parambest2a <- c()  
#     M2a <- Myx_vars
#     count = 1    
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, yv, xv, chy, M2a[count: (count+j-1)], zv)
#       bestm <- [[1]]
#       indexbestm <- [[2]] 
#       Bf2a[j] <- bestm
#       parambest2a[j] <- indexbestm
#       count = count + j;
#     }
#     
#     Bf2b <- c()
#     parambest2b <- c()  
#     M2b <- Myz_allvars
#     count = 1    
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chy, M2b[count: (count+j-1)], zv)
#       bestm <- [[1]]
#       indexbestm <- [[2]] 
#       Bf2b[j] <- bestm
#       parambest2b[j] <- indexbestm
#       count = count + j;
#     }
    
    M3 <- Mz_allvars
    count = 1   
    for (j in 1:nmodelterms)
    {
      tmp = polyfitbayes(indnr, xv, yv, chz, M3[count:(count+j-1)], zv)
      bestm <- tmp[[1]]
      indexbestm <- tmp[[2]]
      
      Bf3[j] <- bestm
      parambest3[j] <- indexbestm
      count = count + j;
    } 

    ## comment out lines 247-273 for comparing Bayes Factor of 
    ## two-variable and three-variable-models (see descriptions above)
#     Bf3a <- c()
#     parambest3a <- c() 
#     M3a <- Mzx_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chz, M3a[count: (count+j-1)], zv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf3a[j] <- bestm
#       parambest3a[j] <- indexbestm
#       count = count + j;
#     }
#  
#     Bf3b <- c()
#     parambest3b <- c() 
#     M3b <- Mzy_vars
#     count = 1  
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chz, M3b[count: (count+j-1)], zv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf3b[j] <- bestm
#       parambest3b[j] <- indexbestm
#       count = count + j;
#     }

    ## comment out lines 278,279,281,282,284,285 for comparing Bayes Factor of 
    ## two-variable and three-variable-models (see descriptions above)
    print(c("Bayes factors of the best models per number of modelterms for dx:", Bf1), quote=F)
#     print(Bf1a)
#     print(Bf1b)
    print(c("Bayes factors of the best models per number of modelterms for dy:", Bf2), quote=F)
#     print(Bf2a)
#     print(Bf2b)
    print(c("Bayes factors of the best models per number of modelterms for dz:", Bf3), quote=F)
#     print(Bf3a)
#     print(Bf3b)
  }

  ############################# indnr == 3 ends, indnr == 4 begins #################################
  
  if (indnr == 4)
  {
    nmodelterms = paramnr
    
    # calling selectterms function
    tmp = selectterms(indnr, paramnr, SEx, SEy, SEz, SEv) 
    Mx <- tmp[[1]]
    My <- tmp[[2]]
    Mz <- tmp[[3]]
    Mv <- tmp[[4]]
#     Mxy <- tmp[[5]]  ## the user might want to uncomment linwa 300-323
#     Mxz <- tmp[[7]]  ## if he/she wants to compare Bayes Factors of models with
#     Myx <- tmp[[6]]  ## two, three and four indicators 
#     Myz <- tmp[[8]]  ## if this is these and related lines are uncommented
#     Mzx <- tmp[[9]]  ## it is necessary to uncomment related parts in selectterms.R (follow the comment instructions in selectterms.R)
#     Mzy <- tmp[[10]]
#     Mxv <- tmp[[11]]
#     Myv <- tmp[[12]]
#     Mzv <- tmp[[13]]
#     Mvx <- tmp[[14]]
#     Mvy <- tmp[[15]]
#     Mvz <- tmp[[16]]
#     Mxyz <- tmp[[17]]
#     Mxyv <- tmp[[18]]
#     Mxzv <- tmp[[19]]
#     Myxz <- tmp[[20]]
#     Myxv <- tmp[[21]]
#     Myzv <- tmp[[22]] 
#     Mzxy <- tmp[[23]]
#     Mzxv <- tmp[[24]]
#     Mzyv <- tmp[[25]]
#     Mvxy <- tmp[[26]]
#     Mvxz <- tmp[[27]]
#     Mvyz <- tmp[[28]]
    
    Mx_allvars = c()
    My_allvars = c()
    Mz_allvars = c()
    Mv_allvars = c()
    ## uncomment lines 331-354 for comparing Bayes Factor of 
    ## two-variable-, three-variable- and four-variable-models (see descriptions above)
#     Mxy_vars = c()
#     Mxz_vars = c()
#     Myx_vars = c()
#     Myz_vars = c()
#     Mzx_vars = c()
#     Mzy_vars = c()
#     Mxv_vars = c()
#     Myv_vars = c()
#     Mzv_vars = c()
#     Mvx_vars = c()
#     Mvy_vars = c()
#     Mvz_vars = c()
#     Mxyz_vars = c() 
#     Mxyv_vars = c()
#     Mxzv_vars = c() 
#     Myxz_vars = c() 
#     Myxv_vars = c() 
#     Myzv_vars = c() 
#     Mzxy_vars = c() 
#     Mzxv_vars = c() 
#     Mzyv_vars = c() 
#     Mvxy_vars = c() 
#     Mvxz_vars = c() 
#     Mvyz_vars = c() 
    
    count = 1
    for (ii in 1:nmodelterms)
    {
      M <- combs(1:97, ii) 
      Mx_allvars <- c(Mx_allvars, M[Mx[ii],])
      My_allvars <- c(My_allvars, M[My[ii],])
      Mz_allvars <- c(Mz_allvars, M[Mz[ii],])
      Mv_allvars <- c(Mv_allvars, M[Mv[ii],])
      ## uncomment lines 365-388 for comparing Bayes Factor of 
      ## two-variable-, three-variable- and four-variable-models (see descriptions above)
#       Mxy_vars <- c(Mxy_vars, M[Mxy[ii],])
#       Mxz_vars <- c(Mxz_vars, M[Mxz[ii],])
#       Myx_vars <- c(Myx_vars, M[Myx[ii],])
#       Myz_vars <- c(Myz_vars, M[Myz[ii],])
#       Mzx_vars <- c(Mzx_vars, M[Mzx[ii],])
#       Mzy_vars <- c(Mzy_vars, M[Mzy[ii],])
#       Mxv_vars <- c(Mxv_vars, M[Mxv[ii],])
#       Myv_vars <- c(Myv_vars, M[Myv[ii],])
#       Mzv_vars <- c(Mzv_vars, M[Mzv[ii],])
#       Mvx_vars <- c(Mvx_vars, M[Mvx[ii],])
#       Mvy_vars <- c(Mvy_vars, M[Mvy[ii],])
#       Mvz_vars <- c(Mvz_vars, M[Mvz[ii],])
#       Mxyz_vars <- (Mxyz_vars, M[Mxyz[ii],]) 
#       Mxyv_vars <- (Mxyv_vars, M[Mxyv[ii],]) 
#       Mxzv_vars <- (Mxzv_vars, M[Mxzv[ii],])  
#       Myxz_vars <- (Myxz_vars, M[Myxz[ii],])  
#       Myxv_vars <- (Myxv_vars, M[Myxv[ii],])  
#       Myzv_vars <- (Myzv_vars, M[Myzv[ii],])  
#       Mzxy_vars <- (Mzxy_vars, M[Mzxy[ii],])  
#       Mzxv_vars <- (Mzxv_vars, M[Mzxv[ii],])  
#       Mzyv_vars <- (Mzyv_vars, M[Mzyv[ii],])  
#       Mvxy_vars <- (Mvxy_vars, M[Mvxy[ii],])  
#       Mvxz_vars <- (Mvxz_vars, M[Mvxz[ii],]) 
#       Mvyz_vars <- (Mvyz_vars, M[Mvyz[ii],])  
      count = count + ii
    }
    
    print(c("Selected model terms (dx):"), quote=F)
    print(Mx_allvars)
    print(c("Selected model terms (dy):"), quote=F)
    print(My_allvars)
    print(c("Selected model terms (dz):"), quote=F)
    print(Mz_allvars)
    print(c("Selected model terms (dv):"), quote=F)
    print(Mv_allvars)
    ## uncomment lines 398-421 for comparing Bayes Factor of 
    ## two-variable-, three-variable- and four-variable-models (see descriptions above)
#     print(Mxy_vars)
#     print(Mxz_vars)
#     print(Myx_vars)
#     print(Myz_vars)
#     print(Mzx_vars)
#     print(Mzy_vars)
#       print(Mxv_vars)
#       print(Myv_vars)
#       print(Mzv_vars)
#       print(Mvx_vars)
#       print(Mvy_vars)
#       print(Mvz_vars)
#       print(Mxyz_vars)
#       print(Mxyv_vars)
#       print(Mxzv_vars) 
#       print(Myxz_vars) 
#       print(Myxv_vars) 
#       print(Myzv_vars) 
#       print(Mzxy_vars) 
#       print(Mzxv_vars) 
#       print(Mzyv_vars) 
#       print(Mvxy_vars) 
#       print(Mvxz_vars) 
#       print(Mvyz_vars) 
    
    # creating empty vectors to be filled in the next step with Bayes factors 
    # and best parameters calling polyfitbayes function
    Bf1 <- c()
    parambest1 <- c()
    Bf2 <- c()
    parambest2 <- c()
    Bf3 <- c()
    parambest3 <- c()
    Bf4 <- c()
    parambest4 <- c()
    
    M1 <- Mx_allvars
    count = 1    
    for (j in 1:nmodelterms)
    {
      tmp = polyfitbayes(indnr, xv, yv, chx, M1[count:(count+j-1)], zv, vv)
      bestm <- tmp[[1]]
      indexbestm <- tmp[[2]]
      
      Bf1[j] <- bestm
      parambest1[j] <- indexbestm
      count = count + j;
    }  

    ## uncomment lines 450-534 for comparing Bayes Factor of 
    ## two-variable-, three-variable- and four-variable-models (see descriptions above)
#     Bf1a <- c()
#     parambest1a <- c() 
#     M1a <- Mxy_vars
#     count = 1  
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chx, M1a[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
    
#       Bf1a[j] <- bestm
#       parambest1a[j] <- indexbestm
#       count = count + j;
#     } 
#     
#     Bf1b <- c()
#     parambest1b <- c() 
#     M1b <- Mxz_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chx, M1b[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]   
#       Bf1b[j] <- bestm
#       parambest1b[j] <- indexbestm
#       count = count + j;
#     } 
#     
#     Bf1c <- c()
#     parambest1c <- c() 
#     M1c <- Mxv_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chx, M1c[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]] 
#       Bf1c[j] <- bestm
#       parambest1c[j] <- indexbestm
#       count = count + j;
#     }
# 
#     Bf1d <- c()
#     parambest1d <- c() 
#     M1d <- Mxyz_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chx, M1d[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]] 
#       Bf1d[j] <- bestm
#       parambest1d[j] <- indexbestm
#       count = count + j;
#     }
#
#     Bf1e <- c()
#     parambest1e <- c() 
#     M1e <- Mxyv_vars
#     count = 1
#     
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chx, M1e[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]] 
#       Bf1e[j] <- bestm
#       parambest1e[j] <- indexbestm
#       count = count + j;
#     }
#  
#     Bf1f <- c()
#     parambest1f <- c() 
#     M1f <- Mxzv_vars
#     count = 1 
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chx, M1f[count: (count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]] 
#       Bf1f[j] <- bestm
#       parambest1f[j] <- indexbestm
#       count = count + j;
#     }
    
    M2 <- My_allvars
    count = 1   
    for (j in 1:nmodelterms)
    {
      tmp = polyfitbayes(indnr, xv, yv, chy, M2[count:(count+j-1)], zv, vv)
      bestm <- tmp[[1]]
      indexbestm <- tmp[[2]] 
      
      Bf2[j] <- bestm
      parambest2[j] <- indexbestm
      count = count + j;
    } 

    ## uncomment lines 551-633 for comparing Bayes Factor of 
    ## two-variable-, three-variable- and four-variable-models (see descriptions above)
#     Bf2a <- c()
#     parambest2a <- c() 
#     M2a <- Myx_vars
#     count = 1 
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chy, M2a[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]] 
#       Bf2a[j] <- bestm
#       parambest2a[j] <- indexbestm
#       count = count + j;
#     }
# 
#     Bf2b <- c()
#     parambest2b <- c() 
#     M2b <- Myz_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chy, M2b[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]] 
#       Bf2b[j] <- bestm
#       parambest2b[j] <- indexbestm
#       count = count + j;
#     }
# 
#     Bf2c <- c()
#     parambest2c <- c() 
#     M2c <- Myv_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chy, M2c[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf2c[j] <- bestm
#       parambest2c[j] <- indexbestm
#       count = count + j;
#     }  
# 
#     Bf2d <- c()
#     parambest2d <- c() 
#     M2d <- Myxz_vars
#     count = 1  
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chy, M2d[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf2d[j] <- bestm
#       parambest2d[j] <- indexbestm
#       count = count + j;
#     } 
# 
#     Bf2e <- c()
#     parambest2e <- c() 
#     M2e <- Myxv_vars
#     count = 1
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chy, M2e[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf2e[j] <- bestm
#       parambest2e[j] <- indexbestm
#       count = count + j;
#     }
# 
#     Bf2f <- c()
#     parambest2f <- c() 
#     M2f <- Myzv_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chy, M2f[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf2f[j] <- bestm
#       parambest2f[j] <- indexbestm
#       count = count + j;
#     }
    
    M3 <- Mz_allvars
    count = 1    
    for (j in 1:nmodelterms)
    {
      tmp = polyfitbayes(indnr, xv, yv, chz, M3[count:(count+j-1)], zv, vv)
      bestm <- tmp[[1]]
      indexbestm <- tmp[[2]]
      
      Bf3[j] <- bestm
      parambest3[j] <- indexbestm
      count = count + j;
    } 

    ## uncomment lines 650-732 for comparing Bayes Factor of 
    ## two-variable-, three-variable- and four-variable-models (see descriptions above)
#     Bf3a <- c()
#     parambest3a <- c() 
#     M3a <- Mzx_vars
#     count = 1 
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chz, M3a[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf3a[j] <- bestm
#       parambest3a[j] <- indexbestm
#       count = count + j;
#     }
#
#     Bf3b <- c()
#     parambest3b <- c() 
#     M3b <- Mzy_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chz, M3b[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf3b[j] <- bestm
#       parambest3b[j] <- indexbestm
#       count = count + j;
#     }
#
#     Bf3c <- c()
#     parambest3c <- c() 
#     M3c <- Mzv_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chz, M3c[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf3c[j] <- bestm
#       parambest3c[j] <- indexbestm
#       count = count + j;
#     }
#
#     Bf3d <- c()
#     parambest3d <- c() 
#     M3d <- Mzxy_vars
#     count = 1  
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chz, M3d[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf3d[j] <- bestm
#       parambest3d[j] <- indexbestm
#       count = count + j;
#     }
#
#     Bf3e <- c()
#     parambest3e <- c() 
#     M3e <- Mzxv_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chz, M3e[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf3e[j] <- bestm
#       parambest3e[j] <- indexbestm
#       count = count + j;
#     }
# 
#     Bf3f <- c()
#     parambest3f <- c() 
#     M3f <- Mzyv_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chz, M3f[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf3f[j] <- bestm
#       parambest3f[j] <- indexbestm
#       count = count + j;
#     }
    
    M4 <- Mv_allvars
    count = 1   
    for (j in 1:nmodelterms)
    {
      tmp = polyfitbayes(indnr, xv, yv, chv, M4[count:(count+j-1)], zv, vv)
      bestm <- tmp[[1]]
      indexbestm <- tmp[[2]]
      
      Bf4[j] <- bestm
      parambest4[j] <- indexbestm
      count = count + j;
    }
    
    ## comment out lines 749-831 for comparing Bayes Factor of 
    ## two-variable-, three-variable- and four-variable-models (see descriptions above)
#     Bf4a <- c()
#     parambest4a <- c() 
#     M4a <- Mvx_vars
#     count = 1    
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chv, M4a[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf4a[j] <- bestm
#       parambest4a[j] <- indexbestm
#       count = count + j;
#     }
# 
#     Bf4b <- c()
#     parambest4b <- c() 
#     M4b <- Mvy_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chv, M4b[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf4b[j] <- bestm
#       parambest4b[j] <- indexbestm
#       count = count + j;
#     }
# 
#     Bf4c <- c()
#     parambest4c <- c() 
#     M4c <- Mvz_vars
#     count = 1    
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chv, M4c[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf4c[j] <- bestm
#       parambest4c[j] <- indexbestm
#       count = count + j;
#     }
# 
#     Bf4d <- c()
#     parambest4d <- c() 
#     M4d <- Mvxy_vars
#     count = 1    
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chv, M4d[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf4d[j] <- bestm
#       parambest4d[j] <- indexbestm
#       count = count + j;
#     }
# 
#     Bf4e <- c()
#     parambest4e <- c() 
#     M4e <- Mvxz_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chv, M4e[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf4e[j] <- bestm
#       parambest4e[j] <- indexbestm
#       count = count + j;
#     }
#
#     Bf4f <- c()
#     parambest4f <- c() 
#     M4f <- Mvyz_vars
#     count = 1   
#     for (j in 1:nmodelterms)
#     {
#       tmp = polyfitbayes(indnr, xv, yv, chv, M4f[count:(count+j-1)], zv, vv)
#       bestm <- [[1]]
#       indexbestm <- [[2]]
#       Bf4f[j] <- bestm
#       parambest4f[j] <- indexbestm
#       count = count + j;
#     }

    ## uncomment lines 773-778, 780-785, 787-792, 794-799 for comparing Bayes Factor of 
    ## two-variable-, three-variable- and four-variable-models (see descriptions above)
    print(c("Bayes factors of the best models per number of modelterms for dx:", Bf1), quote=F)
#     print(Bf1a)
#     print(Bf1b)
#     print(Bf1c)
#     print(Bf1d)
#     print(Bf1e)
#     print(Bf1f)
    print(c("Bayes factors of the best models per number of modelterms for dy:", Bf2), quote=F)
#     print(Bf2a)
#     print(Bf2b)
#     print(Bf2c)
#     print(Bf2d)
#     print(Bf2e)
#     print(Bf2f)
    print(c("Bayes factors of the best models per number of modelterms for dz:", Bf3), quote=F)
#     print(Bf3a)
#     print(Bf3b)
#     print(Bf3c)
#     print(Bf3d)
#     print(Bf3e)
#     print(Bf3f)
    print(c("Bayes factors of the best models per number of modelterms for dv:", Bf4), quote=F)
#     print(Bf4a)
#     print(Bf4b)
#     print(Bf4c)
#     print(Bf4d)
#     print(Bf4e)
#     print(Bf4f)
  }  
}