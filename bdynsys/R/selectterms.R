# selectterms selects/indexes terms which are part of the models to be fitted in the Bayesian fitting procedure
# it takes the SE argument from bayesfac, SE containts the sum square of errors of the fitted models

selectterms <- function(indnr, paramnr, SEx, SEy, SEz, SEv)
{ 
  if (indnr == 2)
  {
    nmodelterms = paramnr
    
    # defining whith polynomial terms contain which indicator as predictor
    Xterms <- rbind(c(2, 4, 10, 11, 14, 16)) 
    Yterms <- rbind(c(3, 5, 12, 13, 15, 17))
    XYterms <- rbind(c(6, 7, 8, 9))
    
    # creating empty vectors to be fillen in the next step with indexes of selected models
    Mx_all <- c()
    My_all <- c()
    
    wantedarray = 1:17   
    for (modterm in 1:nmodelterms)
    {
      comb <- combs(wantedarray, modterm)
      indexset = 1:choose(17, modterm)
      ind = which.min(SEx[modterm, indexset])
      Mx_all[modterm] <- indexset[ind]
    }
    
    wantedarray = 1:17   
    for (modterm in 1:nmodelterms)
    {
      comb <- combs(wantedarray, modterm)
      indexset = 1:choose(17, modterm)
      ind = which.min(SEy[modterm, indexset])
      My_all[modterm] <- indexset[ind]
    }
    return(list(Mx_all, My_all))
  }
  
  
  ############################# indnr == 2 ends, indnr == 3 begins #################################
  
  if (indnr == 3)
  {
    nmodelterms = paramnr
    
    # defining whith polynomial terms contain which indicator as predictor
    Xterms <- rbind(c(2, 5, 28, 29, 34, 37))
    Yterms <- rbind(c(3, 6, 30, 31, 35, 38))
    Zterms <- rbind(c(4, 7, 32, 33, 36, 39))
    XYterms <- rbind(c(8, 11, 14, 15))
    XZterms <- rbind(c(10, 13, 16, 17))
    YZterms <- rbind(c(9, 12, 18, 19))
    XYZterms <- rbind(c(20, 21, 22, 23, 24, 25, 26, 27)) 
    
    # creating empty vectors to be fillen in the next step with indexes of selected models
    Mx_all <- c()
    My_all <- c()
    Mz_all <- c()
    
    # X variable
    wantedarray = 1:39
    for(modterm in 1:nmodelterms)
    {
      comb <- combs(wantedarray, modterm)
      indexset = 1:choose(39, modterm)
      ind = which.min(SEx[modterm, indexset])
      Mx_all[modterm] <- indexset[ind]
    }

## the user should uncomment lines 76-106, if he/she wants to compare
## Bayes Factors of models with two and tree indicators in bayesfac.R
#     wantedarray = sort(c(1, Xterms, Yterms, XYterms)) 
#     Mxy <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 39)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEx[modterm, indexset])
#       Mxy[modterm] <- indexset[ind]      
#     } 
    
#     wantedarray = sort(c(1, Xterms, Zterms, XZterms)) 
#     Mxz <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 39)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEx[modterm, indexset])
#       Mxz[modterm] <- indexset[ind]      
#     }
    
    # Y varialbe
    wantedarray = 1:39
    for(modterm in 1:nmodelterms)
    {
      comb <- combs(wantedarray, modterm)
      indexset = 1:choose(39, modterm)
      ind = which.min(SEy[modterm, indexset])
      My_all[modterm] <- indexset[ind]
    }

## the user should comment out lines 120-150, if he/she wants to compare
## Bayes Factors of models with two and tree indicators in bayesfac.R
#     wantedarray = sort(c(1, Xterms, Yterms, XYterms))  
#     Myx <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 39)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEy[modterm, indexset])
#       Myx[modterm] <- indexset[ind]      
#     } 
#     
#     wantedarray = sort(c(1, Yterms, Zterms, YZterms))   
#     Myz <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 39)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEy[modterm, indexset])
#       Myz[modterm] <- indexset[ind]      
#     }
    
    # Z variable
    wantedarray = 1:39
    for(modterm in 1:nmodelterms)
    {
      comb <- combs(wantedarray, modterm)
      indexset = 1:choose(39, modterm)
      ind = which.min(SEz[modterm, indexset])
      Mz_all[modterm] <- indexset[ind]
    } 
    
## the user should  comment out lines 164-194, if he/she wants to compare
## Bayes Factors of models with two and tree indicators in bayesfac.R    
#     wantedarray = sort(c(1, Xterms, Zterms, XZterms))  
#     Mzx <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 39)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEz[modterm, indexset])
#       Mzx[modterm] <- indexset[ind]      
#     } 
#     
#     wantedarray = sort(c(1, Yterms, Zterms, YZterms)) 
#     Mzy <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 39)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEz[modterm, indexset])
#       Mzy[modterm] <- indexset[ind]      
#     }   
    return(list(Mx_all, My_all, Mz_all))
  }
  
  ############################# indnr == 3 ends, indnr == 4 begins #################################
  
  if (indnr == 4)
  {
    nmodelterms = paramnr
    
    # defining whith polynomial terms contain which indicator as predictor
    Xterms <- rbind(c(2, 6, 82, 83, 90, 94))
    Yterms <- rbind(c(3, 7, 84, 85, 91, 95))
    Zterms <- rbind(c(4, 8, 86, 87, 92, 96))
    Vterms <- rbind(c(5, 9, 88, 89, 93, 97))
    XYterms <- rbind(c(10, 16, 22, 23))
    XZterms <- rbind(c(12, 18, 24, 25))
    XVterms <- rbind(c(13, 19, 28, 29))
    YZterms <- rbind(c(11, 17, 26, 27))
    YVterms <- rbind(c(14, 20, 30, 31))
    ZVterms <- rbind(c(15, 21, 32, 33))
    XYZterms <- rbind(c(34, 35, 36, 46, 47, 48, 58, 62))
    XYVterms <- rbind(c(37, 40, 42, 49, 54, 56, 59, 63))
    XZVterms <- rbind(c(38, 41, 44, 51, 52, 55, 60, 64))
    YZVterms <- rbind(c(39, 43, 45, 50, 53, 57, 61, 65))
    XYZVterms <- rbind(c(66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81))
    
    # creating empty vectors to be fillen in the next step with indexes of selected models
    Mx_all <- c()
    My_all <- c()
    Mz_all <- c()
    Mv_all <- c()
    
    # for X variable 
    wantedarray = 1:97
    for(modterm in 1:nmodelterms)
    {
      comb <- combs(wantedarray, modterm)
      indexset = 1:choose(97, modterm)
      ind = which.min(SEx[modterm, indexset])
      Mx_all[modterm] <- indexset[ind]
    }
    
## the user should uncomment lines 239-333, if he/she wants to compare
## Bayes Factors of models with two, tree and four indicators in bayesfac.R    
#     wantedarray = sort(c(1, Xterms, Yterms, XYterms)) 
#     Mxy <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEx[modterm, indexset])
#       Mxy[modterm] <- indexset[ind]      
#     } 
#     
#     wantedarray = sort(c(1, Xterms, Zterms, XZterms)) 
#     Mxz <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEx[modterm, indexset])
#       Mxz[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Vterms, XVterms))   
#     Mxv <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEx[modterm, indexset])
#       Mxv[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Yterms, Zterm, XYterms, XZterms, YZterms, XYZterms))
#     Mxyz <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEx[modterm, indexset])
#       Mxyz[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Yterms, Vterm, XYterms, XVterms, YVterms, XYVterms)) 
#     Mxyv <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, findindexM(indexM, 97))
#       } 
#       ind = which.min(SEx[modterm, indexset])
#       Mxyv[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Vterms, Zterm, XVterms, XZterms, ZVterms, XZVterms)) 
#     Mxzv <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEx[modterm, indexset])
#       Mxzv[modterm] <- indexset[ind]      
#     }
    
    #for Y variable
    wantedarray = 1:97
    for(modterm in 1:nmodelterms)
    {
      comb <- combs(wantedarray, modterm)
      indexset = 1:choose(97, modterm)
      ind = which.min(SEy[modterm, indexset])
      My_all[modterm] <- indexset[ind]
    }
    
## the user should uncomment lines 347-441, if he/she wants to compare
## Bayes Factors of models with two, tree and four indicators in bayesfac.R      
#     wantedarray = sort(c(1, Xterms, Yterms, XYterms)) 
#     Myx <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEy[modterm, indexset])
#       Myx[modterm] <- indexset[ind]      
#     } 
#     
#     wantedarray = sort(c(1, Yterms, Zterms, YZterms)) 
#     Myz <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEy[modterm, indexset])
#       Myz[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Yterms, Vterms, YVterms))   
#     Myv <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEy[modterm, indexset])
#       Myv[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Yterms, Zterm, XYterms, XZterms, YZterms, XYZterms)) 
#     Myxz <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEy[modterm, indexset])
#       Myxz[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Yterms, Vterm, XYterms, XVterms, YVterms, XYVterms)) 
#     Myxv <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEy[modterm, indexset])
#       Myxv[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Yterms, Vterms, Zterm, YVterms, YZterms, ZVterms, YZVterms)) 
#     Myzv <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEy[modterm, indexset])
#       Myzv[modterm] <- indexset[ind]      
#     }
    
    # for Z variable
    wantedarray = 1:97
    for(modterm in 1:nmodelterms)
    {
      comb <- combs(wantedarray, modterm)
      indexset = 1:choose(97, modterm)
      ind = which.min(SEz[modterm, indexset])
      Mz_all[modterm] <- indexset[ind]
    }

## the user should uncomment lines 455-549, if he/she wants to compare
## Bayes Factors of models with two, tree and four indicators in bayesfac.R  
#     wantedarray = sort(c(1, Xterms, Zterms, XZterms)) 
#     Mzx <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEz[modterm, indexset])
#       Mzx[modterm] <- indexset[ind]      
#     } 
#     
#     wantedarray = sort(c(1, Yterms, Zterms, YZterms))
#     Mzy <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEz[modterm, indexset])
#       Mzy[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Zterms, Vterms, ZVterms)) 
#     Mzv <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEz[modterm, indexset])
#       Mzv[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Yterms, Zterm, XYterms, XZterms, YZterms, XYZterms))
#     Mzxy <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEz[modterm, indexset])
#       Mzxy[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Zterms, Vterm, XZterms, XVterms, ZVterms, XZVterms)) 
#     Mzxv <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEz[modterm, indexset])
#       Mzxv[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Yterms, Vterms, Zterm, YVterms, YZterms, ZVterms, YZVterms)) 
#     Mzyv <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEz[modterm, indexset])
#       Mzyv[modterm] <- indexset[ind]      
#     }
    
    # for V variable
    wantedarray = 1:97
    for(modterm in 1:nmodelterms)
    {
      comb <- combs(wantedarray, modterm)
      indexset = 1:choose(97, modterm)
      ind = which.min(SEv[modterm, indexset])
      Mv_all[modterm] <- indexset[ind]
    }

## the user should uncomment lines 563-657, if he/she wants to compare
## Bayes Factors of models with two, tree and four indicators in bayesfac.R  
#     wantedarray = sort(c(1, Xterms, Vterms, XVterms))  
#     Mvx <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEv[modterm, indexset])
#       Mvx[modterm] <- indexset[ind]      
#     } 
#     
#     wantedarray = sort(c(1, Vterms, Yterms, YVterms)) 
#     Mvy <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEv[modterm, indexset])
#       Mvy[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Zterms, Vterms, ZVterms)) 
#     Mvz <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEv[modterm, indexset])
#       Mvz[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Yterms, Vterm, XYterms, XVterms, YVterms, XYVterms))
#     Mvxy <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEv[modterm, indexset])
#       Mvxy[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Xterms, Zterms, Vterm, XZterms, XVterms, ZVterms, XZVterms)) 
#     Mvxz <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEv[modterm, indexset])
#       Mvxz[modterm] <- indexset[ind]      
#     }
#     
#     wantedarray = sort(c(1, Yterms, Vterms, Zterm, YVterms, YZterms, ZVterms, YZVterms))
#     Mvyz <- c()
#     for (modterm = 1:nmodelterms)
#     {
#       comb <- combs(wantedarray, modterm)
#       indexset = c()
#       for (t in 1:nrow(comb))
#       {
#         indexM <- comb[t, 1:modterm]
#         m <- findindexM(indexM, 97)
#         indexset <- c(indexset, m)
#       } 
#       ind = which.min(SEv[modterm, indexset])
#       Mvyz[modterm] <- indexset[ind]      
#     } 
    return(list(Mx_all, My_all, Mz_all, Mv_all))
  }
  
  ## the user should use the below return command, if he/she wants to compare
  ## Bayes Factors of models with two and three indicators in bayesfac.R  
  #   return(list(Mx_all, My_all, Mz_all, Mxy, Myx, Mxz, Myz, Mzx, Mzy)) 
  
  ## the user should use the below return command if he/she wants to compare
  ## Bayes Factors of models with two, tree and four indicators in bayesfac.R  
#   return(list(Mx_all, My_all, Mz_all, Mv_all, Mxy, Myx, Mxz, Myz, Mzx, Mzy,  
#               Mxv, Myv, Mzv, Mvx, Mvy, Mvz, Mxyz, Mxyv, Mxzv, Myxz, Myxv, Myzv, 
#               Mzxy, Mzxv, Mzyv, Mvxy, Mvxz, Mvyz))  
}