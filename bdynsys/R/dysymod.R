# Main script for the modeling procedure

dysymod <- function(indnr, paramnr, var1, var2, chVar1, chVar2, mvar1, mvar2, var3, chVar3, mvar3, var4, chVar4, mvar4)
{  
  if (indnr == 2)
  {
    # defining how many polynomial terms to be included, recommended 5-6
    # and how many models are to be returned per number of term, per default best 3 models for each number of terms
    nterms = 17;
    nmodelterms = paramnr;
    nmodels = 3; 
    
    # definition of polynomial terms
    term = vector('list', nterms)
    term[[1]] = ''
    term[[2]] = ('/x') 
    term[[3]] = ('/y')
    term[[4]] = ('x')
    term[[5]] = ('y')
    term[[6]] = ('/(x*y)')
    term[[7]] = ('x/y')
    term[[8]] = ('y/x')
    term[[9]] = ('x*y')
    term[[10]] = ('x^2')
    term[[11]] = ('/x^2')
    term[[12]] = ('y^2')
    term[[13]] = ('/y^2')
    term[[14]] = ('x^3')
    term[[15]] = ('y^3')
    term[[16]] = ('/x^3')
    term[[17]] = ('/y^3')
    
    # scaling terms with means
    scaling = vector('list', nterms)
    scaling[[1]] = 1
    scaling[[2]] = mvar1
    scaling[[3]] = mvar2
    scaling[[4]] = 1/mvar1
    scaling[[5]] = 1/mvar2
    scaling[[6]] = mvar1*mvar2
    scaling[[7]] = mvar2/mvar1
    scaling[[8]] = mvar1/mvar2
    scaling[[9]] = 1/(mvar1*mvar2)
    scaling[[10]] = 1/(mvar1*mvar1)
    scaling[[11]] = mvar1*mvar1
    scaling[[12]] = 1/(mvar2*mvar2)
    scaling[[13]] = mvar2*mvar2
    scaling[[14]] = 1/(mvar1*mvar1*mvar1)
    scaling[[15]] = 1/(mvar2*mvar2*mvar2)
    scaling[[16]] = mvar1*mvar1*mvar1
    scaling[[17]] = mvar2*mvar2*mvar2
    scaling <- unlist(scaling)
    
    # creating array were the selected models will be stored
    selmod1 <- array(0, dim=c(nmodelterms, nmodels, nmodelterms))
    selmod2 <- array(0, dim=c(nmodelterms, nmodels, nmodelterms))
    
    # Computing the mean square errors for all models and storing them in
    # SEtesty and SEtestx, calls the bestfit function
    SEtestx <- (bestfitmod(indnr, paramnr, var1, var2, chVar1))
    SEtesty <- (bestfitmod(indnr, paramnr, var1, var2, chVar2))
    # please uncomment if the user would like to check the SEtestx/y visually
    # save(SEtestx, SEtesty, file = "SE.RData")
    
    # Sorting the models according to their mean square errors and choosing 
    # the best three models for each number of modelterms
    for (ii in 1:nmodelterms)
    {
      # creating all possible model combinations
      M <- combs(1:nterms, ii)   
      # sorting SEtestx/y and getting the indexes of the sorted SEtestx/y values
      idx1 = order(SEtestx[ii,], na.last = TRUE, decreasing = FALSE)    
      idx2 = order(SEtesty[ii,], na.last = TRUE, decreasing = FALSE) 
      
      for (jj in 1:nmodels) 
      {        
        selmod1[ii, jj, 1:ii] <- M[idx1[jj],]
        selmod2[ii, jj, 1:ii] <- M[idx2[jj],]
      }    
    }
    
    # displays the three best models for each nmodelterms(paramnr) for x/var1 (indicator1)
    # calling polyfitreg function
    for (i in 1:nmodelterms)
    {
      for (j in 1:nmodels)
      {
        modsel = matrix(0, i, 1)
        modsel[1:i] <- selmod1[i, j, 1:i]

        tmp = polyfitreg(indnr, var1, var2, chVar1, modsel) 
        B <- tmp[[1]]
        meansqerr <- tmp[[2]]
        
        rmsx <- meansqerr*(mvar1*mvar1) 
        
        # Betas
        B_out_x <- array(list(), dim=c(nmodelterms, nmodels)) 
        B_out_x[[i, j]] <- matrix(0, nterms, 1)
        # uncommend if the user would like to check the unscaled Betas
        #print(B)
        B_out_x[[i, j]][modsel] <- B 
        B_out_x[[i, j]] <- B_out_x[[i, j]]*(mvar1*scaling)
        B <- B*mvar1*scaling[modsel] 
        
        # printing results
        out_text_x <- array(list(), dim=c(i, j))
        out_text_x[[i, j]] <- cbind('dx', '=')
        
        for (k in 1:i)
          out_text_x[[i, j]] <- cbind(out_text_x[[i, j]], '+', format(B[k], digits = 2), term[[modsel[k]]])    

        print(c('Model #', as.character(j), 'of size', as.character(i), 'is:'), quote=FALSE) 
        write.table(format(out_text_x[[i, j]], justify="left", width=1), row.names=FALSE, col.names=FALSE, quote=FALSE)
      } 
      print('--------------------------------------------') 
    }  
    #   # displays the three best models for each nmodelterms for y/var2 (indicator2)
    #   # calling polyfitreg function
    for (i in 1:nmodelterms)
    {
      for (j in 1:nmodels)
      {
        modsel = matrix(0, i, 1)
        modsel[1:i] <- selmod2[i, j, 1:i]

        tmp = polyfitreg(indnr, var1, var2, chVar2, modsel) 
        B <- tmp[[1]]
        meansqerr <- tmp[[2]]
        
        rmsy = meansqerr*(mvar2*mvar2) 
        
        # Betas
        B_out_y <- array(list(), dim=c(nmodelterms, nmodels)) 
        B_out_y[[i, j]] <- matrix(0, nterms, 1)
        #print(B)
        B_out_y[[i, j]][modsel] <- B 
        B_out_y[[i, j]] <- B_out_y[[i, j]]*(mvar2*scaling)
        B <- B*mvar2*scaling[modsel] 
        
        # printing results
        out_text_y <- array(list(), dim=c(i, j))
        out_text_y[[i, j]] <- cbind('dy', '=')
        
        for (k in 1:i)
          out_text_y[[i, j]] <- cbind(out_text_y[[i, j]], '+', format(B[k], digits = 2), term[[modsel[k]]])
        
        print(c('Model #', as.character(j), 'of size', as.character(i), 'is:'), quote=FALSE)
        write.table(format(out_text_y[[i, j]], justify="left", width=1), row.names=FALSE, col.names=FALSE, quote=FALSE)
      } 
      print('--------------------------------------------')
    } 
    return(list(B_out_x=B_out_x, B_out_y=B_out_y, out_text_x=out_text_x, out_text_y=out_text_y, 
                scaling=scaling, SEtestx=SEtestx, SEtesty=SEtesty))
  }
  
  ############################# indnr == 2 ends, indnr == 3 begins #################################
  
  if (indnr == 3)
  {
    # defining how many polynomial terms to be included, recommended 5-6
    # and how many models are to be returned per number of term, per default best 3 models for each number of terms
    nterms = 39;
    nmodelterms = paramnr;
    nmodels = 3;
    
    # definition of polynomial terms
    term = vector('list', nterms)
    term[[1]] = ''
    term[[2]] = ('/x') 
    term[[3]] = ('/y')
    term[[4]] = ('/z')
    term[[5]] = ('x')
    term[[6]] = ('y')
    term[[7]] = ('z')
    term[[8]] = ('/(x*y)')
    term[[9]] = ('/(y*z)')
    term[[10]] = ('/(x*z)')
    term[[11]] = ('x*y')
    term[[12]] = ('y*z')
    term[[13]] = ('x*z')
    term[[14]] = ('x/y')
    term[[15]] = ('y/x')
    term[[16]] = ('x/z')
    term[[17]] = ('z/x')
    term[[18]] = ('y/z')
    term[[19]] = ('z/y')
    term[[20]] = ('x/(y*z)')
    term[[21]] = ('y/(x*z)')
    term[[22]] = ('z/(x*y)')
    term[[23]] = ('(x*y)/z')
    term[[24]] = ('(y*z)/x')
    term[[25]] = ('(z*x)/y')
    term[[26]] = ('x*y*z')
    term[[27]] = ('1/(x*y*z)')
    term[[28]] = ('x^2')
    term[[29]] = ('/x^2')
    term[[30]] = ('y^2')
    term[[31]] = ('/y^2')
    term[[32]] = ('z^2')
    term[[33]] = ('/z^2')
    term[[34]] = ('x^3')
    term[[35]] = ('y^3')
    term[[36]] = ('z^3')
    term[[37]] = ('/x^3')
    term[[38]] = ('/y^3')
    term[[39]] = ('/z^3')
    
    # scaling terms with means
    scaling = vector('list', nterms)
    scaling[[1]] = 1
    scaling[[2]] = mvar1
    scaling[[3]] = mvar2
    scaling[[4]] = mvar3
    scaling[[5]] = 1/mvar1
    scaling[[6]] = 1/mvar2
    scaling[[7]] = 1/mvar3
    scaling[[8]] = mvar1*mvar2
    scaling[[9]] = mvar2*mvar3
    scaling[[10]] = mvar1*mvar3
    scaling[[11]] = 1/(mvar1*mvar2)
    scaling[[12]] = 1/(mvar2*mvar3)
    scaling[[13]] = 1/(mvar1*mvar3)
    scaling[[14]] = mvar2/mvar1
    scaling[[15]] = mvar1/mvar2
    scaling[[16]] = mvar3/mvar1
    scaling[[17]] = mvar1/mvar3
    scaling[[18]] = mvar3/mvar2
    scaling[[19]] = mvar2/mvar3
    scaling[[20]] = (mvar2*mvar3)/mvar1
    scaling[[21]] = (mvar1*mvar3)/mvar2
    scaling[[22]] = (mvar1*mvar2)/mvar3
    scaling[[23]] = mvar3/(mvar1*mvar2)
    scaling[[24]] = mvar1/(mvar2*mvar3)
    scaling[[25]] = mvar2/(mvar1*mvar3)
    scaling[[26]] = 1/(mvar1*mvar2*mvar3)
    scaling[[27]] = mvar1*mvar2*mvar3
    scaling[[28]] = mvar1^(-2)
    scaling[[29]] = mvar1^2
    scaling[[30]] = mvar2^(-2)
    scaling[[31]] = mvar2^2
    scaling[[32]] = mvar3^(-2)
    scaling[[33]] = mvar3^2
    scaling[[34]] = 1/(mvar1*mvar1*mvar1)
    scaling[[35]] = 1/(mvar2*mvar2*mvar2)
    scaling[[36]] = 1/(mvar3*mvar3*mvar3)
    scaling[[37]] = mvar1*mvar1*mvar1
    scaling[[38]] = mvar2*mvar2*mvar2
    scaling[[39]] = mvar3*mvar3*mvar3
    scaling <- unlist(scaling)
    
    # creating arraya were the selected models will be stored
    selmod1 <- array(0, dim=c(nmodelterms, nmodels, nmodelterms))
    selmod2 <- array(0, dim=c(nmodelterms, nmodels, nmodelterms))
    selmod3 <- array(0, dim=c(nmodelterms, nmodels, nmodelterms))
    
    # Computing the mean square errors for all models and storing them in
    # SEtestx, SEtesty, SEtestz, calls the bestfit function
    SEtestx <- (bestfitmod(indnr, paramnr, var1, var2, chVar1, var3))
    SEtesty <- (bestfitmod(indnr, paramnr, var1, var2, chVar2, var3))
    SEtestz <- (bestfitmod(indnr, paramnr, var1, var2, chVar3, var3))
    # please uncomment if the user would like to check the SEtestx/y/z visually
    # save(SEtestx, SEtesty, SEtestz, file = "SE.RData")
    
    # Sorting the models according to their mean square errors and choosing 
    # the best three models for each number of modelterms
    for (ii in 1:nmodelterms)
    {
      # creating all possible model combinations
      M <- combs(1:nterms, ii)   
      # sorting SEtestx/y/z and getting the indexes of the sorted SEtestx/y/z values
      idx1 = order(SEtestx[ii,], na.last = TRUE, decreasing = FALSE)    
      idx2 = order(SEtesty[ii,], na.last = TRUE, decreasing = FALSE)
      idx3 = order(SEtestz[ii,], na.last = TRUE, decreasing = FALSE)
      
      for (jj in 1:nmodels) 
      {        
        selmod1[ii, jj, 1:ii] <- M[idx1[jj],]
        selmod2[ii, jj, 1:ii] <- M[idx2[jj],]
        selmod3[ii, jj, 1:ii] <- M[idx3[jj],]
      }    
    }
    
    # displays the three best models for each nmodelterms(paramnr) for x/var1 (indicator1)
    # calling polyfitreg function
    for (i in 1:nmodelterms)
    {
      for (j in 1:nmodels)
      {
        modsel = matrix(0, i, 1)
        modsel[1:i] <- selmod1[i, j, 1:i]
        
        tmp = polyfitreg(indnr, var1, var2, chVar1, modsel, var3) 
        B <- tmp[[1]]
        meansqerr <- tmp[[2]]
        
        rmsx <- meansqerr*(mvar1*mvar1) 
        
        # Betas
        B_out_x <- array(list(), dim=c(nmodelterms, nmodels)) 
        B_out_x[[i, j]] <- matrix(0, nterms, 1)
        # uncommend if the user would like to check the unscaled Betas
        #print(B)
        B_out_x[[i, j]][modsel] <- B 
        B_out_x[[i, j]] <- B_out_x[[i, j]]*(mvar1*scaling)
        B <- B*mvar1*scaling[modsel] 
        
        # printing results
        out_text_x <- array(list(), dim=c(i, j))
        out_text_x[[i, j]] <- cbind('dx', '=')
        
        for (k in 1:i)
          out_text_x[[i, j]] <- cbind(out_text_x[[i, j]], '+', format(B[k], digits = 2), term[[modsel[k]]])    
        
        print(c('Model #', as.character(j), 'of size', as.character(i), 'is:'), quote=FALSE) 
        write.table(format(out_text_x[[i, j]], justify="left", width=1), row.names=FALSE, col.names=FALSE, quote=FALSE)
      } 
      print('--------------------------------------------') 
    }  
    #   # displays the three best models for each nmodelterms for y/var2 (indicator2)
    #   # calling polyfitreg function
    for (i in 1:nmodelterms)
    {
      for (j in 1:nmodels)
      {
        modsel = matrix(0, i, 1)
        modsel[1:i] <- selmod2[i, j, 1:i]
        
        tmp = polyfitreg(indnr, var1, var2, chVar2, modsel, var3) 
        B <- tmp[[1]]
        meansqerr <- tmp[[2]]
        
        rmsy = meansqerr*(mvar2*mvar2) 
        
        # Betas
        B_out_y <- array(list(), dim=c(nmodelterms, nmodels)) 
        B_out_y[[i, j]] <- matrix(0, nterms, 1)
        #print(B)
        B_out_y[[i, j]][modsel] <- B 
        B_out_y[[i, j]] <- B_out_y[[i, j]]*(mvar2*scaling)
        B <- B*mvar2*scaling[modsel] 
        
        # printing results
        out_text_y <- array(list(), dim=c(i, j))
        out_text_y[[i, j]] <- cbind('dy', '=')
        
        for (k in 1:i)
          out_text_y[[i, j]] <- cbind(out_text_y[[i, j]], '+', format(B[k], digits = 2), term[[modsel[k]]])
        
        print(c('Model #', as.character(j), 'of size', as.character(i), 'is:'), quote=FALSE)
        write.table(format(out_text_y[[i, j]], justify="left", width=1), row.names=FALSE, col.names=FALSE, quote=FALSE)
      } 
      print('--------------------------------------------')
    } 
    #   # displays the three best models for each nmodelterms for z/var3 (indicator3)
    #   # calling polyfitreg function
    for (i in 1:nmodelterms)
    {
      for (j in 1:nmodels)
      {
        modsel = matrix(0, i, 1)
        modsel[1:i] <- selmod3[i, j, 1:i]
        
        tmp = polyfitreg(indnr, var1, var2, chVar3, modsel, var3) 
        B <- tmp[[1]]
        meansqerr <- tmp[[2]]
        
        rmsz = meansqerr*(mvar3*mvar3) 
        
        # Betas
        B_out_z <- array(list(), dim=c(nmodelterms, nmodels)) 
        B_out_z[[i, j]] <- matrix(0, nterms, 1)
        #print(B)
        B_out_z[[i, j]][modsel] <- B 
        B_out_z[[i, j]] <- B_out_z[[i, j]]*(mvar3*scaling)
        B <- B*mvar3*scaling[modsel] 
        
        # printing results
        out_text_z <- array(list(), dim=c(i, j))
        out_text_z[[i, j]] <- cbind('dz', '=')
        
        for (k in 1:i)
          out_text_z[[i, j]] <- cbind(out_text_z[[i, j]], '+', format(B[k], digits = 2), term[[modsel[k]]])
        
        print(c('Model #', as.character(j), 'of size', as.character(i), 'is:'), quote=FALSE)
        write.table(format(out_text_z[[i, j]], justify="left", width=1), row.names=FALSE, col.names=FALSE, quote=FALSE)
      } 
      print('--------------------------------------------')
    }
    return(list(B_out_x=B_out_x, B_out_y=B_out_y, B_out_z=B_out_z, out_text_x=out_text_x, out_text_y=out_text_y, 
                out_text_z=out_text_z, scaling=scaling, SEtestx=SEtestx, SEtesty=SEtesty, SEtestz=SEtestz))
  }
  
  ############################# indnr == 3 ends, indnr == 4 begins #################################
  
  if (indnr == 4)
  {
    # defining how many polynomial terms to be included, recommended 4
    # and how many models are to be returned per number of term, per default best 3 models for each number of terms
    nterms = 97;
    nmodelterms = paramnr;
    nmodels = 3;
    
    # definition of polynomial terms
    term = vector('list', nterms)
    term[[1]] = ''
    term[[2]] = ('/x') 
    term[[3]] = ('/y')
    term[[4]] = ('/z')
    term[[5]] = ('/v')
    term[[6]] = ('x')
    term[[7]] = ('y')
    term[[8]] = ('z')
    term[[9]] = ('v')
    term[[10]] = ('/(x*y)')
    term[[11]] = ('/(y*z)')
    term[[12]] = ('/(x*z)')
    term[[13]] = ('/(x*v)')
    term[[14]] = ('/(y*v)')
    term[[15]] = ('/(z*v)')
    term[[16]] = ('x*y')
    term[[17]] = ('y*z')
    term[[18]] = ('x*z')
    term[[19]] = ('x*v')
    term[[20]] = ('y*v')
    term[[21]] = ('z*v')
    term[[22]] = ('x/y')
    term[[23]] = ('y/x')
    term[[24]] = ('x/z')
    term[[25]] = ('z/x')
    term[[26]] = ('y/z')
    term[[27]] = ('z/y')
    term[[28]] = ('x/v')
    term[[29]] = ('v/x')
    term[[30]] = ('y/v')
    term[[31]] = ('v/y')
    term[[32]] = ('z/v')
    term[[33]] = ('v/z')
    term[[34]] = ('x/(y*z)')
    term[[35]] = ('y/(x*z)')
    term[[36]] = ('z/(x*y)')
    term[[37]] = ('v/(x*y)')
    term[[38]] = ('v/(x*z)')
    term[[39]] = ('v/(y*z)')
    term[[40]] = ('x/(y*v)')
    term[[41]] = ('x/(z*v)')
    term[[42]] = ('y/(x*v)')
    term[[43]] = ('y/(z*v)')
    term[[44]] = ('z/(x*v)')
    term[[45]] = ('z/(y*v)')
    term[[46]] = ('(x*y)/z')
    term[[47]] = ('(y*z)/x')
    term[[48]] = ('(z*x)/y')
    term[[49]] = ('(x*y)/v')
    term[[50]] = ('(y*z)/v')
    term[[51]] = ('(z*x)/v')
    term[[52]] = ('(x*v)/z')
    term[[53]] = ('(y*v)/z')
    term[[54]] = ('(y*v)/x')
    term[[55]] = ('(z*v)/x')
    term[[56]] = ('(v*x)/y')
    term[[57]] = ('(v*z)/y')
    term[[58]] = ('x*y*z')
    term[[59]] = ('x*y*v')
    term[[60]] = ('x*v*z')
    term[[61]] = ('v*y*z')
    term[[62]] = ('1/(x*y*z)')
    term[[63]] = ('1/(x*y*v)')
    term[[64]] = ('1/(x*v*z)')
    term[[65]] = ('1/(v*y*z)')
    term[[66]] = ('x/(v*y*z)')
    term[[67]] = ('y/(x*v*z)')
    term[[68]] = ('z/(x*y*v)')
    term[[69]] = ('v/(x*y*z)')
    term[[70]] = ('(x*y*z)/v')
    term[[71]] = ('(x*y*v)/z')
    term[[72]] = ('(x*v*z)/y')
    term[[73]] = ('(v*y*z)/x')
    term[[74]] = ('(x*y)/(v*z)')
    term[[75]] = ('(x*z)/(v*y)')
    term[[76]] = ('(x*v)/(y*z)')
    term[[77]] = ('(y*z)/(v*x)')
    term[[78]] = ('(y*v)/(z*x)')
    term[[79]] = ('(z*v)/(x*y)')
    term[[80]] = ('1/(x*y*z*v)')
    term[[81]] = ('x*y*z*v')
    term[[82]] = ('x^2')
    term[[83]] = ('/x^2')
    term[[84]] = ('y^2')
    term[[85]] = ('/y^2')
    term[[86]] = ('z^2')
    term[[87]] = ('/z^2')
    term[[88]] = ('v^2')
    term[[89]] = ('/v^2')
    term[[90]] = ('x^3')
    term[[91]] = ('y^3')
    term[[92]] = ('z^3')
    term[[93]] = ('v^3')
    term[[94]] = ('/x^3')
    term[[95]] = ('/y^3')
    term[[96]] = ('/z^3')
    term[[97]] = ('/v^3')
    
    # scaling terms with means
    scaling = vector('list', nterms)
    scaling[[1]] = 1
    scaling[[2]] = mvar1
    scaling[[3]] = mvar2
    scaling[[4]] = mvar3
    scaling[[5]] = mvar4
    scaling[[6]] = 1/mvar1
    scaling[[7]] = 1/mvar2
    scaling[[8]] = 1/mvar3
    scaling[[9]] = 1/mvar4
    scaling[[10]] = mvar1*mvar2
    scaling[[11]] = mvar2*mvar3
    scaling[[12]] = mvar1*mvar3
    scaling[[13]] = mvar1*mvar4
    scaling[[14]] = mvar2*mvar4
    scaling[[15]] = mvar3*mvar4
    scaling[[16]] = 1/(mvar1*mvar2)
    scaling[[17]] = 1/(mvar2*mvar3)
    scaling[[18]] = 1/(mvar1*mvar3)
    scaling[[19]] = 1/(mvar1*mvar4)
    scaling[[20]] = 1/(mvar2*mvar4)
    scaling[[21]] = 1/(mvar3*mvar4)
    scaling[[22]] = mvar2/mvar1
    scaling[[23]] = mvar1/mvar2
    scaling[[24]] = mvar3/mvar1
    scaling[[25]] = mvar1/mvar2
    scaling[[26]] = mvar3/mvar2
    scaling[[27]] = mvar2/mvar3
    scaling[[28]] = mvar4/mvar1
    scaling[[29]] = mvar1/mvar4
    scaling[[30]] = mvar4/mvar2
    scaling[[31]] = mvar2/mvar4
    scaling[[32]] = mvar4/mvar3
    scaling[[33]] = mvar3/mvar4
    scaling[[34]] = (mvar2*mvar3)/mvar1
    scaling[[35]] = (mvar1*mvar3)/mvar2
    scaling[[36]] = (mvar1*mvar2)/mvar3
    scaling[[37]] = (mvar1*mvar2)/mvar4
    scaling[[38]] = (mvar1*mvar3)/mvar4
    scaling[[39]] = (mvar2*mvar3)/mvar4
    scaling[[40]] = (mvar2*mvar4)/mvar1
    scaling[[41]] = (mvar3*mvar4)/mvar1
    scaling[[42]] = (mvar1*mvar4)/mvar2
    scaling[[43]] = (mvar3*mvar4)/mvar2
    scaling[[44]] = (mvar1*mvar4)/mvar3
    scaling[[45]] = (mvar2*mvar4)/mvar3
    scaling[[46]] = mvar3/(mvar1*mvar2)
    scaling[[47]] = mvar1/(mvar2*mvar3)
    scaling[[48]] = mvar2/(mvar1*mvar3)
    scaling[[49]] = mvar4/(mvar1*mvar2)
    scaling[[50]] = mvar4/(mvar2*mvar3)
    scaling[[51]] = mvar4/(mvar1*mvar3)
    scaling[[52]] = mvar3/(mvar1*mvar4)
    scaling[[53]] = mvar3/(mvar2*mvar4)
    scaling[[54]] = mvar1/(mvar2*mvar4)
    scaling[[55]] = mvar1/(mvar3*mvar4)
    scaling[[56]] = mvar2/(mvar1*mvar4)
    scaling[[57]] = mvar2/(mvar3*mvar4)
    scaling[[58]] = 1/(mvar1*mvar2*mvar3)
    scaling[[59]] = 1/(mvar1*mvar2*mvar4)
    scaling[[60]] = 1/(mvar1*mvar4*mvar3)
    scaling[[61]] = 1/(mvar4*mvar2*mvar3)
    scaling[[62]] = mvar1*mvar2*mvar3
    scaling[[63]] = mvar1*mvar2*mvar4
    scaling[[64]] = mvar1*mvar4*mvar3
    scaling[[65]] = mvar4*mvar2*mvar3
    scaling[[66]] = (mvar4*mvar2*mvar3)/mvar1
    scaling[[67]] = (mvar1*mvar4*mvar3)/mvar2
    scaling[[68]] = (mvar1*mvar2*mvar4)/mvar3
    scaling[[69]] = (mvar1*mvar2*mvar3)/mvar4
    scaling[[70]] = mvar4/(mvar1*mvar2*mvar3)
    scaling[[71]] = mvar3/(mvar1*mvar2*mvar4)
    scaling[[72]] = mvar2/(mvar1*mvar4*mvar3)
    scaling[[73]] = mvar1/(mvar4*mvar2*mvar3)
    scaling[[74]] = (mvar4*mvar3)/(mvar1*mvar2)
    scaling[[75]] = (mvar4*mvar2)/(mvar1*mvar3)
    scaling[[76]] = (mvar2*mvar3)/(mvar1*mvar4)
    scaling[[77]] = (mvar4*mvar1)/(mvar2*mvar3)
    scaling[[78]] = (mvar3*mvar1)/(mvar2*mvar4)
    scaling[[79]] = (mvar1*mvar2)/(mvar3*mvar4)
    scaling[[80]] = mvar1*mvar2*mvar3*mvar4
    scaling[[81]] = 1/(mvar1*mvar2*mvar3*mvar4)
    scaling[[82]] = mvar1^(-2)
    scaling[[83]] = mvar1^2
    scaling[[84]] = mvar2^(-2)
    scaling[[85]] = mvar2^2
    scaling[[86]] = mvar3^(-2)
    scaling[[87]] = mvar3^2
    scaling[[88]] = mvar4^(-2)
    scaling[[89]] = mvar4^2
    scaling[[90]] = 1/(mvar1*mvar1*mvar1)
    scaling[[91]] = 1/(mvar2*mvar2*mvar2)
    scaling[[92]] = 1/(mvar3*mvar3*mvar3)
    scaling[[93]] = 1/(mvar4*mvar4*mvar4)
    scaling[[94]] = mvar1*mvar1*mvar1
    scaling[[95]] = mvar2*mvar2*mvar2
    scaling[[96]] = mvar3*mvar3*mvar3
    scaling[[97]] = mvar4*mvar4*mvar4
    scaling <- unlist(scaling)
    
    # creating arraya were the selected models will be stored
    selmod1 <- array(0, dim=c(nmodelterms, nmodels, nmodelterms))
    selmod2 <- array(0, dim=c(nmodelterms, nmodels, nmodelterms))
    selmod3 <- array(0, dim=c(nmodelterms, nmodels, nmodelterms))
    selmod4 <- array(0, dim=c(nmodelterms, nmodels, nmodelterms))
    
    # Computing the mean square errors for all models and storing them in
    # SEtestx, SEtesty, SEtestz, SEtestv calls the bestfit function
    SEtestx <- (bestfitmod(indnr, paramnr, var1, var2, chVar1, var3, var4))
    SEtesty <- (bestfitmod(indnr, paramnr, var1, var2, chVar2, var3, var4))
    SEtestz <- (bestfitmod(indnr, paramnr, var1, var2, chVar3, var3, var4))
    SEtestv <- (bestfitmod(indnr, paramnr, var1, var2, chVar4, var3, var4))
    # please uncomment if the user would like to check the SEtestx/y/z/v visually
    # save(SEtestx, SEtesty, SEtestz, SEtestv, file = "SE.RData")
    
    # Sorting the models according to their mean square errors and choosing 
    # the best three models for each number of modelterms
    for (ii in 1:nmodelterms)
    {
      # creating all possible model combinations
      M <- combs(1:nterms, ii)   
      # sorting SEtestx/y/z and getting the indexes of the sorted SEtestx/y/z values
      idx1 = order(SEtestx[ii,], na.last = TRUE, decreasing = FALSE)    
      idx2 = order(SEtesty[ii,], na.last = TRUE, decreasing = FALSE)
      idx3 = order(SEtestz[ii,], na.last = TRUE, decreasing = FALSE)
      idx4 = order(SEtestv[ii,], na.last = TRUE, decreasing = FALSE)
      
      for (jj in 1:nmodels) 
      {        
        selmod1[ii, jj, 1:ii] <- M[idx1[jj],]
        selmod2[ii, jj, 1:ii] <- M[idx2[jj],]
        selmod3[ii, jj, 1:ii] <- M[idx3[jj],]
        selmod4[ii, jj, 1:ii] <- M[idx4[jj],]
      }    
    }
    
    # displays the three best models for each nmodelterms(paramnr) for x/var1 (indicator1)
    # calling polyfitreg function
    for (i in 1:nmodelterms)
    {
      for (j in 1:nmodels)
      {
        modsel = matrix(0, i, 1)
        modsel[1:i] <- selmod1[i, j, 1:i]
        
        tmp = polyfitreg(indnr, var1, var2, chVar1, modsel, var3, var4) 
        B <- tmp[[1]]
        meansqerr <- tmp[[2]]
        
        rmsx <- meansqerr*(mvar1*mvar1) 
        
        # Betas
        B_out_x <- array(list(), dim=c(nmodelterms, nmodels)) 
        B_out_x[[i, j]] <- matrix(0, nterms, 1)
        # uncommend if the user would like to check the unscaled Betas
        #print(B)
        B_out_x[[i, j]][modsel] <- B 
        B_out_x[[i, j]] <- B_out_x[[i, j]]*(mvar1*scaling)
        B <- B*mvar1*scaling[modsel] 
        
        # printing results
        out_text_x <- array(list(), dim=c(i, j))
        out_text_x[[i, j]] <- cbind('dx', '=')
        
        for (k in 1:i)
          out_text_x[[i, j]] <- cbind(out_text_x[[i, j]], '+', format(B[k], digits = 2), term[[modsel[k]]])    
        
        print(c('Model #', as.character(j), 'of size', as.character(i), 'is:'), quote=FALSE) 
        write.table(format(out_text_x[[i, j]], justify="left", width=1), row.names=FALSE, col.names=FALSE, quote=FALSE)
      } 
      print('--------------------------------------------') 
    }  
    #   # displays the three best models for each nmodelterms for y/var2 (indicator2)
    #   # calling polyfitreg function
    for (i in 1:nmodelterms)
    {
      for (j in 1:nmodels)
      {
        modsel = matrix(0, i, 1)
        modsel[1:i] <- selmod2[i, j, 1:i]
        
        tmp = polyfitreg(indnr, var1, var2, chVar2, modsel, var3, var4) 
        B <- tmp[[1]]
        meansqerr <- tmp[[2]]
        
        rmsy = meansqerr*(mvar2*mvar2) 
        
        # Betas
        B_out_y <- array(list(), dim=c(nmodelterms, nmodels)) 
        B_out_y[[i, j]] <- matrix(0, nterms, 1)
        #print(B)
        B_out_y[[i, j]][modsel] <- B 
        B_out_y[[i, j]] <- B_out_y[[i, j]]*(mvar2*scaling)
        B <- B*mvar2*scaling[modsel] 
        
        # printing results
        out_text_y <- array(list(), dim=c(i, j))
        out_text_y[[i, j]] <- cbind('dy', '=')
        
        for (k in 1:i)
          out_text_y[[i, j]] <- cbind(out_text_y[[i, j]], '+', format(B[k], digits = 2), term[[modsel[k]]])
        
        print(c('Model #', as.character(j), 'of size', as.character(i), 'is:'), quote=FALSE)
        write.table(format(out_text_y[[i, j]], justify="left", width=1), row.names=FALSE, col.names=FALSE, quote=FALSE)
      } 
      print('--------------------------------------------')
    } 
    #   # displays the three best models for each nmodelterms for z/var3 (indicator3)
    #   # calling polyfitreg function
    for (i in 1:nmodelterms)
    {
      for (j in 1:nmodels)
      {
        modsel = matrix(0, i, 1)
        modsel[1:i] <- selmod3[i, j, 1:i]
        
        tmp = polyfitreg(indnr, var1, var2, chVar3, modsel, var3, var4) 
        B <- tmp[[1]]
        meansqerr <- tmp[[2]]
        
        rmsz = meansqerr*(mvar3*mvar3) 
        
        # Betas
        B_out_z <- array(list(), dim=c(nmodelterms, nmodels)) 
        B_out_z[[i, j]] <- matrix(0, nterms, 1)
        #print(B)
        B_out_z[[i, j]][modsel] <- B 
        B_out_z[[i, j]] <- B_out_z[[i, j]]*(mvar3*scaling)
        B <- B*mvar3*scaling[modsel] 
        
        # printing results
        out_text_z <- array(list(), dim=c(i, j))
        out_text_z[[i, j]] <- cbind('dz', '=')
        
        for (k in 1:i)
          out_text_z[[i, j]] <- cbind(out_text_z[[i, j]], '+', format(B[k], digits = 2), term[[modsel[k]]])
        
        print(c('Model #', as.character(j), 'of size', as.character(i), 'is:'), quote=FALSE)
        write.table(format(out_text_z[[i, j]], justify="left", width=1), row.names=FALSE, col.names=FALSE, quote=FALSE)
      } 
      print('--------------------------------------------')
    }
    #   # displays the three best models for each nmodelterms for v/var4 (indicator4)
    #   # calling polyfitreg function
    for (i in 1:nmodelterms)
    {
      for (j in 1:nmodels)
      {
        modsel = matrix(0, i, 1)
        modsel[1:i] <- selmod4[i, j, 1:i]
        
        tmp = polyfitreg(indnr, var1, var2, chVar4, modsel, var3, var4) 
        B <- tmp[[1]]
        meansqerr <- tmp[[2]]
        
        rmsv = meansqerr*(mvar4*mvar4) 
        
        # Betas
        B_out_v <- array(list(), dim=c(nmodelterms, nmodels)) 
        B_out_v[[i, j]] <- matrix(0, nterms, 1)
        #print(B)
        B_out_v[[i, j]][modsel] <- B 
        B_out_v[[i, j]] <- B_out_v[[i, j]]*(mvar4*scaling)
        B <- B*mvar4*scaling[modsel] 
        
        # printing results
        out_text_v <- array(list(), dim=c(i, j))
        out_text_v[[i, j]] <- cbind('dv', '=')
        
        for (k in 1:i)
          out_text_v[[i, j]] <- cbind(out_text_v[[i, j]], '+', format(B[k], digits = 2), term[[modsel[k]]])
        
        print(c('Model #', as.character(j), 'of size', as.character(i), 'is:'), quote=FALSE)
        write.table(format(out_text_v[[i, j]], justify="left", width=1), row.names=FALSE, col.names=FALSE, quote=FALSE)
      } 
      print('--------------------------------------------')
    }
    return(list(B_out_x=B_out_x, B_out_y=B_out_y, B_out_z=B_out_z, B_out_v=B_out_v, out_text_x=out_text_x, 
                out_text_y=out_text_y, out_text_z=out_text_z, out_text_v=out_text_v, scaling=scaling, 
                SEtestx=SEtestx, SEtesty=SEtesty, SEtestz=SEtestz, SEtestv=SEtestv))
  }
}