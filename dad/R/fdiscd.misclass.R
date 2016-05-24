fdiscd.misclass <-
function(x, classe, crit = 1, gaussiand = TRUE, kern = NULL, windowh = NULL)  
{
  # Number of variables:
  p <- ncol(x)-1

  # If necessary: change the last variable of x (the groups) into a factor:
  if (!is.factor(x[, p+1]))  x[, p+1] <- as.factor(x[, p+1])
  
  group<-as.factor(x[, p+1]);
  nb.groups <- nlevels(group)
  group.name<-levels(group);

  # Change the variables of classe (the groups and classes) into factor:
  if (!is.factor(classe[, 1]))  classe[, 1] <- as.factor(classe[, 1])
  if (!is.factor(classe[, 2]))  classe[, 2] <- as.factor(classe[, 2])
  nb.classe <- nlevels(classe[, 2])
  classe.name<-levels(classe[, 2]);

  # Check that the levels of x[, p+1] are exactly the levels of classe[, 1]
  if (!identical(levels(classe[, 1]),levels(x[, p+1])))
   {stop("The names of the groups in the two data frames (x and classe arguments are different).")
   }
   
  # Add the classes to x, as the (p+2)th column:
  x <- data.frame(x, classe = factor(NA, levels = levels(classe[, 2])))
  for (n.x in levels(classe[, 1]))
    x[x[, p+1] == n.x, "classe"] <- classe[classe[, 1] == n.x, 2]
  
  # Means and variances/covariances per group:
  moy.x <- by(as.data.frame(x[, 1:p]), x[, p+1], colMeans)
  var.x <- by(as.data.frame(x[, 1:p]), x[, p+1], var)

  # Means and variances/covariances per classe:
  moy.cl <- by(as.data.frame(x[, 1:p]), x[, p+2], colMeans)
  var.cl <- by(as.data.frame(x[, 1:p]), x[, p+2], var)

  # If the densities are estimated by the kernel method:
  if (!gaussiand)
    {if (is.null(kern))  kern = "gauss" }  else
    kern = ""

  # Which method will be used to estimate the densities and compute the inner products?
  choix = character(4)
  # Choice of the dimension
  # "m" : multivariate ; "u" : univariate
  if (p > 1)              
    {choix[1] = "m"
    }  else
    {choix[1] = "u"
    }
  # Choice of the distribution type
  # "g" : gaussian distributions; "n" : non gaussian distributions 
  if (gaussiand)           
    {choix[2] = "g"                    
    }  else
    {choix[2] = "n"
    }
  # Choice of the kernel 
  # "g" : gaussian kernel; "." : not applicable
  # This option offers a limited choice as the only available kernel is the
  # gaussian kernel 
   if (kern == "gauss")
    {choix[3] = "g"
    }  else
    {choix[3] = "."
    }
  # Choice of the window 
  # The window is given by the user in the "windowh" parameter as 
  # "n" : positive number (common to all densities) with which the variance 
  #       matrices are multiplied 
  # "a" : NULL, that is the matrice variance of each density is multiplied by the 
  #       AMISE window (internally computed by the "bandwidth.parameter" function).
  # "." : not applicable
  if (gaussiand)        
    {choix[4] = "."
    }  else
    {if (is.null(windowh))
      {choix[4] = "a"
      }  else
      {if (p == 1)
         {if (length(windowh) == 1)
            {choix[4] = "n"
            } else
            {choix[4] = "l"
            }
         }  else
         {if (is.list(windowh))
            {choix[4] = "l"
            }  else
            {choix[4] = "n"
            }
         } 
      }
    }                                  
  choix = paste(choix, collapse = "")  
  
  # Criterion 1 (affinity to classe representative):
  # Computation of the inner-products:
if (crit == 1)  {
  # square norms of the groups:
  Wf <- rep(NA, nb.groups)
  names(Wf) <- levels(x[, p+1])
  Wfg <- matrix(NA, nrow=nb.groups, ncol=nb.classe, dimnames = list(levels(x[, p+1]), levels(x[, p+2])))
  Wg <- rep(NA, nb.classe)
  names(Wg) <- levels(x[, p+2])
  switch(choix,
    # Case: multivariate, gaussian distribution and estimated parameters
    mg.. =
        {for (i in levels(group))
           {Wf[i] <- l2d.gp(moy.x[[i]],var.x[[i]],moy.x[[i]],var.x[[i]])
            for (j in levels(classe[, 2]))
              Wfg[i,j]<-l2d.gp(moy.x[[i]],var.x[[i]],moy.cl[[j]],var.cl[[j]])
           }
         for (j in levels(classe[, 2]))
            Wg[j] <- l2d.gp(moy.cl[[j]],var.cl[[j]],moy.cl[[j]],var.cl[[j]])
         },
    # Case univariate, gaussian distributions with parametres internally estimed 
    ug.. =  
        {for (i in levels(group))
           {Wf[i] <- l2d.gp.u(moy.x[[i]],var.x[[i]],moy.x[[i]],var.x[[i]])
            for (j in levels(classe[, 2]))
              Wfg[i,j]<-l2d.gp.u(moy.x[[i]],var.x[[i]],moy.cl[[j]],var.cl[[j]])
           }
         for (j in levels(classe[, 2]))
            Wg[j] <- l2d.gp.u(moy.cl[[j]],var.cl[[j]],moy.cl[[j]],var.cl[[j]])
         },
    # Case: multivariate, non Gaussian distribution, density estimated using 
    # Gaussian kernel and AMISE window 
    mnga =
        {nbL.f <- by(x[,1:p], INDICES=x[, p+1], FUN=nrow)
        wL.f <- bandwidth.parameter(p,nbL.f)
        nbL.g <- by(x[,1:p], INDICES=x[, p+2], FUN=nrow)
        wL.g <- bandwidth.parameter(p,nbL.g)
        # Multiplication of the variance by the window parameter
        varLwL.f <- var.x
        for (i in 1:nb.groups)
          {varLwL.f[[i]] <- var.x[[i]]*(wL.f[[i]]^2)}
        varLwL.g <- var.cl
        for (j in 1:nb.classe)
          {varLwL.g[[j]] <- var.cl[[j]]*(wL.g[[j]]^2)}
        for (i in 1:nb.groups)
          {n.x <- group.name[i]
          ind.i <- which(x[, p+1]==n.x)
          xi<-x[ind.i,1:p]
          Wf[n.x] <- l2d.kgw(xi,varLwL.f[[i]],xi,varLwL.f[[i]])
          for(j in 1:nb.classe)
           {n.cl <- classe.name[j]
           #Test if the determinant of 
            #var.x[[i]]*(wL[[i]]^2)+var.cl[[j]]*(wL[[j]]^2) is different from zero
            if (abs(det(varLwL.f[[i]]+varLwL.g[[j]])) < .Machine$double.eps)
              stop("The matrices are not invertible")
            ind.j <- which(x[, p+2]==n.cl)
            xj<-x[ind.j,1:p]
            Wfg[n.x, n.cl]<-l2d.kgw(xi,varLwL.f[[i]],xj,varLwL.g[[j]])
           }
          }
        for(j in 1:nb.classe)
          {n.cl <- classe.name[j]
           ind.j <- which(x[, p+2]==n.cl)
           xj<-x[ind.j,1:p]
           Wg[n.cl] <- l2d.kgw(xj,varLwL.g[[j]],xj,varLwL.g[[j]])}
        },
    # Case univariate, non gaussian distributions estimated by gaussian kernel
    # method, and AMISE windows 
    unga =
        {nbL.f <- by(as.data.frame(x[,1:p]),INDICES=group,FUN=NROW);
        wL.f <- bandwidth.parameter(p,nbL.f)
        nbL.g <- by(x[,1:p], INDICES=x[, p+2], FUN=NROW)
        wL.g <- bandwidth.parameter(p,nbL.g)
        # Multiplication of the variance by the window parameter
        varLwL.f<-var.x
        for (i in 1:nb.groups)
          {varLwL.f[[i]]<-var.x[[i]]*(wL.f[[i]]^2)}
        varLwL.g <- var.cl
        for (j in 1:nb.classe)
          {varLwL.g[[j]] <- var.cl[[j]]*(wL.g[[j]]^2)}
        for (i in 1:nb.groups)
          {n.x <- group.name[i]
          ind.i <- which(x[, p+1]==n.x)
          xi<-x[ind.i,1:p]
          Wf[n.x] <- l2d.kgw.u(xi,varLwL.f[[i]],xi,varLwL.f[[i]])
          for(j in 1:nb.classe)
            {n.cl <- classe.name[j]
            # Test if the variances are different from zero
            if (varLwL.f[[i]]+varLwL.g[[j]] < .Machine$double.eps)
              {stop("One variance or more is equal to zero")
              }
            ind.j <- which(x[, p+2]==n.cl)
            xj<-x[ind.j,1:p]
            Wfg[n.x,n.cl]<-l2d.kgw.u(xi,varLwL.f[[i]],xj,varLwL.g[[j]])
            }
          };
        for(j in 1:nb.classe)
          {n.cl <- classe.name[j]
          ind.j <- which(x[, p+2]==n.cl)
          xj<-x[ind.j,1:p]
          Wg[n.cl] <- l2d.kgw.u(xj,varLwL.g[[j]],xj,varLwL.g[[j]])}
        },
    # Case: multivariate, non gaussian distributions estimed by gaussian kernel
    # method, and bandwith parameter, common to all densities, given by the user
    mngn =
        {nbL.f <- by(x[,1:p], INDICES=x[, p+1], FUN=nrow)
        nbL.g <- by(x[,1:p], INDICES=x[, p+2], FUN=nrow)
        # Multiplication of the variance by the window parameter
        varLwL.f <- var.x
        for (i in 1:nb.groups)
          {varLwL.f[[i]]<-var.x[[i]]*(windowh^2)}
        varLwL.g <- var.cl
        for (j in 1:nb.classe)
          {varLwL.g[[j]] <- var.cl[[j]]*(windowh^2)}
        for (i in 1:nb.groups)
          {n.x <- group.name[i]
           ind.i <- which(x[, p+1]==n.x)
           xi<-x[ind.i,1:p]
           Wf[n.x] <- l2d.kgw(xi,varLwL.f[[i]],xi,varLwL.f[[i]])
           for(j in 1:nb.classe)
           {n.cl <- classe.name[j]
            # Test if the determinant of
            #var.x[[i]]*(wL[[i]]^2)+var.cl[[j]]*(wL[[j]]^2) is different from zero
            if (abs(det(varLwL.f[[i]]+varLwL.g[[j]])) < .Machine$double.eps)
              stop("The matrices are not invertible")
            ind.j <- which(x[, p+2]==n.cl)
            xj<-x[ind.j,1:p]
            Wfg[n.x,n.cl]<-l2d.kgw(xi,varLwL.f[[i]],xj,varLwL.g[[j]])
           }
          }
          for(j in 1:nb.classe)
           {n.cl <- classe.name[j]
            ind.j <- which(x[, p+2]==n.cl)
            xj<-x[ind.j,1:p]
            Wg[n.cl] <- l2d.kgw(xj,varLwL.g[[j]],xj,varLwL.g[[j]])}
        },
    # Case univariate, non gaussian distributions estimed by gaussian kernel
    # method, and bandwith parameter, common to all densities, given by the user    
    ungn =
        {nbL.f <- by(x[,1:p], INDICES=x[, p+1], FUN=NROW)
        nbL.g <- by(x[,1:p], INDICES=x[, p+2], FUN=NROW)
        # Multiplication of the variance by the window parameter
        varLwL.f <- var.x
        for (i in 1:nb.groups)
          {varLwL.f[[i]]<-var.x[[i]]*(windowh^2)}
        varLwL.g <- var.cl
        for (j in 1:nb.classe)
          {varLwL.g[[j]] <- var.cl[[j]]*(windowh^2)}
        for (i in 1:nb.groups)
          {n.x <- group.name[i]
           ind.i <- which(x[, p+1]==n.x)
           xi<-x[ind.i,1:p]
           Wf[n.x] <- l2d.kgw.u(xi,varLwL.f[[i]],xi,varLwL.f[[i]])
          for(j in 1:nb.classe)
           {n.cl <- classe.name[j]
           # Test if the variances are different from zero
            if (varLwL.f[[i]]+varLwL.g[[j]] < .Machine$double.eps)
              {stop("One variance or more is equal to zero")
              }
            ind.j <- which(x[, p+2]==n.cl)
            xj<-x[ind.j,1:p]
            Wfg[n.x,n.cl]<-l2d.kgw.u(xi,varLwL.f[[i]],xj,varLwL.g[[j]])
           }
          }
          for(j in 1:nb.classe)
           {n.cl <- classe.name[j]
            ind.j <- which(x[, p+2]==n.cl)
            xj<-x[ind.j,1:p]
            Wg[n.cl] <- l2d.kgw.u(xj,varLwL.g[[j]],xj,varLwL.g[[j]])}
        }
    )
}

  # Criterion 2 or 3 (affinity to class center):
  # Computation of the inner-products (W matrix):
if (crit %in% 2:3)
 {W<-matrix(NA, ncol=nb.groups, nrow=nb.groups, dimnames = list(levels(x[, p+1]), levels(x[, p+1])))
  switch(choix,
    # Case: multivariate, gaussian distribution and estimated parameters
    mg.. =
        {for (i in 1:nb.groups)
           {for (j in 1:i)
              {# Test if det(var.x[[i]]+var.x[[j]]) is different from zero
              if (abs(det(var.x[[i]]+var.x[[j]])) < .Machine$double.eps)
                {stop(paste(c("Singular covariance matrix (group ", unique(names(var.x)[c(i, j)]), ")"), collapse = ""))
                }
                W[i,j]<-l2d.gp(moy.x[[i]],var.x[[i]],moy.x[[j]],var.x[[j]])
              }
           }
         },
    # Case univariate, gaussian distributions with parametres internally estimed 
    ug.. =  
        {for (i in 1:nb.groups)
          {for (j in 1:i)
            {#Test if var.x[[i]]+var.x[[j]] is different from zero
            if (abs(var.x[[i]]+var.x[[j]]) < .Machine$double.eps)
              {stop(paste(c("Singular covariance matrix (group ", unique(names(var.x)[c(i, j)]), ")"), collapse = ""))
              }
            W[i,j]<-l2d.gp.u(moy.x[[i]],var.x[[i]],moy.x[[j]],var.x[[j]])
            }
          }
        },
    # Case: multivariate, non Gaussian distribution, density estimated using 
    # Gaussian kernel and AMISE window 
    mnga =
        {nbL<-by(x[,1:p], INDICES=x[, p+1], FUN=nrow);
        wL<-bandwidth.parameter(p,nbL)
        # Multiplication of the variance by the window parameter
        varLwL<-var.x
        for (i in 1:nb.groups)
          {varLwL[[i]]<-var.x[[i]]*(wL[[i]]^2)}
        for (i in 1:nb.groups)
          {ind.i <- which(x[, p+1]==group.name[i])
          xi<-x[ind.i,1:p]
          for(j in 1:i)
           {#Test if the determinant of 
            #var.x[[i]]*(wL[[i]]^2)+var.x[[j]]*(wL[[j]]^2) is different from zero
            if (abs(det(varLwL[[i]]+varLwL[[j]])) < .Machine$double.eps)
              stop("The matrices are not invertible")
            ind.j <- which(x[, p+1]==group.name[j])
            xj<-x[ind.j,1:p]
            W[i,j]<-l2d.kgw(xi,varLwL[[i]],xj,varLwL[[j]])
           }
          }
        },
    # Case univariate, non gaussian distributions estimated by gaussian kernel
    # method, and AMISE windows 
    unga =
        {nbL<-by(as.data.frame(x[,1:p]),INDICES=group,FUN=NROW);
        wL<-bandwidth.parameter(p,nbL)
        # Multiplication of the variance by the window parameter
        varLwL<-var.x
        for (i in 1:nb.groups)
          {varLwL[[i]]<-var.x[[i]]*(wL[[i]]^2)
          }
        for (i in 1:nb.groups)
          {ind.i <- which(x[, p+1]==group.name[i])
          xi<-x[ind.i,1:p]
          for(j in 1:i)
            {# Test if the variances are different from zero
            if (varLwL[[i]]+varLwL[[j]] < .Machine$double.eps)
              {stop("One variance or more is equal to zero")
              }
            ind.j <- which(x[, p+1]==group.name[j])
            xj<-x[ind.j,1:p]
            W[i,j]<-l2d.kgw.u(xi,varLwL[[i]],xj,varLwL[[j]])
            }
          };
        },
    # Case: multivariate, non gaussian distributions estimed by gaussian kernel
    # method, and bandwith parameter, common to all densities, given by the user
    mngn =
        {nbL<-by(x[,1:p],INDICES=group,FUN=nrow);
        # Multiplication of the variance by the window parameter
        varLwL<-var.x
        for (i in 1:nb.groups)
          {varLwL[[i]]<-var.x[[i]]*(windowh^2)}
        for (i in 1:nb.groups)
          {ind.i <- which(x[, p+1]==group.name[i])
          xi<-x[ind.i,1:p]
          for(j in 1:i)
           {ind.j <- which(x[, p+1]==group.name[j])
           xj<-x[ind.j,1:p]
           W[i,j]<-l2d.kgw(xi,varLwL[[i]],xj,varLwL[[j]])
           }
          }
        },
    # Case univariate, non gaussian distributions estimed by gaussian kernel
    # method, and bandwith parameter, common to all densities, given by the user    
    ungn =
        {nbL<-by(as.data.frame(x[,1:p]),INDICES=group,FUN=NROW);
        # Multiplication of the variance by the window
        varLwL<-var.x
        for (i in 1:nb.groups)
          {varLwL[[i]]<-var.x[[i]]*(windowh^2)}
        for (i in 1:nb.groups)
          {ind.i <- which(x[, p+1]==group.name[i])
          xi<-x[ind.i,1:p]
          for(j in 1:i)
           {ind.j <- which(x[, p+1]==group.name[j])
            xj<-x[ind.j,1:p]
            W[i,j]<-l2d.kgw.u(xi,varLwL[[i]],xj,varLwL[[j]]) }};
        }
    )
    for (i in 1:(nb.groups-1))
      {for (j in (i+1):nb.groups)                             
        {W[i,j]<-W[j,i] }};
    # End of the computation of the matrix W
 }
  
  # Computation of the distances between groups and classes:
  distances <- matrix(nrow = nrow(classe), ncol = nlevels(classe[, 2]),
                      dimnames = list(as.character(classe[, 1]), levels(classe[, 2])))
  switch(crit,
      "1" = {  # Criterion 1. Affinity to classe representative:
        for (n.x in rownames(distances))  if (sum(x[, p+1] == n.x) > 0)  {
          for (n.cl in colnames(distances))  {
            Wf.tmp <- Wf[n.x]
            if (!(n.x %in% classe[classe[, 2] == n.cl, 1]))
              {Wfg.tmp <- Wfg[n.x, n.cl]
              Wg.tmp <- Wg[n.cl]
              }  else
              {ind.x <- which(x[, p+1] == n.x)
              ind.cl <- which((x[, p+2] == n.cl)&(x[, p+1] != n.x))
              moy.cl.tmp <- colMeans(as.data.frame(x[ind.cl, 1:p]))
              var.cl.tmp <- var(as.data.frame(x[ind.cl, 1:p]))
              switch(choix,
                # Case: multivariate, gaussian distribution and estimated parameters
                mg.. =
                   {Wfg.tmp<-l2d.gp(moy.x[[n.x]],var.x[[n.x]],moy.cl.tmp,var.cl.tmp)
                    Wg.tmp <- l2d.gp(moy.cl.tmp,var.cl.tmp,moy.cl.tmp,var.cl.tmp)
                    },
                # Case univariate, gaussian distributions with parametres internally estimed 
                ug.. =  
                    {Wfg.tmp<-l2d.gp.u(moy.x[[n.x]],var.x[[n.x]],moy.cl.tmp,var.cl.tmp)
                     Wg.tmp <- l2d.gp.u(moy.cl.tmp,var.cl.tmp,moy.cl.tmp,var.cl.tmp)
                     },
                # Case: multivariate, non Gaussian distribution, density estimated using 
                # Gaussian kernel and AMISE window 
                mnga =
                    {nbL <- length(ind.cl)
                    wL <- bandwidth.parameter(p,nbL)
                    # Multiplication of the variance by the window parameter
                    varLwL.tmp <- var.cl.tmp*(wL^2)
                    Wfg.tmp<-l2d.kgw(x[ind.x, 1:p],varLwL.f[[n.x]],x[ind.cl, 1:p],varLwL.tmp)
                    Wg.tmp <- l2d.kgw(x[ind.cl, 1:p],varLwL.tmp,x[ind.cl, 1:p],varLwL.tmp)
                    },
                # Case univariate, non gaussian distributions estimated by gaussian kernel
                # method, and AMISE windows 
                unga =
                    {nbL <- length(ind.cl)
                    wL <- bandwidth.parameter(p,nbL)
                    # Multiplication of the variance by the window parameter
                    varLwL.tmp <- var.cl.tmp*(wL^2)
                    Wfg.tmp <- l2d.kgw.u(x[ind.x, 1:p],varLwL.f[[n.x]],x[ind.cl, 1:p],varLwL.tmp) 
                    Wg.tmp <- l2d.kgw.u(x[ind.cl, 1:p],varLwL.tmp, x[ind.cl, 1:p],varLwL.tmp)
                    },
                # Case: multivariate, non gaussian distributions estimed by gaussian kernel
                # method, and bandwith parameter, common to all densities, given by the user
                mngn =
                    {nbL <- length(ind.cl)
                    # Multiplication of the variance by the window parameter
                    varLwL.tmp <- var.cl.tmp*(windowh^2)
                    Wfg.tmp <- l2d.kgw(x[ind.x, 1:p],varLwL.f[[n.x]],x[ind.cl, 1:p],varLwL.tmp) 
                    Wg.tmp <- l2d.kgw(x[ind.cl, 1:p],varLwL.tmp, x[ind.cl, 1:p],varLwL.tmp)
                    },
                # Case univariate, non gaussian distributions estimed by gaussian kernel
                # method, and bandwith parameter, common to all densities, given by the user    
                ungn =
                    {nbL <- length(ind.cl)
                    # Multiplication of the variance by the window parameter
                    varLwL.tmp <- var.cl.tmp*(windowh^2)
                    Wfg.tmp <- l2d.kgw.u(x[ind.x, 1:p],varLwL.f[[n.x]],x[ind.cl, 1:p],varLwL.tmp)
                    Wg.tmp <- l2d.kgw.u(x[ind.cl, 1:p],varLwL.tmp, x[ind.cl, 1:p],varLwL.tmp)
                    }
                )
              }
            # Computation of the group-classe distances:
            distances[n.x, n.cl] <- sqrt(Wf.tmp - 2*Wfg.tmp + Wg.tmp)
          }
        }
      },
      "2" = {  # Criterion 2. Affinity to classe center:
        for (n.x in rownames(distances))  {
          if (sum(x[, p+1] == n.x) > 0) {
            d2ff <- W[n.x, n.x]
            for (n.cl in colnames(distances))  {
              ind.cl <- which((x[, p+2] == n.cl)&(x[, p+1] != n.x))
              Tj <- length(unique(x[ind.cl, p+1]))
              x.cl <- unique(x[ind.cl, p+1])
              d2fg <- sum(W[n.x, x.cl])/Tj
              d2gg <- sum(W[x.cl, x.cl])/(Tj^2)
              distances[n.x, n.cl] <- sqrt(d2ff - 2*d2fg + d2gg)
            }
          }
        }
      },
      "3" = {  # Criterion 3. Affinity to classe center (weighted):
        for (n.x in rownames(distances))  {
          if (sum(x[, p+1] == n.x) > 0) {
            d2ff <- W[n.x, n.x]
            for (n.cl in colnames(distances))  {
              ind.cl <- which((x[, p+2] == n.cl)&(x[, p+1] != n.x))    
              nr <- table(as.factor(as.character(x[ind.cl, p+1])))
              nj <- length(ind.cl)
              x.cl <- unique(x[ind.cl, p+1])
              # Parts of W matrix which will be used in the calculus of scalar products:
              Wij <- W[n.x, x.cl]
              Wij <- Wij[names(nr)]
              Wjj <- W[x.cl, x.cl]
              Wjj <- Wjj[names(nr), names(nr)]
              d2fg <- sum(nr*Wij)/nj
              d2gg <- as.numeric(rbind(nr) %*% Wjj %*% cbind(nr))/(nj^2)
              distances[n.x, n.cl] <- sqrt(d2ff - 2*d2fg/nj + d2gg/(nj^2))
            }
          }
        }
      }
    )

  group.calc <- unlist(apply(distances, 1, which.min))
  names(group.calc) <- unlist(lapply(strsplit(names(group.calc), ".", fixed = TRUE), head, 1))
  
  misclass <- classe
  colnames(misclass)[2] <- "prior.class"
  rownames(misclass) <- classe[, 1]
  lev.class <- levels(misclass[, 2])
  misclass <- data.frame(misclass, predicted.class = factor(NA, levels = lev.class))
  misclass[names(group.calc), "predicted.class"] <-  lev.class[group.calc]
  misclass <- data.frame(misclass, misclassed = (misclass[, 2] != misclass[, 3]))

  # Allocations per classe:
  alloc.per.classe <- table(misclass[, 2:3])
  
  # Changing of distances into proximities
  inv.distances <- 1/(distances)
  sum.inv.distances <- rowSums(inv.distances)
  prox <- inv.distances/sum.inv.distances
  
  # Misclassification percents per classe
  classe.ratio <- diag(alloc.per.classe)/rowSums(alloc.per.classe)
  classe.ratio <- classe.ratio*100
  misclass.per.classe <- 100 - classe.ratio
  
  results <- list(classification = misclass,
              confusion.mat = alloc.per.classe,
              misalloc.per.classe = misclass.per.classe,
              misclassed = sum(misclass$misclassed, na.rm = TRUE)/sum(!is.na(misclass$misclassed)),
              distances = distances,
              proximities = round(prox, 3)*100)
  
  class(results) <- "fdiscd.misclass"
  
  return(results)
}
