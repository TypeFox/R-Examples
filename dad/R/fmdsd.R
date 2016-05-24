fmdsd <-
function(x, gaussiand=TRUE, windowh=NULL, kern = NULL, normed=FALSE,
      data.centered=FALSE, data.scaled=FALSE, 
      common.variance=FALSE, nb.factors=3, nb.values=10, sub.title="",
			plot.eigen=TRUE, plot.score=FALSE, nscore=1:3, filename=NULL)
{
#---------------
# Preliminaries
#---------------
# p denotes the dimension of the data
p<-ncol(x)-1;

# The initial data is preserved in x0
# (if the data are centered or reduced, x will contain them)
x0<-x

# Rename the last column of x as group
last.column.name=colnames(x)[ncol(x)]
colnames(x)[ncol(x)] <- "group"
group<-as.factor(x$group);
nb.groups<-length(levels(group));
group.name<-levels(group);

# Control and error message
# on data
if (!prod(apply(as.data.frame(x[,1:p]), 2, is.numeric)))
  stop("The variables must be numerical!")
if (max(is.na(x)) == 1)
  stop("There are NAs in x")
# on the window or window parameter
if (!is.null(windowh))
  {if (is.numeric(windowh))
    {if (length(windowh) > 1)
      {stop("windowh must be either a numeric value, either a list of matrix")
      }
    if (windowh < 1e-12)
      {stop("windowh must be strictly positive!")
      }
    } else 
    {if (is.list(windowh))
      {if (is.null(names(windowh)))
         {stop("the elements of the windowh list must be named")
         } else 
         {if (min(names(windowh)==group.name)<1)
            {stop("the names of the windowh list must be the group names")
            }
         }
      if (p >1)
         {if (min(unlist(lapply(windowh, det))) < 1e-12)
            {stop("All elements of windowh must be positive matrices!")
            }
         } else               
         {if (min(unlist(windowh)) < 1e-12)
            {stop("All elements of windowh must be strictly positive!")
            }
         }
      }
    }
  }

# For now, the only choice of kernel is the Gaussian kernel.
if (!gaussiand)
  {if (is.null(kern))  kern = "gauss" }  else
  kern = ""

# Control and error message
if (min(table(group)) <= 1)
  stop("There should be more than one observation in each group")

# Control of nb.factors and nb.values
if (nb.groups < nb.values)
 {nb.values <- nb.groups ;
 }
if (nb.groups < nb.factors)
 {nb.factors <- nb.groups
 }

# Mean per group
meanL<-by(as.data.frame(x[,1:p]),INDICES=group,FUN=colMeans);
meanL0<-meanL
if(data.centered)
  {# Centering data
  for (i in 1:nb.groups)
     {meanL[[i]]<-numeric(p)
     x[x$group==group.name[i],1:p]=scale(x[x$group==group.name[i],1:p],scale=F)
     }
  }

# Variance per group
varL<-by(as.data.frame(x[,1:p]),INDICES=group,FUN=var);
varL0<-varL

# Correlation matrix or correlation coefficient per group
corL=varL
if (p>1)
   {corL<-by(as.data.frame(x[,1:p]),INDICES=group,FUN=cor);
   } else
   {for (i in 1:nb.groups)
     {corL[[i]]<-1
     }
   }

if(data.scaled)
 {varL<-corL
   for (i in 1:nb.groups)
     {x[x$group==group.name[i],1:p]=scale(x[x$group==group.name[i],1:p])
     }
 }

if(common.variance)
 {for (i in 1:nb.groups)
   {varL[[i]]<-var(x[,1:p])
   }
 }
 
# Skewness et kurtosis coefficients per group
#require(e1071)
skewnessL <- by(as.data.frame(x[, 1:p]), INDICES=group, FUN=apply, 2, skewness)
kurtosisL <- by(as.data.frame(x[, 1:p]), INDICES=group, FUN=apply, 2, kurtosis)

#---------------
# Calculus of the inner products matrix W
#---------------
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
# "l" : list of (definite positive) matrices
# "n" : positive number (common to all densities) with which the variance 
#       matrices are multiplied 
# "a" : NULL, that is the matrice variance of each density is multiplied by the 
#       AMISE window (internally computed by the "" function).
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

# Calculus of the inner products between densities
switch(choix,
  # Case: multivariate, gaussian distribution and estimated parameters
  mg.. =
      {W = mat.ip.l2d.gp(meanL, varL)
       },
  # Case univariate, gaussian distributions with parametres internally estimed 
  ug.. =  
      {W = mat.ip.l2d.gp.u(meanL, varL)
      },
  # Case: multivariate, non Gaussian distribution, density estimated using 
  # Gaussian kernel and AMISE window 
  mnga =
      {nbL<-by(x[,1:p],INDICES=group,FUN=nrow);
      wL<-bandwidth.parameter(p,nbL)
      # Multiplication of the variance by the window parameter
      varLwL<-varL
      for (i in 1:nb.groups)
        {varLwL[[i]]<-varL[[i]]*(wL[[i]]^2)}
      W = mat.ip.l2d.kgw(x, varLwL)
      },
  # Case univariate, non gaussian distributions estimated by gaussian kernel
  # method, and AMISE windows 
  unga =
      {nbL<-by(as.data.frame(x[,1:p]),INDICES=group,FUN=NROW);
      wL<-bandwidth.parameter(p,nbL)
      # Multiplication of the variance by the window parameter
      varLwL<-varL
      for (i in 1:nb.groups)
        {varLwL[[i]]<-varL[[i]]*(wL[[i]]^2)
        }
      W = mat.ip.l2d.kgw.u(x, varLwL)
      },
  # Case: multivariate, non gaussian distributions estimed by gaussian kernel
  # method, and bandwith parameter, common to all densities, given by the user
  mngn =
      {nbL<-by(x[,1:p],INDICES=group,FUN=nrow);
      # Multiplication of the variance by the window parameter
      varLwL<-varL
      for (i in 1:nb.groups)
        {varLwL[[i]]<-varL[[i]]*(windowh^2)}
      W = mat.ip.l2d.kgw(x, varLwL)
      },
  # Case univariate, non gaussian distributions estimed by gaussian kernel
  # method, and bandwith parameter, common to all densities, given by the user    
  ungn =
      {nbL<-by(as.data.frame(x[,1:p]),INDICES=group,FUN=NROW);
      # Multiplication of the variance by the window
      varLwL<-varL
      for (i in 1:nb.groups)
        {varLwL[[i]]<-varL[[i]]*(windowh^2)}
      W = mat.ip.l2d.kgw.u(x, varLwL)
      },
  # Case: multivariate, non gaussian distributions estimated by gaussian kernel
  # method, and windows given as a list of matrices
  mngl =
      {W = mat.ip.l2d.kgw(x, windowh)
      },
  
    # Case univariate, non gaussian distributions estimated by gaussian kernel
    # method, and windows given as a list of numbers
  ungl =
      {W = mat.ip.l2d.kgw.u(x, windowh)
      }
  )

norme<-vector("numeric",nb.groups);
for (i in 1:nb.groups)
 {norme[i]<-sqrt(W[i,i])};
if(normed)
 {# Calculus of the matrix W of the normed pca
  for (i in 1:nb.groups)
   {for (j in 1:i)
     {W[i,j]<-W[i,j]/(norme[i]*norme[j])}};
  for (i in 1:(nb.groups-1))
   {for (j in (i+1):nb.groups)
     {W[i,j]<-W[j,i]}};
 }
 
matdist <- diag(0, nb.groups, nb.groups)
dimnames(matdist) <- list(levels(group), levels(group))
if (!normed)
 {for (i in 2:nb.groups)  for (j in 1:(i-1))
   {matdist[i, j] <- matdist[j, i] <- sqrt(W[i, i] + W[j, j] - 2*W[i, j])
   }
 }   else
 {for (i in 2:nb.groups)  for (j in 1:(i-1))
   {matdist[i, j] <- matdist[j, i] <- sqrt(2 - 2*W[i, j])
   }
 }
matdist <- as.dist(matdist)

add.cst <- FALSE

cmds <- cmdscale(matdist, k = nb.factors, eig = TRUE, add = add.cst, x.ret = TRUE)

# Eigenvalues to display
ep <- cmds$eig
epaff <- ep[1:nb.values]

# Scores to display
coor <- cmds$points

# Change of the names of the grouping variable in the returned results
names(attributes(meanL)$dimnames) <- last.column.name  
names(attributes(varL)$dimnames) <- last.column.name
names(attributes(corL)$dimnames) <- last.column.name
names(attributes(skewnessL)$dimnames) <- last.column.name
names(attributes(kurtosisL)$dimnames) <- last.column.name

# Creation of the list of results in an object of class 'fpcad' :
results <- list( call=match.call(),
  group=last.column.name,
  variables=colnames(x)[1:p],
  inertia=data.frame(eigenvalue=epaff,
          inertia=round(1000*epaff/sum(abs(ep)))/10),
#  contributions=data.frame(group.name,PC=cont),
#  qualities=data.frame(group.name,PC=qual),
  scores=data.frame(group.name,PC=coor),
  norm=data.frame(group.name,norm=norme),
  means=meanL,
  variances=varL,
  correlations=corL,
  skewness=skewnessL,
  kurtosis=kurtosisL);
# Change of the name of the grouping variable in the data frames in results
#colnames(results$contributions)[1]=last.column.name
#colnames(results$qualities)[1]=last.column.name
colnames(results$scores)[1]=last.column.name
#colnames(results$norm)[1]=last.column.name
class(results) <- "fmdsd"

# Save this list in a file (if a filename is given)
if (! is.null(filename))
  {save(results,file = filename)
  }
 
# Barplot of the inertia
if (plot.eigen)
  {inertia=results$inertia$inertia
  barplot(inertia, main="Inertia",
          names.arg=1:length(inertia),
          cex.names=1)
  }
 
# Plotting the scores on the principal planes (1,2),(1,3) et (2,3)
# Plane (a,b)
if (plot.score)
  {plot(results, nscore=nscore, sub.title=sub.title)
  }

# Returning the list of results
return(results)
}
