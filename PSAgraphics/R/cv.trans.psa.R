cv.trans.psa <- function(covariates, fcol = NULL){

#cv.trans.psa takes the covariate data frame and replaces each categorical 
#covariate with n >=3 levels with n new binary covariate columns, one for each level.
# covariates: dataframe of covariates
# fcol: columns containing categorical covariates to be exploded.  If 0 then routine identifies
# categorical columns.
X <- covariates

#J: col2 will be vector with columns that are factors, either user defined or found internally.

n.rows <- dim(X)[1]
col2 <- fcol 
if (is.null(fcol))
  { xclass <- sapply(X, is.numeric) 
    for (i in 1:dim(X)[2]){if (!xclass[i]) {col2[i] <- i} else{col2[i] <- 0}}
  }
  
#J: get rid of possible 0s in col2
col2 <- col2[col2[]!=0] 
if(sum(col2)==0){stop("No Categorical/Factor Columns Identified")}


#J: Binary factors become 0:1.  n-level factors (n >= 3) become n distinct columns, 
#coded 0:1 for each level in the obvious way.  

fac.num <- length(col2)
d <- matrix(0, nrow = dim(X)[1], ncol = fac.num) 
fac.size <- NULL
fac.codes <- NULL
         
for (i in 1:fac.num) 
  {
# First find the levels of each factor, create the names for each of the new columns,
# and code (in the matrix 'd') the factors numerically if they weren't already so coded. 
   fac.i <- X[, col2[i]]
   fac.levels <- matrix(sort(unique(fac.i))) 
   num.levels <- length(fac.levels)
   colnames(fac.levels) <- dimnames(X)[[2]][col2[i]] 
   rownames(fac.levels) <- paste(dimnames(X)[[2]][col2[i]], 1:num.levels)
   fac.codes <- rbind(fac.codes, fac.levels)
   colnames(fac.codes) <- c("Levels")

# fac.size will be 1 for a binary factor, number of levels for other factors.  
   if (num.levels == 2){fac.size[i] = 1} else {fac.size[i] = num.levels} 
   for (j in 1:n.rows) 
     {for (k in 1:num.levels) 
       {if (fac.levels[k] == fac.i[j]) 
          d[j,i] = k-1
       }
     }
  }

# dd has one column for each binary factor, and one column for each *level* of other factors.  Same for dlabel
dd = matrix(0, nrow = n.rows, ncol = sum(fac.size)) 
dlabel = matrix(0, nrow = 1, ncol = sum(fac.size)) 
ko = 0 
for (p in 1:fac.num) 
  {h1 = ko + 1 
   ko = ko + fac.size[p] 
   h2 = ko 
   
    if (fac.size[p] > 2) 
     {
       fac.p <- X[, col2[p]]
       mmdp <- matrix(0, n.rows, length(unique(d[,p]))) 
       dd[ ,h1:h2]  <- ifelse(d[,p] == (col(mmdp)-1), 1, 0) 
       dlabel[ ,h1:h2] =paste(dimnames(X)[[2]][col2[p]], sort(unique(fac.p)), sep='_' )
     } 
    if (fac.size[p] == 1) 
     {
       fac.p <- X[, col2[p]]      
       dd[,h2] = d[,p] 
       dlabel[,h2] =paste(dimnames(X)[[2]][col2[p]], sort(unique(fac.p))[2], sep='_' )
     } 
  
#  if (fac.size[p] > 2) 
#     {mmdp <- matrix(0, n.rows, length(unique(d[,p]))) 
#      dd[ ,h1:h2]  <- ifelse(d[,p] == (col(mmdp)-1), 1, 0) 
#      dlabel[ ,h1:h2] = paste(dimnames(X)[[2]][col2[p]], 1:fac.size[p], sep='_')
#     } 
#   if (fac.size[p] == 1) 
#     {dd[,h2] = d[,p] 
#      dlabel[,h2] = dimnames(X)[[2]][col2[p]]
#     } 
     
     
   X[,col2[p]] = d[,p]
  } 

colnames(dd) = dlabel[1,] 
X2 = cbind(X,dd) 
X = X2[ ,-col2] 
  
  
out <- list(X)
names(out) <- c("covariates.transformed")
return(out)

}
