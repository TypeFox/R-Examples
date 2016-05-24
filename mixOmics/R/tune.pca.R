# Copyright (C) 2011
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence in Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# Leigh Coonan, Student, University of Quuensland, Australia
# Fangzhou Yao, Student, University of Queensland, Australia
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


tune.pca <-
function(X, 
         ncomp = NULL, 
         center = TRUE, 	# sets the mean of the data to zero, ensures that the first PC describes the direction of the maximum variance
         scale = FALSE, 	# variance is unit accross different units
         max.iter = 500,      
         tol = 1e-09)
     
{
    X = as.matrix(X) 	
    X = scale(X, center = center, scale = scale)
    sc = attr(X, "scaled:scale")
    if (any(sc == 0)) {
            stop("cannot rescale a constant/zero column to unit variance.")###???? is this where this should be?
          }
    ncomp. = ncomp

## added value for ncomp if NULL
    if (is.null(ncomp.)) {
        ncomp. = min(nrow(X), ncol(X))
        }


## added warning
    if (ncomp. > min(ncol(X), nrow(X))) {
        stop("Use smaller 'ncomp'")
        }

# Find eigenvalues
    is.na.X = is.na(X)
    na.X = FALSE
    if (any(is.na.X)) na.X = TRUE

    {
    if(na.X){
        result = nipals(X, ncomp = min(ncol(X), nrow(X)), reconst = FALSE, max.iter = 500, tol = 1e-09)
        result$sdev =(result$eig/sqrt(max(1, nrow(X) - 1)))^2
        }
    if(!na.X){
        result = pcasvd(X, ncomp = ncomp., center = center, scale = scale)
        result$sdev = result$sdev^2
        }
    }
    result$ncomp = ncomp.
    
#  list eigenvalues, prop. of explained varience and cumulative proportion of explained variance
        prop.var = as.vector(result$sdev/sum(result$sdev))
        cum.var = as.vector(cumsum(prop.var))
        result$sdev = as.vector(result$sdev)
        names(result$sdev) = paste("PC", 1:length(result$sdev), sep = "")
        names(prop.var) = paste("PC", 1:length(prop.var), sep = "")
        names(cum.var) = paste("PC", 1:length(cum.var), sep = "")

# reduce the output if too many components
# note: if NA values, we have an estimation of the variance using NIPALS
	if(result$ncomp > 10){
        cat("Eigenvalues for the first 10 principal components:", "\n")
        print(result$sdev[1:10])
        cat("\n")
	if(!na.X) {
	        cat("Proportion of explained variance for the first 10 principal components:", "\n")
	        print(prop.var[1:10])
	        cat("\n")
	        cat("Cumulative proportion explained variance for the first 10 principal components:", "\n")
	        print(cum.var[1:10])
	        cat("\n")
		} else {
		cat("Estimated proportion of explained variance for the first 10 principal components:", "\n")
        	print(prop.var[1:10])
        	cat("\n")
        	cat("Estimated cumulative proportion explained variance for the first 10 principal components:", "\n")
        	print(cum.var[1:10])
        	cat("\n")
		} # end NA.X
	cat("(for all principal components, see object$var, object$prop.var and object$cum.var)")
        cat("\n")
	}
	else{
        cat("Eigenvalues for the first ", result$ncomp, "principal components:", "\n")
        print(result$sdev[1:result$ncomp])
        cat("\n")
		if(!na.X) {
        	cat("Proportion of explained variance for the first ", result$ncomp, "principal components:", "\n")
        	print(prop.var[1:result$ncomp])
        	cat("\n")
        	cat("Cumulative proportion explained variance for the first ", result$ncomp, "principal components:", "\n")
        	print(cum.var[1:result$ncomp])
        	cat("\n")
		}
		else{
		cat("Estimated proportion of explained variance for the first ", result$ncomp, "principal components:", "\n")
        	print(prop.var[1:result$ncomp])
        	cat("\n")
        	cat("Estimated cumulative proportion explained variance for the first ", result$ncomp, "principal components:", "\n")
        	print(cum.var[1:result$ncomp])
        	cat("\n")
		} # end NA.X
	}

# Plot the principal components and explained variance
# note: if NA values, we have an estimation of the variance using NIPALS
    if(!na.X) {
        barplot(prop.var[1:result$ncomp], names.arg = 1:result$ncomp, xlab = "Principal Components", 
                ylab = "Proportion of Explained Variance")
    } 
	else {
        barplot(prop.var[1:result$ncomp], names.arg = 1:result$ncomp, xlab = "Principal Components", 
                ylab = "Estimated Proportion of Explained Variance")
	}

	r = list(var = result$sdev, prop.var = prop.var, cum.var = cum.var)
	return(invisible(r))
}
