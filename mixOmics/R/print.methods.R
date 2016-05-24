# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# Leigh Coonan, Student, University of Queensland, Australia
# Fangzhou Yao, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia and
# Shangai University of Finance and Economics, Shanghai, P.R. China
# Jeff Coquery, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia and
# Sup Biotech, Paris, France

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


#------------------ print method for pls ------------------#
print.pls <-
function(x, ...) 
{

    mode = paste("'", x$mode, "'", sep = "")
	
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
	
    cat(" PLS with a", mode, "mode with", x$ncomp, "PLS components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")

    cat(" No variable selection. \n\n")
	
    cat(" Available components: \n", 
        "-------------------- \n")
	
	cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}

#------------------ print method for plsda ------------------#
print.plsda <-
function(x, ...) 
{
	
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
	
    cat(" PLS-DA (regression mode) with", x$ncomp, "PLS-DA components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")

    cat(" No variable selection. \n\n")
	
    cat(" Available components: \n", 
        "-------------------- \n")
	
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}


#----------------- print method for spls ------------------#
print.spls <-
function(x, ...)
{

    mode = paste("'", x$mode, "'", sep = "")
    keepX = paste("[", x$keepX, "]", sep = "")
    keepY = paste("[", x$keepY, "]", sep = "")
	
	cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
	
    cat(" sPLS with a", mode, "mode with", x$ncomp, "sPLS components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")
	
    cat(" Selection of", keepX, "variables on each of the sPLS components on the X data set. \n")
    cat(" Selection of", keepY, "variables on each of the sPLS components on the Y data set. \n\n")

    cat(" Available components: \n", 
        "-------------------- \n")
	
	cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}


#----------------- print method for splsda ------------------#
print.splsda <-
function(x, ...)
{

    keepX = paste("[", x$keepX, "]", sep = "")
    keepY = paste("[", x$keepY, "]", sep = "")
	
	cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")	
    cat(" sPLS-DA (regression mode) with", x$ncomp, "sPLS-DA components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")
	
    cat(" Selection of", keepX, "variables on each of the sPLS-DA components on the X data set. \n")
    cat(" No Y variables can be selected. \n\n")

    cat(" Available components: \n", 
        "-------------------- \n")
	
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}


#------------------ print method for rcc ------------------#
print.rcc <-
function(x, ...)
{

    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")

    cat(" rCCA with", x$ncomp, "components and regularization parameters", x$lambda[1], "and", x$lambda[2], "for the X and Y data. \n")
    cat(" You entered data X of dimensions :", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions :", nrow(x$Y), ncol(x$Y), "\n\n")
	
    cat(" Available components: \n", 
        "-------------------- \n")
	
	cat(" canonical correlations: see object$cor \n")
	cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}


#------- print for summary with (s)PLS object or rcc ---------#
print.summary <-
function(x, ...)
{

    print.gap = 4
    what = x$what
    digits = x$digits

	#--------------------- output pls/spls ---------------------#
    if(x$method == "pls" | x$method == "spls"){

        if (x$method == "pls") {
            cat(" PLS mode:", x$mode)
			cat("\n Number of variates considered:", x$ncomp, "\n")
        }	
        else {
            cat(" sPLS mode:", x$mode)
			cat("\n Number of variates considered:", x$ncomp)
			cat("\n Number of X-variables selected on each of the sPLS components:", x$keepX)
			cat("\n Number of Y-variables selected on each of the sPLS components:", x$keepY, "\n")
        }			

        #---------- affichage communaute ----------#
        if (any(what == "all") || any(what == "communalities")) { 
            cat("\n\n Communalities Analysis:\n",
                "----------------------")
				
            cat("\n X-Variables vs their own Variates: see object$CM.X$own \n")
            cat("\n X-Variables vs the opposite Variates: see object$CM.X$opp \n")
            cat("\n Y-Variables vs their own Variates: see object$CM.Y$opp \n")
            cat("\n Y-Variables vs the opposite Variates: see object$CM.Y$opp \n")
        }

        #--------- affichage redondance -----------#
        if (any(what == "all") || any(what == "redundancy")) {
            cat("\n\n Redundancy Analysis:\n",
                "-------------------\n")
				
            cat("\n X-Variables vs their own Variates: see object$Rd.X$own \n")
            cat("\n X-Variables vs the opposite Variates: see object$Rd.X$opp \n")
            cat("\n Y-Variables vs their own Variates: see object$Rd.Y$opp \n")
            cat("\n Y-Variables vs the opposite Variates: see object$Rd.Y$opp \n")
        }

        #---------- tableau VIP ---------#
        if (any(what == "all") || any(what == "VIP")) {
            cat("\n\n", "Variable Importance in the Projection (VIP): see object$VIP \n",
                        "------------------------------------------- \n\n")
        }

    }  #end if pls

    # ---------------------- output rcc ------------------------#
    if(x$method == "rcc" ){
        print.gap = 4
        if (any(what == "all")) {
            cat(" Number of canonical variates considered:", x$ncomp, "\n")
            cat("\n Canonical correlations:",
                "\n ----------------------\n")
            print(round(x$can.cor, digits = digits), print.gap = print.gap)
        }

        #-- affichage communaute --#
        if (any(what == "all") || any(what == "communalities")) { 
            cat("\n\n Canonical Communalities Analysis:\n",
                "--------------------------------")

            cat("\n X-Variables vs their own Canonical Variates: see object$Cm.X$own \n")
            cat("\n X-Variables vs the opposite Canonical Variates: see object$Cm.X$opp \n")
            cat("\n Y-Variables vs their own Canonical Variates: see object$Cm.Y$own \n")
            cat("\n Y-Variables vs the opposite Canonical Variates: see object$Cm.Y$opp \n")
        }

        #--------- affichage redondance -----------#
        if (any(what == "all") || any(what == "redundancy")) {
            cat("\n\n Redundancy Analysis:\n",
                "-------------------\n")
				
            cat("\n X-Variables vs their own Variates: see object$Rd.X$own \n")
            cat("\n X-Variables vs the opposite Variates: see object$Rd.X$opp \n")
            cat("\n Y-Variables vs their own Variates: see object$Rd.Y$opp \n")
            cat("\n Y-Variables vs the opposite Variates: see object$Rd.Y$opp \n")
        }

    }  #end rcc
}


# ------------------------ print for pca --------------------------------
print.pca <- function(x, ...) {
    
 
    per.var = as.vector(x$sdev/sum(x$sdev))
    cum.var=as.vector(cumsum(per.var))
    x$sdev=as.vector(x$sdev)
    names(x$sdev) = paste("PC", 1:length(x$sdev), sep = "")
    names(per.var) = paste("PC", 1:length(per.var), sep = "")
    names(cum.var) = paste("PC", 1:length(cum.var), sep = "")
    

    cat("Eigenvalues for the first ", x$ncomp, "principal components:", "\n")
    print(x$sdev[1:x$ncomp])
    cat("\n")
    
    
    if(!x$NA.X) {
        cat("Proportion of explained variance for the first ", x$ncomp, "principal components:", "\n")
        print(per.var[1:x$ncomp])
        cat("\n")

        cat("Cumulative proportion explained variance for the first ", x$ncomp, "principal components:", "\n")
        print(cum.var[1:x$ncomp])
        cat("\n")

        cat(" Other available components: \n", "-------------------- \n")
        cat(" loading vectors: see object$rotation \n")
    }
}

# ------------------------ print for spca -------------------------
print.spca <-
function(x, ...) 
{
	
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
	
    cat(" sparse pCA with", x$ncomp, "principal components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")

    cat(" Selection of", x$keepX, "variables on each of the principal components on the X data set. \n")

    cat(" Available components: \n", 
        "-------------------- \n")
	
    cat(" loading vectors: see object$rotation \n")
    cat(" principal components: see object$x \n")
    cat(" cumulative explained variance: see object$varX \n")
    cat(" variable names: see object$names \n")

}

# ------------------------ print for ipca -------------------------
print.ipca <-
function(x, ...) 
{
	
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
	
    cat(" IPCA with", x$ncomp, "independent components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")

    cat(" Available components: \n", 
        "-------------------- \n")
	
    cat(" unmixing matrix: see object$unmixing \n")
    cat(" independent principal components: see object$x \n")
    cat(" mxing matrix: see object$mixing \n")
    cat(" kurtosis: see object$kurtosis \n")
    cat(" variable names: see object$names \n")
    cat(" independent loading vectors: see object$loadings \n")
}

# ------------------------ print for sipca -------------------------
print.sipca <-
function(x, ...) 
{
	
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
	
    cat(" Sparse IPCA with", x$ncomp, "independent components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
	
	cat(" Selection of", x$keepX, "variables on each of the principal components on the X data set. \n")

    cat(" Available components: \n", 
        "-------------------- \n")
	
    cat(" unmixing matrix: see object$unmixing \n")
	cat(" independent principal components: see object$x \n")
    cat(" mxing matrix: see object$mixing \n")
    cat(" kurtosis: see object$kurtosis \n")
    cat(" variable names: see object$names \n")
    cat(" independent loading vectors: see object$loadings \n")

}

# ------------------------ print for rgcca -------------------------
print.rgcca <- function(x, ...) {
    
    cat("\nCall:\n", deparse(x$class, width.cutoff = 500), "\n\n")
    
    # components
    for(k in 1:length(x$blocks)){
        cat(" rGCCA with", x$ncomp[[k]], "components on block", k, "named", x$names$blocks[k], "\n")
    }
    cat("\n")
    
    # dimension
    for(k in 1 : length(x$blocks)){
        cat(" Dimension of block", k, 'is ', dim(x$blocks[[k]]), "\n")
    }
    cat("\n")
    cat(" Available components: \n", "-------------------- \n")
    
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}


# ------------------------ print for sgcca -------------------------
print.sgcca<- function(x, ...){
    
    cat("\nCall:\n", deparse(x$class, width.cutoff = 500), "\n\n")
    
    # components
    for(k in 1 : length(x$blocks)){
        cat(" sGCCA with", x$ncomp[[k]], "components on block", k, "named", x$names$blocks[k], "\n")
    }
    cat("\n")
    
    # dimension
    for(k in 1 : length(x$blocks)){
        cat(" Dimension of block", k, 'is ', dim(x$blocks[[k]]), "\n")
    }
    cat("\n")
    
    # selected variables
    list.select = list()
    for(k in 1:length(x$blocks)){
        list.select[[k]] = apply(x$loadings[[k]], 2, function(x){sum(x!=0)})
        cat(" Selection of", list.select[[k]], "variables on each of the sGCCA components on the block", k, "\n")
    }
    cat("\n")
    cat(" Available components: \n", "-------------------- \n")
    
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}

