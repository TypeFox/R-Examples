# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
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


plsda <-
function(X, 
         Y, 
         ncomp = 2,
         max.iter = 500,		 
         tol = 1e-06,
         near.zero.var = TRUE)
{
    
    #  Testing the input Y
    if (is.null(dim(Y)))
    {
        Y = as.factor(Y)	
        ind.mat = unmap(as.numeric(Y))					
    }else{
        stop("'Y' should be a factor or a class vector.")						
    }		
 
    result = pls(X, ind.mat, ncomp = ncomp, mode = "regression",
                 max.iter = max.iter, tol = tol,near.zero.var=near.zero.var)
     
    cl = match.call()
    cl[[1]] = as.name('plsda')
    result$call = cl
    	
    result$ind.mat = ind.mat
    result$names$Y = levels(Y)
    
    class(result) = "plsda"
    return(invisible(result))	
}
