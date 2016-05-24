# Copyright (C) 2009
# Kim-Anh Le Cao,
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
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


# object: a pls/spls object
# comp: to display the variables selected on dimension 'comp'
# names.X, names.Y: set to true means that the X and Y data frames have row names (see example below with srbct)

selectVar <-
function(...) UseMethod("selectVar")


get.name.and.value=function(x,comp)
{
    name.var=names(sort(abs(x[,comp]), decreasing = T)[1:sum(x[,comp]!=0)])
    value.var=x[name.var,comp]
    return(list(name = name.var, value = data.frame(value.var)))
}


# ------------------ for all object  --------------------
selectVar.pls  = selectVar.plsda  =
selectVar.sgcca = selectVar.rgcca = 
selectVar.pca =selectVar.sipca=
function(object, comp =1, block=NULL, ...)
{

    # check arguments
    # -----------------
    if(length(comp) > 1)
    stop("Expecting one single value for 'comp'")
    if(is.null(block)){
        if(any(comp > object$ncomp))
        stop("'comp' is greater than the number of components in the fitted model")
    }else{
        if(any(comp > object$ncomp[block]))
        stop("'comp' is greater than the number of components in the fitted model for the block you specified")
    }
    
    if(is.null(block))
    {
        null.block=TRUE
        block=1:length(object$loadings)
    }else{null.block=FALSE}
    
    # main function: get the names and values of the non zero loadings
    # -----------------
    out=lapply(object$loadings[block],get.name.and.value,comp=comp)
    
    
    # outputs
    # ----------
    #if all blocks are considered by default (null.block=TRUE) and it's a DA analysis, then we don't show Y
    if(null.block)
    {
        if(length(grep("plsda",class(object)))>0)
        {
            out=out[-2] #remove Y
            out=out[[1]]
        }
        if(length(grep("sgccda",class(object)))>0)
        {
            out=out[-object$indY] #remove Y
        }
    }
    
    if(length(grep("pca",class(object)))>0 | length(grep("sipca",class(object)))>0)
    {
        out=out[[1]]
    }
    
    #we add comp as an output
    out$comp=comp
    
    return(out)
}

