# Copyright (C) 2015
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
# Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France

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


# ---------------------------------------------------------

# to do later
# create a S3 function?



plotContrib = function(object,
                        contrib="max",  # choose between 'max" or "min"
                        method = "mean", # choose between 'mean" or "median"
                        block=1, #single value
                        comp = 1,
                        show.ties = TRUE,
                        col.ties="white",
                        ndisplay = NULL,
                        cex.name = 0.7,
                        cex.legend = 0.8,
                        name.var=NULL,
                        complete.name.var=FALSE,
                        legend = TRUE,
                        legend.color = NULL,
                        main = NULL,
                        legend.title = 'Outcome',
                        plot = TRUE
) {
    
    # -------------------
    # input checks
    # ------------------
    
    class.object=c("plsda","splsda","sgccda")
    if(!any(class(object)%in%class.object)) stop("plotContrib is only available for objects of the following class: ", paste(class.object,collapse=" "))
    class.object = class(object)
    
    # method
    # ----
    if(length(method) !=1 || !method%in%c("mean","median")){
        method = 'median'
        warning('Argument method: we expect either mean or median, set to median by default')
    }
    # cex
    # --
    if(cex.name <= 0){cex.name = 0.7}
    if(cex.legend <= 0){cex.name = 0.8}
    
    ncomp = object$ncomp
    
    if(length(block)>1)
    stop("Incorrect value for 'block', a single value (numeric or character) is expected", call. = FALSE)

    ##selectvar
    selected.var=selectVar(object,comp=comp,block=block) # gives name and values of the blocks in 'block'
    name.selected.var=selected.var[[1]]$name
    value.selected.var=selected.var[[1]]$value

    # ndisplay
    # ------
    # if null set by default to all variables from selectVar
    if(is.null(ndisplay)){
        ndisplay = length(name.selected.var)
    }else if(ndisplay > length(name.selected.var)){
        message("'ndisplay' value is larger than the number of selected variables! It has been reseted to ", length(name.selected.var))
        ndisplay=length(name.selected.var)
    }
    
    name.selected.var=name.selected.var[1:ndisplay]
    value.selected.var=value.selected.var[1:ndisplay,]

    #comp
    # ----
    ncomp = object$ncomp[block]
    if(any(max(comp) > ncomp)){
        stop(paste("Argument comp should be less or equal to", ncomp))
    }
    
    names.block=as.character(names(selected.var)[1]) #it should be one block and ncomp, so we take the first one
    #what follows can be combined if pls-objects get a "block" output
    if(any(class(object)%in%"plsda"))
    {
        X=object[names.block][[1]]
    }
    if(any(class(object)%in%"sgccda"))
    {
        X=object$block[names.block][[1]]
    }
    
    #name.var
    ind.match=match(name.selected.var,colnames(X)) # look at the position of the selected variables in the original data X
    if(!is.null(name.var))
    {
        if(length(name.var)!=ncol(X))
        stop("name.var should be a vector or length ", ncol(X))
        colnames.X=as.character(name.var[ind.match]) # get the
    }else{
        colnames.X=as.character(colnames(X))[ind.match]
    }
    X=X[,name.selected.var] #reduce the problem to ndisplay
    
    #completing colnames.X by the original names of the variables when missing
    if(complete.name.var==TRUE)
    {
        ind=which(colnames.X=="")
        if(length(ind)>0) colnames.X[ind]=colnames(X)[ind]
    }

    if(any(class(object)%in%"sgccda"))
    {
        Y = object$Y
    }else{
        Y = factor(map(object$ind.mat), labels= object$names$Y)
    }
    
    #title
    #-----
    if(!is.null(main) & !is.character(main)){
        warning('main needs to be of type character')
    }
    
    #legend.color
    #-----
    if(!is.null(legend.color) & (length(legend.color) != nlevels(Y))){
        warning('legend.color must be the same length than the number of group, by default set to default colors')
        legend.color = color.mixo(1:10)  # by default set to the colors in color.mixo (10 colors)
    }
    if(is.null(legend.color)) legend.color = color.mixo(1:10) # by default set to the colors in color.mixo (10 colors)
    
    if(col.ties%in%legend.color[1:nlevels(Y)])
    {
        stop("'col.ties' should not be in 'legend.color'")
    }
    
    # --------------------
    # end check inputs
    # ---------------------
    
    # ==================================================
    # First, calculate the contribution of the loading
    # =================================================
    
    # Start: Initialisation
    which.comp = method.group = list()
    which.contrib = data.frame(matrix(FALSE,
    ncol = nlevels(Y) + 2, nrow = ndisplay,
    dimnames = list(name.selected.var, c(paste0("Contrib.", levels(Y)), "Contrib", "GroupContrib"))))
    #    dimnames = list(name.selected.var, c(paste0("Contrib.", levels(Y)), "Contrib", "GroupContrib"))))
    # End: Initialisation
    
    # calculate the max.method per group for each variable, and identifies which group has the max max.method
    for(k in 1:ncol(X)){
        if (method == 'mean'){
            method.group[[k]] = tapply(X[, k], Y, mean,na.rm=TRUE)
        } else if (method == 'median'){
            method.group[[k]] = tapply(X[, k], Y, median,na.rm=TRUE)
        }
        # determine which group has the highest mean/median
        if(contrib=="max")
        {
            which.contrib[k, 1:nlevels(Y)] = (method.group[[k]]) == max((method.group[[k]]))
        }else if(contrib=="min")
        {
            which.contrib[k, 1:nlevels(Y)] = (method.group[[k]]) == min((method.group[[k]]))
        }
        
        
    }
    
    
    # if ties, we set the color to white
    which.contrib$Contrib = apply(which.contrib, 1, function(x){if (length(which(x)) > 1){
        return(col.ties)
    } else { # otherwise we use legend color provided
        ##return(color.mixo(1:nlevels(Y))[which(x)])
        return(legend.color[1:nlevels(Y)][which(x)])
    }})
    
    # we also add an output column indicating the group that is max
    which.contrib$GroupContrib = apply(which.contrib[, 1:(nlevels(Y))], 1, function(x){if (length(which(x)) > 1){
        return("tie")
    } else {
        return(levels(Y)[which(x)])
    }})
    
    method.group = do.call(rbind, method.group)
    
    
    #   # run a kruskall wallis test and adjust the pvalues for multiple testing using Benjamini Hochberg
    # this has been removed as the type of stat test is particular to the data
    #   kruskall.pval = apply(X[, name.selected.var, drop = FALSE], 2, function(x){kruskal.test(x ~ Y)$p.value})
    #   kruskall.pval.adj = p.adjust(kruskall.pval, method = 'BH')
    
    contrib = data.frame(method.group, which.contrib, importance = value.selected.var)
    # End contribution calculation
    
    # =====================================
    # then determine the colors/groups matching max contribution
    # =======================================
    
    # when working with sparse counts in particular and using the median to measure contribution
    # ties to determine the contribution of a variable may happen, in that case remove them, otherwise they are showns as blank
    if (show.ties == FALSE){
        contrib = contrib[!contrib$Contrib %in% col.ties, ]
        colnames.X=rownames(contrib)
    }
    
    #contrib = contrib[1 : min(nrow(contrib), ndisplay), , drop = FALSE]
    
    #display barplot with names of variables
    
    # ==================================
    # represent the barplot
    # ==================================
    #added condition if all we need is the contribution stats
    if(plot){
        # in case the user does not want to show the legend, then margins can be reduced
        if(legend){
            layout(matrix(c(1, 2), 1, 2, byrow = TRUE), c(0.7, 0.3), TRUE)
        }else{
            layout(matrix(c(1, 2), 1, 2, byrow = TRUE), c(0.7, 0.1), TRUE)
        }
        # barplot with contributions
        par(mar = c(5, min(9, max(sapply(colnames.X, nchar))), 4, 0))
        mp = barplot(contrib$importance, horiz = T, las = 1, col = contrib$Contrib, axisnames = TRUE, names.arg = colnames.X, #names.arg = row.names(contrib),
        cex.names = cex.name, cex.axis = 0.7, beside = TRUE,border=NA)
        if(is.null(main)){
            if(class.object[1] %in% c("sgcca","sgccda","spls"))
            title(paste0('Contribution on comp ', comp, "\nBlock '", names.block,"'"))
            else
            title(paste('Contribution on comp', comp))
        }else{
            title(paste(main))
        }
        
        # legend
        if(legend){
            par(mar = c(5, 0, 4, 3) + 0.1)
            plot(1,1, type = "n", axes = FALSE, ann = FALSE)
            legend(0.8, 1, col = legend.color[1:nlevels(Y)], legend = levels(Y), pch = 19, 
            title = paste(legend.title), 
            cex = cex.legend)
        }
        par(mar=c(5.1, 4.1, 4.1, 2.1))
        par(mfrow=c(1,1))
        
    } # end if plot
    
    # ===================================
    # return the contribution matrix
    # ===================================
    return(invisible(list(contrib = contrib)))
}
