# Copyright (C) 2014 
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Amrit Singh, University of British Columbia, Vancouver.
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
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


predict.pls <- predict.spls <-
function(object, newdata,  ...)
{

    #-- validation des arguments --#
    if (missing(newdata))
    stop("No new data available.")
     
    X = object$X
    Y = object$Y
    q = ncol(Y)
    p = ncol(X)
     
    if (length(dim(newdata)) == 2) {
        if (ncol(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ", p, 
            " or a vector of length = ", p, ".")
    }
     
    if (length(dim(newdata)) == 0) {
        if (length(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ", p, 
            " or a vector of length = ", p, ".")
        dim(newdata) = c(1, p) 
    }
     
    #-- initialisation des matrices --#	
    ncomp = object$ncomp
    a = object$loadings$X
    b = object$loadings$Y
    c = object$mat.c
     
    means.X = attr(X, "scaled:center")
    means.Y = attr(Y, "scaled:center")
    sigma.X = attr(X, "scaled:scale")
    sigma.Y = attr(Y, "scaled:scale")
     
    newdata = as.matrix(newdata)
    ones = matrix(rep(1, nrow(newdata)), ncol = 1)
    ##- coeff de regression 
    B.hat = array(0, dim = c(p, q, ncomp))
    ##- prediction
    Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
    Y.hat2 = array(0, dim = c(nrow(newdata), q, ncomp))
    ##- variates
    t.pred = array(0, dim = c(nrow(newdata), ncomp))
    
    variates.X = object$variates$X
    betay = list()
    
    #-- prediction --#
    for(h in 1:ncomp){
        
        dd= coefficients(lm(Y~variates.X[,1:h,drop=FALSE])) #regression of Y on variates.global.X => =loadings.global.Y at a scale factor
        if(q==1){betay[[h]]=(dd[-1])}
        if(q>=2){betay[[h]]=(dd[-1,])}
        
        W = a[, 1:h,drop=FALSE] %*% solve(t(c[, 1:h,drop=FALSE]) %*% a[, 1:h,drop=FALSE])
        B = W %*% drop(betay[[h]])

        Y.temp=scale(newdata,center=means.X,scale=sigma.X) %*% as.matrix(B) #so far: gives a prediction of Y centered and scaled
        Y.temp2=scale(Y.temp,center=FALSE,scale=1/sigma.Y) #so far: gives a prediction of Y centered, with the right scaling
        Y.temp3=scale(Y.temp2,center=-means.Y,scale=FALSE) #so far: gives a prediction of Y with the right centering and scaling
        
        Y.hat[, , h] = Y.temp3 # we add the variance and the mean of Y used in object to predict
        t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[, h]
        B.hat[, , h] = B
    }  #end h
     
    #-- valeurs sortantes --#
    rownames(t.pred) = rownames(newdata)
    colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
    rownames(Y.hat) = rownames(newdata)
    colnames(Y.hat) = colnames(Y)
     
    return(invisible(list(predict = Y.hat, variates = t.pred, B.hat = B.hat,betay=betay)))
}


# -------------------------- for plsda and splsda ---------------------------------------
predict.plsda <- predict.splsda <-
function(object, newdata, 
         method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), ...)  
{
	#-- validation des arguments --#
    if (missing(newdata))
    stop("No new data available.")
     
    X = object$X
    Y = object$Y 
    Yprim = object$ind.mat   
    q = ncol(Yprim)          
    p = ncol(X)
	
    if (length(dim(newdata)) == 2) {
        if (ncol(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ", p, 
            " or a vector of length = ", p, ".")
    }
     
    if (length(dim(newdata)) == 0) {
        if (length(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ", p, 
            " or a vector of length = ", p, ".")
        dim(newdata) = c(1, p) 
    }
     
    #-- initialisation des matrices --#	
    ncomp = object$ncomp
    a = object$loadings$X
    b = object$loadings$Y
    c = object$mat.c
     
    means.X = attr(X, "scaled:center")
    means.Y = attr(Y, "scaled:center")
    sigma.X = attr(X, "scaled:scale")
    sigma.Y = attr(Y, "scaled:scale")
     
    newdata = as.matrix(newdata)
    ones = matrix(rep(1, nrow(newdata)), ncol = 1)
    ##- coeff de regression 
    B.hat = array(0, dim = c(p, q, ncomp))
    ##- prediction
    Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
    ##- variates
    t.pred = array(0, dim = c(nrow(newdata), ncomp))
    variates.X = object$variates$X
    betay = list()
    
    #-- prediction --#
    for(h in 1:ncomp){
        dd= coefficients(lm(Y~variates.X[,1:h,drop=FALSE])) #regression of Y on variates.global.X => =loadings.global.Y at a scale factor
        if(q==1){betay[[h]]=(dd[-1])}
        if(q>=2){betay[[h]]=(dd[-1,])}

        W = a[, 1:h,drop=FALSE] %*% solve(t(c[, 1:h,drop=FALSE]) %*% a[, 1:h,drop=FALSE])
        B = W %*% drop(betay[[h]])

        Y.temp=scale(newdata,center=means.X,scale=sigma.X) %*% as.matrix(B) #so far: gives a prediction of Y centered and scaled
        Y.temp2=scale(Y.temp,center=FALSE,scale=1/sigma.Y) #so far: gives a prediction of Y centered, with the right scaling
        Y.temp3=scale(Y.temp2,center=-means.Y,scale=FALSE) #so far: gives a prediction of Y with the right centering and scaling

        Y.hat[, , h] = Y.temp3 # we add the variance and the mean of Y used in object to predict
        t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[, h]
        B.hat[, , h] = B
    }  #end h
    
    G = matrix(0, nrow = q, ncol = ncomp)
    cls = list()
    
    for (i in 1:q) {
        if(ncomp > 1) {

            G[i, ] = apply(object$variates$X[Yprim[, i] == 1, , drop = FALSE], 2, mean)
        }
        else {
            G[i, ] = mean(object$variates$X[Yprim[, i] == 1, ])
        }
    }	
		
	# ----    max distance -----------------
	
    if (any(method == "all") || any(method == "max.dist")) {
         
	    function.pred = function(x){
            nr = nrow(x)
            tmp = vector("numeric", nr)
            for(j in 1:nr){
                tmp[j] = (which(x[j, ] == max(x[j, ]))[1])
            }
            return(tmp)
        }
        cls$max.dist = matrix(apply(Y.hat, 3, function.pred), ncol = ncomp)
        colnames(cls$max.dist) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
    }

	# ----    centroids distance -----------------

	if (any(method == "all") || any(method == "centroids.dist")) {
     
        cl = matrix(nrow = nrow(newdata), ncol = ncomp)
         
        centroids.fun = function(x, G, h) {
            q = nrow(G)
            x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
            if (h > 1) {
                d = apply((x - G[, 1:h])^2, 1, sum)
            }
            else {
                d = (x - G[, 1])^2
            }
            cl.id = which.min(d)
        }
	    	
        for (h in 1:ncomp) {
            cl.id = apply(matrix(t.pred[, 1:h], ncol = h), 1, centroids.fun, G = G, h = h)
            cl[, h] = cl.id		
        }
        colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
        cls$centroids.dist = cl
    }	

	# ----    mahalanobis distance -----------------
	
    if (any(method == "all") || any(method == "mahalanobis.dist")) {
     
        cl = matrix(nrow = nrow(newdata), ncol = ncomp)
         
        Sr.fun = function(x, G, Yprim, h) {
            q = nrow(G)
            Xe = Yprim %*% G[, 1:h]
            Xr = object$variates$X[, 1:h] - Xe
            Sr = t(Xr) %*% Xr / nrow(Y)
            Sr.inv = solve(Sr)
            x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
            if (h > 1) {
                mat = (x - G[, 1:h]) %*% Sr.inv %*% t(x - G[, 1:h])
                d = apply(mat^2, 1, sum)
            }
            else {
                d = drop(Sr.inv) * (x - G[, 1])^2
            }
            cl.id = which.min(d)
        }
	     	
        for (h in 1:ncomp) {
            cl.id = apply(matrix(t.pred[, 1:h], ncol = h), 1, Sr.fun, G = G, Yprim = Yprim, h = h)
            cl[, h] = cl.id		
        }
        colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
        cls$mahalanobis.dist = cl
    }
	
    #-- valeurs sortantes --#
    if (any(method == "all")) method = "all"
    rownames(t.pred) = rownames(newdata)
    colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
    rownames(Y.hat) = rownames(newdata)
    colnames(Y.hat) = colnames(Y)
    colnames(G) = paste("dim", c(1:ncomp), sep = " ")
     
    return(invisible(list(predict = Y.hat, variates = t.pred, B.hat = B.hat, 
		                  centroids = G, method = method, class = cls)))
}
