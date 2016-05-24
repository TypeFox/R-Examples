# file OrdinalLogisticBiplot/R/PlotClusters.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

PlotClusters<- function(A, Groups = ones(c(nrow(A), 1)),
                         colors = NULL, chulls = TRUE,centers = TRUE,
                         ConfidentLevel=0.95) {

    if (!is.factor(Groups)){ 
      stop("The group variable must be a factor")
    }
    
    levellab = levels(Groups)
    g = length(levels(Groups))
    if (is.null(colors)) {
        if (g > 1) {
          palette(rainbow(g))
          colors = palette()
        }else{
          rep("red", g)
        }
    }else{
       if (length(colors) == 1){
         colors = rep(colors[1], g)
       }else if(length(colors) < g){
               colors = rep(colors[1], g)
             }else if(length(colors) == nrow(A)){
                     colors = c(1:g)
                   }else{
                     colors = colors[1:g]
                   }    
    }
    
    Sizes = zeros(c(g, 1))
    Means = zeros(c(g, 2))
    X = list()
    for (i in 1:g) {
      print(i)
      X[[i]] = A[which(Groups == levellab[i]), ]
      if(!is.null(dim(X[[i]]))){
          Sizes[i] = dim(X[[i]])[1]
          Means[i, ] = apply(X[[i]], 2, mean)
      }else{
          Sizes[i] = 1
          Means[i, ] = X[[i]]
      }      
    }
    
    if(!is.null(ConfidentLevel)){
      if((ConfidentLevel > 0) & (ConfidentLevel < 1)){
         distances = matrix(0,nrow(A),1)
          for(i in 1:nrow(A)){
            distances[i] = sqrt((A[i,1]-Means[Groups[i],1])^2 + (A[i,2]-Means[Groups[i],2])^2)
          }
          cutPoint = quantile(distances,ConfidentLevel)
          inConfidentLevel = matrix(0,nrow(A),1)
          inConfidentLevel = as.numeric(distances < cutPoint)
          ACL = A[which(inConfidentLevel == 1), ]
          GroupsCL = Groups[which(as.factor(inConfidentLevel) == 1)] 
          SizesCI = zeros(c(g, 1))
          MeansCI = zeros(c(g, 2))
          XCI = list()
          for (i in 1:g) {
            XCI[[i]] = ACL[which(GroupsCL == levellab[i]), ]
            SizesCI[i] = dim(XCI[[i]])[1]
            MeansCI[i, ] = apply(XCI[[i]], 2, mean)
          }
          X = XCI
          Groups = GroupsCL
          Means = MeansCI
      }
    }
    
    
    if(chulls){
      for (i in 1:g) {
        hpts = chull(X[[i]])
        hpts <- c(hpts, hpts[1])
        if(!is.null(nrow(X[[i]]))){
            lines(X[[i]][hpts, ], col = colors[i])
        }
      }
    }
    if(centers){
      for (i in 1:g) {
        points(Means[i,1],Means[i,2],pch=19,col=colors[i],cex=1)
        text(Means[i,1],Means[i,2],levels(Groups)[i],pos=2,offset=0.2,cex=0.5,col="black")
      }
    }

}
