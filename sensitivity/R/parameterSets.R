##
## Copyright (c) Joseph Guillaume <josephguillaume@gmail.com>
##
## Based on design of parameterSets function by Felix Andrews in hydromad package
##  hydromad.catchment.org
##

parameterSets<-function(par.ranges,samples,method=c("sobol","innergrid","grid")){
  method=match.arg(method)
  if(is.null(names(par.ranges))) names(par.ranges)=make.names(par.ranges)
  
  switch(method,
         "sobol"={
           ## Sample a sobol sequence
           if (!requireNamespace("randtoolbox", quietly = TRUE)){ 
              stop('The package randtoolbox is missing, but is required to create 
                  a sample with method="sobol"')
           }
           if (requireNamespace("randtoolbox", quietly = TRUE)){ 
             pts <- randtoolbox::sobol(samples,length(par.ranges))
           }
           ## Scale
           for(i in 1:length(par.ranges)) 
             pts[,i]<-pts[,i]*(diff(par.ranges[[i]]))+par.ranges[[i]][1]
           return(pts)
         },
         "innergrid"={
           if(length(samples)==1) samples<-rep(samples,length(par.ranges))
           offsets=sapply(par.ranges,diff)/samples/2
           points=lapply(1:length(par.ranges),
                         function(i) seq(par.ranges[[i]][1]+offsets[i],
                                         par.ranges[[i]][2]-offsets[i],
                                         length.out=samples[[1]]))
           names(points)<-names(par.ranges)
           return(as.matrix(do.call(expand.grid,points)))
         },
         "grid"={
           if(length(samples)==1) samples<-rep(samples,length(par.ranges))
           points=lapply(par.ranges,
                         function(r) seq(r[1],
                                         r[2],
                                         length.out=samples[[1]]))
           names(points)<-names(par.ranges)
           return(as.matrix(do.call(expand.grid,points)))
         }
  )
  
}