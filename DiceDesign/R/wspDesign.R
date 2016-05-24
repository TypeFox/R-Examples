#####ALGORITHM of WOOTON SERGENT PHAN-TAN-LUU (WSP)#####


####wspDesign#####

#---------------------------------------------------------------------------|
#args :  design     : matrix (n,d) containing n design points               |
#                     d is the dimension                                    |
#        dmin       : minimal bound for mindist final design                |    
#out :                a list containing the input arguments plus            |
#                   the space filling design according to mindist criterion |
#---------------------------------------------------------------------------|


wspDesign <- function(design,dmin){           #  dmin :must be given

  m <- design
  n <-nrow(m)                     # Number of points in m.
  d <-ncol(m)                     # Dimension of m.
  A<-rbind(m,rep(0.5,d))          # add center point
  D<-as.matrix(dist(A))
  diag(D)<-rep(Inf,n+1)            # Matrix of distances between all pairs of point (both design points and center)

  base<-as.numeric(which.min(D[n+1,1:n])) # the nearest design point to center point
  D=D[1:n,-(n+1)] 
                                 # Matrix of distances between all pairs of design points
 
  


  points<-base                    # vectors containing base points and suppressed points as well
  x<-base                         # vector containing only base points
  while (length(points)<n){      
              
              co<-as.numeric(which(D[base,]<dmin))
              if(length(co)>0)
               {D[,co]<-rep(Inf,n)            
                D[co,]<-rep(Inf,n)
                points<-c(points,co)}   #suppressed points
                    
                
              new_base<-which.min(D[base,])  #new base point
              points<-c(points,base)
              D[base,]<-D[,base]<-rep(Inf,n)
              base<-new_base
              x<-c(x,base)}
              
   Mres <- m[x,]  # Final matrix    #final space-filling design
  
  return(list(InitialDesign=design,dmin=dmin,design=Mres))
  
   return(Mres) 
}


