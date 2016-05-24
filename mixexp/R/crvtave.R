### function to get extreme vertices and centroids ###
crvtave <-function(ndm,conmx)
 {  # Inputs:
    # ndm   - highest order of centroids requested
    # conmx - matrix of constraints
    #
rc<-dim(conmx) #dimension of constraint matrix
ncon2<-rc[1]   #number of constraints
nvrr<-rc[2]-1  # number of mixture variables
rtheta2<-conmx
dim(rtheta2)<-ncon2*rc[2] #vector of columns of conmx
 if(0 > ndm | ndm > (nvrr-1)) 
   stop(" The maximum order of centroid requested must be between ", 0," and ", (nvrr-1),"\n")  
eflag<-Eflags(ndm,nvrr,ncon2,rtheta2) 
 if(eflag[1]<0)
   stop(" There are inconsistent constraints -- revise constraint matrix and rerun","\n")
 if(eflag[2]<0)
   stop(" Too many vertices, this function only works when (#vertices + #centroids )<=1000","\n")
 if(eflag[3]<0)
   stop(" ncm must be between 0 and NCON when calling ALLNR","\n")
 if(eflag[4]<0)
   stop(" Too many centroids, this function only works when (#vertices + #centroids )<=1000","\n") 
nvrtr<-Nrows(ndm,nvrr,ncon2,rtheta2)
v<-Vertcen(ndm,nvrr,ncon2,rtheta2)
v<-v[1:(nvrtr*(nvrr+1))]
vtcn<-matrix(v,ncol=(nvrr+1))
# make column labels for output matrix
dimnames(vtcn)<-list(NULL,c(paste("x",1:nvrr,sep=""),"dimen"))
vtcn
 } 
### End of function to get extreme vertices and centroids ###    
