#####MINIMAL SPANNING TREE (ALGORITHM of PRIM(1957))#####

####mstCriteria#####

#---------------------------------------------------------------------------|
#args :  design     : the design                                            |
#        plot2d     : default="FALSE",  show 2D projections of the mst      |
#out     l     : list with edges and both mean and standard value of the    |
#                lengthes                                                   |
#---------------------------------------------------------------------------|



mstCriteria <- function(design,plot2d="FALSE")

{
  m <- design
  n<-nrow(m)             #number of points 
  D<-as.matrix(dist(m))  #distances matrix two points two
  diag(D)<-rep(Inf,n)    
  x<-rep(0,n-1)          #Initialization of the vector such that each coordinate is an edge length 
  Peak<-1                #First peak of the tree
  p<-Peak
  Tree<-list(c())
  Tree<-rep(Tree,n)      #Initialization of a void list
 
  
#-----------------------------------------------------------|
  if(plot2d)
  {plot(m,pch=19,col="blue",xlab="",ylab="")}        
#-----------------------------------------------------------|


  for( i in 1:(n-1))               
                     
    { 

       if(length(Peak)==1)
              {
                 new_p<-which.min(as.numeric(D[p,]))
                 x[i]<-D[p,new_p]
                 D[new_p,Peak]=D[Peak,new_p]<-Inf
                 Peak<-c(Peak,new_p)
               }


       else{
                 A<-D[Peak,]
                 w<-which.min(A)
                 reste<-w%%nrow(A)
                 if(reste!=0) 

                    {p<-Peak[reste]
                    new_p<-w %/% nrow(A)+1}

               else
                    {p<-Peak[nrow(A)]
                    new_p<-w %/% nrow(A)}
       
       x[i]<-D[p,new_p]
       D[new_p,Peak]<-D[Peak,new_p]<-Inf
       Peak<-c(Peak,new_p)
 
            }    

       Tree[[p]]<-c(Tree[[p]],new_p)
       Tree[[new_p]]<-c(Tree[[new_p]],p)                             
       if (plot2d=="TRUE")
          {lines(rbind(m[p,],m[new_p,]),col='red')}
    }
       return(list(tree=Tree,stats=c(mean(x),sqrt(var(x)))))
}

