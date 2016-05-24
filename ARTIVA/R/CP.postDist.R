CP.postDist <-
function(CPsamples, burn_in=NULL, segMinLength=2){

  if(is.null(burn_in))burn_in=dim(CPsamples)[1]/4
  
 
    # number of ChangePoints
    CPnumber=apply(CPsamples!=0,1,sum)-2
    CPnumberPostDist=array(0,max(CPnumber)+1)
    iterationsNumber=length(CPnumber)
    for(i in 0:max(CPnumber)){
      CPnumberPostDist[i+1]=sum(CPnumber==i)/iterationsNumber
    }
    
   # CP position posterior distribution
    CPpositionPostDist=array(0, max(CPsamples))
    for(i in (CPsamples[1,1]+1):(max(CPsamples[1,]-1))){
      CPpositionPostDist[i]=sum(CPsamples==i)/iterationsNumber
    }

    # CP position estimation 

  nbCP=which(CPnumberPostDist==max(CPnumberPostDist))-1
  CPpos=CPposition(nbCP,CPpositionPostDist,CPsamples,segMinLength)
 
    return(list(CPnumber=CPnumberPostDist,CPposition=CPpositionPostDist,estimatedCPpos=CPpos, estimatedCPnumber=nbCP ))
}
