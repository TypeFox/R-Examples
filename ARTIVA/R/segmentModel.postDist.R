segmentModel.postDist <-
function(CPnumberPostDist,CPpositionPostDist,TFnumber,CPsamples,coefSamples,segMinLength, edgesThreshold=0.5, TFnames=NULL, CPpos=NULL){

  if((dim(coefSamples)[2] %% (TFnumber+1))!=0)stop("Parameter TFnumber does not correspond to the coefSamples dimension")
  edgesPostDist=NULL
  coeffMean=NULL
  
  if(is.null(CPpos)){
    nbCP=which(CPnumberPostDist==max(CPnumberPostDist))-1
  
    if(nbCP>0) {
      ##CPpos=c(CPsamples[1,1],sort(order(CPpositionPostDist,decreasing=TRUE)[1:nbCP]),max(CPsamples[1,]))
    
      totalPos=length(CPpositionPostDist)
      nbCPtoFind=nbCP
      CPtemp=c(CPsamples[1,1],totalPos)
      orderCPtemp=order(CPpositionPostDist,decreasing=TRUE)
    
      while(nbCPtoFind>0){
        while(sum(abs(CPtemp-orderCPtemp[1])<segMinLength)>0){
          orderCPtemp=orderCPtemp[-1]
        }
        CPtemp=c(CPtemp,orderCPtemp[1])
        orderCPtemp=orderCPtemp[-1]
        nbCPtoFind=nbCPtoFind-1
      }
      CPpos=sort(CPtemp)
      
    }else{
      CPpos=c(CPsamples[1,1],max(CPsamples[1,]))	
    }
  }
     
  ## for each segment:
  for(seg in 1:(length(CPpos)-1)){
    start = CPpos[seg]
    end = CPpos[seg+1]
    ## search for index of iterations where phase is encountered
    segRow = which(apply(CPsamples, 1, hasSeg, start, end))
    
    if(length(segRow)==0){
      print(paste(" !!!! Problem segment in" , seg, "starting from",start, "to", end,": no such segment in the CPsamples !!!!"))
      print("Please try function runARTIVA again with more iterations")
      edgesPostDist=rbind(edgesPostDist,array(NA,TFnumber))
      coeffMean=rbind(coeffMean,array(NA,TFnumber))
    }else{
      ## positions of start at each iteration
      segCol = apply(CPsamples[segRow,], 1, segPos, start)
      
      ## compute edges distribution
      edgesSamples=t(t(coefSamples!=0)*(1:(TFnumber+1)))

      ## find models of corresponding phase
      ## (use of unique index of matrices - by columns)
      ## term1 = vertical vector containing the column number of the begining of the selected phase, with length=length(segCol)
      ## term2 = horizontal vector of 1, with length= TFnumber
      ## term1%*%term2 = matrix of size [TFnumber x  length(segCol)]
      ## term3 = matrix of size [TFnumber x  length(segCol)] with 1st column: 1, 2nd column : 2, etc...
      edgesCols=matrix((segCol-1)*(TFnumber+1),length(segCol),1)%*%matrix(1,1,TFnumber)+t(matrix(1:TFnumber,TFnumber,length(segCol)))

      ## list of the present edges 
      allEdges = edgesSamples[c((edgesCols-1)*nrow(edgesSamples)+segRow)]

      ## proportion of the edges ([-1] to remove 0 which indicates the absence of an edge)
      postDist=round(table(allEdges)[-1]/length(segRow),4)
      tmp=array(0, TFnumber)
      tmp[which(1:TFnumber %in% names(postDist))]=postDist
      
      edgesPostDist=rbind(edgesPostDist,tmp)

      ## edges coefficients estimation
      allCoef=coefSamples[c((edgesCols-1)*nrow(coefSamples)+segRow)]
      edgesSelection=which(tmp>edgesThreshold)
      tmpCoef=array(0, TFnumber)
 
      for (e in edgesSelection){
        tmpCoef[e]=mean(allCoef[which(allEdges==e)])
      }
      coeffMean=rbind(coeffMean, round(tmpCoef,5))
      ##
    }
  }
  
  row.names(edgesPostDist)=paste("seg", c(1:(length(CPpos)-1)), sep="")
  row.names(coeffMean)=paste("seg", c(1:(length(CPpos)-1)), sep="")
 
  
  return(list(CPpos=CPpos,edgesPostDist=edgesPostDist, edgesCoeff=coeffMean))
  
}
