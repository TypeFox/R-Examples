getNodes=function(lowerf,upperf,sampling,error,relative){
 tol=10^-7

 
 if (abs(upperf-lowerf)<tol)
      return(list())
 # put the lower/upper frequency in a node
 lowerNode=getNodesAux(lowerf,sampling,error,"lower",abs(upperf-lowerf),relative);
 upperNode=getNodesAux(upperf,sampling,error,"upper",abs(upperf-lowerf),relative);
 #cat(lowerNode,"\n")
 #cat(upperNode,"\n")
 return(mergeNodes(lowerNode,upperNode,sampling,error,relative))
 
}
 
 mergeNodes=function(ln,un,sampling,error,relative){

 # lowerNode==upperNode
 if((ln[[1]]==un[[1]])&&(ln[[2]]==un[[2]])){
 finalNodes=list(ln[[1]],ln);
 return(finalNodes)
 }
 #upperNode is a child of lowerNode
 if(isChild(un,ln)){
    #new lower node verifying the cover conditions
    newLnode=c(ln[[1]]+1,ln[[2]]*2)
    finalNodes=mergeNodes(newLnode,un,sampling,error,relative);
    return(finalNodes)
 }
 #lowerNode is a child of upperNode
 if(isChild(ln,un)){
    #new upper node verifying the cover conditions
    newUnode=c(un[[1]]+1,un[[2]]*2+1)
    finalNodes=mergeNodes(ln,newUnode,sampling,error,relative);
    return(finalNodes)
 }
 #general case: there is overlap.there may be an uncovered band
 uncoverLowF=sampling/(2^(ln[[1]]+1))*(ln[[2]]+1);
 uncoverHighF=sampling/(2^(un[[1]]+1))*un[[2]];
 mediumNodes=getNodes(uncoverLowF,uncoverHighF,sampling,error,relative);
 if (length(mediumNodes)==0){#there was not an uncovered band
  finalNodes=list(max(ln[[1]],un[[1]]),ln,un)
 }else{

   finalNodes=c(list(max(ln[[1]],un[[1]]),ln),(mediumNodes[2:length(mediumNodes)]),list(un))
 }
   return(finalNodes)
}

#node1 is child of node2
isChild=function(node1,node2){
#if node1 is child of node2-->the level of the node1 must be greater than the level of the node2
if (node1[[1]]>node2[[1]]){
  c1=getC(node1[[1]],node1[[2]])
  c2=getC(node2[[1]],node2[[2]])
  #the parent of node1 at the level of node2 is codified with the first node2[[1]](level of node2) numbers of c1
  parentc1=c1[1:length(c2)]
  auxBool=prod(c2==parentc1)
  ischild=(auxBool==1)
}else{ #node1 is not child of node2
  ischild=FALSE
}
return (ischild)
}