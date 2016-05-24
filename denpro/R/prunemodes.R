prunemodes<-function(lst,modenum=1,num=NULL,exmalim=NULL,maxnum=NULL)
{
# prunes from a level set tree "lst" the modes with "num" 
# smallest excess masses 
# or the modes with smaller excess mass than "exmalim"

if (is.null(num)){
    curmodenum<-moodilkm(lst$parent)$lkm
    num<-curmodenum-modenum
}

go.on<-TRUE
nn<-1
while (go.on){

  len<-length(lst$parent)
  child.frekve<-matrix(0,len,1)
  for (i in 1:len){
     if (lst$parent[i]>0) 
     child.frekve[lst$parent[i]]<-child.frekve[lst$parent[i]]+1
  }

  ml<-moodilkm(lst$parent)
  mode.list<-ml$modloc
  roots.of.modes<-matrix(0,length(mode.list),1)
  for (aa in 1:length(mode.list)){
      node<-mode.list[aa]
      while ((lst$parent[node]>0) && (child.frekve[lst$parent[node]]==1)){ 
          node<-lst$parent[node]
      }
      roots.of.modes[aa]<-node
  }

  em<-excmas(lst)
  or<-order(em[roots.of.modes])
  smallest<-ml$modloc[or[1]]
  if (nn==1) exma.of.modes<-em[roots.of.modes]

  node<-smallest
  emsmallest<-em[node]

  if ((is.null(exmalim)) || ((!is.null(exmalim)) && (emsmallest<=exmalim))){

     rem.list<-c(node)
     while ((lst$parent[node]>0) && (child.frekve[lst$parent[node]]==1)){ 
           node<-lst$parent[node]
           rem.list<-c(rem.list,node)
     }

     for (kk in 1:length(rem.list)){
        remo<-rem.list[kk]
        for (ll in 1:length(lst$parent)){
            if (lst$parent[ll]>remo) lst$parent[ll]<-lst$parent[ll]-1
        }
        lst$parent<-lst$parent[-remo]
     }
     lst$level<-lst$level[-rem.list]
     lst$volume<-lst$volume[-rem.list]
     lst$center<-lst$center[,-rem.list]
     lst$distcenter<-lst$distcenter[,-rem.list]
     lst$proba<-lst$proba[-rem.list]
     lst$infopointer<-lst$infopointer[-rem.list]
  }
  else if ((!is.null(exmalim)) && (emsmallest>exmalim)) go.on<-FALSE

  nn<-nn+1
  if ((nn>num) && (is.null(exmalim))) go.on<-FALSE
  if ((!is.null(maxnum)) && (nn>maxnum)) go.on<-FALSE 
}

lst$exma.of.modes<-exma.of.modes

return(lst=lst)
}


