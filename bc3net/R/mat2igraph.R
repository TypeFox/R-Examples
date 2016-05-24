.mat2igraph=function(bnet,weighted=TRUE){

       nam <- rownames(bnet)      
       edges <- which(bnet>0, arr.ind=TRUE)

       if(nrow(edges)>0){
         
         if(weighted==TRUE){
            weight=sapply(1:nrow(edges),function(i){
              a=edges[i,1]
              b=edges[i,2]
              w=bnet[a,b]
              w
            })
         }else{
           weight=rep(1,nrow(edges))
         }
        
         edges=data.frame(edgesA=nam[edges[,1]],edgesB=nam[edges[,2]],weight)
         bnet=graph.data.frame(edges[,1:2],directed=FALSE)
         E(bnet)$weight=edges[,3]
      
         add=nam[!nam%in%V(bnet)$name]
         if(length(add)>0){
            bnet=add.vertices(bnet,nv=length(add))
            V(bnet)$name[is.na(V(bnet)$name)]=add
         }
      }else{
         bnet=graph(c(),n=length(nam), directed=FALSE)
         V(bnet)$name=nam
      }  

       
      return(bnet)
}  
