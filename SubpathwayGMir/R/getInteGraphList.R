getInteGraphList <-
function(graphList,relations){
    add.vertex.shape("triangle", clip=vertex.shapes("circle")$clip,
                 plot=mytriangle)
    InteGraphList<-graphList
    if(length(graphList)!=0){
	for(t in 1:length(graphList)){
        if(vcount(graphList[[t]])>0){
		V(graphList[[t]])$name<-V(graphList[[t]])$names
	    addedpair<-data.frame();
        for(j in 1:vcount(graphList[[t]])){
		   nodes<-unlist(strsplit(V(graphList[[t]])$name[j]," "))
	       node_symbol<-toupper(getSymbolFromGene(nodes))
	   	   matched_id<-which(toupper(as.character(relations[,2]))%in%node_symbol)
		   if(length(matched_id)>0)
		   {
		   hit<-relations[matched_id,]
		   addedpair <- rbind(addedpair,data.frame(j,hit))
		   }
		 }
	   addedpair<-as.matrix(addedpair)
       if(nrow(addedpair)>0){
	   addedmirna<-unique(addedpair[,2])
	   InteGraphList[[t]]<-add.vertices(InteGraphList[[t]],length(addedmirna),id=addedmirna,name=addedmirna,
	                            names=addedmirna,type="mirna",reaction="target",graphics_type="triangle",graphics_width="72",
                                graphics_height="34",graphics_name=addedmirna,graphics_fgcolor="dimgray",
								graphics_bgcolor="#FF00FF",graphics_coords="unknow",graphics_x="40",
								graphics_y="30",link="unknow"
                                )
	   for(i in 1:nrow(addedpair)){
	   InteGraphList[[t]]<-add.edges(InteGraphList[[t]],as.numeric(c(match(addedpair[i,2],V(InteGraphList[[t]])$names),addedpair[i,1])),directed=is.directed(graphList[[t]]))
	   }
       }
	  }
	 }
   }
  
   return(InteGraphList)
}