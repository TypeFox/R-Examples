getInteGraphList<-function(graphList,relations){
   # integrate miRNAs into pathways
    InteGraphList<-graphList
    if(length(graphList)!=0){
	for(t in 1:length(graphList)){
        if(vcount(graphList[[t]])>0){
		V(graphList[[t]])$name<-V(graphList[[t]])$names
	    addedpair<-data.frame();
        # convert KO  to gene symbol in patwhays
        for(j in 1:vcount(graphList[[t]])){
		   
           nodes<-unlist(strsplit(V(graphList[[t]])$name[j]," "))
	       #if(graphList[[t]]$org=="ko"){
	      # node_symbol<-toupper(getSymbolFromKO(nodes))
	       #}else{
	       node_symbol<-toupper(getSymbolFromGene(nodes))
	       #}
	   
          # match the symbol for pathway genes with the targets symbol for mirna.
	       matched_id<-which(toupper(as.character(relations[,2]))%in%node_symbol)
	       if(length(matched_id)>0){
		   hit<-relations[matched_id,]
	       matched_mir<-unique(as.character(hit[,1]))
	       matched_len<-length(matched_mir);matched_gene<-c()
	       for(m in 1:matched_len){
	       matched_gene<-c(matched_gene,paste(unique(as.character(hit[which(hit[,1]%in%matched_mir[m]),2])),collapse=";"))
	       }
	       value<-unname(as.matrix(data.frame(j,paste(nodes,collapse=";"),matched_mir,matched_gene)))
		   addedpair<-rbind(addedpair,value)
		  }
		 #
       }
   ###################################
       addedpair<-as.matrix(addedpair)
	   if(dim(addedpair)[1]>1){
	   addedmirna<-unique(addedpair[,3])
	   for(tt in 1:length(addedmirna)){
	   InteGraphList[[t]]<-add.vertices(InteGraphList[[t]],length(addedmirna[tt]),id=addedmirna[tt],name=addedmirna[[tt]],
	                            names=addedmirna[tt],type="lncRNA",reaction="target",graphics_type="sphere",graphics_width="46",
                                graphics_height="17",graphics_name=addedmirna[tt],graphics_fgcolor="#000000",
								graphics_bgcolor="#FF2AA9",graphics_coords="unknow",graphics_x=as.character(runif(1,min(as.numeric(V(InteGraphList[[t]])$graphics_x)),
								max(as.numeric(V(InteGraphList[[t]])$graphics_x)))),
								graphics_y=as.character(runif(1,min(as.numeric(V(InteGraphList[[t]])$graphics_y)),max(as.numeric(V(InteGraphList[[t]])$graphics_y)))),link="unknow"
                                )
								}
	   for(i in 1:dim(addedpair)[1]){
	   InteGraphList[[t]]<-add.edges(InteGraphList[[t]],as.numeric(c(match(addedpair[i,3],V(InteGraphList[[t]])$names),addedpair[i,1])),directed=is.directed(graphList[[t]]))
  	   }
       }
	  }
	 }
   }
  
   return(InteGraphList)
}
