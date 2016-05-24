getEdgeLabel <-
function(graph){
     edge.name<-E(graph)$subtype_name
     edge.value<-E(graph)$subtype_value
     #edge.label<-E(graph)$subtype_value
     edge.label<-rep("",len=length(edge.name))
     for(i in seq(edge.name)){
         edge_i<-unlist(strsplit(edge.name[i],";"))
        if("phosphorylation" %in% edge_i){
             edge.label[i]<-paste("+p",edge.label[i],sep=" ")
        }
        if("dephosphorylation" %in% edge_i){
             edge.label[i]<-paste("-p",edge.label[i],sep=" ")
        }
        if("glycosylation"  %in% edge_i){
             edge.label[i]<-paste("+g",edge.label[i],sep=" ")
        }
        if("ubiquitination"  %in% edge_i){
             edge.label[i]<-paste("+u",edge.label[i],sep=" ")
        }
        if("methylation"  %in% edge_i){
             edge.label[i]<-paste("+m",edge.label[i],sep=" ")
        }
        if("missing interaction"  %in% edge_i){
             edge.label[i]<-paste("/",edge.label[i],sep=" ")
        }
        if("dissociation"  %in% edge_i){
             edge.label[i]<-paste("|",edge.label[i],sep=" ")
        }
        if("binding/association"  %in% edge_i){
             edge.label[i]<-paste("---",edge.label[i],sep=" ")
         }
        if("repression"  %in% edge_i){
             edge.label[i]<-paste("-e-|",edge.label[i],sep=" ")
        }
        if("expression"  %in% edge_i){
             edge.label[i]<-paste("-e->",edge.label[i],sep=" ")
        }
        if("inhibition"  %in% edge_i){
             edge.label[i]<-paste("--|",edge.label[i],sep=" ")
        }
        if("activation"  %in% edge_i){
             edge.label[i]<-paste("-->",edge.label[i],sep=" ")
        }
        if("indirect effect"  %in% edge_i){
             edge.label[i]<-paste("..>",edge.label[i],sep=" ")
        }
        if("state change"  %in% edge_i){
             edge.label[i]<-paste("...",edge.label[i],sep=" ")
        }
        if("compound" %in% edge_i){
             compound<-V(graph)[V(graph)$id==edge.value[i]]$graphics_name
	         if(length(compound)==1){
                 edge.label[i]<-paste(compound,edge.label[i],sep=" ")
	         }    
        }           
    }
    return (edge.label)
}
