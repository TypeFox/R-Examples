##' @title GO enrichment annotation
##' 
##' @description GO enrichment annotation of the functional molecules or networks.
##' 
##' @param data An igraph object or vector of vertex names.
##' @param onto GO categories, three possible values are \code{MF} for GO function, \code{BP} for GO process, and \code{CC} for GO componet
##' @param pvalue Significant level. Default value is \code{0.05}.
##' @return Numeric vector, pvalue with the length as the size of data.
##' @references The Gene Ontology Consortium (January 2008). The Gene Ontology project in 2008. Nucleic Acids Res. 36 (Database issue): D440?C4.
##' @export
##' @examples 
##' entrez<-data.frame(c("121549","51160","83878","11338","196477","9319","608","7015"))
##' net<-construction(input=entrez,hierarchy=0,species="human",db="Biogrid",ID.type="Entrez Gene")
##' res<-enrichment.annotation(net,pvalue=0.05,onto="CC")

enrichment.annotation<-function(data,onto=c("MF","BP","CC"),pvalue=0.05)
{
	onto<-match.arg(onto)
	if(is.igraph(data)){
    gene<-V(data)$name
	}else if(is.character(data)){
    gene<-data
	}
  gene <- gene[!is.na(gene)]
	if(onto=="MF"){
    GO_function<-NULL
		data("GO_function", envir = environment())
		res<-enrichment0(gene,GOset=GO_function[,c("GeneID","GO_ID","GO_term")],pvalue)
	}else if(onto=="BP"){
    GO_process<-NULL
		data("GO_process", envir = environment())
		res<-enrichment0(gene,GOset=GO_process[,c("GeneID","GO_ID","GO_term")],pvalue)
	}else if(onto=="CC"){
    GO_component<-NULL
		data("GO_component", envir = environment())
		res<-enrichment0(gene,GOset=GO_component[,c("GeneID","GO_ID","GO_term")],pvalue)
	}
	return(res)
}

##  calculate the pvalue
enrichment0<-function(gene,GOset,pvalue)
{
	n0<-length(gene)
	n1<-nrow(GOset)
	gene2go<-GOset[!is.na(match(GOset$GeneID,gene)),]
	rownames(gene2go)<-NULL
	GOids<-table(gene2go$GO_ID)
	GOidset<-table(GOset[!is.na(match(GOset$GO_ID,names(GOids))),"GO_ID"])
	p_value<-1-phyper(GOids-1,GOidset,n1-GOidset,n0)
	p_value<-sort(p_value[p_value<=pvalue])
	x<- data.frame(matrix(0,length(names(p_value)),6))
	colnames(x)<-c("GO_ID","GO_term","GeneID","p.value")
	x$GO_ID<-names(p_value)
	x$GO_term<-gene2go$GO_term[match(names(p_value),gene2go$GO_ID,nomatch=0)]
	id<-lapply(names(p_value),function(GO){
		which(!is.na(match(gene2go$GO_ID,GO)))})
	x$Evidence<-unlist(lapply(id,function(i){
		paste(unique(gene2go$Evidence[i]),collapse=",")}))
	x$GeneID<-unlist(lapply(id,function(i){
		paste(unique(gene2go$GeneID[i]),collapse=",")}))
	x$p.value<-p_value
	x$PubMed<-unlist(lapply(id,function(i){
		paste(unique(gene2go$PubMed[i]),collapse=",")}))
	return(x)
}
