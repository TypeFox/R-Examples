# deg: proteins list you want to entrich in Entrez Genes id
# cutoff: is for the significance of the enrichment (ex 0.05)
#category: which set of ontology ("BP","MF","CC")
#FRDcorrection: 1/0 (you do FDR or correction, for the p value. False Discovery Rate)
#DelRedundency: 1/0 wether you want to delete the redundant between significant GOTERMS.
#GOTERMS: is a binary tree, so there can be redundancy in the children


enrichmentGO <- function( deg, cutoff, category, FDRcorrection, DelRedundency ){

  if( (category!="BP") & (category!="CC") & (category!="MF") ){
    print('Error: please input the correct category BP, CC or MF ')
    return()
  }

  if( (FDRcorrection!=0) & (FDRcorrection!=1) ){
    print('Error: please input FDRcorrection, 0 for no FDR correction and 1 for BH correction')
    return()
  }

  if( (DelRedundency!=0) & (DelRedundency!=1) ){
    print('Error: please input DelRedundency, 0 for keep the GO redundency and 1 for delete the GO redundency')
    return()
  }


#Change folder path with respect to where you save the folder
gene2go = utils::read.table("./go.gene2go.general.human", header=F, sep="\t")
gene2go = base::as.matrix(gene2go)
#change path with respect to the folder you are working on
term = utils::read.table("./go.term", header=F, sep="\t")
term2term = base::as.matrix( utils::read.table("./go.graph_path", header=F, sep="\t") )
term2term = term2term[term2term[,1]!=term2term[,2]  ,  ]


gene_in_go = unique( gene2go[,1] )
deg_in_go = intersect( gene_in_go, deg )
if( length(deg_in_go)==0 ){
  print('None of the genes annotate in GO')
  return()
}

nondeg_in_go = setdiff( gene_in_go, deg_in_go )

goid_all = unique( gene2go[,2] )

if(category=="BP"){
  goid_sub = intersect(goid_all, term[ term[,3]=="biological_process" ,1])
}else if(category=="CC"){
  goid_sub = intersect(goid_all, term[ term[,3]=="cellular_component" ,1])
}else{
  goid_sub = intersect(goid_all, term[ term[,3]=="molecular_function" ,1])
}



gene2go_sub = gene2go[ is.element( gene2go[,2],goid_sub ) , ]
num_gene_in_terms = table( factor(gene2go_sub[,2],levels=goid_sub) )
num_gene_in_terms = base::as.matrix(  cbind( as.numeric( names(num_gene_in_terms)) , num_gene_in_terms ) )
num_gene_in_terms = num_gene_in_terms[ order(num_gene_in_terms[,1] ), ]
#num_gene_in_terms = tapply( gene2go_sub[,1], as.factor(gene2go_sub[,2]), length )

deg2go_sub = gene2go_sub[ is.element( gene2go_sub[,1],deg_in_go ) , ]
if( length( dim(deg2go_sub) )==0 ){
  num_deg_in_terms = table( factor(deg2go_sub[2],levels=goid_sub ) )
}else{
  num_deg_in_terms = table( factor(deg2go_sub[,2],levels=goid_sub ) )
}
num_deg_in_terms = base::as.matrix(  cbind( as.numeric( names(num_deg_in_terms)) , num_deg_in_terms ) )
num_deg_in_terms = num_deg_in_terms[ order(num_deg_in_terms[,1] ), ]
#num_deg_in_terms = tapply(  deg2go_sub[,1], as.factor(deg2go_sub[,2]), length )

# num_deg_in_terms = num_deg_in_terms[ order(names(num_deg_in_terms)) ]
# num_gene_in_terms = num_gene_in_terms[ order(names(num_gene_in_terms)) ]

num_deg_in_go = rep( length(deg_in_go), length(goid_sub) )
num_nondeg_in_go = rep( length(nondeg_in_go), length(goid_sub) )

num_matrix = cbind(num_deg_in_terms, num_deg_in_go, num_nondeg_in_go, num_gene_in_terms[,2])

pvalue = apply( num_matrix, 1, function(x) stats::phyper(x[2]-1, x[3], x[4], x[5],lower.tail=F) )


if(FDRcorrection==1){
  qvalue = stats::p.adjust(pvalue,"BH")

  pvalue = qvalue
}

pvalue_matrix = cbind( num_matrix, pvalue  )

sig_goid = pvalue_matrix[ pvalue_matrix[,6] <= cutoff, 1 ]

if( length(sig_goid)==0   ){
  print('No significant GO')
  return()
}

if(DelRedundency==1){
  sig_goid_father = term2term[ is.element(term2term[,2],sig_goid )  , 1 ]
  sig_goid_child = setdiff( sig_goid, sig_goid_father )

  sig_goid = sig_goid_child
}

sig_term = term[ is.element(term[,1],sig_goid)  ,]
sig_go = pvalue_matrix[ is.element( pvalue_matrix[,1], sig_goid ) , ]

if( length(sig_goid)>1 ){
  sig_go = sig_go[ order(sig_go[,1]) ,]
  sig_term = sig_term[ order(sig_term[,1]) ,]
  sig_term_info = cbind(sig_go, sig_term[, c(4,2,5,3)] )
}else{
  sig_term_info = data.frame( c(sig_go, sig_term[c(4,2,5,3)] ) )
}

sig_term_info = sig_term_info[ order(sig_term_info[,6])  ,]
sig_term_info[,4] = sig_term_info[,3] + sig_term_info[,4]
colnames(sig_term_info)=c("GOtermIndex","#OverlappingGenes","#InputGeneList","#TotalGeneInGO","#GeneInGOterm","Pvalues","GOtermAcc","GOtermDescription","TreeDepthOfGOterm","SubOntology")
rownames(sig_term_info)=NULL
return( sig_term_info )

}









