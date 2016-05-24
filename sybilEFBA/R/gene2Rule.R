# function to get state of expression of GPR rules from gene expression data

gene2Rule <-function(model,geneExpr,selected_rxns=NULL){
 # GenExpr        a data frame: geneID, expr_val
 # Return
 # ruleExpr: rxn_id,expr_val
 # N.B.: GPR rules are in form Sum-of-products (AND to OR) 
 
 # rules in brief:
 #1-Complexes: min, 2-isoenzymes: sum
 #3-multifunctioning: divide by count
 
#TO DO:
# execlude blocked rxns when evaluating gene2Rule
#nrxn=react_num(model)
	if(length(selected_rxns)==0){
			ugpr=as.vector(unique(gpr(model)))
	}else{
			ugpr=as.vector(unique(gpr(model)[selected_rxns]))		
	}
	ugpr=ugpr[ugpr!=""]
	gprExpr=NULL;
	for (v_rule in ugpr){
	   rl=gsub("\\)"," ) ",v_rule)# 
	   rl=gsub("\\("," ( ",rl)# 
	   
	   pr=lapply(strsplit(unlist(strsplit(rl," or "))," and "),function(x) gsub("[() ]","",x))
	   # rules ar in for Sum-of-Product and only two levels ( any rule can be written in this form if not containing NOT)
			expr_val=0;
			for( p in 1:length(pr)){
				gene_ind=match(pr[[p]],geneExpr$GeneID)#cope with repetitions of genes in complexes which(geneExpr$geneID %in% pr[[p]])
				if(length(gene_ind)<length(pr[[p]])){
					warning(sprintf("Rule %s containing gene names not in geneID list, term no: %d term: %s ",v_rule,p,pr[[p]][1]))
				}else{
						expr_val=expr_val + min(geneExpr[gene_ind,"expr_val"])#mean
				}
			}

		cnt=sum(gpr(model)==v_rule)#account for multi-functioning genes
		gprExpr=rbind(gprExpr,cbind(rxn_id=react_id(model)[gpr(model)==v_rule],expr_val=expr_val/cnt,gpr=v_rule,cnt=rep(cnt,cnt)))
		#print(sprintf("Rule: %s cnt: %d expr: %f",v_rule,cnt,expr_val/cnt)) 
	}
print(is(gprExpr))
return(gprExpr)
}


#######################################################################################################
 
 
 
  