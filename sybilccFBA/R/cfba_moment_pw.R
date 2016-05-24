cfba_moment_pw <- function (model,mod2=NULL, Kcat,MW=NULL,selected_rxns=NULL,verboseMode=2,
objVal=NULL,RHS=NULL,solver=SYBIL_SETTINGS("SOLVER"),medval=NULL){
#MOMENT pairwise OR like MATLAB implementation
#model: modelorg object
#mod2: modelorg with only irreversible reactions.
#        It can be sent  to save time of recalculating it with each call.
#Kcat: kcat values in unit 1/S. Contains three slots: reaction name,direction(dirxn),value(val)
#MW: list of molecular weights of all genes, using readfaa.r, in units g/mol
#selected_rxns: optional parameter to apply to subset of reactions
#verboseMode: defines the level of output messages
#RHS: the budget C, for EColi 0.27
#objVal:when not null the problem will be to find the minimum budget that give the specified obective value(biomass)
#solver: cplexAPI or glpkAPI
#medval: median of Kcat values , used for missing values

#return problem object containing moment model with the given kcat and molecular weights.

#require(sybil)
geneCol=NULL;#column in problem coresponding to given gene

if(length(selected_rxns)==0)
    {selected_rxns=react_id(mod2)} #exclude s0001 in 25 rxns +4 rxns

	if(length(mod2)==0){
		mod2=mod2irrev(model)
	}
prob <- sysBiolAlg(mod2, algorithm = "fba",solver=solver)
n=react_num(mod2)
m=met_num(mod2)

Kcat[,2]<-as.numeric(Kcat[,2])*60*60 ;# convert to /hr
MW[,2]=as.numeric(MW[,2]); #g/mol instead of g/mmol
#medval=3600*22.6#median(Kcat[,2])#
if(length(medval)==0) medval=median(Kcat[,2])
#print(medval)

colid=n+1
rowind=m+1
aux_id=1
#Add variables for all genes
	for( g in allGenes(mod2)){
		if(g!="s0001"){
					geneCol=rbind(geneCol,cbind(gene=g,Col=colid))
					addCols(lp = problem(prob),1)
					changeColsBnds(lp = problem(prob),colid,lb=0,ub=1000)
					colid=colid+1;			
				}}
	
	for (r in selected_rxns){
	    vind=which(react_id(mod2)==r)
		
		if(substring(r,nchar(r)-1)=='_b'){v_rxn=substring(r,1,nchar(r)-2);v_dirxn=-1
		}else{if(substring(r,nchar(r)-1)=='_f'){v_rxn=substring(r,1,nchar(r)-2);v_dirxn=1
		}else{v_rxn=r;v_dirxn=1}
		}
		if(v_rxn %in% Kcat[Kcat[,"dirxn"]==v_dirxn,1]){
			kval=as.numeric(Kcat[Kcat[,"dirxn"]==v_dirxn & Kcat[,1]==v_rxn,"val"]);
		}else{kval=medval;}
	    
		rl=gpr(mod2)[react_id(mod2)==r]
		if(rl=="s0001" || rl=="") next;
# remove s0001 gene from rules in iAF model
		if (rl=='( b0875  or  s0001 )') rl='b0875';
		if (rl=='( b1377  or  b0929  or  b2215  or  b0241  or  s0001  or  b3875  or  b1319  or  b0957 )')
			 rl='( b1377  or  b0929  or  b2215  or  b0241 or  b3875  or  b1319  or  b0957 )';
		if (rl=='( s0001  or  b0451 )') rl='(b0451)';
		if (rl=='( s0001  or  b3927 )') rl='( b3927 )';
		#ijo model
	if(rl=="(b3927 or s0001)") rl="b3927";
	if(rl=="(s0001 or b0957 or b3875 or b2215 or b0241 or b1319 or b1377 or b0929)" )
	     rl="(b0957 or b3875 or b2215 or b0241 or b1319 or b1377 or b0929)"
    if(rl=="(s0001 or b0875)") rl="b0875";
	if(rl=="(b0451 or s0001)") rl="b0451";
	
  print(r)
	rl=gsub("\\)"," ) ",rl)# 
	rl=gsub("\\("," ( ",rl)# 
		# to be sure that 'and' will be a whole word and not part of any gene name
	pr=lapply(strsplit(unlist(strsplit(rl," or "))," and "),function(x) gsub("[() ]","",x))
	#print(pr)
	#print(length(pr[[1]]))
	if( length(pr)==1) {# no OR (only one term) 
		if(length(pr[[1]])==1){
#########################################1-one Gene----------------
	#			print("one Gene")
				gene=pr[[1]]
			#get the corresponding col of the gene
			colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
	#print(gene);print(colind);	
			addRowsToProb(lp = problem(prob),
				  i = rowind ,              type = "U",
				  lb = 0,              ub = 0,
				  cind = list(c(vind,colind)),
				  nzval = list(c(1,-1*kval))
				  ,rnames = r
				 )
				 rowind=rowind+1;
    }else{#straight ANDs
#########################################2-straight ANDs----------------
	#	print("straight ANDs")
		#add aux variable
		addCols(lp = problem(prob),1)
		aux_id=colid
		colid=colid+1;			
		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(vind,aux_id)),
              nzval = list(c(1,-1*kval))
			  ,rnames = r
			 )
			 rowind=rowind+1; 
		for(gene in pr[[1]]){			
		  	colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
			addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(aux_id,colind)),
              nzval = list(c(1,-1))
			#  ,rnames = paster
			 )
			 rowind=rowind+1;   
		}
		}
		}else{#OR/AND
#########################################3-straight OR----------------
		#add row for the rxn
	#	print("OR/AND")
        row_vals=rep(0,colid+length(pr))
		#row_vals[vind]=1
		 for( p in 1:length(pr)){#ORed gene list
           if( length(pr[[p]])==1 ){
				gene=pr[[p]]
				colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
		#changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=colind,val=-1*kval)
				row_vals[colind]=-1*kval
			}else{#AND
		  #add auxiliary var then add rows
#########################################4-mixed OR and AND (sum of products)----------------		  
	#	  print("aux ands")
			addCols(lp = problem(prob),1)
			  aux_id=colid
			colid=colid+1;			
          
			#changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=aux_id,val=-1*kval)	
			row_vals[aux_id]=-1*kval
			for(g in pr[[p]]){
				colind=as.numeric(geneCol[geneCol[,"gene"]==g,"Col"])
				
				 addRowsToProb(lp = problem(prob),
				  i = rowind ,              type = "U",
				  lb = 0,              ub = 0,
				  cind = list(c(aux_id,colind)),          nzval = list(c(1,-1))
				  ,rnames = paste("Aux_",aux_id,sep="")
				 )
							 rowind=rowind+1;
				}
			}
     #     print(p)
		  }		
		  
			
			if(length(pr)==2)
				{row_vals[vind]=1;
				addRowsToProb(lp = problem(prob),
				  i = rowind ,              type = "U",
				  lb = 0,              ub = 0,
				cind = list(which(row_vals!=0)),        nzval = list(row_vals[which(row_vals!=0)])
					,rnames = r
				 )		  
				 rowind=rowind+1;  
			}else{#add rows for pairs
				 #Add aux variable,define row: 2 genes -> 1 var ; 1 gene+ var ->var2, last: var+gene+rxn 
				 #Example R_AMPTASEPG:i in (1960,4337,4336) 3 rows
		#		print("multiple OR")			
			genes=which(row_vals!=0)
		#		 print(genes)
		#		 print(row_vals[genes])
			addCols(lp = problem(prob),1)
			  aux_id=colid
			colid=colid+1;		
		
          	 addRowsToProb(lp = problem(prob),
				  i = rowind ,              type = "U",
				  lb = 0,              ub = 0,
				cind = list(c(genes[1],genes[2],aux_id)),        nzval = list(c(-1,-1,1))
					,rnames = paste("AuxOR_",r,sep="")
				)
				rowind=rowind+1;  
				 o_aux_id=aux_id
		#		 print(o_aux_id)
				 if(length(genes)>3){
				 for(i in (3:(length(genes)-1))){															
					addCols(lp = problem(prob),1)
					n_aux_id=colid
					colid=colid+1;			
					addRowsToProb(lp = problem(prob),
						  i = rowind ,              type = "U",
						  lb = 0,              ub = 0,
						cind = list(c(genes[i],o_aux_id,n_aux_id)),        nzval = list(c(-1,-1,1))
							,rnames = r
						)
						rowind=rowind+1;  
						 o_aux_id=n_aux_id
				 }}
				 	addRowsToProb(lp = problem(prob),
						  i = rowind ,              type = "U",
						  lb = 0,              ub = 0,
						cind = list(c(genes[length(genes)],o_aux_id,vind)),        nzval = list(c(-1*kval,-1*kval,1))
							,rnames = r
						)
						rowind=rowind+1;  
				
				 }
			 }
			 }
###########################
if(length(MW )>0 ){
   lnz=NULL
   lcind=NULL
   #print(colid)
   for(v in c((n+1):(colid-1))){
		# print(v);
		changeColsBnds(lp = problem(prob),v,lb=0,ub=1000)
	}
   for(g in geneCol[,"gene"]) {
		colind=as.numeric(geneCol[geneCol[,"gene"]==g,"Col"])
		#changeColsBnds(lp = problem(prob),colind,lb=0,ub=1000)
		if(g %in% MW[,1]){# consider only genes with computed MW.
			lnz=cbind(lnz,MW[MW[,1]==g,2])
	#		print(MW[MW[,1]==g,2])
			lcind=cbind(lcind,colind)
			}
		}

if(length(objVal)>0){	#test value of upper bound
	objc=getObjCoefs(problem(prob),j=c(1:n))
	addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "L",
              lb = objVal,              ub = 0,
              cind = c(1:n),          nzval = objc
			  ,rnames = paste("ObjC_",rowind,sep="")
			 )
						 rowind=rowind+1;  
	changeObjCoefs(lp = problem(prob),j=c(1:n),obj_coef=rep(0,n))
	#print("old obj..")
	#set new obj function
	changeObjCoefs(lp = problem(prob),j=as.numeric(lcind),obj_coef=lnz)
	setObjDir(lp = problem(prob),"min")
}else{	# Add solvant constraint	
   addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = RHS,#C denotes the total weight of proteins,
              cind = list(lcind),          nzval = list(lnz)
			  ,rnames = "MW"
			 )
	 }
}#MW
	
	##############################################	
	if (verboseMode > 2) {                      
				            fname=format(Sys.time(), "momentPW_%Y%m%d_%H%M.lp");
					       print(sprintf("Writing the problem to: %s/%s...",getwd(),fname));
		                	writeProb(lp=problem(prob), fname)
							if(solver=="cplexAPI" && verboseMode > 3){
								Nnz=NULL
								nr=getNumRows(problem(prob))
								lp=problem(prob)
								for(r in c(1:nr)){
								  print(r)
									 rv=cplexAPI::getRowsCPLEX(env = lp@oobj@env, lp = lp@oobj@lp, begin = r-1, end = r-1)
									 #write rv,matind,val
									 for(cl in 1:length(rv$matind)){
										   Nnz=rbind(Nnz,cbind(i=r,j=rv$matind[cl],val=rv$matval[cl]))
										}
								}
								write.csv(file="nnzpw.csv",Nnz)
							}
                }
	############################
	sol=optimizeProb(prob)
	
	return(list(sol,geneCol,prob))
	}
	