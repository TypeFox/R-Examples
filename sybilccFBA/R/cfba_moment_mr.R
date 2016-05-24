#given  Kcat vector (in 1/Sec), and molecular weights
#1-Add gi vars
#2-add gpr constraints: add f/b for reversible rxns react_rev??
#3-add Density constraint
#... ----------
#mr:multiple reactions
cfba_moment_mr <- function (model,mod2=NULL, Kcat,MW=NULL,selected_rxns=NULL,verboseMode=2,objVal=NULL,
RHS=NULL,solver=SYBIL_SETTINGS("SOLVER"),medval=NULL){
#Multiple OR are not done pairwise but evaluated in one line in contrast to MOMENT original implementation.
#addCols
#Kcat measurement for set of rxns
#MW measuremnt for gene using readfaa.r

Kcat[,2]<-as.numeric(Kcat[,2])*60*60 ;# convert to 1/hr
MW[,2]=as.numeric(MW[,2]); #g/mol instead of g/mmol
if(length(mod2)==0){
		mod2=mod2irrev(model)
	}
if(length(medval)==0) medval=median(Kcat[,2])

if(length(selected_rxns)==0)
    {selected_rxns=react_id(model)[gpr(model)!="" & gpr(model)!="s0001"]} #exclude s0001 in 25 rxns +4 rxns

prob <- sysBiolAlg(mod2, algorithm = "fba",solver=solver)
n=react_num(mod2)
m=met_num(mod2)

colid=n+1
rowind=m+1
aux_id=1
geneCol=NULL;#column in problem corresponding to given gene

# add all variables once
#Add variables for all genes
	for( g in allGenes(model)){
		if(g!="s0001"){
		            cnt=length(grep(g,gpr(model)))
		            geneCol=rbind(geneCol,cbind(gene=g,rxn=NA,Col=colid,cnt))
					addCols(lp = problem(prob),1)
					changeColsBnds(lp = problem(prob),colid,lb=0,ub=1000)
					colid=colid+1;
					if(cnt>1){
						for(r in grep(g,gpr(model)) ){
							geneCol=rbind(geneCol,cbind(gene=g,rxn=react_id(model)[r],Col=colid,cnt))
							addCols(lp = problem(prob),1)
							changeColsBnds(lp = problem(prob),colid,lb=0,ub=1000)
							colid=colid+1;			
							}
							}
				}}

#ng=length(allGenes(model))
for (r in selected_rxns){
	rl=gpr(model)[react_id(model)==r]
	#rl=gsub("s0001",)
	# remove s0001 gene from rules in iAF model as done by original MOMENT paper simulations
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
	
	if (verboseMode > 2) print(r);
	if(r %in% Kcat[Kcat[,"dirxn"]==1,1]){
		kval=as.numeric(Kcat[Kcat[,"dirxn"]==1 & Kcat[,1]==r,"val"]);
		}else{kval=medval;}
	
	if(r %in% Kcat[Kcat[,"dirxn"]==-1,1]){
		rkval=as.numeric(Kcat[Kcat[,"dirxn"]==-1 & Kcat[,1]==r,"val"]);
		}else{rkval=medval;}
	
	if(!react_rev(model)[react_id(model)==r]){ #if ! reversible
			vind=which(react_id(mod2)==r)
			vind_b=0
			lnz=1
		}else{# Add row with median value to rev direction
			#vind=c(which(react_id(mod2)==paste(r,"_f",sep="")),which(react_id(mod2)==paste(r,"_b",sep="")))		
			#lnz=c(1,1)
			vind=which(react_id(mod2)==paste(r,"_f",sep=""))
			vind_b=which(react_id(mod2)==paste(r,"_b",sep=""))		
			lnz=1
		}
	rl=gsub("\\)"," ) ",rl)# 
	rl=gsub("\\("," ( ",rl)# 
		# to be sure that 'and' will be a whole word and not part of any gene name
	pr=lapply(strsplit(unlist(strsplit(rl," or "))," and "),function(x) gsub("[() ]","",x))
	if (verboseMode > 2){
		print(pr)
		print(length(pr[[1]]))
	}
	if( length(pr)==1) {# no OR (only one term) 
	#add row vi<=gi*kcat
		if(length(pr[[1]])==1){
#########################################1-one Gene----------------
		    gene=pr[[1]]
		#get the corresponding col of the gene
		#colind=as.numeric(geneCol[geneCol[,"gene"]==gene & is.na(geneCol[,"rxn"]),"Col"])
		generow=geneCol[geneCol[,"gene"]==gene & is.na(geneCol[,"rxn"]),]
		if(generow[4]==1){#if cnt==1
			colind=as.numeric(generow[3])
		}else{
			colind=as.numeric(geneCol[geneCol[,"gene"]==gene & geneCol[,"rxn"]==r & !is.na(geneCol[,"rxn"]),"Col"])
		}
		if (verboseMode > 2) print(list(vind,colind));
	
		#AddRow Vind - Kcat*colind <=0
		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(vind,colind)),
              nzval = list(c(lnz,-1*kval))
			  ,rnames = r
			 )
			 rowind=rowind+1;
          if(vind_b>0){
		  		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(vind_b,colind)),
              nzval = list(c(lnz,-1*rkval))
			  ,rnames = paste(r,"_b",sep="")
			 )
			 rowind=rowind+1;            
		  }			 
		}else{#straight ANDs
#########################################----2-straight ANDs----------------
		print("straight ANDs")
		#add aux variable
		addCols(lp = problem(prob),1)
		aux_id=colid
		colid=colid+1;			
		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(vind,aux_id)),
              nzval = list(c(lnz,-1*kval))
			  ,rnames = r
			 )
			 rowind=rowind+1;            
		if(react_rev(model)[react_id(model)==r]){
    		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(vind_b,aux_id)),
              nzval = list(c(lnz,-1*rkval))
			  ,rnames = paste(r,"_b",sep="")
			 )
			 rowind=rowind+1;            
		}
		for(gene in pr[[1]]){			
				#colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
				generow=geneCol[geneCol[,"gene"]==gene & is.na(geneCol[,"rxn"]),]
				if(generow[4]==1){#if cnt==1
					colind=as.numeric(generow[3])
				}else{
					colind=as.numeric(geneCol[geneCol[,"gene"]==gene & geneCol[,"rxn"]==r & !is.na(geneCol[,"rxn"]),"Col"])
				}
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
#########################################------3-straight OR----------------
		#add row for the rxn
		   row_vals=rep(0,colid+length(pr))
			row_vals[vind]=lnz
			row_vals_b=rep(0,colid+length(pr))
			row_vals_b[vind_b]=lnz
		 for( p in 1:length(pr)){#ORed gene list
           if( length(pr[[p]])==1 ){
				gene=pr[[p]]
				#colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
				generow=geneCol[geneCol[,"gene"]==gene & is.na(geneCol[,"rxn"]),]
				if(generow[4]==1){#if cnt==1
					colind=as.numeric(generow[3])
				}else{
					colind=as.numeric(geneCol[geneCol[,"gene"]==gene & geneCol[,"rxn"]==r & !is.na(geneCol[,"rxn"]),"Col"])
				}

          	    #changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=colind,val=-1*kval)
				row_vals[colind]=-1*kval
				if(react_rev(model)[react_id(model)==r]){
				 #changeMatrixRow(lp=problem(prob),i=rxn_row_id_b,j=colind,val=-1*rkval)
				  row_vals_b[colind]=-1*rkval
				}
		  }else{#AND
		  #add auxiliary var then add rows
#########################################--------4-mixed OR and AND (sum of products)----------------		  
		  print("aux ands")
			addCols(lp = problem(prob),1)
			  aux_id=colid
			colid=colid+1;			
          
			#changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=aux_id,val=-1*kval)	
			row_vals[aux_id]=-1*kval
			if(react_rev(model)[react_id(model)==r]){			
			   # changeMatrixRow(lp=problem(prob),i=rxn_row_id_b,j=aux_id,val=-1*rkval)
				row_vals_b[aux_id]=-1*rkval			   
			}
			for(g in pr[[p]]){
                #colind=as.numeric(geneCol[geneCol[,"gene"]==g,"Col"])
        		generow=geneCol[geneCol[,"gene"]==g & is.na(geneCol[,"rxn"]),]
				if(generow[4]==1){#if cnt==1
					colind=as.numeric(generow[3])
				}else{
					colind=as.numeric(geneCol[geneCol[,"gene"]==g & geneCol[,"rxn"]==r & !is.na(geneCol[,"rxn"]),"Col"])
				}

		     addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(aux_id,colind)),          nzval = list(c(1,-1))
			  ,rnames = paste("Aux_",aux_id,sep="")
			 )
						 rowind=rowind+1;  
			}
          if (verboseMode > 2) print(p);
		  }		  
		  #changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=which(row_vals!=0),val=row_vals[which(row_vals!=0)])	
		}  
			addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
            cind = list(which(row_vals!=0)),        nzval = list(row_vals[which(row_vals!=0)])
			    ,rnames = r
			 )
			#rxn_row_id=rowind
			rowind=rowind+1;  
			if(react_rev(model)[react_id(model)==r]){
             		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(which(row_vals_b!=0)),        nzval = list(row_vals_b[which(row_vals_b!=0)])
			  ,rnames = paste(r,"_b",sep="")
			 )
			#rxn_row_id_b=rowind
			rowind=rowind+1;  
	       }
		
				 #rowind=rowind+1;            
	}
}	#for r
#Add multiple reactions constraints
for(g in  which(is.na(geneCol[,"rxn"]) & as.numeric(geneCol[,"cnt"])>1)){ 
	cnt=as.numeric(geneCol[g,"cnt"])
	stCol=as.numeric(geneCol[g,"Col"])
	gene=geneCol[g,"gene"]
    addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(stCol:(stCol+cnt)),        nzval = list(c(-1,rep(1,cnt)))
			  ,rnames = gene
			 )
    rowind=rowind+1; 
}

############change column bounds ##########
#n +1 ,colid
#changeColsBnds(lp = problem(prob),c(n+1:colid-1),lb=rep(0,colid-n),ub=rep(1000,colid-n))
###############II. Add MW constraint------------------
## Sum(gi. MWi)<= C[gdwpr/gDW]

if(length(MW )>0 ){
   lnz=NULL
   lcind=NULL
   if (verboseMode > 2) print(colid);
   for(v in c((n+1):(colid-1))){
   if (verboseMode > 2) print(v);
   changeColsBnds(lp = problem(prob),v,lb=0,ub=1000)
		}
   for(g in which(is.na(geneCol[,"rxn"]))) {
		colind=as.numeric(geneCol[g,"Col"])
		#changeColsBnds(lp = problem(prob),colind,lb=0,ub=1000)
		if(geneCol[g,"gene"] %in% MW[,1]){# consider only genes with computed MW.
			lnz=cbind(lnz,MW[MW[,1]==geneCol[g,"gene"],2])
			if (verboseMode > 2) print(MW[MW[,1]==geneCol[g,"gene"],2]);
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
}
			 #######################################################
	if (verboseMode > 2) {                      
				            fname=format(Sys.time(), "moment_%Y%m%d_%H%M.lp");
					       print(sprintf("Writing the problem to: %s/%s...",getwd(),fname));
		                	writeProb(lp=problem(prob), fname)
							if(solver=="cplexAPI"){
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
								write.csv(file="nnz.csv",Nnz)
							}
                }
	############################
	
	sol=optimizeProb(prob)
	
	return(list(sol,geneCol,prob))
}