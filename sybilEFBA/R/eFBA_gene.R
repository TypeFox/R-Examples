################################################
# Function: eFBA
#
# Performs an expression based flux balance analysis
# 
# Take into account the expression status of genes
#  first step find FB as the max biomass using simpleFBA
# use FB=cTx as a new constraint in a MILP with a new ojective function
#      Mininmize: Sum(expr(g)<>flux(g)  // log2 expr level can be used as scaling
#   where expr(g)=1 if g expressed value>T1, else 0
#          flux(g)=0 if (all reactions catalyzed by g) have no flux(Threshold T2 and sum of all), 0 else.
############
#Multifunctional enzymes is making a problem: no solution exists  4/2/2015
# for multifunctional enzymes add a constraint g=g1 or g2 or g3 ....
eFBA_gene <- function (model, expressionData,
                  Tf =0.01,#SYBIL_SETTINGS("TOLERANCE"),  # threshold on flux
				  pct_objective=100,
 		           lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
			   solver = SYBIL_SETTINGS("SOLVER"),
			    method = SYBIL_SETTINGS("METHOD"),
					 #solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM")
			  solverParm=data.frame(CPX_PARAM_EPRHS=1e-6),
			  testgpr=NULL,
                  verboseMode = 2, ...){
#PARAMETERS:
#===========
# model                 Sybil model structure (class modelorg)
# expressionData        a dataframe (gene expression data): gene ID, and expression level
# Tf                    Threshold value for flux(g)

# RETURN
# optsol class
# Second obj function optimal value

## getRxnEqn: returns the equation for a given rxn, used for output
getRxnEqn=function(rxn=1){
	m=S(model)[,rxn]
	# input ==> output
	mcf=ifelse(abs(m)>1,paste("(",abs(m),")",sep=""),"")
	eqn=paste(gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m<0)],sep=""))),
		ifelse(react_rev(model)[rxn],"<==>","-->"),
	    gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m>0)],sep=""))) ,sep=" ")
	
return(eqn)
}
##--------------------------------------------------------------------------##
 # check prerequisites 
    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
    }
    
#require(sybilCPLEX)
##-------------------------------------Prepare Problem object --------------------##
# get OptObj instance
#lp <- sybil::sysBiolAlg(model, algorithm = "fba", ...)

			 
##-------------------------------------------------------------------------------##


	# # Run FBA
          orig_sol = sybil::optimizeProb(model,solver=solver);# solver and method already in lp
		FB =  lp_obj(orig_sol)*pct_objective/100# $obj;  ##objvalue sol.f
	    if ( length(checkSolStat(lp_stat(orig_sol),solver))!=0 ){## checkSolStat
			print('No feasible solution - for initial FBA problem\n');
			break;
	    }
	  if (verboseMode > 3) {
	  	print(orig_sol);
	  	}
	   #if (length(max_growth_rate)>0 ) return (orig_sol);
##-------------------------------------------------------------------------------##
# get flux status Flux(g)
# 
#	exprStatus=ifelse(expressionData$level<Te,0,1);
#   genflux as integer variables yi


#formulate new problem eFBA problem
    #1- add new n variables integer vars ,
    
    #2- add new constraint for FB
    #3- set new obj function if expr(gi)=0 then add yi to obj else add -yi. 
    #4- Add new n constraints(identifying yi's): 2 for irreversible and 4 for reversible rxns
	#5- Add constraints for multifunctional enzymes

        nc     <- react_num(model)
        nr     <- met_num(model)
        ng    <- length(allGenes(model))
        geneNames=allGenes(model)
        if(length( testgpr)==0){
                  gpr_rxn_ind=(gpr(model)!="" ) # only rxns having GPR rule
         }else { gpr_rxn_ind=testgpr;}
		
		print(sprintf("%s: Calculate GPR state ....",format(Sys.time(), "%d-%m-%Y %X")))
            x <- logical(length(allGenes(model)));
		   	gi=expressionData[,1] %in% allGenes(model);##expressionData$Locus
		   	gj=match(expressionData[,1],allGenes(model) )  ;
		   	x <- rep(NA,length(x));  # missing genes will be ignored and not included in objective function or as 0 state.
			x[gj]=ifelse(expressionData[,2][gi]==0,FALSE,ifelse(expressionData[,2][gi]==1,TRUE,NA)); #NA undefined
			rxnStatus =rep(0,nc);# 0:will not be included in objective, -1: ON(Max), 1 OFF (Min)
		   	for (i in 1:nc){
		    		# test if gene is important or not in identifying rule state
		       		if(gprRules(model)[i]!="") {
						tmp=eval(parse(text = gprRules(model)[i]));
						rxnStatus[i]=ifelse(is.na(tmp),0,ifelse(tmp,-1,1));#coefficient of rxn in objective
						#=0 when no rule exists 
					}
			}
			print("state after GPR state evaluation:");
			print(table(rxnStatus))
		rxnStatus[!gpr_rxn_ind]=0;# exclude rxns not in testgpr
		print("state after testgpr:");
			print(table(rxnStatus))
		# test flux variabilty of rxns before adding constraints
		# if a reaction fixed or having a range different from gene state it will be ignored
			print(sprintf("%s: Calculate FVA ....",format(Sys.time(), "%d-%m-%Y %X")))
            rxnfva=(rxnStatus!=0)
			fv <- fluxVar(model,react=which(rxnfva),solver=solver)
			fv_min=lp_obj(fv)[1:sum(rxnfva)];fv_max=lp_obj(fv)[(sum(rxnfva)+1):(2*sum(rxnfva))]
		
		#print(length(fv
		rxnStatus[react_pos(react(fv))][abs(fv_min-fv_max)<Tf]=0;# exclude rxns with fixed values
		#state fixed: exclude rxns that must be ON (state can't be changed
		rxnStatus[react_pos(react(fv))][(fv_min>Tf) | (fv_max < -Tf)]=0;# cutoff
		print(sprintf("Number of rxns in FVA: %d, Number of rxns found fixed: %d,
		\n rxns that are always ON: %d ",sum(rxnfva),sum(abs(fv_min-fv_max)<Tf),sum(((fv_min>Tf) | (fv_max < -Tf)) & abs(fv_min-fv_max)>=Tf) ))
    	print("state after FVA:");
		print(table(rxnStatus))
		gpr_rxn_ind[rxnStatus==0]=FALSE;
		rxnstname=ifelse(rxnStatus==0,"No rule/Unk",ifelse(rxnStatus==-1,"ON","OFF"));
		print(table(rxnstname))
		
	     
        #  the problem: minimize:
        #  
        #            |      |     
        #         S  |  0   |  = b
        #            |      |     
        #       ------------------
        #            |      |
        #        c^T |  0   | = FB
        #            |      |
        #       ------------------
        #          identifying constraints of fi
		#       ------------------
		#       Linear constraints of GPR
        #       ------------------
        #  lb   wt_lb|  0   |
        #  ub   wt_ub|  1   |
        #            |      |
        #  obj    0  |  2*ei-1    where ei =0 gene i NOT expressed/ 1 otherwise

        nGpr=sum(gpr_rxn_ind);
		rules=gpr(model)[gpr_rxn_ind]
        is_irrev=(gpr_rxn_ind & lowbnd(model)>=0);
	    nIrrev=sum(is_irrev);
		
		is_rev=(gpr_rxn_ind & lowbnd(model)<0);
        nRev=sum(is_rev);
		INF=max(uppbnd(model))+1;
		#--------------------identify variable for genes:05/02/2015--------------------------
		colid=nc+nGpr+nRev+1;
		geneCol=NULL;#column in problemcoresponding to given gene

		#Add variables for all genes
		for( g in allGenes(model)){
				cnt=length(grep(g,gpr(model)[gpr_rxn_ind]))#include only rules of reactions
				geneCol=rbind(geneCol,cbind(gene=g,rxn=NA,Col=colid,cnt,varname=paste("g",g,sep="_"),inputStat=x[grep(g,allGenes(model))]))
				colid=colid+1;
				if(cnt>1){
					for(r in grep(g,gpr(model)[gpr_rxn_ind]) ){
						geneCol=rbind(geneCol,cbind(gene=g,rxn=(react_id(model)[gpr_rxn_ind])[r],Col=colid,cnt,varname=paste("g",g,r,sep="_"),inputStat=x[grep(g,allGenes(model))]))
						colid=colid+1;			
					}
				}
		}
		nvg=colid-(nc+nGpr+nRev+1);
		# ---------------------------------------------
        # constraint matrix
        # ---------------------------------------------

       # the initial matrix dimensions  variables:   nc:rates, nc:flux(r)(0,1) , ng
	   # integer variables:fi(nGpr),yi(nRev),gene state ng
	   # rows nIrrev:2, nRev:4  ,Gpr:nGpr         ---02/02/2015--------
       LHS <- Matrix::Matrix(0, nrow = nr+nIrrev*2+nRev*4+nGpr+1, ncol = (nc+nGpr+nRev+nvg), sparse = TRUE)

       # rows for the initial S
       LHS[1:nr,1:nc] <- S(model)

       # fix the value of the objective function
       LHS[(nr+1),1:nc] <- obj_coef(model)

       # rows for the identifying constraint of flux variables
       #diag is wrong : it cant be seq, only the first 810 will be used regardless of gpr existence
	# Irreversible reactions 2 constraints
	if(nIrrev>0){
    #  Constraint 1: vi-M*fi<=Tf
       ii=matrix(c((nr+2)   :(nr+nIrrev+1)  ,((1:nc)[is_irrev] )),ncol=2)
       LHS[ii]<- 1
		
		if(nIrrev>1){
			diag(LHS[(nr+2)   :(nr+nIrrev+1)  ,((nc+1):(nc+nIrrev))   ]) <- -INF   # M*fi
		}else{# diag function fails when it is one row
			LHS[(nr+2)   :(nr+nIrrev+1)  ,((nc+1):(nc+nIrrev))  ] <- -INF
		}
		
		#  Constraint 2: Tf*fi-vi<=0
       ii=matrix(c((nr+nIrrev+2)   :(nr+2*nIrrev+1)  ,((1:nc)[is_irrev] )),ncol=2)
       LHS[ii]<- -1
		
		if(nIrrev>1){
			diag(LHS[(nr+nIrrev+2)   :(nr+2*nIrrev+1)  ,((nc+1):(nc+nIrrev))   ]) <- Tf   # Tf*fi
		}else{# diag function fails when it is one row
			LHS[(nr+nIrrev+2)   :(nr+2*nIrrev+1)  ,((nc+1):(nc+nIrrev))   ] <- Tf
		}
	}#Irrev
	# Irreversible reactions 4 constraints, one additional integer variable
	if(nRev>0){
    #  Constraint 1: vi-M*fi<=Tf
       ii=matrix(c((nr+2*nIrrev+2)   :(nr+2*nIrrev+nRev+1)  ,((1:nc)[is_rev] )),ncol=2)
       LHS[ii]<- 1
	    
		if(nRev>1){
			diag(LHS[(nr+2*nIrrev+2)   :(nr+2*nIrrev+nRev+1)  ,((nc+nIrrev+1):(nc+nIrrev+nRev))   ]) <- -INF   # M*fi
		}else{# diag function fails when it is one row
			LHS[(nr+2*nIrrev+2)   :(nr+2*nIrrev+nRev+1)  ,((nc+nIrrev+1):(nc+nIrrev+nRev))   ] <- -INF
		}
		#---------------------------------------------------------#
		#  Constraint 2: Tf*fi-vi-My<=0
        ii=matrix(c((nr+2*nIrrev+nRev+2)   :(nr+2*nIrrev+2*nRev+1)  ,((1:nc)[is_rev] )),ncol=2)
        LHS[ii]<- -1
		if(nRev>1){
			diag(LHS[(nr+2*nIrrev+nRev+2)   :(nr+2*nIrrev+2*nRev+1)  ,((nc+nIrrev+1):(nc+nIrrev+nRev))   ]) <- Tf   # Tf*fi
			diag(LHS[(nr+2*nIrrev+nRev+2)   :(nr+2*nIrrev+2*nRev+1)  ,((nc+nIrrev+nRev+1):(nc+nIrrev+2*nRev))   ]) <- -INF   # -M*yi
		}else{# diag function fails when it is one row
			LHS[(nr+2*nIrrev+nRev+2)   :(nr+2*nIrrev+2*nRev+1)  ,((nc+nIrrev+1):(nc+nIrrev+nRev))   ] <- Tf   # Tf*fi
			LHS[(nr+2*nIrrev+nRev+2)   :(nr+2*nIrrev+2*nRev+1)  ,((nc+nIrrev+nRev+1):(nc+nIrrev+2*nRev))   ] <- -INF   # -M*yi
		}
       # ---------------------------------------------
    #  Constraint 3: -vi-M*fi<=Tf
       ii=matrix(c((nr+2*nIrrev+2*nRev+2)   :(nr+2*nIrrev+3*nRev+1)  ,((1:nc)[is_rev] )),ncol=2)
       LHS[ii]<- -1
	    
		if(nRev>1){
			diag(LHS[(nr+2*nIrrev+2*nRev+2)   :(nr+2*nIrrev+3*nRev+1)  ,((nc+nIrrev+1):(nc+nIrrev+nRev))   ]) <- -INF   # M*fi
		}else{# diag function fails when it is one row
			LHS[(nr+2*nIrrev+2*nRev+2)   :(nr+2*nIrrev+3*nRev+1)  ,((nc+nIrrev+1):(nc+nIrrev+nRev))   ] <- -INF
		}
		#---------------------------------------------------------#
		#  Constraint 4: Tf*fi+vi-M(1-y)<=0
        ii=matrix(c((nr+2*nIrrev+3*nRev+2)   :(nr+2*nIrrev+4*nRev+1)  ,((1:nc)[is_rev] )),ncol=2)
        LHS[ii]<- 1
		if(nRev>1){
			diag(LHS[(nr+2*nIrrev+3*nRev+2)   :(nr+2*nIrrev+4*nRev+1)  ,((nc+nIrrev+1):(nc+nIrrev+nRev))   ]) <- Tf   # Tf*fi
			diag(LHS[(nr+2*nIrrev+3*nRev+2)   :(nr+2*nIrrev+4*nRev+1)  ,((nc+nIrrev+nRev+1):(nc+nIrrev+2*nRev))   ]) <- INF   # M*yi
		}else{# diag function fails when it is one row
			LHS[(nr+2*nIrrev+3*nRev+2)   :(nr+2*nIrrev+4*nRev+1)  ,((nc+nIrrev+1):(nc+nIrrev+nRev))   ] <- Tf   # Tf*fi
			LHS[(nr+2*nIrrev+3*nRev+2)   :(nr+2*nIrrev+4*nRev+1)  ,((nc+nIrrev+nRev+1):(nc+nIrrev+2*nRev))   ] <- INF   # M*yi
		}
	}#nRev
	#---------------------------------------------------------#
	         # lower and upper bounds for COLUMNS(variables) and ROWS(constraints)
    # ---------------------------------------------
     
     	lower  <- c(lowbnd(model), rep(0, nIrrev+2*nRev), rep(0, nvg))
     	upper  <- c(uppbnd(model), rep(1, nIrrev+2*nRev), rep(1, nvg))
     	#RHS 1:nr+1 : as FBA model , FB,      ?? why -2*ub!!
		#rhs(model) : removed from model, 30/3/2013
     	#rlower <- c(rhs(model), FB, -2*uppbnd(model)[gpr_rxn_ind],-2*uppbnd(model)[gpr_rxn_ind],rep(0, nGpr))  # ,rep(0,ngpr)
    	#rupper <- c(rhs(model), FB, rep(0, 3*nGpr))  #,rep(0,ngpr)
        rlower <- c(rep(0,nr), FB, rep(-INF,nIrrev), rep(-INF,nIrrev), rep(-2*INF,nRev), rep(-2*INF,nRev), rep(-2*INF,nRev), rep(-INF,nRev),rep(0,nGpr))  # ,rep(0,ngpr)
    	rupper <- c(rep(0,nr), FB, rep(Tf,nIrrev),   rep(0,nIrrev),    rep(Tf,nRev),     rep(0,nRev),      rep(Tf,nRev),     rep(INF,nRev),rep(0,nGpr))  #,rep(0,ngpr)
        
		gprtype=rep("E",nGpr);
    #-----------------------------------------------------------------------------------------#
	# constraints of gpr: NOT is not considered
	message("start gpr....");
	 lastRow=dim(LHS)[1];

	 row_i=nr+2*nIrrev+4*nRev+1;
     rxn_map=c(which(is_irrev),which(is_rev))# mapping from rxn to fi:at first irreversible rxns then reversible
	for(i in 1:nGpr){  # i : is the rxn number
		 #rl=rules[i];
		 rl=gpr(model)[rxn_map[i]];
# 				# search for brackets : add aux variable for each bracket  in the same way as above
				# Consider only SUM of PRODUCT (AND to [OR]), without NOT, one level
		row_i=row_i+1;
		if (verboseMode > 2) {
			print(c(i,rxn_map[i],rl))
		}
		# may make a problem if gene name contains 'or' like YOR123
		#replace )or with ) or & or( with or (
		rl=gsub("\\)"," ) ",rl)# 
		rl=gsub("\\("," ( ",rl)# 

		pr=lapply(strsplit(unlist(strsplit(rl," or "))," and "),function(x) gsub("[() ]","",x))
		if( length(pr)==1) {# no OR (only one term) 
			#LHS[row_i,nc+nGpr+nRev+which(geneNames %in% pr[[1]])]=1;
			    for(gene in pr[[1]]){
					#get the corresponding col of the gene
					generow=geneCol[geneCol[,"gene"]==gene & is.na(geneCol[,"rxn"]),]
					if(generow[4]==1){#if cnt==1
						colind=as.numeric(generow[3])
					}else{
						colind=as.numeric(geneCol[geneCol[,"gene"]==gene & geneCol[,"rxn"]==react_id(model)[rxn_map[i]] & !is.na(geneCol[,"rxn"]),"Col"])
					}
					LHS[row_i,colind]=1;
				}
			LHS[row_i,nc+i]=-length(pr[[1]]);  # this is for the rxn
			#rlower[row_i]=-0.1; 
			rupper[row_i]=length(pr[[1]] )-1;#+0.1
			if( length(pr[[1]]) >1)  gprtype[i]="R"; # Range
		 } else {# first the sum row
		     LHS[row_i,nc+i]=length(pr);  # this is for the rxn
		     #rlower[row_i]=-0.1; # put to zero 
		     rupper[row_i]=length(pr)-1;	
		     gprtype[i]="R";
		     for( p in 1:length(pr)){
                     	 if( length(pr[[p]])==1 ){ 
			                #LHS[row_i,nc+nGpr+nRev+which( geneNames %in% pr[[p]] ) ]=-1;
							gene=pr[[p]]
							#get the corresponding col of the gene							
							generow=geneCol[geneCol[,"gene"]==gene & is.na(geneCol[,"rxn"]),]
							#print(generow)
							#print(react_id(model)[rxn_map[i]])
							if(generow[4]==1){#if cnt==1
								colind=as.numeric(generow[3])
							}else{
								colind=as.numeric(geneCol[geneCol[,"gene"]==gene & geneCol[,"rxn"]==react_id(model)[rxn_map[i]] & !is.na(geneCol[,"rxn"]),"Col"])
							}	
							#print(colind)
							LHS[row_i,colind]=-1;					
			             }else{       # add auxiliary variable : add row for it then column update new row and SUM row
							  crow=Matrix::Matrix(0, nrow = 1, ncol =dim(LHS)[2])
							  LHS=rBind(LHS,crow)
							  rlower=c(rlower,0); 				  rupper=c(rupper,0);
							  lastRow=lastRow+1;
							  ccol=Matrix::Matrix(0, nrow = lastRow, ncol =1)
							  LHS=cBind(LHS,ccol) # new column
							  lower=c(lower,0);	  upper=c(upper,1);
							  LHS[lastRow,dim(LHS)[2]]=-length(pr[[p]]);  #Coef of Aux
							  #LHS[lastRow,nc+nGpr+nRev+which(geneNames %in% pr[[p]])]=1;
								for(gene in pr[[p]]){
									#get the corresponding col of the gene
									generow=geneCol[geneCol[,"gene"]==gene & is.na(geneCol[,"rxn"]),]
									if(generow[4]==1){#if cnt==1
										colind=as.numeric(generow[3])
									}else{
										colind=as.numeric(geneCol[geneCol[,"gene"]==gene & geneCol[,"rxn"]==react_id(model)[rxn_map[i]] & !is.na(geneCol[,"rxn"]),"Col"])
									}
									LHS[lastRow,colind]=1;
								}
							  rupper[lastRow]=length(pr[[p]])-1;#RHS of Aux constraint
							  LHS[row_i,dim(LHS)[2]]=-1;  # add the Aux variable to the Sum row
							}
			}#for p
		}#if len
       }# for i
######################
# Constraints of multifunctional enzymes
######################
  aux_cnt=dim(LHS)[1]-(nr+2*nIrrev+4*nRev+1);
  mrg=geneCol[is.na(geneCol[,"rxn"]) & geneCol[,"cnt"]>1 & !is.na(geneCol[,"inputStat"]),];
  if(length(mrg[,1])>0){
	  for(gr in (1:length(mrg[,1]))){
	  #Add constraint x=x1 or x2 or ... xn, 0<= n*x-x1-x2...<=n-1
		  crow=Matrix::Matrix(0, nrow = 1, ncol =dim(LHS)[2])
		  LHS=rBind(LHS,crow)
		  rlower=c(rlower,0); 				  rupper=c(rupper,as.numeric(mrg[gr,"cnt"])-1);
		  lastRow=dim(LHS)[1]#lastRow+1;
		  LHS[lastRow,as.numeric(mrg[gr,"Col"])]=as.numeric(mrg[gr,"cnt"])
		  sc=as.numeric(geneCol[!is.na(geneCol[,"rxn"]) & geneCol[,"gene"]==mrg[gr,"gene"],"Col"]);
		  LHS[lastRow,sc]=-1
	  }
   }
	mr_cnt=dim(LHS)[1]-(nr+2*nIrrev+4*nRev+1+aux_cnt)
###############################
       message("gpr ... ... OK")
       # ---------------------------------------------
       # objective function
       # ---------------------------------------------
       # Notes: 1-use gene status
       # Min(ifelse(expr(rxn(i),-yi,yi)     math. |x-y|=x-y  when x>=y and y-x when y>=x 
       # then the difference will be Sum(rxn(i))+obj 
       #
       
        
    num_constr=dim(LHS)[1];
	num_var=dim(LHS)[2];
	
	#print(paste("no cons:",num_constr,"no var:",num_var));
	#x:expression state true,false,na
    #cobj <- c(rep(0, nc+nGpr+nRev),ifelse(is.na(x),0,ifelse(x,-1,1)) ,rep(0, num_var-(nc+nGpr+nRev+ng))) # NA->0
	cobj <- rep(0,num_var);
	cobj[as.numeric(geneCol[is.na(geneCol[,"rxn"]),"Col"])] <- ifelse(is.na(x),0,ifelse(x,-1,1));
	
	cNames=paste(c(rep("x",nc),rep("Irv",nIrrev),rep("Rev",nRev),rep("y",nRev),rep("",nvg)),
		 		               c(1:nc,which(is_irrev),which(is_rev),which(is_rev),geneCol[,"varname"]),sep="" ) ;
	if((num_var-(nc+nGpr+nRev+nvg))>0){cNames=c(cNames,paste(rep("Aux",(num_var-(nc+nGpr+nRev+nvg))),1:(num_var-(nc+nGpr+nRev+nvg)),sep=""))}
	#browser();#Q to exit
	# ---------------------------------------------
        # build problem object
    # ---------------------------------------------
	switch(solver,
            # ----------------------- #
            "glpkAPI" = {
                out <- vector(mode = "list", length = 4)
                prob <- glpkAPI::initProbGLPK();# new problem
               
                out[[1]] <-  glpkAPI::addRowsGLPK(prob, nrows=num_constr)
		outj <- glpkAPI::addColsGLPK(prob, ncols=num_var)
		#setColNameGLPK(prob,c(1:(num_var)),paste(c(rep("x",nc),rep("rxn",nGpr),rep("g",ng),rep("Aux",(num_var-(nc+nGpr+ng)))),
		#               c(1:nc,which(gpr_rxn_ind),1:ng,1:(num_var-(nc+nGpr+ng))),sep="" ) );
		mapply(glpkAPI::setColNameGLPK, j = c(1:(num_var)), cname = cNames, MoreArgs = list(lp = prob));

		glpkAPI::setObjDirGLPK(prob, glpkAPI::GLP_MIN);
                ## note: when FX or LO value taken from lb and when UP take from UP and when R
               rtype <- c(rep(glpkAPI::GLP_FX, nr),glpkAPI::GLP_LO, rep(glpkAPI::GLP_UP, 2*nIrrev+4*nRev), ifelse(gprtype=="E",glpkAPI::GLP_FX,glpkAPI::GLP_DB),
			   rep(glpkAPI::GLP_DB,aux_cnt),rep(glpkAPI::GLP_DB,mr_cnt))
	       print(table(rtype));
	        # set the right hand side Sv = b

	       out[[4]] <- glpkAPI::setRowsBndsGLPK(prob, c(1:num_constr), lb=rlower, ub=rupper,type=rtype )
        
	        # add upper and lower bounds: ai <= vi <= bi
                cc <- glpkAPI::setColsBndsObjCoefsGLPK(prob, c(1:num_var),   lower,   upper,  cobj   )
                #print(cc);			
        	#nzLHS=nzijr(LHS);
            #cc <- loadMatrixGLPK(prob,   nzLHS$ne,   nzLHS$ia,  nzLHS$ja,  nzLHS$ar     )
			TMPmat <- as(LHS, "TsparseMatrix")
            cc <- glpkAPI::loadMatrixGLPK(prob,length(TMPmat@x),ia  = TMPmat@i + 1,ja  = TMPmat@j + 1,ra  = TMPmat@x)
                
            ctype <- c(rep(glpkAPI::GLP_CV, nc), rep(glpkAPI::GLP_BV,num_var-nc)); # number of integer variable
		    glpkAPI::setColsKindGLPK(prob,c(1:(num_var)),ctype);

                if (verboseMode > 2) {                      
				        fname=format(Sys.time(), "glpk_eFBA_%Y%m%d_%H%M.lp");
					print(sprintf("writing problem to: %s/%s...",getwd(),fname));
		                	glpkAPI::writeLPGLPK(prob,fname);
		                	print("Solving...");
                }
                
                ## Solve
                glpkAPI::setMIPParmGLPK(glpkAPI::PRESOLVE,glpkAPI::GLP_ON);
		lp_ok=glpkAPI::solveMIPGLPK(prob);
		 if (verboseMode > 2) {
		 	print(glpkAPI::return_codeGLPK(lp_ok));
		 	}
		lp_stat=glpkAPI::mipStatusGLPK(prob);
		 if (verboseMode > 2) {
		 	print(glpkAPI::status_codeGLPK(lp_stat));
		 	}
		lp_obj=glpkAPI::mipObjValGLPK(prob);
		colVal=glpkAPI::mipColsValGLPK(prob);
		newFlux=colVal
		colst=cbind(cNames,val=colVal,lower,upper,ctype);
	       ## --- 
	        #newFlux=floor(newFlux/Tf)*Tf#sybil:::.floorValues(newFlux,tol=Tf);
	        sol_geneStat=ifelse(newFlux[(nc+nGpr+nRev+1) : (nc+nGpr+nRev+nvg)]==1,"ON","OFF");
	        newFlux=newFlux[1:nc];
	        newStat=ifelse(abs(newFlux)>Tf,1,0);

            },
            # ----------------------- ----------------------- ------------------------- ----------------------- ----------------------#
           "cplexAPI" = {
                 out <- vector(mode = "list", length = 3)
				 prob <- cplexAPI::openProbCPLEX()
				 out <- cplexAPI::setIntParmCPLEX(prob$env, cplexAPI::CPX_PARAM_SCRIND, cplexAPI::CPX_OFF)		                
                
                # when R: rngval+rhs   to rhs   (rngval<0), E,L,R,U?
                rtype <- c(rep("E",nr),"G",  rep("L", 2*nIrrev+4*nRev),gprtype,rep("R",aux_cnt),rep("R",mr_cnt))
                prob$lp<- cplexAPI::initProbCPLEX(prob$env)
				cplexAPI::chgProbNameCPLEX(prob$env, prob$lp, "eFBA gene cplex");
                 out[[1]] <- cplexAPI::newRowsCPLEX(prob$env, prob$lp,
                                         nrows=num_constr, rhs=rupper, sense=rtype,rngval=rlower-rupper)

                 out[[2]] <- cplexAPI::newColsCPLEX(prob$env, prob$lp,
                                         num_var, obj=cobj, lb=lower, ub=upper,cnames=cNames)
          
				TMPmat <- as(LHS, "TsparseMatrix")
				out[[3]] <- cplexAPI::chgCoefListCPLEX(prob$env, prob$lp,#lp@oobj@env, lp@oobj@lp,
                                   nnz = length(TMPmat@x),
                                   ia  = TMPmat@i,
                                   ja  = TMPmat@j,
                                   ra  = TMPmat@x);				
	      
                if(num_var>nc){ ctype<-c(rep('C', nc), rep('B',num_var - nc));# integer variables
				}else{ctype<-rep('C', nc);}
				status = cplexAPI::copyColTypeCPLEX (prob$env, prob$lp, ctype);
				check <- cplexAPI::setObjDirCPLEX(prob$env, prob$lp, cplexAPI::CPX_MIN);
		
	         if (verboseMode > 2) {                      
                   fname=format(Sys.time(), "Cplex_eFBA_gene_%Y%m%d_%H%M.lp");
	 			    print(sprintf("writing problem to: %s/%s ...",getwd(),fname));
		 			cplexAPI::writeProbCPLEX(prob$env, prob$lp,fname);
		       }
              #   --------------------------------------------------------------
              print("Solving...");
	       lp_ok     <- cplexAPI::mipoptCPLEX(prob$env, prob$lp);
	       print(lp_ok);
	       sol=cplexAPI::solutionCPLEX(prob$env, prob$lp);
                if (verboseMode > 2) {
                	print(sol)
                	}
              if(is(sol)[1]=="cpxerr"){
					print(sol);
					stop("Execution Terminated! error in MILP")
              }
              lp_obj=sol$objval;
               lp_stat   <- cplexAPI::getStatCPLEX(prob$env, prob$lp)
              
               colst=cbind(cNames,val=sol$x,lower,upper,ctype);
                 
               #have flux and gene ON /flux and gene OFF?
               newFlux=floor(sol$x/Tf)*Tf#sybil:::.floorValues(sol$x,tol=Tf);
               newFlux=newFlux[1:nc];
               sol_geneStat=ifelse(sol$x[(nc+nGpr+nRev+1):(nc+nGpr+nRev+nvg)]==1,"ON","OFF");
               newStat=ifelse(abs(newFlux)>Tf,1,0);
               # when using binary variables as indicators: .floorValues(.ceilValues(newFlux[(nc+1):(2*nc)],tol=Tf/10),tol=Tf);
           },
           # ----------------------- #
           {
               wrong_type_msg(solver)
            }
        )
        
 	       print(sprintf("Number of Rxns with gene rule ON=%d ",length(rxnStatus[rxnStatus==-1]) ));

	       print(sprintf("Rxns with gene OFF=%d",length(rxnStatus[rxnStatus==1]) ));  
	       print(sprintf("Total Nr of selected rules=%d ",length(rxnStatus[rxnStatus==0]) ));
	       print(sprintf("Total number of rxns: %d",length(rxnStatus) ));
	       print(sprintf("The difference from objective: %.0f ", (lp_obj+length(rxnStatus[rxnStatus==-1])) ));#may contain a problem

 		origFlux=getFluxDist(orig_sol);
		
 		#excReact = findExchReact(model)[1];# 1 is position
		#excReactPos=react_pos(excReact$exchange);
        excReact = findExchReact(model);
        excReactPos=(react_id(model) %in% react_id(excReact));#excReact$exchange

        is_excReact=((c(1:length(react_id(model)))) %in% excReactPos);
               
                origStat=ifelse(abs(origFlux)>Tf,1,0);
				print(sprintf("Difference from FBA fluxes: %d",sum(abs(origStat-newStat))));
				
               rxnstname=ifelse(rxnStatus==0,"No rule/Unk",ifelse(rxnStatus==-1,"ON","OFF"));
               rxn_st <- cbind(gpr=gpr(model),react=react_id(model),reactName=react_name(model),is_excReact=is_excReact,
                                 expr=rxnstname,lb=lowbnd(model),
                                       #  eqns=sapply(c(1:length(react_id(model))),function(x)getRxnEqn(x)) ,
                                          origFlux,origStat, newFlux,newStat,difference=abs(origStat-newStat));
 
               #gene_st <- cbind(geneLocus=allGenes(model),inputStat=x,SolStat=sol_geneStat)
			   gene_st <- cbind(varname=geneCol[,"varname"],inputStat=geneCol[,"inputStat"],SolStat=sol_geneStat)
 
  remove(prob);
#----------------------------------------  ***  --------------------------------#
	optsol <- list(        ok = lp_ok,
			       obj = FB,
			       stat = lp_stat,
			       fluxes = newFlux,
			       origfluxes=origFlux,
			       rxn=rxn_st,
			       new_objFn=lp_obj,
			       gene_stat=sol_geneStat,
			       colst=colst
			  )
	return(optsol);
}
#               write.csv(rxn_st,"rst.csv");
#              write.csv(gene_st,"gene.csv");
               
#               output_rep=table(list(expr=rxnstname,orig=origStat,newsol=newStat));
#               write.csv(output_rep,"eFBA_report.csv");
