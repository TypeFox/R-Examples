findFluxGeneExpr=function(model,fluxes,threshold=1e-6,
                 lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
			   solver = SYBIL_SETTINGS("SOLVER"),
			    method = SYBIL_SETTINGS("METHOD"),
                  solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM"),
				  verboseMode = 2){
# given the flux rate use gpr to formulate a LP to give solution of  min required gene expression
# gpr_rxn
#-------------------------------------------------     ========   ---------------------------------------------------#

# use only genes as variables and substitute Flux(rxn) with the state of the flux
	ng    <- length(allGenes(model))

	gpr_rxn_ind=(gpr(model)!="")
	gpr_rxn=sum(gpr_rxn_ind);
	rules=gpr(model)[gpr_rxn_ind]
        geneNames=allGenes(model)
 
	nc=ng;
	nr=gpr_rxn;

       #LHS <- as.matrix.csr(0, nrow = nr, ncol = nc)
		LHS <- Matrix::Matrix(0, nrow = nr, ncol = nc, sparse = TRUE)

       # ---------------------------------------------
            # lower and upper bounds for COLUMNS(variables) and ROWS(constraints)
       # ---------------------------------------------
     
     	lower  <-  rep(0, nc)
     	upper  <-  rep(1, nc)
     	#RHS 1:nr+1 : as FBA model , FB,      ?? why -2*ub!!
     	rlower <- rep(0, nr) 
	    rupper <- rep(0, nr)
        gprtype=rep("E",nr);
        
        gprFluxStat=ifelse(abs(fluxes[gpr_rxn_ind])>threshold,1,0);
        lastRow=nr;
		

	for(i in 1:gpr_rxn){  # i : is the rxn number
		rl=rules[i];
 				# search for brackets : add aux variable for each bracket  in the same way as above
				# Consider only SUM of PRODUCT (AND to [OR])
		#pr=lapply(strsplit(unlist(strsplit(rules[i],"or")),"and"),function(x) gsub("[() ]","",x))
		rl=gsub("\\)"," ) ",rl)# 
		rl=gsub("\\("," ( ",rl)# 
		# to be sure that 'and' will be a whole word and not part of any gene name
		pr=lapply(strsplit(unlist(strsplit(rl," or "))," and "),function(x) gsub("[() ]","",x))
		if( length(pr)==1) {# no OR (only one term) 
			LHS[i,which(geneNames %in% pr[[1]])]=1;
		        rupper[i]=length(pr[[1]]);  
		        if(gprFluxStat[i]==1){# gene is required
		                	   rlower[i]=length( pr[[1]] );		                	    
		        }else{
		                       rlower[i]=0;    # soft constraint: becaue when it is tight no solution found
		                	   gprtype[i]="R";
		        }     
		 } else {# first the sum row
			      rupper[i]=length(pr);
			      rlower[i]=gprFluxStat[i];  # soft constraint in case of missing flux
			      gprtype[i]="R";
		     for( p in 1:length(pr)){
                     	 if( length(pr[[p]])==1 ){ 
			              LHS[i,which( geneNames %in% pr[[p]] ) ]=1;
			}else{       # add auxiliary variable : add row for it then column update new row and SUM row
				  #LHS=rbind(LHS,as.matrix.csr(0, nrow = 1, ncol =dim(LHS)[2]))
				  crow=Matrix::Matrix(0, nrow = 1, ncol =dim(LHS)[2])
				  LHS=rBind(LHS,crow)
				  rlower=c(rlower,0); 				  rupper=c(rupper,0);
				  lastRow=lastRow+1;
				  #LHS=cbind(LHS,as.matrix.csr(0, nrow = lastRow, ncol =1)) # new column
				  ccol=Matrix::Matrix(0, nrow = lastRow, ncol =1)
				  LHS=cBind(LHS,ccol)
				  lower=c(lower,0);	  upper=c(upper,1);

				  LHS[lastRow,dim(LHS)[2]]=-length(pr[[p]]);  
				  LHS[lastRow,which(geneNames %in% pr[[p]])]=1;
				  rupper[lastRow]=length(pr[[p]])-1;

				  LHS[i,dim(LHS)[2]]=1;  # add the Aux variable to Sum row
			}
		}#for p
		}#if len
       }# for i
       print(rupper)
       print(rlower)
       print(gprFluxStat)
       num_constr=dim(LHS)[1];
	num_var=dim(LHS)[2];
	aux_cnt=num_constr-ng;
	print(paste("no cons:",num_constr,"no var:",num_var));

#  objective Min number of genes ON
       	cobj <- rep(0, num_var)
       	cobj[1:ng]=1;
	#--------
	cNames=paste(c(rep("g",ng),rep("Aux",(num_var-ng))),
		 		               c(1:ng,1:(num_var-ng)),sep="" ) ;
# formulate and solve
	 switch(solver,
            # ----------------------- #
            "glpkAPI" = {
                out <- vector(mode = "list", length = 4)
                prob <- glpkAPI::initProbGLPK();# new problem
               
                out[[1]] <-  glpkAPI::addRowsGLPK(prob, nrows=num_constr)
				outj <- glpkAPI::addColsGLPK(prob, ncols=num_var)
				#setColNameGLPK(prob,c(1:(num_var)),paste(c(rep("g",ng),rep("Aux",(num_var-ng))),
				#			   c(1:ng,1:(num_var-ng)),sep="" ) );
				mapply(glpkAPI::setColNameGLPK, j = c(1:(num_var)), cname = cNames, MoreArgs = list(lp = prob));

				glpkAPI::setObjDirGLPK(prob, glpkAPI::GLP_MIN);
                ## note: when FX or LO value taken from lb and when UP take from UP and when R
                rtype <- c(ifelse(gprtype=="E",glpkAPI::GLP_FX,glpkAPI::GLP_DB),rep(glpkAPI::GLP_DB,aux_cnt))
	       
	        # set the right hand side Sv = b

				out[[4]] <- glpkAPI::setRowsBndsGLPK(prob, c(1:num_constr), lb=rlower, ub=rupper,type=rtype )
        
	        # add upper and lower bounds: ai <= vi <= bi
                cc <- glpkAPI::setColsBndsObjCoefsGLPK(prob, c(1:num_var),   lower,   upper,  cobj   )
               # print(cc);			
        	#nzLHS=nzijr(LHS);
				TMPmat <- as(LHS, "TsparseMatrix")
                #cc <- loadMatrixGLPK(prob,   nzLHS$ne,   nzLHS$ia,  nzLHS$ja,  nzLHS$ar     )
				cc <- glpkAPI::loadMatrixGLPK(prob,length(TMPmat@x),ia  = TMPmat@i + 1,ja  = TMPmat@j + 1,ra  = TMPmat@x)

                ctype <-  rep(glpkAPI::GLP_BV,num_var); # number of integer variable
				glpkAPI::setColsKindGLPK(prob,c(1:(num_var)),ctype);

                #writeLPGLPK(prob,"FluxGeneGLPK.lp");
                #print("problem saved gplk");
                 if (verboseMode > 2) {                      
				            fname=format(Sys.time(), "findFluxGeneStat_glpk_%Y%m%d_%H%M.lp");
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
		newFlux=glpkAPI::mipColsValGLPK(prob);
	       ## --- 	       
               sol_geneStat=ifelse(newFlux[1:ng]==1,"ON","OFF");
                
            },
            # ----------------------- ----------------------- ------------------------- ----------------------- ----------------------#
           "cplexAPI" = {
                out <- vector(mode = "list", length = 3)
				prob <- cplexAPI::openProbCPLEX()
				out <- cplexAPI::setIntParmCPLEX(prob$env, cplexAPI::CPX_PARAM_SCRIND, cplexAPI::CPX_OFF)
		                
                
				# when R: rngval+rhs   to rhs   (rngval<0)
                rtype <- c(gprtype,rep("R",aux_cnt))
                #lp@oobj$lp<- initProbCPLEX(lp@oobj$env)
				#chgProbNameCPLEX(lp@oobj$env, lp@oobj$lp, "eFBA");
				prob$lp<- cplexAPI::initProbCPLEX(prob$env)
				cplexAPI::chgProbNameCPLEX(prob$env, prob$lp, "fluxGeneExpr cplex");
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
            
            ctype<- rep('B',num_constr);# integer variables
		    status = cplexAPI::copyColTypeCPLEX (prob$env, prob$lp, ctype);
		    check <- cplexAPI::setObjDirCPLEX(prob$env, prob$lp, cplexAPI::CPX_MIN);
		
	         if (verboseMode > 2) {  
				#print("writing problem... ");				
				fname=format(Sys.time(), "findFluxGeneStat_clex_%Y%m%d_%H%M.lp");
				print(sprintf("writing problem to: %s/%s...",getwd(),fname));
				cplexAPI::writeProbCPLEX(prob$env, prob$lp,fname);
			}
              #   --------------------------------------------------------------
           print("Solving...");
	       lp_ok     <- cplexAPI::mipoptCPLEX(prob$env, prob$lp);
	       print(lp_ok);
	       sol=cplexAPI::solutionCPLEX(prob$env, prob$lp);
               print(sol)
               lp_obj=sol$objval;
               lp_stat   <- getSolStat(prob$lp);
              
               #have flux and gene ON /flux and gene OFF?
               
               sol_geneStat=ifelse(sol$x[1:ng]==1,"ON","OFF");
				#=ifelse(newFlux[(nc+gpr_rxn+1) : (nc+gpr_rxn+ng)]==1,"ON","OFF");
			   # when using binary variables as indicators: .floorValues(.ceilValues(newFlux[(nc+1):(2*nc)],tol=Tf/10),tol=Tf);
           },
           # ----------------------- #
           {
               wrong_type_msg(solver)
            }
        )
        
return(list(Locus=allGenes(model),State=sol_geneStat))
}

# for testing: f1=eFBA_gene(allOff), f2=FBA: compare (count ON)