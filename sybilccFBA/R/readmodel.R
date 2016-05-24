# single quote make a problem
# mat=read.table("C:\\Users\\Abdelmonem\\Dropbox\\Sybil\\CCFBA\\MOMENT\\MOMENT\\mat.csv", sep=";",header=T,quote="") 
# mets=read.table("C:\\Users\\Abdelmonem\\Dropbox\\Sybil\\CCFBA\\MOMENT\\MOMENT\\mets.csv", sep=";",header=T,quote="") 
# rxns=read.table("C:\\Users\\Abdelmonem\\Dropbox\\Sybil\\CCFBA\\MOMENT\\MOMENT\\rxns.csv", sep=";",header=T,quote="") 
# rbnds=read.table("C:\\Users\\Abdelmonem\\Dropbox\\Sybil\\CCFBA\\MOMENT\\MOMENT\\rbnds.csv", sep=";",header=T,quote="")
# cbnds=read.table("C:\\Users\\Abdelmonem\\Dropbox\\Sybil\\CCFBA\\MOMENT\\MOMENT\\cbnds.csv", sep=";",header=T,quote="")
# load(Moment_mdl)
#mp  <- system.file(package = "sybilccFBA", "extdata")
 readmodel <- function(mat,mets,rxns,rbnds,cbnds,solver="glpkAPI"){
nr=6705
nc=4991
nf=3234
nm=1674# 6 metabolites
   
   LHS <- Matrix::Matrix(0, nrow = nr, ncol = nc  )
   
   for (i in c(1:length(mat[,1]))){
     LHS[mat[i,1],mat[i,2]]=mat[i,3]
   }
 #map matrix
 cobj=c(rxns[,"ocf"],rep(0,nc-nf))
 
######
if(solver=="cplexAPI"){
				prob <- cplexAPI::openProbCPLEX()
				out <- cplexAPI::setIntParmCPLEX(prob$env, cplexAPI::CPX_PARAM_SCRIND, cplexAPI::CPX_OFF)
		                
                cplexAPI::chgProbNameCPLEX(prob$env, prob$lp, "Moment cplex");
                #E,L,"G"
                rtype <- c(rep("E",nm),rep("L",nr-nm))
                cplexAPI::setObjDirCPLEX(prob$env, prob$lp, cplexAPI::CPX_MAX);
				rupper=c(rbnds[1:nm,1],rbnds[c((nm+1):nr),2])
                rupper[nr]=0.27 # set to limit in paper 48% of protein is metabolic
				cplexAPI::newRowsCPLEX(prob$env, prob$lp,
                                         nrows=nr, rhs=rupper, sense=rtype)
				upper=cbnds[,2]
				lower=cbnds[,1]
	upper[2609]=0#ac
	upper[2729]=1000#glc	
upper[2835]=0#pyruvate	
upper[2705]=0 #pr_index= find(strcmp(model.rxns ,'REV_EX_fru(e)'))
upper[2774]=0#pr_index= find(strcmp(model.rxns ,'REV_EX_lac_L(e)'))
#pr_index= find(strcmp(model.grRules ,'(b0003)'))
#2723
#pr_index= find(strcmp(model.rxns ,'REV_EX_galt(e)'))
                cplexAPI::newColsCPLEX(prob$env, prob$lp,
                                         nc, obj=cobj, lb=lower, ub=upper)#,cnames=cNames)
				
				print(sprintf("%s : step 2: nzijr....",format(Sys.time(), "%d-%m-%Y %X")))
                        # constraint matrix
				TMPmat <- as(LHS, "TsparseMatrix")
				cplexAPI::chgCoefListCPLEX(prob$env, prob$lp,#lp@oobj@env, lp@oobj@lp,
                                   nnz = length(TMPmat@x),
                                   ia  = TMPmat@i,
                                   ja  = TMPmat@j,
                                   ra  = TMPmat@x);
							
						#                nzLHS=nzijr();#LHS    
#print(sprintf("%s : step 3: chngCoef....",format(Sys.time(), "%d-%m-%Y %X")))				
                      fname=format(Sys.time(), "Cplex_moment_%Y%m%d_%H%M.lp");
					  print(sprintf("Writing problem to file: %s/%s  ...",getwd(),fname));
					  cplexAPI::writeProbCPLEX(prob$env, prob$lp, fname);

               lp_ok     <- cplexAPI::lpoptCPLEX(prob$env, prob$lp);
	           print(lp_ok);
				sol=cplexAPI::solutionCPLEX(prob$env, prob$lp);
				print(sprintf("GLC upt=%f, AC=%f Pyr=%f fruc=%f Lac=%f",sol$x[2729],sol$x[2609],sol$x[2835],sol$x[2705],sol$x[2774]))
				colst=sol$x;
	       }else{
		          prob <- glpkAPI::initProbGLPK();# new problem
               
                 glpkAPI::addRowsGLPK(prob, nrows=nr)
				outj <- glpkAPI::addColsGLPK(prob, ncols=nc)
				#setColNameGLPK(prob,c(1:(nc+nIrrev+2*nRev )),cNames );
				glpkAPI::setObjDirGLPK(prob, glpkAPI::GLP_MAX);
                ## note: when FX or LO value taken from lb and when UP take from UP
               #rtype <- c(rep(glpkAPI::GLP_FX, nr),glpkAPI::GLP_LO, rep(glpkAPI::GLP_UP, nIrrev+2*nRev  ))# objective >=
	             rtype <- c(rep(glpkAPI::GLP_FX,nm),rep(glpkAPI::GLP_UP,nr-nm))
                
	        # set the right hand side Sv = b
				
                rlower=c(rbnds[1:nm,1],rbnds[(nm+1):nr,2])#use ub for UP type
				rupper=rbnds[1:nr,2]
				rupper[nr]=0.27
				#out[[4]] <- glpkAPI::setRowsBndsGLPK(prob, c(1:nr), lb=rlower, ub=rupper,type=rtype )
				glpkAPI::setRowsBndsGLPK(prob, c(1:nr), lb=rlower, ub=rupper,type=rtype )
	        # add upper and lower bounds: ai <= vi <= bi
				upper=cbnds[,2]
				lower=cbnds[,1]
	# upper[2609]=0#ac
	# upper[2729]=0#glc	
	# upper[2835]=0#pyruvate	
# upper[2705]=0 #pr_index= find(strcmp(model.rxns ,'REV_EX_fru(e)'))
# upper[2774]=0#pr_index= find(strcmp(model.rxns ,'REV_EX_lac_L(e)'))
# #pr_index= find(strcmp(model.grRules ,'(b0003)'))
# upper[2723]=0#pr_index= find(strcmp(model.rxns ,'REV_EX_galt(e)'))
# upper[2718]=0#galactose
# upper[2783]=0#maltose
# upper[2709]=0#fum
# upper[2741]=0#glycerol
# upper[2782]=1000#which(rxns$react_id=="REV_EX_mal_L(e)")

#moment KO
# upper[1479]=lower[1479]=0#LACZ
# upper[1489]=lower[1489]=0#
# upper[3023]=lower[3023]=0#pr_index= find(strcmp(model.rxns ,'REV_LCTStpp'))

	cc <- glpkAPI::setColsBndsObjCoefsGLPK(prob, c(1:nc),   lower,   upper,  cobj   )
                
				TMPmat <- as(LHS, "TsparseMatrix")
                cc <- glpkAPI::loadMatrixGLPK(prob,length(TMPmat@x),ia  = TMPmat@i + 1,ja  = TMPmat@j + 1,ra  = TMPmat@x)
               
               #if (verboseMode > 2) {                      
					fname=format(Sys.time(), "glpk_eFBA_%Y%m%d_%H%M.lp");
					print(sprintf("Writing problem to file: %s/%s ...",getwd(),fname));
                	glpkAPI::writeLPGLPK(prob,fname);
                	print(format(Sys.time(), "Testing time : %Y%m%d %X Solving..."));
                #}
				
				## Solve
            #glpkAPI::setMIPParmGLPK(glpkAPI::PRESOLVE,glpkAPI::GLP_ON);
			lp_ok=glpkAPI::solveSimplexGLPK(prob);
			glpkAPI::return_codeGLPK(lp_ok);
			lp_stat=glpkAPI::getSolStatGLPK(prob);
			glpkAPI::status_codeGLPK(lp_stat);
			lp_obj=glpkAPI::getObjValGLPK(prob);
			colst=glpkAPI::getColsPrimGLPK(prob);
			newFlux=colst
			print(sprintf("GLC upt=%f, AC=%f Pyr=%f fruc=%f Lac=%f galt=%f",colst[2729],colst[2609],colst[2835],colst[2705],colst[2774]
			,colst[2723]))
	 
#	  write.csv(file="iAFmomentKO.csv",cbind(colst[1:3234],rxns))
#	  write.csv(file="gi.csv",colst[3235:4495])
	  #fd1flx=getRevFlux(model,mod2,sol[[1]]$fluxes[1:react_num(mod2)])
		   }
		   
return(colst);		   
# aa=LHS[6705,]
# bb=aa[aa!=0]
# sum(bb * sol$x[3247:4506])

}