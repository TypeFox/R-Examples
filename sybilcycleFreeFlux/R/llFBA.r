## llFBA
llFBA <- function(model
	 ,lpdir = SYBIL_SETTINGS("OPT_DIRECTION")
	 ,solver = SYBIL_SETTINGS("SOLVER")
         ,method = SYBIL_SETTINGS("METHOD")
         ,solverParm=data.frame(CPX_PARAM_EPRHS=1e-7)
		,verboseMode = 2
) {
# find FBA solution without cycles and with minor changes to fluxs that are not in cycles
# Output FBA solution > different from mintotal flux 
# cycles can be enumerated

##  features to be added : Same obj as FBA
##  return optObj
# solve MILP to get loopless FBA

# nzijr<-function(LHS){
	                # A=as.matrix(LHS);
				# eps=1e-10;
				# ni=dim(A)[1]; nj=dim(A)[2];
				# ne=0; ia=NULL; ja=NULL; ra=as.double(NULL);
				# for (i in 1:ni){
				# for(j in 1:nj){
				 # if (abs(A[i,j])>eps){
				   # ne=ne+1;
				   # ia=c(ia,i);
				   # ja=c(ja,j);
				   # ra=c(ra,A[i,j]); 
				 # }
				# }
				# }
		                # nzLHS=list(ne=ne,ia=ia,ja=ja,ar=ra);
# return(nzLHS);
# }

        if (!is(model, "modelorg")) {
            stop("needs an object of class modelorg!")
        }
    
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
        #            |      |
        #         -1 | -ub  | <= 0    ] 
        #            |      |           } -ub*yi <= vi <= ub * yi
        #         +1 | -ub  | <= 0    ]
        #            |      |
        #       ------------------
        #  lb   wt_lb|  0   |
        #  ub   wt_ub|  1   |
 
    ## add constraints loopless
    #excReact = findExchReact(model)[1];# 1 is position
    #excReactPos=react_pos(excReact$exchange);
    excReact = findExchReact(model);
	excReactPos=which(react_id(model) %in% react_id(excReact));#excReact$exchange

	is_excReact=((c(1:length(react_id(model)))) %in% excReactPos);
	active=!(lowbnd(model)==0 & uppbnd(model)==0)
	

	loopicias = (!is_excReact) & active;
    Sn=S(model)[,loopicias ];
    #Ninternal = sparseNull(sparse(Sn));
	print("Calculating NULL matrix..");
	Nint=MASS::Null(t(Sn))  # transpose used here to give the same result as Matlab fun (I don't know the difference)
	# when the same input is used it gives different result
	
	nc     <- react_num(model)
	nr     <- met_num(model)
	ng=dim(Sn)[2] # nonExchng active rxns
	Nintc = dim(Nint)[2]
 	#Nintr  = dim(Nint)[1]
 	
  nRows = nr+4*ng+Nintc
  nCols=nc+2*ng ## vars: v|G|a

##1- A [S  
 print("Building LHS...");
  LHS <- Matrix::Matrix(0, nrow = nRows, ncol = nCols, sparse = TRUE)
        	#as.matrix.csr(0, nrow = nRows, ncol = nCols)
          
          # rows for S
           LHS[1:nr,1:(nc)] <- S(model)
         # Nint' G=0
          LHS[(nr+1):(nr+Nintc),(nc+1):(nc+ng)] <- t(Nint)
         
         # which columns in Sint are in Nint
         # v<=1000 a
           ii=matrix(c((nr+Nintc+1)   :(nr+Nintc+ng)  ,which(loopicias) ),ncol=2)
	   LHS[ii ] <- 1;
	   diag(LHS[((nr+Nintc+1)   :(nr+Nintc+ng) ) , ((nc+ng+1):(nc+2*ng)) ] ) = -1000;
	         
	         #          -v+1000a<=1000  (if v is negative a must be zero)                                     
	 ii=matrix(c((nr+Nintc+ng+1)   :(nr+Nintc+2*ng)  ,which(loopicias) ),ncol=2)
		   LHS[ii ] <- -1;
	   diag(LHS[((nr+Nintc+ng+1)   :(nr+Nintc+2*ng) ) , ((nc+ng+1):(nc+2*ng)) ] ) = 1000;

      	   # G+1001 a<=1000
      	   diag(LHS[((nr+Nintc+2*ng+1)   :(nr+Nintc+3*ng) ) , ((nc+1):(nc+ng)) ] ) = 1;#G
       	   diag(LHS[((nr+Nintc+2*ng+1)   :(nr+Nintc+3*ng) ) , ((nc+ng+1):(nc+2*ng)) ] ) = 1001;#a
       	    # -G - 1001 a <=-1
	   diag(LHS[((nr+Nintc+3*ng+1)   :(nr+Nintc+4*ng) ) , ((nc+1):(nc+ng)) ] ) = -1;#G
       	   diag(LHS[((nr+Nintc+3*ng+1)   :(nr+Nintc+4*ng) ) , ((nc+ng+1):(nc+2*ng)) ] ) = -1001;#a
      # bounds ---***********-------------********-----------
       	   lower  <- c(lowbnd(model), rep(-1000, ng),rep(0, ng))
            upper  <- c(uppbnd(model), rep(1000, ng),rep(1, ng))
            
            rlower <- c(rep(0,nr),rep(0, Nintc) ,rep(-1000, ng) ,rep(-1000, ng),rep(-1000, ng),rep(-2001, ng))
            rupper <- c(rep(0,nr),rep(0, Nintc) ,rep(0, ng) ,rep(1000, ng),rep(1000, ng),rep(-1, ng))

		#write.csv(file="bnd.csv",cbind(rlower=rlower,upp=rupper));
        # ---------------------------------------------
        # objective function
        # ---------------------------------------------

        cobj <- c(obj_coef(model), rep(0, 2*ng))
	

	#rtype <- c(rep(glpkAPI::GLP_FX, nr+Nintc), rep(glpkAPI::GLP_LO, 3*ng))
	#ctype <- c(rep(GLP_CV, nc), rep(GLP_BV,2*ng ));


   cNames=paste(c(rep("x",nc),rep("G",ng),rep("a",ng)),
   		 		 		               c(1:nc,which(loopicias),which(loopicias)),sep="_" ) ;
   		 		 		               
  # ---------------------------------------------
        # build problem object
        # ---------------------------------------------
print("Building the problem...");
        switch(solver,
            # ----------------------- #
            "glpkAPI" = {
                out <- vector(mode = "list", length = 4)
                prob <- glpkAPI::initProbGLPK();# new problem
               
                out[[1]] <-  glpkAPI::addRowsGLPK(prob, nrows=nRows)
		outj <- glpkAPI::addColsGLPK(prob, ncols=nCols )
		glpkAPI::setColNameGLPK(prob,c(1:(nCols )),cNames );
		
		glpkAPI::setObjDirGLPK(prob, glpkAPI::GLP_MAX);
		
                ## note: when FX or LO value taken from lb and when UP take from UP
                rtype <- c(rep(glpkAPI::GLP_FX, nr+Nintc), rep(glpkAPI::GLP_UP, 4*ng))
	       	ctype <- c(rep(glpkAPI::GLP_CV, nc), rep(glpkAPI::GLP_CV, ng),rep(glpkAPI::GLP_BV,ng ));

	        # set the right hand side Sv = b

	       out[[4]] <- glpkAPI::setRowsBndsGLPK(prob, c(1:(nCols)), lb=rlower, ub=rupper,type=rtype )
        
	        # add upper and lower bounds: ai <= vi <= bi
                cc <- glpkAPI::setColsBndsObjCoefsGLPK(prob, c(1:(nCols )),   lower,   upper,  cobj   )
                if (verboseMode > 2) { print(cc);	}		
        	#nzLHS=nzijr(LHS);
        	#print(nzLHS);
             #   cc <- glpkAPI::loadMatrixGLPK(prob,   nzLHS$ne,   nzLHS$ia,  nzLHS$ja,  nzLHS$ar     )
                 TMPmat <- as(LHS, "TsparseMatrix")
                cc <- glpkAPI::loadMatrixGLPK(prob,length(TMPmat@x),ia  = TMPmat@i + 1,
				ja  = TMPmat@j + 1,ra  = TMPmat@x)
                             
			  #print(cc);
                #ctype <- c(rep(GLP_CV, nc), rep(GLP_BV,nIrrev+2*nRev ));
		glpkAPI::setColsKindGLPK(prob,c(1:(nCols )),ctype);

               if (verboseMode > 2) {                      
		        fname=format(Sys.time(), "glpk_llFBA_%Y%m%d_%H%M.lp");
			print(sprintf("Writing problem to file: %s/%s ...",getwd(),fname));
                	glpkAPI::writeLPGLPK(prob,fname);
                	print("Solving...");
                }
                ## Solve
                	glpkAPI::setMIPParmGLPK(glpkAPI::PRESOLVE,glpkAPI::GLP_ON);
		lp_ok=	glpkAPI::solveMIPGLPK(prob);
			glpkAPI::return_codeGLPK(lp_ok);
		lp_stat=	glpkAPI::mipStatusGLPK(prob);
			glpkAPI::status_codeGLPK(lp_stat);
		lp_obj=	glpkAPI::mipObjValGLPK(prob);
		colst=	glpkAPI::mipColsValGLPK(prob);
		newFlux=colst
	       ## --- 
	            

            },
            # ----------------------- #
           "cplexAPI" = {
                out <- vector(mode = "list", length = 3)
		prob <- openProbCPLEX()
		setIntParmCPLEX(prob$env, CPX_PARAM_SCRIND, CPX_OFF)
		                
                chgProbNameCPLEX(prob$env, prob$lp, "llFBA cplex");
                
                rtype <- c(rep("E", nr+Nintc), rep("L", 4*ng))
                #rtype <- c(rep("E",nr+1), rep("L", 2*nIrrev+4*nRev  ))
                

                out[[1]] <- newRowsCPLEX(prob$env, prob$lp,
                                         nrows=nRows, rhs=rupper, sense=rtype)

                out[[2]] <- newColsCPLEX(prob$env, prob$lp,
                                         nCols , obj=cobj, lb=lower, ub=upper,cnames=cNames)

                print("Calc nzLHS");
				#nzLHS=nzijr(LHS);        
                # out[[3]] <- chgCoefListCPLEX(prob$env, prob$lp,
                                             # nzLHS$ne,
                                             # nzLHS$ia  - 1,
                                            # nzLHS$ja - 1,
                                             # nzLHS$ar)
                 TMPmat <- as(LHS, "TsparseMatrix")
				cplexAPI::chgCoefListCPLEX(prob$env, prob$lp,#lp@oobj@env, lp@oobj@lp,
                                   nnz = length(TMPmat@x),
                                   ia  = TMPmat@i,
                                   ja  = TMPmat@j,
                                   ra  = TMPmat@x);
				
  		ctype<-c(rep('C', nc+ng), rep('B',ng ));
		status = copyColTypeCPLEX (prob$env, prob$lp, ctype);
		check <- setObjDirCPLEX(prob$env, prob$lp, CPX_MAX);# should use same as lpdir
		
                if (verboseMode > 2) {                      
                      fname=format(Sys.time(), "Cplex_llFBA_%Y%m%d_%H%M.lp");
			 print(sprintf("Writing problem to file: %s/%s  ...",getwd(),fname));
			writeProbCPLEX(prob$env, prob$lp, fname);
		      
		       print("Solving...");
		       }
	       lp_ok     <- mipoptCPLEX(prob$env, prob$lp);
	       print(lp_ok);
	       sol=solutionCPLEX(prob$env, prob$lp);
	       
               if (verboseMode > 3) {print(sol);}	
               
               lp_obj=sol$objval;
               lp_stat   <- getStatCPLEX(prob$env, prob$lp)
              
                #colst=sol$x;
                
               #newFlux=.floorValues(sol$x,tol=Tf);
               newFlux=sol$x;#newFlux[1:nc];
               #newStat=ifelse(abs(newFlux)>Tf,1,0);
               # when using binary variables as indicators: .floorValues(.ceilValues(newFlux[(nc+1):(2*nc)],tol=Tf/10),tol=Tf);
           },
		   "sybilGUROBI" = {
                out <- vector(mode = "list", length = 3)
		#prob <- openProbCPLEX()
		prob<-initProb()
	#	setIntParmCPLEX(prob$env, CPX_PARAM_SCRIND, CPX_OFF)
		                
               # chgProbNameCPLEX(prob$env, prob$lp, "llFBA cplex");
                
                rtype <- c(rep("E", nr+Nintc), rep("L", 4*ng))
                #rtype <- c(rep("E",nr+1), rep("L", 2*nIrrev+4*nRev  ))
                

                #out[[1]] <- newRowsCPLEX(prob$env, prob$lp,
                 #                        nrows=nRows, rhs=rupper, sense=rtype)
				 #out[[1]] <- addRowsToProb(prob,nRows,type=rtype,lb=,ub=rupper,)

                out[[2]] <- newColsCPLEX(prob$env, prob$lp,
                                         nCols , obj=cobj, lb=lower, ub=upper,cnames=cNames)

                print("Calc nzLHS");
				#nzLHS=nzijr(LHS);        
                # out[[3]] <- chgCoefListCPLEX(prob$env, prob$lp,
                                             # nzLHS$ne,
                                             # nzLHS$ia  - 1,
                                            # nzLHS$ja - 1,
                                             # nzLHS$ar)
                 TMPmat <- as(LHS, "TsparseMatrix")
				cplexAPI::chgCoefListCPLEX(prob$env, prob$lp,#lp@oobj@env, lp@oobj@lp,
                                   nnz = length(TMPmat@x),
                                   ia  = TMPmat@i,
                                   ja  = TMPmat@j,
                                   ra  = TMPmat@x);
				
  		ctype<-c(rep('C', nc+ng), rep('B',ng ));
		status = copyColTypeCPLEX (prob$env, prob$lp, ctype);
		check <- setObjDirCPLEX(prob$env, prob$lp, CPX_MAX);# should use same as lpdir
		
                if (verboseMode > 2) {                      
                      fname=format(Sys.time(), "Cplex_llFBA_%Y%m%d_%H%M.lp");
			 print(sprintf("Writing problem to file: %s/%s  ...",getwd(),fname));
			writeProbCPLEX(prob$env, prob$lp, fname);
		      
		       print("Solving...");
		       }
	       lp_ok     <- mipoptCPLEX(prob$env, prob$lp);
	       print(lp_ok);
	       sol=solutionCPLEX(prob$env, prob$lp);
	       
               if (verboseMode > 3) {print(sol);}	
               
               lp_obj=sol$objval;
               lp_stat   <- getStatCPLEX(prob$env, prob$lp)
              
                #colst=sol$x;
                
               #newFlux=.floorValues(sol$x,tol=Tf);
               newFlux=sol$x;#newFlux[1:nc];
               #newStat=ifelse(abs(newFlux)>Tf,1,0);
               # when using binary variables as indicators: .floorValues(.ceilValues(newFlux[(nc+1):(2*nc)],tol=Tf/10),tol=Tf);
           },
           
           # ----------------------- #
           {
               wrong_type_msg(solver)
            }
        )
           		 		 		               
 remove(prob);
# ----------------------------------------  ***  -------------------------------- #
	optsol <- list(        ok = lp_ok,
			           stat = lp_stat,
			            obj=lp_obj,
			       fluxes = newFlux[1:nc],
			       G = newFlux[(nc+1):(nc+ng)],
			       a = newFlux[(nc+ng+1):(nc+2*ng)],
			       cNames=cNames
			  )
	return(optsol);
}


