# Given a wildtype flux (some fluxes can be absent i.e NA), find FBA solution of network given by model,
# such that Sum(|v_i-v_wt|) is minimal (like lMOMA but with option to exclude some reactions)

#find Manhaten-Distance-Closest Flux
findMDCFlux <- function(model
        ,wtflux  # desired flux distribution when wtflux is NA don't include its constraint 
	,objVal = NA #of objval
	, pct_objective=100
	 ,lpdir = SYBIL_SETTINGS("OPT_DIRECTION")
	 ,solver = SYBIL_SETTINGS("SOLVER")
         ,method = SYBIL_SETTINGS("METHOD")
         ,solverParm=data.frame(CPX_PARAM_EPRHS=1e-6)
		,verboseMode = 2) {

        if (!is(model, "modelorg")) {
            stop("needs an object of class modelorg!")
        }
      
      if(is.na(objVal)){
      	sol=optimizeProb(model,solver =solver,method=method,solverParm=solverParm);
      	objVal=lp_obj(sol)*pct_objective/100;
      }
      
        out <- FALSE

        nc     <- react_num(model)
        nr     <- met_num(model)
        nd    <- sum(!is.na(wtflux))
        absMAX <- SYBIL_SETTINGS("MAXIMUM")
        
        nRows=nr+2*nd+1
		nCols=nc+2*nd
	
        #  the problem: minimize 
        #  
        #         S  | 0 |  0   |  = 0
        #         1   | 1 |  0   | >= v_wt 
        #         -1  | 0 |  1   | >= -v_wt 
        #       -------------------------
        #       c  |  0   |  0   | >= old_obj*pct
        #       delta-,delta+  >=0
        #          
        #        obj    0  |  1   |  1   |

        # ---------------------------------------------
        # constraint matrix
        # ---------------------------------------------
	
        # the initial matrix dimensions
     	#LHS <- Matrix::Matrix(0, nrow = nr+2*nc, ncol = 4*nc, sparse = TRUE)
     	LHS <- Matrix::Matrix(0, nrow = nRows, ncol = nCols, sparse = TRUE)
     	#as.matrix.csr(0, nrow = nRows, ncol = nCols)
       
       # rows for S
        LHS[1:nr,1:(nc)] <- S(model)

        # rows for the delta match matrix
			#         1   | 1 |  0   | >= v_wt 

	#diag(LHS[(nr+1)   :(nr+nd),1       :(nc)]) <- 1
	ii=matrix(c((nr+1)   :(nr+nd), which(!is.na(wtflux))),ncol=2)
	LHS[ii ]<-1
	if(nd==1) {LHS[(nr+1)   :(nr+nd),(nc+1)       :(nc+nd)] <- 1;# dia function generates an error when applied to 1 row
	}else{diag(LHS[(nr+1)   :(nr+nd),(nc+1)       :(nc+nd)]) <- 1}
	
	#   -1  | 0 |  1   | >= -v_wt
	#diag(LHS[(nr+nd+1)   :(nr+2*nd),1 :nc ]) <- -1	
	ii=matrix(c((nr+nd+1)   :(nr+2*nd),which(!is.na(wtflux))),ncol=2)
	LHS[ii  ] <- -1
	if(nd==1) {LHS[(nr+nd+1)   :(nr+2*nd),(nc+nd+1) :(nc+2*nd) ] <- 1
    }else{diag(LHS[(nr+nd+1)   :(nr+2*nd),(nc+nd+1) :(nc+2*nd) ]) <- 1}		
       		#          c  |  0   |  0   | >= old_obj*pct
        LHS[(nr+2*nd+1),1       :nc    ] <- obj_coef(model)
        
        if(verboseMode>2) print(sprintf("nrows:%d, ncols=%d, nd=%d,nc=%d,nr=%d",nRows,nCols,nd,nc,nr));
        
        # ---------------------------------------------
        # lower and upper bounds
        # ---------------------------------------------

            # Here, we keep the wild type flux distribution fixed.
            lower  <- c(lowbnd(model), rep(0, 2*nd))
            upper  <- c(uppbnd(model), rep(absMAX, 2*nd))
            
            rlower <- c(rep(0, nr), wtflux[!is.na(wtflux)],-wtflux[!is.na(wtflux)],objVal)
            rupper <- c(rep(0, nr), rep(absMAX, 2*nd+1))

		#write.csv(file="bnd.csv",cbind(rlower=rlower,upp=rupper));
        # ---------------------------------------------
        # objective function
        # ---------------------------------------------

        cobj <- c(rep(0, nc), rep(1, 2*nd))
        
        # ---------------------------------------------
        
        cNames=paste(c(rep("x",nc),rep("dp",nd),rep("dn",nd)),
       		 		 		               c(1:nc,which(!is.na(wtflux)),which(!is.na(wtflux))),sep="_" ) ;
               # ---------------------------------------------
       

        # ---------------------------------------------
        # build problem object
        # ---------------------------------------------
        
        switch(solver,
            # ----------------------- #
            "glpkAPI" = {
                out <- vector(mode = "list", length = 5)
                 prob <- glpkAPI::initProbGLPK();# new problem
              
                    rtype <- c(rep(glpkAPI::GLP_FX, nr), rep(glpkAPI::GLP_LO, 2*nd))
                    if (lpdir == "max") {
                        rtype <- c(rtype, glpkAPI::GLP_LO)
                    }
                    else {
                        rtype <- c(rtype, glpkAPI::GLP_UP)
                    }
                  

                #nzLHS    <- nonZeroElements(LHS)
                TMPmat <- as(LHS, "TsparseMatrix")
                
                #out[[1]] <- addRowsCols(prob, nRows, nCols)
                out[[1]] <-  glpkAPI::addRowsGLPK(prob, nrows=nRows)
		outj <- glpkAPI::addColsGLPK(prob, ncols=nCols)
                #glpkAPI::setColNameGLPK(prob,c(1:(nCols)),paste(c(rep("x",nc),rep("dn",nd),rep("dp",nd)),
				 #              c(1:nc,which(!is.na(wtflux)),which(!is.na(wtflux))),sep="" ) );
				 mapply(glpkAPI::setColNameGLPK, j = c(1:nCols), cname = cNames, MoreArgs = list(lp = prob));

				glpkAPI::setObjDirGLPK(prob, glpkAPI::GLP_MIN);
 
                out[[2]] <- glpkAPI::loadMatrixGLPK(prob,
                                       length(TMPmat@x),
				      TMPmat@i + 1,
				      TMPmat@j + 1,
                                       TMPmat@x
                                      )
            
                # add upper and lower bounds: ai < vi < bi
                out[[3]] <- glpkAPI::setColsBndsObjCoefsGLPK(prob,
                                                c(1:nCols),
                                                lower,
                                                upper,
                                                cobj
                                               )
            
                # set the right hand side Sv = b
                out[[4]] <- glpkAPI::setRowsBndsGLPK(prob,
                                            c(1:nRows),
                                            rlower,
                                            rupper,
                                            rtype
                                           )
	 	  #set precision              
                  parm <- sapply(dimnames(solverParm)[[2]],
                               function(x) eval(parse(text = x))) 
                val  <- solverParm[1,]
                if (method == "interior") {
                    glpkAPI::setInteriorParmGLPK(parm, val)
                    out[[5]] <- TRUE
                }
                else {
                    glpkAPI::setSimplexParmGLPK(parm, val)
                    out[[5]] <- TRUE
                }
	  
             if (verboseMode > 2) {                      
	    		        fname=format(Sys.time(), "glpk_ManhatDist_%Y%m%d_%H%M.lp");
	    			print(sprintf("write problem: %s/%s",getwd(),fname));
	                    	glpkAPI::writeLPGLPK(prob,fname);
	                    	print("Solving...");
                }
                
                #----Solving----------------------
                lp_ok <- glpkAPI::solveSimplexGLPK(prob)
                lp_obj   <- glpkAPI::getObjValGLPK(prob)
		lp_stat   <- glpkAPI::getSolStatGLPK(prob)
		if (is.na(lp_stat)) {
		    lp_stat <- lp_ok
		}
		
		lp_fluxes <- glpkAPI::getColsPrimGLPK(prob)
               # ----------------------- #
            },# end switch case glpk
                 # ----------------------- #
            "cplexAPI" = {
                out <- vector(mode = "list", length = 4)
		prob <- cplexAPI::openProbCPLEX()
                out <- cplexAPI::setIntParmCPLEX(prob$env, cplexAPI::CPX_PARAM_SCRIND, cplexAPI::CPX_OFF)# why?
                
                cplexAPI::chgProbNameCPLEX(prob$env, prob$lp, "ManhatenDist cplex");
                    
                 cplexAPI::setObjDirCPLEX(prob$env, prob$lp, cplexAPI::CPX_MIN)
                 rtype <- c(rep("E", nr), rep("G", 2*nd),"G");# Biomass >=
                    
                                
                #nzLHS    <- nonZeroElements(LHS)
  		TMPmat <- as(LHS, "TsparseMatrix")

                out[[1]] <- cplexAPI::newRowsCPLEX(prob$env, prob$lp,
                                         nRows, rlower, rtype)
                                         
                out[[2]] <- cplexAPI::newColsCPLEX(prob$env, prob$lp,
                                         nCols, cobj, lower, upper)
			# TMPmat
                out[[3]] <- cplexAPI::chgCoefListCPLEX(prob$env, prob$lp,
                                            length(TMPmat@x),
                                       TMPmat@i ,
                                       TMPmat@j ,
                                       TMPmat@x)
                         
                         #set precision 
                         parm <- sapply(dimnames(solverParm)[[2]],
                                    function(x) eval(parse(text = x)))
                 out[[4]]  <- cplexAPI::setDblParmCPLEX(prob$env, parm, solverParm);
                 
                if (verboseMode > 2) {                      
		                      fname=format(Sys.time(), "Cplex_ManhatDist_%Y%m%d_%H%M.lp");
					 print(sprintf("write problem: %s/%s",getwd(),fname));
					cplexAPI::writeProbCPLEX(prob$env, prob$lp,fname);
				      
				       print("Solving...");
		       }
		 
		 
		 #----Solving---------------------- 
		lp_ok    <- cplexAPI::lpoptCPLEX(prob$env, prob$lp)
		lp_obj   <- cplexAPI::getObjValCPLEX(prob$env, prob$lp)
		lp_stat   <- cplexAPI::getStatCPLEX(prob$env, prob$lp)
		if (is.na(lp_stat)) {
		    lp_stat <- lp_ok
		}
		lp_fluxes <- cplexAPI::getProbVarCPLEX(prob$env, prob$lp, 0, nCols - 1)
               # ----------------------- #
            },
            
            
            # -----------default: switch solver ------------ #
            {
                wrong_type_msg(solver)
            }
        )  # end switch solver

#--------------Output-------------------------

        optsol <- list(ok = lp_ok,
	                       obj = lp_obj,
	                       stat = lp_stat,
	                       fluxes = lp_fluxes,
	                       mdcflx=lp_fluxes[1:nc],
	                       wtflx=wtflux
	                  )
	    
    return(optsol)
        
    }
