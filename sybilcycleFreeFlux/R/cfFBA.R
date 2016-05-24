cfFBA <- function(model
,wtflux  # initial flux distribution when wtflux is NA don't include its constraint 
	,objVal = NA #min objval
	,fixExchRxn=TRUE
	,excReactPos=NULL
	 ,lpdir = SYBIL_SETTINGS("OPT_DIRECTION")
	 ,solver = SYBIL_SETTINGS("SOLVER")
         ,method = SYBIL_SETTINGS("METHOD")
        # ,solverParm=data.frame(CPX_PARAM_EPRHS=1e-7)
		,verboseMode = 2
######### ADDED BY GABRIEL ###############
		,retOptSol = TRUE
##########################################
) {
# fix the direction of rxn, are not allowed to increase 
# minimize abs(flux)
# change lower and upper bounds only and obj fun
# forcing the same direction is a problem in some loops reversible rxns flux is c-1000
#N.B:23/10/2012, Keep exchange reactions fixed

        if (is(model, "modelorg")) {
            #stop("needs an object of class modelorg!")
            mod=sysBiolAlg(model,solver =solver,method=method)#,solverParm=solverParm);
        }  
        else  if(is(model, "sysBiolAlg")) {
	 	mod=model;
	 	}
	 else{
	 	stop("lrFBA:needs an object of class modelorg or sysBiolAlg!")
        }
         
         prob=problem(mod);
         
      if(is.na(objVal)){
      	sol=optimizeProb(mod);
      	objVal=sol$obj;#lp_obj(sol);
      }      
#####
nCols=nc(mod)#=react_num(model);
lb=getColsLowBnds(prob,c(1:nCols))#lowbnd(model);
ub=getColsUppBnds(prob,c(1:nCols))#uppbnd(model);
ocf=rep(0,nCols);
# lb=max(min(wtflx,0),lb)    ub=min(max(0,wt),ub)
for (i in 1:nCols ){
	if(!is.na(wtflux[i])){
		if(wtflux[i]<0){ # flux is -ve
			lb[i]=wtflux[i]; # old lb must be -ve also and |lb|>=|wtflx| 
			ub[i]=min(0,ub[i]); # old ub may be -ve, should not be relaxed
		}
		else #if(wtflux[i]>=0)
		{
			lb[i]=max(0,lb[i]);# made a problem in ATPM, should not be relaxed to 0
			ub[i]=wtflux[i]; 
		}
		ocf[i]=sign(wtflux[i]);
		#print(sprintf("",i,lb[i],wtflux[i]))
	}
}        

# force obective constraint in model
mod_ocfs=getObjCoefs(prob,c(1:nCols))
bmrxn=(mod_ocfs==1);
        
#  obj_coef(model)=ocf;
#   lowbnd(model)=lb;
#   uppbnd(model)=ub;

if(fixExchRxn){  # Fixing exchange rxn from changing may preserve some rxns with biological evidence from changing
	if(is.null(excReactPos) && is(model, "modelorg")){
		excReact = findExchReact(model);
		excReactPos=which(react_id(model) %in% react_id(excReact));#excReact$exchange
	}
 	lb[excReactPos]=wtflux[excReactPos]
	ub[excReactPos]=wtflux[excReactPos]
}
# Fix biomass
lb[bmrxn]=objVal;
ub[bmrxn]=objVal;
   
   setObjDir(prob,"min");
 sol=optimizeProb(mod,react=c(1:nCols),lb=lb,ub=ub,obj_coef =ocf)
   #### testing: problem in no of rows
# 
 #--------------Output-------------------------

#        optsol <- list(ok = sol$ok,
#	                       obj = sol$obj,
#	                       stat =sol$stat,
#	                       fluxes = sol$fluxes,
#	                       wtflx=wtflux,
#	                       lb=lb,
#	                       ub=ub,
#	                       ocf=ocf
#	                  )	    



######### ADDED BY GABRIEL ###############
if (isTRUE(retOptSol)) {
            optsol <- new("optsol_optimizeProb",
                          mod_id       = mod_id(model),
                          mod_key      = mod_key(model),
                          solver       = solver(prob),
                          method       = method(prob),
                          algorithm    = algorithm(mod),
                          num_of_prob  = 1L,
                          lp_dir       = factor(getObjDir(prob)),
                          lp_num_rows  = nr(mod),
                          lp_num_cols  = nc(mod),
                          lp_ok        = as.integer(sol[["ok"]]),
                          lp_obj       = sol[["obj"]],
                          lp_stat      = as.integer(sol[["stat"]]),
                          obj_coef     = ocf,
                          obj_func     = printObjFunc(model),
                          fldind       = fldind(mod),
                          fluxdist     = fluxDistribution(fluxes = sol[["fluxes"]],
                                                    nrow = length(sol[["fluxes"]]),
                                                    ncol = 1L),
                          alg_par      = alg_par(mod))
 
}
else{
        optsol <- list(ok = sol$ok,
	                   obj = sol$obj,
	                   stat =sol$stat,
	                   fluxes = sol$fluxes,
	                   wtflx=wtflux,
	                   lb=lb,
	                   ub=ub,
	                   ocf=ocf
	                  )	    

}
##########################################


    return(optsol)        
    }
