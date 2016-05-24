################################################
# Function: dynamicFBA
#
# Performs a dynamic flux balance analysis
# 
# The function dynamicFBA() is inspired by the function
# dynamicFBA() contained in the COBRA Toolbox.
# The algorithm is the same.

dynamicFBA <- function (model,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,exclUptakeRxns,
 		  #lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                  #solver = SYBIL_SETTINGS("SOLVER"),
                  #method = SYBIL_SETTINGS("METHOD"), 
                  #solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM"),
                  fld = FALSE,verboseMode = 2, ...){
#PARAMETERS:
#===========
# model                 Sybil model structure (class modelorg)
# substrateRxns         List of exchange reaction names for substrates
#                       initially in the media that may change (e.g. not
#                       h2o or co2)
# initConcentrations    Initial concentrations of substrates (in the same
#                       structure as substrateRxns)
# initBiomass           Initial biomass (must be non zero)
# timeStep              Time step size
# nSteps                Maximum number of time steps
#
#OPTIONAL PARAMETERS
#===================
# exclUptakeRxns        List of uptake reactions whose substrate
#                       concentrations do not change (Default =
#                       {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'})
# 
#RETURN VALUES:
#=============
# concentrationMatrix   Matrix of extracellular metabolite concentrations
# excRxnNames           Names of exchange reactions for the EC metabolites
# timeVec               Vector of time points
# biomassVec            Vector of biomass values
#
# If no initial concentration is given for a substrate that has an open
# uptake in the model (i.e. model.lb < 0) the concentration is assumed to
# be high enough to not be limiting. If the uptake rate for a nutrient is
# calculated to exceed the maximum uptake rate for that nutrient specified
# in the model and the max uptake rate specified is > 0, the maximum uptake 
# rate specified in the model is used instead of the calculated uptake
# rate.


##--------------------------------------------------------------------------##
 # check prerequisites 
    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
    }
##--------------------------------------------------------------------------##

# Uptake reactions whose substrate concentrations do not change
if (missing(exclUptakeRxns)){
    exclUptakeRxns = c('EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)');
    if (verboseMode > 2){
       print('Default extra cellular uptake reactions will be used: ')
       print(exclUptakeRxns);
    }
}

# Find exchange reactions
excReact = findExchReact(model);
excReactInd=(react_id(model) %in% react_id(excReact));#excReact$exchange
#represent extra cellular reaction with boolean vector.
exclUptakeRxnsInd=is.element(react_id(model) ,exclUptakeRxns);
#Exclude reactions with concentrations that will not be changed 
excReactInd = excReactInd & !exclUptakeRxnsInd;   #excInd & ~ismember(model.rxns,exclUptakeRxns);
#get reaction names
excRxnNames =react_id(model)[excReactInd];                #excRxnNames = model.rxns(excInd);

substrateRxnsInd=(react_id(model) %in% substrateRxns)
# Figure out if substrate reactions are correct: all substrate reactions should be exchange reactions.
missingSub = substrateRxnsInd & !excReactInd;
if (sum(missingSub)!=0){
    print(sum(missingSub));
    print(react_id(model)[missingSub]);
    print('Invalid substrate uptake reaction!');
}
## 	***********************************************************     ##
# Initialize concentrations
#substrateMatchInd = intersect(excRxnNames,substrateRxns);
concentrations=rep(0,length(react_id(model)))#table(excRxnNames);##vector(length=length(excRxnNames),mode="numeric");
#concentrations[1:length(concentrations)]=0;
concentrations[substrateRxnsInd] = initConcentrations;

# Deal with reactions for which there are no initial concentrations
originalBound = -lowbnd(model);# take all to be able to directly update
noInitConcentration = (concentrations==0)&(lowbnd(model)<0)#(concentrations == 0 & originalBound > 0);
concentrations[noInitConcentration] = 1000;

biomass = initBiomass;

# Initialize bounds
uptakeBound =  concentrations/(biomass*timeStep);

# Make sure bounds are not higher than what are specified in the model
aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
uptakeBound[aboveOriginal] = originalBound[aboveOriginal];
lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];

concentrationMatrix = concentrations[excReactInd];
biomassVec = biomass;
timeVec = 0;
##------------------------------------- Prepare Problem object --------------------##
# get OptObj instance
#lpmod <- prepProbObj(model,
#                         nCols      = react_num(model),
#                         nRows      = met_num(model),
#             #            alg        = "FBA",
#                         solver     = solver,
#                         method     = method,
#                         lpdir      = lpdir
#                         #solverParm = solverParm
#             )

lpmod <- sybil::sysBiolAlg(model, algorithm = "fba", ...)

##-----------------------------------------------------------------------------##
 if (verboseMode > 2) print('Step number    Biomass\n');
# Inititialize progress bar ...');
#if (verboseMode == 2)  progr <- .progressBar();

for (stepNo in 1:nSteps){
    
#    if (verboseMode == 2)  progr <- .progressBar(stepNo, nSteps, progr);
      
    # Run FBA
      #sol = simpleFBA(lpmod,fld=TRUE,solver= solver,method= method,);
      sol = sybil::optimizeProb(lpmod);
      mu =  sol$obj;  ##objvalue sol.f
    if ( length(checkSolStat(sol$stat,solver(problem(lpmod))))!=0 ){## checkSolStat
        print('No feasible solution - nutrients exhausted\n');
        break;
    }
   uptakeFlux = sol$fluxes[excReactInd];
   biomass = biomass*exp(mu*timeStep);
    #biomass = biomass*(1+mu*timeStep);
    biomassVec = c(biomassVec,biomass);
    
    # Update concentrations
    concentrations[excReactInd]= concentrations[excReactInd] - uptakeFlux/mu*biomass*(1-exp(mu*timeStep));
    #concentrations = concentrations + uptakeFlux*biomass*timeStep;
    concentrations[concentrations <= 0] = 0;
    concentrationMatrix = c(concentrationMatrix,concentrations[excReactInd]);
      
    # Update bounds for uptake reactions
    uptakeBound[excReactInd] =  concentrations[excReactInd]/(biomass*timeStep);
    # This is to avoid any numerical issues
    uptakeBound[uptakeBound > 1000] = 1000;
    # Figure out if the computed bounds were above the original bounds
    aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
    # Revert to original bounds if the rate was too high
    uptakeBound[aboveOriginal] = originalBound[aboveOriginal];# uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
    uptakeBound=ifelse(abs(uptakeBound) < 1e-9,0,uptakeBound);
## Change lower bounds according to the result of last step
    #lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];  
    uppb_tmp <- getColsUppBnds(problem(lpmod), which(excReactInd));
    changeColsBnds(problem(lpmod),which(excReactInd),lb=-uptakeBound[excReactInd],ub=uppb_tmp);
    
    if (verboseMode > 2) print(paste(stepNo,sep="    ",biomass));
    #waitbar(stepNo/nSteps,h);
    timeVec = c(timeVec,stepNo*timeStep);
}# end loop

## Preparing OUTPUT
#concentrationMatrix,excRxnNames,timeVec,biomassVec
return (optsol_dynamicFBA(solver = solver(problem(lpmod)),
		      method = method(problem(lpmod)),
		      nprob  = stepNo,
		      ncols  = react_num(model),
		      nrows  = met_num(model),
		      fld    = fld,
		      concmat=concentrationMatrix,
		      exRxn=excRxnNames,
                      tmVec=timeVec,  
                      bmVec=biomassVec
                      )
  )
}
