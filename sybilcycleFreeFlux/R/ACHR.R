ACHR <- function (model,W=2000,nPoints=5000,stepsPerPoint=10,
                  solver = SYBIL_SETTINGS("SOLVER"),
                  method = SYBIL_SETTINGS("METHOD")){
#1-warmup
##------------------------------------- Prepare Problem object --------------------##
# get OptObj instance
#lpmod <- prepProbObj(model,
#                         nCols      = react_num(model),
#                         nRows      = met_num(model),
#             #            alg        = "FBA",
#                         solver     = solver,
#                         method     = method
#                         #solverParm = solverParm
#             )
 lpmod <- sybil::sysBiolAlg(model, solver     = solver,
                         method     = method
                         )

##-----------------------------------------------------------------------------##
# first run FVA to get min, max for each rxn
nRxns=react_num(model);
if( W < 2*nRxns){
	stop("Warmup points should be more than double the number of rxns");
}

warmupPts = matrix(rep(0,nRxns*W),ncol=W);

#Generate the points (every iteration generates two points (min,max)
for( i in (1:floor(W/2))){
    # Set the objective function
            if (i <= nRxns){
                ocf = rep(0,nRxns);
                ocf[i] = 1;
            }  else{            # Create random objective function
    		    ocf = runif(nRxns)-0.5;
                }
      changeObjCoefs(problem(lpmod),c(1:nRxns),ocf);
      
    for (minMax in c("min", "max")){
            # Determine the max or min for the rxn	

	setObjDir(problem(lpmod),minMax)

	sol=optimizeProb(lpmod)#,solver= solver,method= method,lpdir=minMax taken fom lpmod
        x = sol$fluxes;
        status = sol$stat;
        if( (status == 5 && solver=="glpkAPI")||(status == 1 && solver=="cplexAPI")){
            validFlag = TRUE;
        } else{
            print ('invalid solution')
            validFlag = FALSE;
            print(status)
        }
        
        # Continue if optimal solution is found
        
        # Move points to within bounds
        x[x > uppbnd(model)] = uppbnd(model)[x > uppbnd(model)];
        x[x < lowbnd(model)] = lowbnd(model)[x < lowbnd(model)];
        
        # Store point
        if (minMax == "min"){
            warmupPts[,2*i-1] = x;
        }else{
            warmupPts[,2*i] = x;
        }
                       
    }
    if (validFlag){
        i = i+1;
    }
}
centerPoint = apply(warmupPts,1,mean);
# Move points in
 warmupPts = warmupPts*.33 + .67*centerPoint;

###########################-----*******--------#######################################################################
###########################-----*******--------#######################################################################

#2-Sample
	N=MASS::Null(t(S(model)))
# Minimum allowed distance to the closest constraint
maxMinTol = 1e-9;
# Ignore directions where u is really small
uTol = 1e-9; 
# Project out of directions that are too close to the boundary
dTol = 1e-14;

totalStepCount = 0;
prevPoint = centerPoint;
nProj=0;

totalCount = nPoints*stepsPerPoint;#nFiles*pointsPerFile*stepsPerPoint;

print('File #\tPoint #\tStep #\tTime\t#Time left\n');

    # Allocate memory for all points
    Points =  matrix(rep(0,nRxns*nPoints),ncol=nPoints);
    
    pointCount = 1;
    while (pointCount <= nPoints){
            
        # Create the random step size vector (one for each step)
        randVector = runif(stepsPerPoint);
      
        stepCount = 1;
        while (stepCount <= stepsPerPoint){
            
            # Pick a random warmup point
            a = ceiling(W*runif(1));
            xa = warmupPts[,a];
            
            # Get a direction from the center point to the warmup point
            u = (xa-centerPoint);
            u = u/norm(as.matrix(u),"F")#u/norm(u)  OR sqrt(sum(a^2));
            
            # Figure out the distances to upper and lower bounds
            distUb = (uppbnd(model) - prevPoint);
            distLb = (prevPoint - lowbnd(model));
            
            # Figure out if we are too close to a boundary
            validDir = ((distUb > dTol) & (distLb > dTol));
            
            # Figure out positive and negative directions
            posDirn = which(u[validDir] > uTol);
            negDirn = which(u[validDir] < -uTol);
            
                   # Figure out all the possible maximum and minimum step sizes
            maxStepTemp = distUb[validDir]/u[validDir];
            minStepTemp = -distLb[validDir]/u[validDir];
            maxStepVec = c(maxStepTemp[posDirn],minStepTemp[negDirn]);
            minStepVec = c(minStepTemp[posDirn],maxStepTemp[negDirn]);
            
            # Figure out the true max & min step sizes
            maxStep = min(maxStepVec);
            minStep = max(minStepVec);
            print(sprintf('step %d: %f  %f',stepCount,minStep,maxStep));
            
            # Find new direction if we're getting too close to a constraint
            if ((abs(minStep) < maxMinTol && abs(maxStep) < maxMinTol) || (minStep > maxStep)){
                print(sprintf('Warning %f %f\n',minStep,maxStep));
                next;
            }
            
            # Pick a rand out of list_of_rands and use it to get a random
            # step distance
            stepDist = randVector[stepCount]*(maxStep-minStep)+minStep;
       
       # Advance to the next point
            curPoint = matrix(prevPoint + stepDist*u,ncol=1);
       
       # Reproject the current point and go to the next step
            if ( (totalStepCount %% stepsPerPoint) == 0){
            	print(sprintf("totalStepCount %d",totalStepCount));
       #     	print(sprintf("dimS:%d %d, dimcP %d,%d",dim(S(model))[1],dim(S(model))[2],dim(curPoint)[1],dim(curPoint)[2]));    
                if (max(abs(S(model)%*%curPoint)) > 1e-9){# if not in feasible region
       #         print(sprintf("dimN:%d %d, dimcP %d,%d",dim(N)[1],dim(N)[2],dim(curPoint)[1],dim(curPoint)[2]));
                  curPoint = N%*%(t(N)%*%curPoint);
                  nProj=nProj+1;
                }
            }
       
       #fprintf('#f\t#f\n',minStep,maxStep);
                curPoint[curPoint > uppbnd(model)] = uppbnd(model)[curPoint > uppbnd(model)];
	        curPoint[curPoint < lowbnd(model)] = lowbnd(model)[curPoint < lowbnd(model)];
    
            prevPoint = curPoint;
            stepCount = stepCount + 1;
            
            # Count the total number of steps
            totalStepCount = totalStepCount + 1;
            
            #recalculate the center point
            centerPoint = ((W+totalStepCount)*centerPoint + curPoint)/(W+totalStepCount+1);
            
        } # Steps per point
        
        # Add the current point to points
        Points[,pointCount] = curPoint;
       
        pointCount = pointCount + 1;
         
    } # Points per cycle
    
    return(list(wpts=warmupPts,totalStepCount=totalStepCount,curPoint=curPoint,Points=Points,W=W,nPoints=nPoints,
    centerPoint=centerPoint,nProj=nProj))
}
