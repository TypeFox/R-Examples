main <-
function(X, Y, initiation, GLOBvar, HYPERvar, OUTvar){
  
  ### assignement of global variables used here ###
  niter = GLOBvar$niter
  burn_in= OUTvar$burn_in
  PSRFactor=OUTvar$PSRFactor
  PSRFactor_nbIterEval=OUTvar$PSRFactor_nbIterEval
  smax = GLOBvar$smax
  q = GLOBvar$q
  n = GLOBvar$n
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  cD = HYPERvar$cD
  alphaD = HYPERvar$alphaD
  betaD = HYPERvar$betaD
  ### end assignement ###
  
  ### assignement of initiation variables used here ###
  # initial state
  E = initiation$initState$E
  Sall = initiation$initState$Sall
  Ball = initiation$initState$Ball
  Sig2all = initiation$initState$Sig2all
  s = initiation$initState$s
  # counters
  EdgesMovesCount = initiation$counters$EdgesMovesCount
  EdgesMovesAcceptation = initiation$counters$EdgesMovesAcceptation
  CPMovesCount = initiation$counters$CPMovesCount
  CPMovesAcceptation = initiation$counters$CPMovesAcceptation
  # storage matrices
  Estock = initiation$listStock$Estock
  Sstock = initiation$listStock$Sstock 
  Bstock = initiation$listStock$Bstock
  Sig2stock=initiation$listStock$Sig2stock    
  PSRF_CP_stock=NULL
  PSRF_edges_stock=NULL
  rcvgce=NULL
  
  # do niter interations
  cvgce=FALSE
  r=2
  if (PSRFactor){
    print(paste("Looking for convergence (up to",niter, "iterations) :"))
  }else{
    print(paste("Running",niter, "iterations:"))
  }
 
  
  while (!cvgce && r <= niter){
    ## Print percentage of accomplished iteration 
    if (!cvgce &&  (r%%(round(niter/10))) == 0) {
      cat(round(10 * r/(round(niter/10))), "% ",sep="")
    }
  

    D = rgamma(1, shape=s+alphaD, rate = 1+betaD)
    rho = computeRho4(s, 0, smax, cD, D)

    ## test if no birth is possible (due to the current vector E and value of segminlength)
    toremove = E
    if(GLOBvar$segMinLength>1) for(i in 1:(GLOBvar$segMinLength-1)) toremove = c(toremove, E-i, E+i)
    ## possible CPs are those not in 'toremove'
    possibleCP = setdiff((1+GLOBvar$dyn):E[length(E)], toremove)
    if(length(possibleCP)==0) {rho = computeRho4(smax, 0, smax, cD, D)}
    
    ## Sample u to choose one of the 4 moves : CP birth, CP death, CP shift, Update phases.
    u1 = runif(1, 0, 1)

    ## Run 1 out of the 4 moves (depending on the value of u)
    if (u1 < rho[1]){
      #print("birth")
      ## CP birth move: return the new model if the move is accepted, the previous model otherwise.
      out = cp.birth(E, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar)
      
    } else {
      if(u1 < rho[2]){
        #print("death")
        ## CP death move: return the new model if the move is accepted, the previous model otherwise.
        out = cp.death(E, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar)

      } else {
        if(u1 < rho[3]){
          #print("shift")
          ## CP shift move: return the new model if the move is accepted, the previous model otherwise.
          out = cp.shift(E, Sall, Ball, Sig2all, X, Y, GLOBvar, HYPERvar)
          
        } else {
          ## Update phases: return the new model if the move is accepted, the previous model otherwise.
          out = segment.update (E, Sall, Ball, Sig2all, X, Y, GLOBvar, HYPERvar, EdgesMovesCount, EdgesMovesAcceptation)
          EdgesMovesCount=out$EdgesMove
          EdgesMovesAcceptation=out$EdgesAccept
        }
      }
    }

    ## Apply changes to the current model
    E = out$E
    Ball = out$Ball
    Sall = (abs(out$Ball)>0)*1
    Sig2all = out$Sig2all
    s = length(E)-2
	
	
    ## Update moves counts
    CPMovesCount[out$move] = CPMovesCount[out$move]+1
    CPMovesAcceptation[out$move] = CPMovesAcceptation[out$move]+out$accept

    ## Stock model and parameters
    Estock[r,1:(s+2)] = E

    ##updated by Sophie 07/07/2011
    if(q==1){
      Sstock[r,1:(s+1)] = Sall[,1:q] 
    }else{
      Sstock[r,1:(s+1)] = Sall[,1:q]  %*% 2^(0:(q-1))
    }

    Bstock[r,1:((s+1)*(q+1))] = c(t(Ball))
    Sig2stock[r,1:length(Sig2all)] = Sig2all

	
    ## added by Sophie 01/10/09

    ##   (dim(unique(Estock[(r/2+1):r,]))[1])>1 <=> at least 2 differing CP configurations were visited
    
    if(PSRFactor & (!cvgce) & r>burn_in &((r-burn_in)%%OUTvar$PSRFactor_freq==0)){
      if( (dim(unique(Estock[(r/2+1):r,]))[1])>1){ 
        
        ##CPstock matrix of size r/2 x n+1 (having 1 in row i and column j of there is a changepoint in position j at iteration i)
        CPstock=matrix(0, r/2, n+1)
        for( i in 1:(r/2)){
          CPstock[i,Estock[r/2+i,]]=1
        }
        CPstock=CPstock[,(Estock[1,1]+1):n]

        ## PSRF computation
        
        PSRF_CP=PSRF(CPstock,OUTvar$PSRFactor_nbSeq)
        
        PSRF_CP_stock=rbind(PSRF_CP_stock,c(r,PSRF_CP))
        
        if(PSRF_CP<OUTvar$PSRF_thres){
          ##Edges
          diff=Estock[(r/2+1):r,-1]-Estock[(r/2+1):r,1:(dim(Estock)[2]-1)]
          diff=as.matrix(diff)
          diff[which(diff==-(n+1))]=0
          
          tmp=NULL
          for (i in 1:dim(diff)[2]){
            tmp=cbind(tmp,matrix(diff[,i],dim(diff)[1],q+1))
          }
          
          tmp2=(Bstock[(r/2+1):r,]!=0)*tmp
          
          ## sum aver each parent /n
          postProbAVG=matrix(0,dim(diff)[1],q+1)
          for (i in 1:dim(diff)[2]){
            postProbAVG=postProbAVG+tmp2[,((q+1)*(i-1)+1):((q+1)*i)]
          }		
          postProbAVG=postProbAVG[,1:q]/(max(Estock[1,])-(Estock[1,1]))
          
          PSRF_edges=PSRF(postProbAVG,OUTvar$PSRFactor_nbSeq)
          PSRF_edges_stock=rbind(PSRF_edges_stock,c(r,PSRF_edges))

          
          if(PSRF_edges < OUTvar$PSRF_thres){
            cat("\n")
            print(paste("**** Convergence obtained (after ",r, " iterations, PSRF threshold = ", OUTvar$PSRF_thres,") ***", sep=""))
            cvgce=TRUE
            niter=r+PSRFactor_nbIterEval-1
            rcvgce=r
            print(paste("RJMCMC iterations for Postdistribution evaluation (",PSRFactor_nbIterEval, " iterations).", sep=""))
          }
          
        }# end if (PSRF_CP<OUTvar$PSRF_thres)
      }#  end if( (dim(unique(Estock[(r/2+1):r,]))[1])>1){
    }#  end if(PSRFactor & (!cvgce) & ...
    
    r=r+1
    
  }# end while r<=
  cat("\n")

  
  if(PSRFactor){
    if(cvgce){
      samples = list(CP=Estock[rcvgce:(r-1),], Edges=Sstock[rcvgce:(r-1),], coeff=Bstock[rcvgce:(r-1),], variance=Sig2stock[rcvgce:(r-1),])
    }else{
      print("Warning convergence was not reached")
      print(paste("CP PSRF=", PSRF_CP))
      if(PSRF_CP<OUTvar$PSRF_thres){
        print(paste("Edge PSRF=", PSRF_edges))
      }
      samples = list(CP=Estock, Edges=Sstock, coeff=Bstock, variance=Sig2stock)
    }
  }else{
    samples = list(CP=Estock, Edges=Sstock, coeff=Bstock, variance=Sig2stock)
  }

  CPMovesAcceptation=round(CPMovesAcceptation/CPMovesCount, 4)*100
  EdgesMovesAcceptation=round(EdgesMovesAcceptation/EdgesMovesCount,4)*100
  
  counters = list(CPMovesCount=CPMovesCount, CPMovesAcceptationPrct=CPMovesAcceptation, EdgesMovesCount=EdgesMovesCount, EdgesMovesAcceptationPrct=EdgesMovesAcceptation, iterations=r-1, rcvgce=rcvgce)

  ## return RJMCMC samples
  return(list(samples=samples, counters=counters))
}
