GxMestimator <- 
function(listdata,fMP,zeroset,closedform,starting,K,coreNumber,usedfullparam,usedmanualinitial,priority,gradientlevel) {

  initiallevel = (names(usedmanualinitial) %in% usedfullparam);

  if(sum(initiallevel) == length(usedfullparam)) {
    if(closedform) { GxMmle = optimizer_closedform(listdata=listdata,fMP=fMP,zeroset=zeroset,starting=starting,coreNumber=coreNumber);
    } else {
      if(priority==1) {
        Output3 = optimizer_numericalintegration(listdata=listdata,fMP=fMP,zeroset=zeroset,starting=starting,K=3,coreNumber=coreNumber);
        starting$initial = Output3$par;
        GxMmle = optimizer_numericalintegration(listdata=listdata,fMP=fMP,zeroset=zeroset,starting=starting,K=K,coreNumber=coreNumber);
      } else { GxMmle = optimizer_numericalintegration(listdata=listdata,fMP=fMP,zeroset=zeroset,starting=starting,K=K,coreNumber=coreNumber);
      }  
    } 
  } else {
    if(closedform) { GxMmle = optimizer_closedform(listdata=listdata,fMP=fMP,zeroset=zeroset,starting=starting,coreNumber=coreNumber);
    } else {       
      if(priority==1) {
        Output3 = optimizer_numericalintegration(listdata=listdata,fMP=fMP,zeroset=zeroset,starting=starting,K=3,coreNumber=coreNumber);
        starting$initial = Output3$par;
        GxMmle = optimizer_numericalintegration(listdata=listdata,fMP=fMP,zeroset=zeroset,starting=starting,K=K,coreNumber=coreNumber);
      } else { 
        Output3 = optimizer_numericalintegration(listdata=listdata,fMP=fMP,zeroset=zeroset,starting=starting,K=3,coreNumber=coreNumber);
        starting$initial = Output3$par;
        if(!is.null(usedmanualinitial)) {
          usedmanualnames = names(usedmanualinitial);
          for (i in 1:length(usedmanualnames)) {
            localind = which(names(Output3$par) %in% usedmanualnames[i]);
            starting$initial[localind] = usedmanualinitial[[i]];
          }
        }
        GxMmle = optimizer_numericalintegration(listdata=listdata,fMP=fMP,zeroset=zeroset,starting=starting,K=K,coreNumber=coreNumber);
      }  
    }     
  }

  if(any(is.na(GxMmle$gradient))) {
    cat('The optimization gets stuck on the boundary for parameter', names(GxMmle$par)[min(which(is.na(GxMmle$gradient)))], '.\n');
    cat('Changing the corresponding manualinitial value may help.\n');
  } else {
    if(max(abs(na.omit(GxMmle$gradient))) > gradientlevel) {
      cat('The greatest gradient for the solution is', max(abs(na.omit(GxMmle$gradient))), '. It is too large to claim convergence.\n'); 
      cat('The optimization might have not reached a minimum value.\n');
      cat('Run function \'checkGxM\' to see whether this is caused by a local singular issue.\n');
      cat('Changing the manualinitial value(s) or setting certain parameter(s) zero may help.\n');
    }
  } 

  return(GxMmle);
}
