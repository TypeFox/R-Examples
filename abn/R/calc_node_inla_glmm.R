###############################################################################
## wrapper for INLA to fit glmm 
## calc_node_glmm.R
##
###############################################################################

calc.node.inla.glmm<-function(child.loc,dag.m.loc,data.df.loc,data.dists.loc,ntrials.loc,exposure.loc,compute.fixed.loc,
                         mean.intercept.loc,prec.intercept.loc,mean.loc,prec.loc,loggam.shape.loc,loggam.inv.scale.loc,verbose.loc){

#print(data.df.loc);
#print(group.var);
group.var<-names(data.df.loc)[length(names(data.df.loc))];## group variable is always the last column

 ## 1. get the formula part of the call - create a string of this
 if(length(which(dag.m.loc[child.loc,-child.loc]==1))==0){## independent node
                        str.eqn.str<-paste(colnames(dag.m.loc)[child.loc],"~1+");
 } else { ## have some covariate
    
    if(dim(dag.m.loc)[1]==2){## special case - 2x2 DAG and so names are not retained when -child.loc
    str.eqn.str<-paste(colnames(dag.m.loc)[child.loc],"~",colnames(dag.m.loc)[-child.loc],"+",sep="");
    } else {str.eqn.str<-paste(colnames(dag.m.loc)[child.loc],"~",paste(names(which(dag.m.loc[child.loc,-child.loc]==1)),collapse="+",sep=""),"+",sep="");}
  }
 #cat(str.eqn.str,"\n");
 #stop("");
 ## new part - add in the function for latent variables
 str.eqn.str<-paste(str.eqn.str,"f(",group.var,",model=\"iid\",hyper=list(theta=list(prior=\"loggamma\",param=c(",
                                           loggam.shape.loc,",",loggam.inv.scale.loc,")))),\n",sep="");
 #print(str.eqn.str);
 #stop("");
 ## 2. data set
 str.data<-"data=data.df.loc,";
 
 ## 3. family
 str.family<-paste("family=\"",data.dists.loc[[child.loc]],"\",",sep="");
 
 ## 4. additional parameter for number of trials (binomial) or exposure (poisson)
 str.extra<-"";
 if(data.dists.loc[[child.loc]]=="binomial"){ str.extra<-paste("Ntrials=ntrials.loc,",sep="");}
 if(data.dists.loc[[child.loc]]=="poisson"){ str.extra<-paste("E=exposure.loc,",sep="");}
 if(data.dists.loc[[child.loc]]=="gaussian"){str.extra<-paste("control.family=list(hyper = list(prec = list(prior=\"loggamma\",param=c(",
                                           loggam.shape.loc,",",loggam.inv.scale.loc,")))),\n",sep="");}

 ## 5. get the full command
 res<-NULL;
 start.str<-"res<-inla(";
 end.str  <-paste("\ncontrol.fixed=list(mean.intercept=",mean.intercept.loc,",\n",
                                   "prec.intercept=",prec.intercept.loc,",\n",
                                   "mean=",          mean.loc,",\n",
                                   "prec=",          prec.loc,",\n", 
                                   "compute=",       compute.fixed.loc,"))\n",sep="");

 r<-NULL; 
 full.command<-paste("r<-try(",start.str,str.eqn.str,str.data,str.family,str.extra,end.str,",silent=TRUE)",sep="");

 ## 6. some debugging - if requested
 if(verbose.loc){cat("commands which are parsed and sent to inla().\n");
             print(full.command);}

 ## 7. now run the actual command - parse and eval - is parsed in current scope and so data.df exists here
 eval(parse(text=full.command));

 if(length(r)==1){ 
       ### INLA failed
       # cat("INLA failed\n");
       return(FALSE);
    } else {

 ## 8. return the results
 if(!compute.fixed.loc){## only want marginal likelihood
          return(res$mlik[1]);## n.b. [1] is so we choose the integrated rather than Gaussian version - debateable which to choose      

                   } else { ## alternatively get *all* the output from inla() - copious
                      return(res);}
}

} ## end of function
