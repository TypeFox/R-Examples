## hzar.deltas<-function(x) {
##   ## Assumes x is a vector or a simple list
##   xS<-sort(unique(x));
##   nS<-length(xS);
##   x0<-xS[1];
##   xN<-xS[nS];
##   dScale<-(2*(xN-x0));
##   dLen<-c(xS[2:nS],xN)-c(x0,xS)[1:nS];
##   return(data.frame(x=xS,dx=as.numeric(dLen/dScale)));
## }

hzar.qScores <- function(x, wt,probs=c(0,0.25,0.5,0.75,1.0)){
  if(is.data.frame(x)){
    #print("A");
    res<-list(q=probs);
    temp.func <- function(index) hzar.qScores(as.numeric(x[,index]),wt=wt,probs=probs);
    if(is.character(names(x))){
      res2<-lapply(names(x), temp.func)
      names(res2)<-names(x);
      return(do.call(data.frame,c(res,res2)));
    }
                                        #print("B");
    return(do.call(data.frame,c(res,lapply(1:(dim(x)[[2]]),
                                           temp.func))));
  }
   #print("C");
  xS<-sort(unique(x));
  nS<-length(xS);
  x0<-xS[1];
  xN<-xS[nS];
  dScale<-(xN-x0);
  dLen<-xS[2:nS]-xS[1:(nS-1)];
  if(is.data.frame(wt))
    wt<-wt[[1]];
  wt<-wt-max(wt);
   print("D");
  
  scores<-c(0,hzar.qScores.getScores(xSeries=xS,
                                 dSeries=dLen/dScale,
                                 raw.x=x,
                                 raw.wt=wt));
  
  cS<-cumsum(scores/sum(scores));
   print("E");
junk.func<-function(prob,cDSeries,xSeries){
    #               print("F");
                  if(length(prob)>1||is.list(prob))
                    return(sapply(prob,junk.func,cDSeries,xSeries));
     #              print("G");
                  if(prob<0) prob<-0;
                  if(prob>1) prob<-1;
                  
                  if(length(which(.Machine$double.eps>(cDSeries-prob)^2))>0)
                    return(mean(xSeries[which(.Machine$double.eps>
                                              (cDSeries-prob)^2)]));
# print("H");
                  i1<-rev(which(cDSeries<prob))[[1]];
                  i2<-which(cDSeries>prob)[[1]];
 #         print("I");                          
                  return(xSeries[[i1]]+
                         (prob-cDSeries[[i1]])*
                         (xSeries [[i2]]-xSeries [[i1]])/
                         (cDSeries[[i2]]-cDSeries[[i1]]));
                  
                }
  return(sapply(probs,junk.func
                ,cS,xS));
}

hzar.qScores.getScores <- function(xSeries,dSeries,raw.x,raw.wt){
  if(length(xSeries)>1e3){
    cutValue<-as.integer(length(xSeries)/2)+1;
    return(c(hzar.qScores.getScores(xSeries[1:cutValue],
                                    dSeries[1:cutValue],
                                    raw.x [raw.x<=xSeries[[cutValue]]],
                                    raw.wt[raw.x<=xSeries[[cutValue]]]),
             hzar.qScores.getScores(xSeries[cutValue:length(xSeries)],
                                    dSeries[cutValue:length(xSeries)],
                                    raw.x [raw.x>=xSeries[[cutValue]]],
                                    raw.wt[raw.x>=xSeries[[cutValue]]])));
  }
  
##  print(xSeries[[1]])
  getScore<-function(x,
                     xS=xSeries,
                     dS=dSeries,
                     rx=raw.x,
                     rwt=raw.wt){
    data=c(rwt[xS[[x]]==rx],rwt[xS[[x+1]]==rx]);
    return(mean(exp(data))*dS[[x]] );
  }


  return(sapply(1:(length(xSeries)-1),getScore));
}
## scores<-c(0,sapply(1:(nS-1),
##                  function(x,xSeries,dSeries,raw.x,raw.wt)
##                  mean(exp(raw.wt[xSeries[[x  ]]==raw.x |
##                             xSeries[[x+1]]==raw.x ])
##                       )*dSeries[[x]],
##                  xSeries=xS,dSeries=dLen/dScale,raw.x=x,raw.wt=wt));
hzar.qScores.dataGroup <- function(dataGroup,probs=c(0.025,0.5,0.975)){
  return(hzar.qScores(as.data.frame(dataGroup$data.mcmc),dataGroup$data.LL$model.LL,probs=probs));
}

hzar.qScores.obsDataGroup <- function(oDG,probs=c(0.025,0.5,0.975)){
  ## Need to do model selection, hard-coding AICc
  temp.AICc<-hzar.AICc.hzar.obsDataGroup(oDG);
  print(selected.model<-rownames(temp.AICc)[temp.AICc[[1]]==min(temp.AICc[[1]])]);
  dG.all<-oDG$data.groups[[selected.model]];
  return(hzar.qScores.dataGroup(dG.all,probs=probs));
}
## junk.getVal <- 
## junk.LL<-do.call(rbind,
##                  lapply(junk.deltas$x,
##                         function(x,data.x,data.LL)
##                         data.frame(x=x,
##                                    xL=mean(exp(data.LL[x==data.x]))),
##                         data.x=junk$center,
##                         data.LL=junk$model.LL))

##  cbind(x=junk.deltas[,"x"],sR=junk.deltas$dx*junk.LL$xL)

##  data.frame(x=junk.deltas[,"x"],sR=junk.deltas$dx*junk.LL$xL)->junk.scores

## junk.scores$score<-junk.scores$sR/sum(junk.scores$sR)

##  junk.scores$x[cumsum(junk.scores$score)>0.025]
