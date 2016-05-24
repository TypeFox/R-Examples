msas <-
function(adjSets){
  adjSets<-adjSets[!is.na(adjSets$openPaths),] #remove unevaluated ones;
  sets<-adjSets[,1:(ncol(adjSets)-2)];
  sufficient<-adjSets$openPaths==0;
  minimal<-numeric();
  for(i in 1:length(sufficient))
  { if(sufficient[i]==FALSE) minimal<-c(minimal, -1)
    else
    { smallerSets<-sum(apply(X=sets[sufficient,], MARGIN=1,
                             function(x) viv(x, unlist(sets[i,]))))-1;
      minimal<-c(minimal, smallerSets);
    }    
  }
  rv<-minimal;
  return(rv);
}