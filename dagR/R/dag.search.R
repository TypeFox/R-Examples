dag.search <-
function (dag, type="brute", allow.unknown = FALSE, trace = FALSE, stop = 0) 
{
  dag$searchType<-NULL;
  dag$searchRes<-NULL;

  if(type=="brute")
   searchRes<-brute.search(dag, allow.unknown=allow.unknown,
                                trace=trace, stop=stop);
 
  rv <- dag;
  rv$searchType<-paste(c(type, ifelse(allow.unknown==FALSE,
                             "unknowns not allowed",
                             "unknowns allowed"),
                             ifelse(stop==0, "no early stopping",
                                             "early stopping possible")),
                       collapse=", ");
  rv$searchRes<-searchRes;
  return(rv);
}

