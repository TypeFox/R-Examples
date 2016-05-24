summary_dagRdag <-
function(dag){

 cat(c("DAG created by dagR version ",
     ifelse(is.null(dag$version),"1.0.1 (or manually)", dag$version),"."),
     sep="");

 cat(c("\n - from X =", dag$names[1], "to Y =", dag$names[length(dag$names)]));

 cat("\n - acyclic:", is.acyclic(dag, maxSecs=5)[[1]]);


 cat(c("\n - ", length(dag$names),   " nodes: X, Y, ",
                length(dag$names)-2, " cov(s) (",
                sum(sapply(X=c(2:(length(dag$names)-1)),
                    FUN<-function(x) is.unknown(x,dag))), 
                                     " unknown/unmeasured)"), sep="");

 if(length(dag$names)>2)
  for(i in 2:(length(dag$names)-1))
  cat(c("\n   - node ", i, " = ", dag$names[i]), sep="");

 if(is.null(dag$symbols)==FALSE)
 { cat(c("\n - alternative node symbols: ", dag$symbols));
 }
 

 cat(c("\n - ", length(dag$arc.type), " arc(s) (",
                length(dag$arc.type[dag$arc.type==1]), " undirected)"),
       sep="");

 cat("\n - currently", ifelse(is.null(dag$adj), "unadjusted",
                       paste(c("adjusted for node(s)", dag$adj,
     "\n   i.e.",      paste(dag$names[dag$adj], collapse=', ')),
     collapse=' ')));

 if(is.null(dag$pathsN)==TRUE)
 { cat("\n - potentially biasing paths not sought");
 } else
 { if(dag$pathsN==0)
   { cat("\n - 0 potentially biasing paths found");
   } else
   {
     if(is.null(dag$path.status)==TRUE)
     { cat("\n -", dag$pathsN,
                   "potentially biasing path(s) found (not evaluated)");
     } else
     { cat("\n - ", dag$pathsN, 
                    " potentially biasing path(s) found (",
                    length(dag$path.status[dag$path.status=="open"]),
                    " open)", sep='');
     }
   }
 }

 if(is.null(dag$searchType)==TRUE)
 { cat("\n - possible adjustment sets not examined");
 } else
 { cat("\n - sufficient adjustment set(s) found:");
   if(nrow(dag$searchRes)==0) 
   { cat("\n   none to check");
   } else
   { dag$searchRes <- dag$searchRes[!is.na(dag$searchRes$openPaths), ];
     suff<-msas(dag$searchRes);
     if(length(suff[suff!=-1])==0)
     { cat("\n   none");
     } else
     { for(i in 1:nrow(dag$searchRes))
       { if(dag$searchRes$openPaths[i]==0)
         { adjset<-na.omit(unlist(dag$searchRes[i,
                                  (1:(ncol(dag$searchRes)-2))]));
           cat(c("\n   ", ifelse(is.na(adjset[1])==FALSE,
                                 paste(adjset, collapse=' '),
                                 "empty"),
                 ifelse(suff[i]==0, " *minimal*", "")));
         }
       }
     }
   }
   cat(c("\n   (search type: ", dag$searchType, ")"), sep='');
 }
 cat("\n");
}
