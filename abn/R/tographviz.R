## tographviz.R --- 
## Author          : Fraser Lewis
## Edited by Marta Pittavino
## Last Modified on: 01/07/2014
################################################################################

tographviz<-function(dag.m,data.df,data.dists,group.var=NULL,outfile,directed=TRUE){

    
  if(!is.null(group.var)){## have group variable so just need to rebuild data.df without this
  if(!(is.character(group.var) && (length(group.var)==1))){stop("name of group variable is not a character?!");}
  if(!length(which(group.var%in%names(data.df)==TRUE))){stop("name of group variable does not match any of those in data.df");}
  group.var.vals<-data.df[,group.var];## get group id data
  data.df<-data.df[,-which(names(data.df)==group.var)];## drop the group variable from original data.frame and overwrite
   
  }
    #some checks
    #check.valid.dag(dag.m=dag.m,data.df=data.df,is.ban.matrix=FALSE,group.var=NULL);
    

    # check dag is in a matrix
    if(!is.matrix(dag.m)){stop("The DAG definition dag.m must be in a matrix");}

    # check data for missing names
    if(is.null(colnames(dag.m)) || is.null(rownames(dag.m))){stop("dag.m must have both row and column names set");}

    # check dimension
    if(dim(dag.m)[1]!=dim(data.df)[2] || dim(dag.m)[2]!=dim(data.df)[2] ){stop("dag.m as dimension inconsistent with data.df - if using grouped data you must supply group.var argument");}

    # check binary
    for(i in 1:dim(dag.m)[1]){for(j in 1:dim(dag.m)[2]){if(dag.m[i,j]!=0 && dag.m[i,j]!=1){stop("dag.m must comprise only 1's or 0's");}}}

    
    ## create header part
    cat(ifelse(directed, "digraph dag {","graph dag {"),"\n\n",file=outfile,append=FALSE);
    # Old version: if(directed){ cat("digraph dag {","\n\n",file=outfile,append=FALSE); }
    # Old version: else{ cat("graph dag {","\n\n",file=outfile,append=FALSE);} 
            for(i in 1:length(colnames(dag.m))){
                       if(data.dists[[i]]=="binomial"){cat(paste("\"",colnames(dag.m)[i],"\"[shape=square];\n",sep=""),file=outfile,append=TRUE);}
                       if(data.dists[[i]]=="gaussian"){cat(paste("\"",colnames(dag.m)[i],"\"[shape=oval];\n",sep=""),file=outfile,append=TRUE);}
                       if(data.dists[[i]]=="poisson"){cat(paste("\"",colnames(dag.m)[i],"\"[shape=diamond];\n",sep=""),file=outfile,append=TRUE);}
            }
            cat("\n\n\n",file=outfile,append=TRUE)
    
 for(i in colnames(dag.m)){##for each variable
             children<-which(dag.m[,i]==1);##get row with children
             if(length(children)>=1){##if have at least one child
             child.nom<-rownames(dag.m)[children];
            # if(directed) {for(j in child.nom){cat("\"",i,"\"","->","\"",j,"\";","\n",sep="",file=outfile,append=TRUE);}}
             #else { for(j in child.nom){cat("\"",i,"\"","--","\"",j,"\";","\n",sep="",file=outfile,append=TRUE);}
           {
             for(j in child.nom){
               cat("\"",i,"\"",ifelse(directed, "->", "--"),"\"",j,"\";","\n",sep="",file=outfile,append=TRUE);
             }
           }
                 }
                }
     ## footer part
            cat("\n}\n",file=outfile,append=TRUE);
}
