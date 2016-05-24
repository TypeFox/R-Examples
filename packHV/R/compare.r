#' Comparing two databases assumed to be identical
#'
#' Compares two data frames assumed to be identical, prints the differences in the console and also returns the results in a data frame
#'
#' @param d1 first data frame
#' @param d2 second data frame
#' @param id character string, primary key of the two data bases
#' @param file.export character string, name of the XLS file exported
#' @return A data frame containing the differences between the two data bases
#' @author Hugo Varet
#' @examples
#' N=100
#' data1=data.frame(id=1:N,a=rnorm(N),
#'                         b=factor(sample(LETTERS[1:5],N,TRUE)),
#'                         c=as.character(sample(LETTERS[1:5],N,TRUE)),
#'                         d=as.Date(32768:(32768+N-1),origin="1900-01-01"))
#' data1$c=as.character(data1$c)
#' data2=data1
#' data2$id[3]=4654
#' data2$a[30]=NA
#' data2$a[31]=45
#' data2$b=as.character(data2$b)
#' data2$d=as.character(data2$d)
#' data2$e=rnorm(N)
#' compare(data1,data2,"id")

# last updated: july 02, 2013

compare=function(d1,d2,id,file.export=NULL){
  # d1: first data frame
  # d2: second data frame
  # id: character string of the primary key of d1 and d2
  call=match.call()
  name_d1=deparse(substitute(d1))
  name_d2=deparse(substitute(d2))
   
  if (!I(id %in% names(d1)) | !I(id %in% names(d2))){stop(paste("The primary key ",id," is not in the two databases\n",sep=""))}  
  if (any(duplicated(d1[,id])) & any(duplicated(d2[,id]))){stop(paste("Duplicates in ",name_d1," and in ",name_d2,"\n",sep=""))}
  if (any(duplicated(d1[,id]))){stop(paste("Duplicates in ",name_d1,"\n",sep=""))}
  if (any(duplicated(d2[,id]))){stop(paste("Duplicates in ",name_d2,"\n",sep=""))}

  cat("Comparing the two databases ",name_d1," and ",name_d2,":\n\n",sep="")
  out_same=TRUE
  out_table=data.frame(matrix(nrow=0,ncol=4,dimnames=list(NULL,c(id,"Covariate",name_d1,name_d2))))
  
  if (length(setdiff(d1[,id],d2[,id]))>0){out_same=FALSE;cat("Individual(s) ",paste(setdiff(d1[,id],d2[,id]),collapse=", ")," in ",name_d1," but not in ",name_d2,"\n",sep="")}
  if (length(setdiff(d2[,id],d1[,id]))>0){out_same=FALSE;cat("Individual(s) ",paste(setdiff(d2[,id],d1[,id]),collapse=", ")," in ",name_d2," but not in ",name_d1,"\n",sep="")}
  if (length(setdiff(names(d1),names(d2)))>0){out_same=FALSE;cat("Variable(s) ",paste(setdiff(names(d1),names(d2)),collapse=", ")," in ",name_d1," but not in ",name_d2,"\n",sep="")}
  if (length(setdiff(names(d2),names(d1)))>0){out_same=FALSE;cat("Variable(s) ",paste(setdiff(names(d2),names(d1)),collapse=", ")," in ",name_d2," but not in ",name_d1,"\n",sep="")}
  mask.ind=intersect(d1[,id],d2[,id])
  mask.var=intersect(names(d1),names(d2))
  d1=d1[d1[,id] %in% mask.ind,mask.var]; d2=d2[d2[,id] %in% mask.ind,mask.var];
  d1=d1[order(d1[,id]),]; d2=d2[order(d2[,id]),];
  
  for (var in setdiff(mask.var,id)){    
    pb_type=FALSE
    types=c(class(d1[,var]),class(d2[,var]))
    if (types[1]!=types[2]){
      pb_type=TRUE
      msg.type=paste(" ",var," is ",types[1]," in ",name_d1," while ",types[2]," in ",name_d2,"\n",sep="")
    }
    
    pb_factor=FALSE
    if (is.factor(d1[,var])){
      if (is.factor(d2[,var])){
        d1[,var]=factor(d1[,var]); d2[,var]=factor(d2[,var])
        if (length(setdiff(levels(d1[,var]),levels(d2[,var])))>0 |
            length(setdiff(levels(d2[,var]),levels(d1[,var])))>0){
          pb_factor=TRUE
          msg.factor=paste(" There are not the same levels of factor in ",name_d1," and ",name_d2,"\n",sep="")
        }
      }
    }
    
    tmp=data.frame(id=d1[,id],v1=d1[,var],v2=d2[,var])
    tmp$v1=as.character(tmp$v1); tmp$v2=as.character(tmp$v2);
    tmp$v1[is.na(tmp$v1)]="NA"; tmp$v2[is.na(tmp$v2)]="NA";
    names(tmp)=c(id,paste(name_d1,var,sep="_"),paste(name_d2,var,sep="_"))
    tmp=tmp[tmp[,2]!=tmp[,3] | is.na(tmp[,2]!=tmp[,3]),]
    
    if (any(c(pb_type,pb_factor,nrow(tmp)>0))){                                 # print differences
      cat("\nVariable ",var,":\n",sep="")
      if (pb_type){cat(msg.type)}
      if (pb_factor){cat(msg.factor)}
      if (nrow(tmp)>0){
        print(tmp,row.names=rep("",nrow(tmp)),justify="left")
        tmp_add=data.frame(id=tmp[,id],Covariate=rep(var,nrow(tmp)),v1=tmp[,2],v2=tmp[,3])
        names(tmp_add)=c(id,"Covariate",name_d1,name_d2)
        out_table=rbind(out_table,tmp_add)
      }
      out_same=FALSE
    }
  }
  cat("\n")
  if (!is.null(file.export)){
    WriteXLS("out_table",file.export,"data bases comparison",Encoding="latin1",AdjWidth=TRUE)
    cat("Export created: ",file.export,"\n",sep="")  
  }
  if (out_same){
    cat("The two databases are identical\n")
  } else{
    cat("There are differences between the two databases\n")
  }
  return(invisible(out_table))
}
