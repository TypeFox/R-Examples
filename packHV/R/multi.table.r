#' Multi cross table
#'
#' Builds a big cross table between several covariates
#'
#' @param data the data frame in which we can find \code{vars}
#' @param vars vector of character string of covariates
#' @return A matrix containing all the contingency tables between the covariates
#' @author Hugo Varet
#' @seealso \code{\link{plot_multi.table}}
#' @examples
#' multi.table(cgd,c("treat","sex","inherit"))

# fonction qui croise toutes les variables qualitatives
multi.table=function(data,vars){
  if (any(!I(vars %in% names(data)))){
    stop(paste("Variable(s) ",paste(vars[!I(vars %in% names(data))],collapse=", ")," not in ",deparse(substitute(data)),"\n",sep=""))
  }
  vars=unique(vars)

  names=NULL
  levels=NULL
  vars.pb=NULL
  for (var in vars){
    if (!(is.factor(data[,var]) | is.character(data[,var]))){
      vars.pb=c(vars.pb,var)
    }
  }
  if (length(vars.pb)>0){
    stop(paste(paste(vars.pb,collapse=", ")," not factor neither character",sep=""))
  }
  
  data=convert_factor(data,vars)
  for (var in vars){names=c(names,var,rep("",nlevels(data[,var])))}
  for (var in vars){levels=c(levels,"",levels(data[,var]))}
  names=names[-length(names)]
  res=matrix("",nrow=length(names),ncol=length(names),dimnames=list(names,names))
  for (var1 in vars){
    for (var2 in vars){
      deb1=which(rownames(res)==var1)
      deb2=which(colnames(res)==var2)
      res[deb1:(deb1+nlevels(data[,var1])-1),deb2:(deb2+nlevels(data[,var2])-1)]=table(data[,var1],data[,var2])
    }
  }
  res2=matrix("",nrow(res)+1,ncol(res)+1)
  colnames(res2)=rownames(res2)=c("",colnames(res))
  res2[-1,-1]=res
  res2[1,]=res2[,1]=levels
  
  res3=matrix("",nrow(res2)+1,ncol(res2)+1)
  colnames(res3)=rownames(res3)=c("",colnames(res2))
  res3[-c(1:2),-c(1:2)]=res2[-1,-1]
  res3[1,]=res3[,1]=c("",res2[1,])
  
  if (length(vars)==2){
    res3=res3[1:(nlevels(data[,vars[1]])+2), c(1:2,c(1:nlevels(data[,vars[2]]))+nlevels(data[,vars[1]])+3)]
  }

  return(noquote(res3))
}

#N=100
#my.data=data.frame(a=factor(sample(1:3,N,TRUE)),b=factor(sample(1:5,N,TRUE)),c=factor(sample(1:4,N,TRUE)))
#multi.table(my.data,c("a","b","c"))
