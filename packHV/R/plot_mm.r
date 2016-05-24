#' Spaghetti plot and plot of the mean at each time
#'
#' Spaghetti plot and plot of the mean at each time
#'
#' @param formula \code{obs~time+(group|id)} or \code{obs~time+(1|id)}
#' @param data data frame in which we can find \code{obs}, \code{time}, \code{group} and \code{id}
#' @param col.spag vector of length \code{nrow(data)} with colors (one for each individual)
#' @param col.mean vector of length \code{length(levels(group))} with colors (one for each group)
#' @param type \code{"spaghettis"}, \code{"mean"} or \code{"both"}
#' @param tick.times boolean, \code{TRUE} to display ticks at each observation time on the x-axis
#' @param xlab character sring, label of the time axis
#' @param ylab character string, label of the y axis
#' @param main character string, main title
#' @param lwd.spag numeric, width of the spaghetti lines, 1 by default
#' @param lwd.mean numeric, width of the mean lines, 4 by default
#' @param ... Other arguments to be passed in \code{\link{plot}}
#' @return None
#' @author Hugo Varet on Anais Charles-Nelson's idea
#' @examples
#' N=10
#' time=rep(1:4,N)
#' obs=1.1*time + rep(0:1,each=2*N) + rnorm(4*N)
#' my.data=data.frame(id=rep(1:N,each=4),time,obs,group=rep(1:2,each=N*2))
#' par(xaxs="i",yaxs="i")
#' plot_mm(obs~time+(group|id),my.data,col.spag=my.data$group,
#'         col.mean=c("blue","red"),type="both",main="Test plot_mm")

plot_mm=function(formula,data,col.spag=1,col.mean=1,type="spaghettis",tick.times=TRUE,
                 xlab=NULL,ylab=NULL,main="",lwd.spag=1,lwd.mean=4,...){
  formumla=deparse(formula)
  name_data=deparse(data)
  formula <- paste(formula, collapse=" ")
  formula <- gsub("[[:space:]]+", " ", formula)
  varnames <- gsub("\\||\\+|~[[:space:]]|\\(|\\)| \\?\\(.*\\)", "", formula)
  varnames <- strsplit(varnames,  "[[:space:]]*(\\+|,|~)[[:space:]]*")
  varnames <- gsub("  ", " ", varnames[[1]])      
  varnames <- strsplit(varnames,  "[[:space:]]")[[1]]
  if (!all(varnames[varnames!="1"] %in% names(data))){
    stop(paste("at least one covariate is not in ",name_data,"\n",sep=""))
  }
  
  id=data[,varnames[4]]; time=data[,varnames[2]]; obs=data[,varnames[1]];
  if (varnames[3]==1){group=rep(1,nrow(data))} else{group=data[,varnames[3]]}
  data=data.frame(id,time,obs,group,col.spag=col.spag)
  data$group=factor(data$group)
  data=data[order(data$time),]; data=data[order(data$id),];
  
  if (!I(type %in% c("spaghettis","mean","both"))){
    stop("type must be equal to \"spaghettis\", \"mean\" or \"both\"")
  }

  if (is.null(xlab)){xlab=varnames[2]}
  if (is.null(ylab)){ylab=varnames[1]}  
  
  plot(data$time,data$obs,col="white",xlab=xlab,ylab=ylab,main=main,bty="n",...)
  if (tick.times){axis(1,at=unique(data$time),labels=FALSE)}
  if (type %in% c("spaghettis","both")){                                        # spaghettis
    for (i in unique(data$id)){
      data.i=data[data$id==i,]
      lines(data.i$time,data.i$obs,col=data.i$col.spag[1],lwd=lwd.spag)
    }
  }
  if (type %in% c("mean","both")){                                              # means
    for (i in 1:length(levels(data$group))){
      data.g=data[data$group==levels(data$group)[i],c("time","obs")]
      data.g=aggregate(data.g,by=list(data.g$time),FUN=function(x){mean(x,na.rm=TRUE)})
      col.g=col.mean[i]
      lines(data.g$time,data.g$obs,col=col.g,lwd=lwd.mean)  
    }               
  }
}

#N=20
#time=rep(1:4,N)
#obs=1.1*time + rep(0:1,each=2*N) + rnorm(4*N)
#my.data=data.frame(id=rep(1:N,each=4),time,obs,group=rep(1:2,each=N*2))
#par(xaxs="i",yaxs="i")
#plot_mm(obs~time+(group|id),my.data,col.spag=my.data$group,col.mean=c("blue","red"),type="both",main="Test plot_mm")
