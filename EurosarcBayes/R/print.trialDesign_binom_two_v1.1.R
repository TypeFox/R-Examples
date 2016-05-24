

##############################################################
# Collect the summary data
.summary_trialDesign_binom_two=function(x){

  if(length(x@reviews)>10){
    n=paste(length(x@reviews),"analyses between",x@reviews[1],"and",x@reviews[length(x@reviews)])
  } else {
    n=paste0(x@reviews,collapse=",")
  }

  dm=dim(x@data)
  alpha=round(1-min(apply(x@data[1:2,3:dm[2]],2,sum)/100),4)

  beta= round(sum(x@data[3:4,2])/100,4)

  expp0=max(x@data[5,3:dm[2]])

  return(c(n,alpha,1-beta,expp0,x@data[5,2],x@post.futility,x@post.efficacy,x@post.toxicity, x@post.no_toxicity))
}


##############################################################
# Print a single summary of design
.print.trialDesign_binom_two=function(x, ...){

  m=data.frame(n=rep("",0),alpha=rep("",0), beta=rep("",0),Exp.P0=rep("",0),Exp.P1=rep("",0),post.futility=rep("",0),post.efficacy=rep("",0),post.toxicity=rep("",0),post.no_toxicity=rep("",0), stringsAsFactors =FALSE)

  m[1,]=.summary_trialDesign_binom_two(x)

  print(m,row.names=FALSE)
}


.show.trialDesign_binom_two=function(object){

  m=data.frame(n=rep("",0),alpha=rep("",0), beta=rep("",0),Exp.P0=rep("",0),Exp.P1=rep("",0),post.futility=rep("",0),post.efficacy=rep("",0),post.toxicity=rep("",0),post.no_toxicity=rep("",0), stringsAsFactors =FALSE)

  m[1,]=.summary_trialDesign_binom_two(object)

  print(m,row.names=FALSE)
}

setMethod(f="print",signature=c(x="trialDesign_binom_two"),definition=.print.trialDesign_binom_two)

setMethod("show","trialDesign_binom_two",definition=.show.trialDesign_binom_two)

##############################################################
# Print a table of a number of design summaries saved in a list
print.list_trialDesign_binom_two=function(x, ...){

  m=data.frame(name=rep("",0),n=rep("",0),alpha=rep("",0),beta=rep("",0),Exp.P0=rep("",0),Exp.P1=rep("",0),post.futility=rep("",0),post.efficacy=rep("",0),post.toxicity=rep("",0),post.no_toxicity=rep("",0), stringsAsFactors =FALSE)

  for(i in names(x)){
    m[i,]=c(i,.summary_trialDesign_binom_two(x[[i]]))
  }
  print(m,row.names=FALSE,...)
}


