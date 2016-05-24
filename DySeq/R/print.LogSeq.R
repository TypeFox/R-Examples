#'print.LogSeq
#'
#'Generates output for LogSeq Objects, see: \code{\link{LogSeq}}
#'
#'
#'@param x a LogSeq object, that should be printed
#'@param ... further arguments passed to or from other methods.
#'@export

print.LogSeq<-function(x, ...){
  
  if(dim(x[[1]])[1]==1){
    warning("Only one sequence was found!
             Aggregation not possible!
             Single case estimates are printed!")
    output<-t(x[[1]])
    rownames(output)<-c("Intercept", "Actor", "Partner", "Interac")
    colnames(output)<-"beta"
    print(output)
    return(output)
  }else{
    if(class(x)!="LogSeq") warning("x should be a LogSeq object!")
    if(length(x[[3]])==1){
      lambdas<-x[[1]]

      intercept<-stats::t.test(lambdas[,1])
      actor<-stats::t.test(lambdas[,2])
      partner<-stats::t.test(lambdas[,3])
      interac<-stats::t.test(lambdas[,4])

      mean_beta<-c(intercept$estimate, actor$estimate, partner$estimate,interac$estimate)
      exp_mean_beta<-exp(mean_beta)
      p.values<-c(intercept$p.value, actor$p.value, partner$p.value,interac$p.value)

      output_seq<-cbind(mean_beta,exp_mean_beta,p.values)
      rownames(output_seq)<-c("intercept","actor","partner","interac")
      output_seq<-as.data.frame(output_seq)
      output<-round(output_seq,3)
      print(output)
      return(output)
    }
    if(length(x[[3]])>1){
      warning("only group one and two will be compared, all other will be ignored!")

      lambdas.x<-x[[1]][x[[3]]==1,]
      lambdas.y<-x[[1]][x[[3]]==2,]

      intercept.x<-stats::t.test(lambdas.x[,1])
      actor.x<-stats::t.test(lambdas.x[,2])
      partner.x<-stats::t.test(lambdas.x[,3])
      interac.x<-stats::t.test(lambdas.x[,4])

      intercept.y<-stats::t.test(lambdas.y[,1])
      actor.y<-stats::t.test(lambdas.y[,2])
      partner.y<-stats::t.test(lambdas.y[,3])
      interac.y<-stats::t.test(lambdas.y[,4])

      intercept.yx<-stats::t.test(lambdas.y[,1],lambdas.x[,1])
      actor.yx<-stats::t.test(lambdas.y[,2],lambdas.x[,2])
      partner.yx<-stats::t.test(lambdas.y[,3],lambdas.x[,3])
      interac.yx<-stats::t.test(lambdas.y[,4],lambdas.x[,4])

      mean_beta.x<-c(intercept.x$estimate, actor.x$estimate, partner.x$estimate,interac.x$estimate)
      p.values.x<-c(intercept.x$p.value, actor.x$p.value, partner.x$p.value,interac.x$p.value)

      mean_beta.y<-c(intercept.y$estimate, actor.y$estimate, partner.y$estimate,interac.y$estimate)
      p.values.y<-c(intercept.y$p.value, actor.y$p.value, partner.y$p.value,interac.y$p.value)

      p.values.yx<-c(intercept.yx$p.value, actor.yx$p.value, partner.yx$p.value,interac.yx$p.value)

      output_seq<-cbind(mean_beta.x,p.values.x, mean_beta.y,p.values.y, p.values.yx)
      rownames(output_seq)<-c("intercept","actor","partner","interac")
      colnames(output_seq)<-c("mean beta g1", "p.value g1", "mean beta g2", "p.value g1", "p.value g1 vs g2")
      output<-round(output_seq,3)
      print(output)
    }
  }
}



