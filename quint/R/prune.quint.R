prune.quint <-
function(tree,pp=1,...){
        object<-tree
       #pp=pruning parameter
       if(names(object$fi[4])=="Difcomponent"){
       stop("Pruning is not possible; The quint object lacks estimates of the biascorrected criterion. Grow again a large tree using the bootstrap procedure." )}
        object$fi[is.na(object$fi[,4]),4]<-0
        object$fi[is.na(object$fi[,5]),5]<-0
       maxrow<-which(object$fi[,4]==max(object$fi[,4]))[1]
       bestrow<-min( which(object$fi[,4]>= (object$fi[maxrow,4]-pp*object$fi[maxrow,6]) ) )
       #cat("Size of pruned tree is",bestrow+1, "leafs","\n")
              con<-object$control
      con$Boot<-FALSE
       con$maxl<-bestrow+1
       besttree<-quint(data=object$data,control=con)
       besttree$fi<-object$fi[1:bestrow,]
       objboot<-list(siboot=object$siboot[1:bestrow,,])
       besttree<-c(besttree,objboot)
        besttree$control$Boot<-object$control$Boot
           class(besttree)<-"quint"  
          return(besttree)}
