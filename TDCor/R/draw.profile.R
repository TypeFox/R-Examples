draw.profile <-
function(dataset,gene,...)

{ 

    Lst <- list(...)

    if( !is.null(Lst$col) )

    { col <- Lst$col}else

    { col<-"black"}



    if( !is.null(Lst$main) )

    { main <- Lst$main}else

    { main<-gene}



    if( !is.null(Lst$type) )

    { type <- Lst$type}else

    { type<-"b"}





  times=c(0,9,seq(12,57,3))

  yl=c(round(min(dataset[gene,]-1)),round(max(dataset[gene,])+1))

  xl=c(0,57)

  i=plot(times,dataset[gene,],xlab="time after stimulus",ylab="transcript accumulation",ylim=yl,xlim=xl,type=type,main=main,col=col,lwd=2)

  return(i)}
