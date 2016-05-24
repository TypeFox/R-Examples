ImagePoly <-
function(poly,dataframe,context.id,names=NULL,show.names=NULL,method="equal",nbr=10,...){
  #method <- "equal"
  #nbr <- 10
  ## Context.id must be a field of both the shapefile and the dataframe

  if (!(context.id %in%names(poly)&context.id%in%names(dataframe)))
    {stop("The field ",context.id , " is not in both the shapefile and the dataframe")}


  polyA <- performAddFields(poly,dataframe,context.id,names)
  names_coded <- as.list(polyA[[2]])

  if (is.null(show.names)){show.names<-names_coded}

  poly <- polyA[[1]]


  P<-mapply(function(x,y){performSinglePlot(poly,x,y,method,nbr)},names_coded,show.names,SIMPLIFY=FALSE,...)


  N <- length(P)
  even <- 2*(N%/%2)==N
  x <- as.integer(N/2)
  
  if(length(P)==1){
    plot(P[[1]],more=FALSE)
  }else{
    if(length(P)==2){
      #for (i in 1:x){plot(P[[i]],split=c(i,1,x,2),more=TRUE)}
      plot(P[[1]],split=c(1,1,2,1),more=TRUE)
      plot(P[[2]],split=c(2,1,2,1),more=FALSE)
    }else{
      if(length(P)==3){
        plot(P[[1]],split=c(1,1,2,2),more=TRUE)
        plot(P[[2]],split=c(2,1,2,2),more=TRUE)
        plot(P[[3]],split=c(1,2,2,2),more=FALSE)
      }else{
        if (!even)
          {
            x <-x+1
            for (i in 1:x){plot(P[[i]],split=c(i,1,x,2),more=TRUE)}
            for (i in (x+1):(length(P)-1)){plot(P[[i]],split=c(i-x,2,x,2),more=TRUE)}
            plot(P[[length(P)]],split=c(length(P)-x,2,x,2),more=FALSE)
          }
        else
          {
            for (i in 1:x){plot(P[[i]],split=c(i,1,x,2),more=TRUE)}
            for (i in (x+1):(length(P)-1)){plot(P[[i]],split=c(i-x,2,x,2),more=TRUE)}
            plot(P[[length(P)]],split=c(length(P)-x,2,x,2),more=FALSE)
          }
      }
    }  
  }
}
