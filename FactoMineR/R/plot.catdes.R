plot.catdes <- function (x,col="deepskyblue",show="all",numchar=10,...){
  l=max(length(x$quanti),length(x$category))		# measure the length of x
  long=rep(0,l)
                                        #print(x$category)
  list.catdes=list(long)
  minimum=0
  maximum=0
  count=0
  for (i in 1:l){
                                        #print(i)
    if (!is.null(x$quanti[[i]])){
      quanti=as.data.frame(x$quanti[[i]])
      quanti.catdes=as.vector(quanti[,1])
      names(quanti.catdes)=rownames(quanti)
    } else quanti.catdes=NULL
    if (!is.null(x$category[[i]])){
      category=as.data.frame(x$category[[i]])
      category.catdes=as.vector(category[,5])
      names(category.catdes)=rownames(category)
    } else category.catdes=NULL
                                        #print(category.catdes)
    
    if (show=="all") catdes.aux=c(quanti.catdes,category.catdes)
    if (show=="quanti") catdes.aux=quanti.catdes					# different options
    if (show=="quali") catdes.aux=category.catdes
    
    if (!is.null(catdes.aux)) {
      count=count+1								#count is a counter of the catdes clusters which are non null.   
      long[i]=length(catdes.aux)
      minimum=min(catdes.aux,minimum)			# find the longest catdes clusters and the smallest.
      maximum=max(catdes.aux,maximum)
    } else long[i]=0
    list.catdes[[i]]=catdes.aux								# list the catdes clusters
  }
  if(count!=0){
    if (count<=4){ 
      numc=count									# design of the graphic window
      numr=1
    } else{
      numc=4
      numr=round(count/4)+1
    }
    dev.new()
    par(las = 3)
    par(mfrow = c(numr, numc))
    for(i in 1:l){
      catdes.aux=list.catdes[[i]]
      if(!is.null(catdes.aux)){
        catdes.aux=sort(catdes.aux,decreasing=FALSE)					#plot the catdes for every cluster in the graphic window
        barplot(catdes.aux, width =c(1,1), col = col, border = "black", 
                ylim = c(minimum-1,maximum+1),xlim=c(0,max(long)+1), 
                main = unique(names(x$category)[i],names(x$quanti)), cex.names = 1, ylab="v.test", 
                names.arg = substr(names(catdes.aux), 1, numchar))
	
      }
    }
  }    
  par(las = 0)
}




