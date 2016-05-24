###################################################
### chunk number 1: 
###################################################
create.disProg <- function(week, observed, state, start=c(2001,1), freq=52, neighbourhood=NULL, populationFrac=NULL,epochAsDate=FALSE){
  namesObs <-colnames(observed)
  
  # check whether observed contains only numbers
  if(!all(sapply(observed, is.numeric))){
    stop("\'observed\' must be a matrix with numbers\n")
  }
   
  #univariate timeseries ?
  if(is.vector(observed)){
    observed <- matrix(observed,ncol=1)
    namesObs <- deparse(quote(observed))
  } else {  # ensure we have a matrix
    observed <- as.matrix(observed)
  }
  
  if(missing(state)){
    state <- 0*observed
  } else if(is.vector(state)){
    state <- matrix(state,ncol=1)
  } else {
    state <- as.matrix(state)
  }
    
  #check number of columns of observed and state
  nAreas <- ncol(observed)
  nObs <- nrow(observed)
  if(ncol(observed) != ncol(state)){
    #if there is only one state-vector for more than one area, repeat it
    if(ncol(state)==1) {
      state <- matrix(rep(state,nAreas),ncol=nAreas,byrow=FALSE)
    } else { 
      cat('wrong dimensions of observed and state \n')
      return(NULL)
    }
  }
  
  #check neighbourhood matrix
  # neighbourhood can be a matrix or an array of dimension c(nAreas,nAreas, nrow(observed))
  if(!is.null(neighbourhood) ) {
    dimNhood <- dim(neighbourhood)
    if(length(dimNhood)==2 & any(dimNhood != nAreas)) {
      cat('wrong dimensions of neighbourhood matrix \n')
      return(NULL)
    } else if (length(dimNhood)==3 & (any(dimNhood[1:2] != nAreas) | (dimNhood[3] != nrow(observed)) )){
      cat('wrong dimensions of neighbourhood matrix \n')
      return(NULL)
    }
  } else {
     # no neighbourhood specified
     neighbourhood <- matrix(NA,nrow=nAreas,ncol=nAreas)
  }
  
  if(is.null(populationFrac)) {
    populationFrac <- matrix(1/ncol(observed),nrow=nObs, ncol=ncol(observed))
  } else {
    # make sure populationFrac is a matrix
    populationFrac <- as.matrix(populationFrac)
    # check dimensions
    if(nrow(populationFrac)!= nObs | ncol(populationFrac)!= nAreas)
      stop("dimensions of \'populationFrac\' and \'observed\' do not match\n")
    # check whether populationFrac contains only numbers
    if(!all(sapply(populationFrac, is.numeric))){
      stop("\'populationFrac\' must be a matrix with real numbers\n")
  }
  }

  #labels for observed and state
  if(is.null(namesObs)){
    namesObs <- paste(deparse(quote(observed)),1:nAreas,sep="")

  }
  
  colnames(observed) <- namesObs
  colnames(state) <- namesObs
  
  res <- list("week"=week, "observed"=observed, "state"=state, "start"=start, "freq"=freq,  "neighbourhood"=neighbourhood, "populationFrac"=populationFrac,"epochAsDate"=epochAsDate)
  class(res) <- "disProg"
  return(res)
}

print.disProg <- function(x, ...) {
  cat( "-- An object of class disProg -- \n" )
  cat( "freq:\t\t", x$freq,"\n" )
  cat( "start:\t\t", x$start,"\n" )
  cat( "dim(observed):\t", dim(x$observed), "\n\n")

  n <- 1
  cat("Head of observed:\n")
  print(head(x$observed,n))

  #cat("\nhead of neighbourhood:\n")
  #print( head(x$neighbourhood,n))
}

###################################################
### chunk number 2: 
###################################################
sumNeighbours <- function(disProgObj){

  observed <- disProgObj$observed
  neighbours <- matrix(nrow=nrow(observed),ncol=ncol(observed))
    
  for(i in 1:ncol(observed)){
    #only one neighbour
    if(sum(disProgObj$neighbourhood[,i])==1)
      neighbours[,i] <- observed[,disProgObj$neighbourhood[,i]==1]
    #more than one neighbour
    else
      neighbours[,i] <- apply(observed[,disProgObj$neighbourhood[,i]==1], MARGIN=1, sum)
  }
  return(neighbours)
}


###################################################
### chunk number 3: 
###################################################
aggregate.disProg <- function(x,...){
  #aggregate observed counts
  observed <- apply(x$observed,MARGIN=1,sum)
  #aggregate states
  state <- apply(x$state,MARGIN=1,sum)
  state[state > 1] <- 1
  
  #create univariate disProg object
  x <- create.disProg(week=x$week, observed=observed, state=state, freq=x$freq,start=x$start)
  return(x)
}


###################################################
### chunk number 4: 
###################################################

plot.disProg.one <- function(x, title = "", xaxis.years=TRUE, quarters=TRUE, startyear = x$start[1], firstweek = x$start[2], ylim=NULL, xlab="time", ylab="No. infected",type="hh",lty=c(1,1),col=c(1,1), outbreak.symbol = list(pch=3, col=3),legend.opts=list(x="top", legend=c("Infected", "Outbreak"),lty=NULL,pch=NULL,col=NULL),...) {


  observed <- x$observed
  state    <- x$state

  # width of the column
  tab <- 0.5

  # left/right help for constructing the columns
  observedxl <- (1:length(observed))-tab
  observedxr <- (1:length(observed))+tab
        
  # control where the highest value is
  max <- max(observed)
        
  #if ylim is not specified
  if(is.null(ylim)){
    ylim <- c(-1/20*max, max)
  }

  #Plot the results using one Large plot call
  matplot(x=cbind(observedxl, observedxr),y=cbind(observed, observed),xlab=xlab,ylab=ylab,
          type=type,lty=lty, col=col, ylim=ylim,axes = !(xaxis.years),...)
  
  #Show the outbreaks
  if (!is.null(outbreak.symbol)) {
    for(i in 1:length(observed)){
      matlines( c(i-tab, i+tab), c(observed[i],observed[i]) )
      if(state[i] == 1)
        matpoints( i, ylim[1], pch=outbreak.symbol$pch, col=outbreak.symbol$col)
    }
  }
        
  title(title)
  cex <- par()$cex.axis
        
  #Label of x-axis 
  if(xaxis.years){        
    # get the number of quarters lying in range for getting the year and quarter order
    obsPerYear <- x$freq
    obsPerQuarter <- x$freq/4
    myat.week <- seq(ceiling((obsPerYear-firstweek+1)/obsPerQuarter) * obsPerQuarter + 1, length(observed)+(floor((obsPerYear-firstweek + 1)/obsPerQuarter) * obsPerQuarter +1), by=obsPerQuarter)
    # get the right year order
    year <- (myat.week - obsPerYear) %/% obsPerYear + startyear
    # function to define the quarter order
    quarterFunc <- function(i) { switch(i+1,"I","II","III","IV")}
    # get the right number and order of quarter labels
    quarter <- sapply( (myat.week-1) %/% obsPerQuarter %% 4, quarterFunc)
    # get the positions for the axis labels
    myat.week <- myat.week - (obsPerYear - firstweek + 1)

    # construct the computed axis labels
    if (quarters) {
      if (cex == 1) {
        mylabels.week <- paste(year,"\n\n",quarter,sep="")
      } else {
        mylabels.week <- paste(year,"\n",quarter,sep="")
      }
    } else {
      mylabels.week <- paste(year,sep="")
    }
      
        
    axis( at=myat.week , labels=mylabels.week , side=1, line = 1 )
    axis( side=2 )
  }
        
  #should there be a legend? 
  if(!is.null(legend.opts) && (class(legend.opts) == "list")) {
    #Fill empty (mandatory) slots in legend.opts list
    if (is.null(legend.opts$lty)) legend.opts$lty = c(lty[1],NA)
    if (is.null(legend.opts$col)) legend.opts$col = c(col[1],outbreak.symbol$col)
    if (is.null(legend.opts$pch)) legend.opts$pch = c(NA,outbreak.symbol$pch)
    if (is.null(legend.opts$x))   legend.opts$x = "top"
    if (is.null(legend.opts$legend)) legend.opts$legend = c("Infected", "Outbreak")

    #Create the legend
    do.call("legend",legend.opts)
  }
  
  invisible()
}

plot.disProg <- function(x, title = "", xaxis.years=TRUE, startyear = x$start[1], firstweek = x$start[2], as.one=TRUE, same.scale=TRUE, ...){
  observed <- x$observed
  state    <- x$state
  
  #univariate timeseries ?
  if(is.vector(observed))
    observed <- matrix(observed,ncol=1)
  if(is.vector(state))
    state <- matrix(state,ncol=1)
  
  nAreas <- ncol(observed)
  max <- max(observed)
  
  #check if x is multivariate or univariate 
  
  #multivariate time series
  if(nAreas > 1){
    #all areas in one plot -- not supported in sts
    if(as.one){
      matplot(observed,type="l",lty=1:nAreas,col=1:nAreas,ylim=c(0, 1.1*max),xlab="time",ylab="No. of Infected", axes=!xaxis.years)
      #If no legend.opts is specified or not set to null
      if ((is.na(pmatch("legend.opts",names(list(...))))) | 
          (!is.na(pmatch("legend.opts",names(list(...)))) & (!is.null(list(...)$legend.opts)))) {
        legend.opts <- list(...)$legend.opts
        if (is.null(legend.opts$x)) legend.opts$x = "topleft"
        if (is.null(legend.opts$legend)) legend.opts$legend = colnames(observed)
        if (is.null(legend.opts$col)) legend.opts$col = 1:nAreas
        if (is.null(legend.opts$lty)) legend.opts$lty = 1:nAreas
        if (is.null(legend.opts$ncol)) legend.opts$ncol = 5
        if (is.null(legend.opts$bty)) legend.opts$bty = "n"

        do.call("legend",legend.opts)
      }

      title(title)
      
      if(xaxis.years){  #todo: move this as output of ONE function
          # get the number of quarters lying in range for getting the year and quarter order
          myat.week <- seq(ceiling((52-firstweek+1)/13) * 13 + 1, length(observed)+(floor((52-firstweek + 1)/13) * 13 +1), by=13)
          # get the right year order
          year <- (myat.week - 52) %/% 52 + startyear
          # function to define the quarter order
          quarterFunc <- function(i) { switch(i+1,"I","II","III","IV")}
          # get the right number and order of quarter labels
          quarter <- sapply( (myat.week-1) %/% 13 %% 4, quarterFunc)
          # get the positions for the axis labels
          myat.week <- myat.week - (52 - firstweek + 1)

          # construct the computed axis labels
          cex <- par()$cex.axis
          if (cex == 1) {
            mylabels.week <- paste(year,"\n\n",quarter,sep="")
          } else {
            mylabels.week <- paste(year,"\n",quarter,sep="")
          }
        
          axis( at=myat.week , labels=mylabels.week , side=1, line = 1 )
          axis( side=2 )
       }
      
    } else {  #plot each area
      #set window size     
      par(mfrow=magic.dim(nAreas),mar=c(2,1,1,1))
      
      if(same.scale)
        ylim <- c(-1/20*max, max)
      else
        ylim <- NULL
      
      #plot areas
      k <- 1:nAreas
      sapply(k, function(k) {
         plot.disProg.one(create.disProg(x$week, observed[,k], state[,k], freq=x$freq,start=x$start), 
                          title = "", startyear = startyear, firstweek = firstweek, 
                          xaxis.years=xaxis.years, ylim=ylim, legend.opts=NULL, ... )   
         mtext(colnames(observed)[k],line=-1.3)     
         })
      #reset graphical params
      par(mfrow=c(1,1), mar=c(5, 4, 4, 2)+0.1)
    }
  } else {  #univariate time series
    plot.disProg.one(x=x, title = title, startyear = startyear, firstweek = firstweek, xaxis.years=xaxis.years, ...)
  }
  invisible()
}


###################################################
### chunk number 5: 
###################################################

plot.survRes.one <- function(x, method=x$control$name, disease=x$control$data, domany=FALSE,ylim=NULL,xaxis.years=TRUE,startyear = 2001, firstweek = 1, xlab="time", ylab="No. infected", main=NULL, type="hhs",lty=c(1,1,2),col=c(1,1,4), outbreak.symbol = list(pch=3, col=3),alarm.symbol=list(pch=24, col=2),legend.opts=list(x="top",legend=c("Infected", "Upperbound", "Alarm", "Outbreak"),lty=NULL,col=NULL,pch=NULL), ...) {

  ################## Handle the NULL arguments ########################################################
  if (is.null(main)) main = paste("Analysis of ", as.character(disease), " using ", as.character(method),sep="") 
  #No titles are drawn when more than one is plotted.
  if (domany) main = ""

  survResObj <- x
  observed <- survResObj$disProgObj$observed[survResObj$control$range]
  state    <- survResObj$disProgObj$state[survResObj$control$range]
        
  #print(list(...))
  # width of the column
  tab <- 0.5

  # left/right help for constructing the columns
  observedxl <- (1:length(observed))-tab 
  observedxr <- (1:length(observed))+tab
  upperboundx <- (1:length(survResObj$upperbound)) #-0.5 
  
  # control where the highest value is
  max <- max(max(observed),max(survResObj$upperbound))
  
  #if ylim is not specified
  #if(is.null(ylim)){
  #  ylim <- c(-1/20*max, max)
  #}

#~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (is.null(ylim)) {
    max <- max(max(observed), max(survResObj$upperbound))
    ylim <- c(-1/20 * max, max)
  } else {
    max <- ylim[2]
  }
  
  #ensure that there is enough space for the alarm/outbreak symbols   
  if(ylim[1]>=0)
    ylim[1] <- -1/20*max
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

  #Generate the matrices to plot
  xstuff <- cbind(observedxl, observedxr, upperboundx) #no adjusting + min(x$control$range) - 1
  ystuff <- cbind(observed, observed, survResObj$upperbound)
        

  #Plot the results using one Large plot call (we do this by modifying
  #the call). 
  matplot(x=xstuff,y=ystuff,xlab=xlab,ylab=ylab,main=main,ylim=ylim,axes = !(xaxis.years),type=type,lty=lty,col=col,...)

  
  if (!is.null(survResObj$aggr)) {
    points(upperboundx+tab,survResObj$aggr,col=1)
  }
  
  for(i in 1:length(observed)){
    matlines( c(i-tab, i+tab), c(observed[i],observed[i]),col=col[1])
    if(survResObj$alarm[i] == 1)
      matpoints( i, -1/40*max, pch=alarm.symbol$pch, col=alarm.symbol$col)
    if(state[i] == 1)
      matpoints( i, -1/20*max, pch=outbreak.symbol$pch, col=outbreak.symbol$col)
  }

  # check where to place the legend. If the left upper side is free place it there
  if (max * 2/3 >= max(
                max(observed[1:floor(1/4 * length(observed))]),
                max(survResObj$upperbound[1:floor(1/4 * length(survResObj$upperbound))])
                )) {
    xlegpos <- 0
  }

  #Label of x-axis 
  if(xaxis.years){        
    # get the number of quarters lying in range for getting the year and quarter order
    myat.week <- seq(ceiling((52-firstweek+1)/13) * 13 + 1, length(observed)+(floor((52-firstweek + 1)/13) * 13 +1), by=13)
    # get the right year order
    year <- (myat.week - 52) %/% 52 + startyear
    # function to define the quarter order
    quarterFunc <- function(i) { switch(i+1,"I","II","III","IV")}
    # get the right number and order of quarter labels
    quarter <- sapply( (myat.week-1) %/% 13 %% 4, quarterFunc)
    # get the positions for the axis labels
    myat.week <- myat.week - (52 - firstweek + 1)

    # construct the computed axis labels
    #cex <- par()$cex.axis
    #if (cex == 1) {
    mylabels.week <- paste(year,"\n\n",quarter,sep="")
    #} else {
    #  mylabels.week <- paste(year,"\n",quarter,sep="")
    #}
    
    axis( at=myat.week , labels=mylabels.week , side=1, line = 1 )
    axis( side=2 )
  }
  
  if(!is.null(legend.opts) && (class(legend.opts) == "list")) {
    #Fill empty (mandatory) slots in legend.opts list
    if (is.null(legend.opts$lty)) legend.opts$lty = c(lty[1],lty[3],NA,NA)
    if (is.null(legend.opts$col)) legend.opts$col = c(col[1],col[3],alarm.symbol$col,outbreak.symbol$col)
    if (is.null(legend.opts$pch)) legend.opts$pch = c(NA,NA,alarm.symbol$pch,outbreak.symbol$pch)
    if (is.null(legend.opts$x))   legend.opts$x = "top"
    if (is.null(legend.opts$legend)) 
      legend.opts$legend = c("Infected", "Upperbound", "Alarm", "Outbreak")

    do.call("legend",legend.opts)
  }
  
  
  invisible()
}


#the main function -- cant we do better than this?
plot.survRes <- function(x, method=x$control$name, disease=x$control$data, xaxis.years=TRUE,startyear = 2001, firstweek = 1, same.scale=TRUE,...) {
  observed <- x$disProgObj$observed
  state <- x$disProgObj$state
  alarm <- x$alarm

  #univariate timeseries ?
  if(is.vector(observed))
    observed <- matrix(observed,ncol=1)
  if(is.vector(state))
    state <- matrix(state,ncol=1)
  if(is.vector(alarm)) 
    alarm <- matrix(alarm,ncol=1)
  nAreas <- ncol(observed)
  max <-  max(max(observed),max(x$upperbound))

  #multivariate time series
  if(nAreas > 1){
    #all areas in one plot 
      #set window size     
      par(mfrow=magic.dim(nAreas),mar=c(2,1,1,1))
      
      if(same.scale) {
        ylim <- c(-1/20*max, max)
      } else {
        ylim <- NULL
      }
      
      #plot areas
      k <- 1:nAreas
      sapply(k, function(k) {
        #Create the survRes
        dP <- create.disProg(x$disProgObj$week, observed[,k], state[,k],start=x$start)
        obj <- list(alarm=alarm[,k],disProgObj=dP,control=x$control,upperbound=x$upperbound[,k])
        class(obj) <- "survRes"
        plot.survRes.one(obj,startyear = startyear, firstweek = firstweek, 
                         xaxis.years=xaxis.years, ylim=ylim, legend.opts=NULL,domany=TRUE,... )   
         mtext(colnames(observed)[k],line=-1.3)     
         })
      #reset graphical params
      par(mfrow=c(1,1), mar=c(5, 4, 4, 2)+0.1)
    }
  else {  #univariate time series
    plot.survRes.one(x=x, startyear = startyear, firstweek = firstweek, xaxis.years=xaxis.years, domany=FALSE,...)
  }
  invisible()
}



###################################################
### chunk number 6: 
###################################################
magic.dim <- function(k){
  if(k==1)
    return(c(1,1))
  
  #factorize k  
  factors <- primeFactors(k)
  
  #find the best factorization of k into two factors
  res <- bestCombination(factors)
    
  #if k is a prime or the difference between the two factors of k is too large
  #rather use the roots of the next square number greater than k 
  
  #up is root of the smallest square number >= k
  up <- ceiling(sqrt(k))
  #low is root of the biggest square number < k
  low <- up -1
  
  if(diff(res) >5){
    # e.g. k=11 is a prime, the next square number is 16 so up=4 and low=3
    # low^2 = 9 < 11 is naturally too small, up^2=16 > 11 so c(4,4) is a solution
    # but low*up = 3*4 = 12 > 11 is also adequate and a better solution
    if((k - low^2) < up)
      res <- c(low,up)
    else
      res <- c(up,up)
  }  
  
  return(sort(res))
  
}


###################################################
### chunk number 7: 
###################################################
primeFactors <- function(x){
  if(x==1)
    return(1)
    
  factors<- numeric(0)  
  i<-1
  
  #start with i=2 and divide x by i (as often as possible) then try division by i+1
  #until all factors are found, i.e. x=1 
  while(i < x){
    i <- i+1   
    
    while((x %% i)==0){
      # each time a new factor i is found, save it and proceed with x = x/i 
      # e.g. k=20: 2 is a factor of x=20, continue with x = 10 = 20/2  
      #            2 is a factor of x=10, continue with x = 5 = 10/2
      #            3 and 4 are no factors of x = 5   
      #            5 is a factor of x = 5, continue with x = 1 
      # result: 20 = c(2, 2, 5)
      factors <- c(factors, i)
      x <- x/i
    }
  }
  return(factors)
}


###################################################
### chunk number 8: 
###################################################
######################################################################
# Given a prime number factorization of a number, e.g. 36
# yields x=c(2,2,3,3)
# and parition x into two groups, such that the product of the numbers
# in group one is as similar as possible to the product
# of the numbers of group two. This is useful in magic.dim
#
# Params:
#  x - the prime number factorization
#
# Returns:
#  c(prod(set1),prod(set2))
######################################################################

bestCombination <- function(x) {
  #Compute the power set of 0:1^length(x), i.e. a binary indicator for
  #variable stating whether to include it in set 1 or not.
  combos <- as.matrix(expand.grid(rep(list(0:1),length(x))))
  mode(combos) <- "logical"
  
  #Small helper function, given a vector of length(x) stating whether
  #to include an element in set1 or not, compute the product
  #of set1 and set2=x\backslash set1
  #set1: all those for which include is TRUE, set2: all those for which
  #include is FALSE
  setsize <- function(include) { c(prod(x[include]),prod(x[!include])) }

  #Compute the product of set1 and set2 for each possible combination
  sizes <- apply(combos,MARGIN=1,FUN=setsize)
  #Calculate the combination, where x is as close to y as possible
  bestConfig <- combos[which.min(abs(diff(sizes))),]
  #Return this setsize of this configuration
  return(setsize(bestConfig))
}




