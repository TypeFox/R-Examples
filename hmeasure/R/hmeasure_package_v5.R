HMeasure <- function(true.class, scores,
                     severity.ratio=NA, threshold=0.5, level=0.95
                     ){


  ####################
  ### INPUT CHECKS ###
  ####################

  # try to catch mistaken order of arguments
  if (is.matrix(true.class) || is.data.frame(true.class)){
	stop(	'True class should be a vector, not a matrix / data frame.  Consider the order of the arguments.'	)
  }

  # no missing values in the labels allowed
  if (any(is.na(true.class))){
    stop('Missing values in class labels are not allowed.')}


  # relabel, and make sure there are only 2 class labels
  true.class <- relabel(true.class) 
  
  # row names can confuse and are otherwise useless - remove them
  rownames(scores) <- NULL
  rownames(true.class) <- NULL
 
  # turn scores into a data frame (if it were not one already)
  if (is.vector(scores)){
    scores <- as.data.frame(scores)
    # message('Scores coerced from vector to data frame')
    
  }

  if (is.matrix(scores)){
    n <- dim(scores)[1]
    k <- dim(scores)[2]
  
  # in the case of a matrix, throw a warning if columns (classifiers) > rows (data)
    if (n < k) {
      warning(gettextf(
        'Consider transposing score matrix: number of classifiers (columns) = %d exceeds number %d of datapoints (rows)', 
        k, n), domain = NA)
    }

    scores <- as.data.frame(scores)
    # message('Scores coerced from matrix to data frame')
    
  }
  

  
  if (dim(scores)[1]!=length(true.class)){
    stop('Label vector provided has different length than respective classifier scores')
  }

  # only look at complete cases in the score data frame
  if (any(is.na(scores))){
    warning(	'Missing entries detected in matrix of scores. Respective entries will be disregarded'	)
  }
  complete.rows <- complete.cases(scores)
  scores <- subset(scores,subset=complete.rows)
  true.class <- subset(true.class,subset=complete.rows)
  rownames(scores) <- NULL
  rownames(true.class) <- NULL

  # now that format is correct, get sample size and number of classifiers  
  n <- dim(scores)[1]
  k <- dim(scores)[2]
     
   
     
  # THRESHOLD - if only one value for the threshold has been provided
  # (e.g., the default of 0.5), use the same for all classifiers
  # else check that the array of thresholds has one per classifier
  
  if (length(threshold) == 1){
  	threshold <- rep(threshold,k)
  } else {  
    if (length(threshold)<k){
    	warning(	'Threshold must either be a single value, or a vector of length equal to the number of classifiers employed. The default value of 0.5 will be used.')
    }
  }

  ############################
  ### INPUT CHECK COMPLETE ###
  ############################

  #############################################
  ### SINGLE CLASSIFIER FUNCTION DEFINITION ###
  #############################################
  # to keep the code tidy, we implement an internal function for a single classifier 
  HMeasure.single <- function(y, s, classifier.name=NULL,
				severity.ratio=severity.ratio,
  				threshold=threshold, level=level	){

    # PROCESSING 
    n <- length(s)

    
    # this is a numeric version of the class labels
    n1 <- sum(y) 
    n0 <- n-n1
    pi0 <- n0/n
    pi1 <- n1/n
  
    # retrieve severity ratio - set to default if absent
    if (is.na(severity.ratio)){
      severity.ratio <- pi1/pi0
    } 



    # order data into increasing scores
    zord <- order(s)
    sc <- s[zord]

    # note: we make no assumptions about the range of s

    # COMPUTE ROC CURVE

    # Calculate raw ROC, replacing any tied sequences by a diagonal
    # Raw ROC starts at F0[1]=0, F1[1]=0, and ends at F0[K1]=1, F1[K1]=1.
    Get.Score.Distributions <- function(y,s,n1,n0){
      # tapply(y,s,sum) counts the instances of each unique score, and ranks them by score
      s1 <- unname(tapply(y, s, sum))/n1
      s1 <- c(0,s1,1-sum(s1)) # make sure to add the points 0,0 and 1,1
      s0 <- unname(tapply(1-y, s, sum))/n0
      s0 <- c(0,s0,1-sum(s0)) # make sure to add the points 0,0 and 1,1
      
      # number of unique scores
      S <- length(s1)
      # what were r0i and r1i in ML paper are now the empirical cdfs
      F1 <- cumsum(s1)
      F0 <- cumsum(s0)
      return(list(F1=F1,F0=F0,s1=s1,s0=s0,S=S))
    }



    out.scores <- Get.Score.Distributions(y=y,s=s,n1=n1,n0=n0)
    AUC <- 1- sum(out.scores$s0 * (out.scores$F1 - 0.5 * out.scores$s1))
    # if the AUC < .5, switch signs and repeat
    switched <- FALSE
    the.criterion <- AUC < 0.5
    if (the.criterion){
      switched <- TRUE      
      s <- 1-s
      out.scores <- Get.Score.Distributions(y,s,n1,n0)
      if (is.null(classifier.name)){
      warning('ROC curve mostly lying under the diagonal. Switching scores.', domain = NA)
	} else {
      warning(gettextf( 'ROC curve of %s mostly lying under the diagonal. Switching scores.',
      	classifier.name), domain = NA)
	}
    }
    F1 <- out.scores$F1
    F0 <- out.scores$F0 
    s0 <- out.scores$s0
    s1 <- out.scores$s1
    S <- out.scores$S


	# get misclassification statistics
	misclass.out <- misclassCounts(as.numeric(s>threshold),true.class)
	misclass.metrics <- misclass.out$metrics
	temp <- misclass.out$conf.matrix
	misclass.conf <- data.frame(
		TP=temp[1,1], FP=temp[2,1],
		TN=temp[2,2], FN=temp[1,2])
      			
    
    # get aggregate statistics:
    AUC <- 1- sum(s0 * (F1 - 0.5 * s1)) # REPLACING TIED SCORES BY A DIAGONAL
    Gini <- 2*AUC - 1
    KS <- max(abs(F0 - F1))
    cost.parameter <- severity.ratio/(1+severity.ratio)
    MER <- min(pi0*(1-F0)+pi1*F1)
    MWL <- 2*min(cost.parameter*pi0*(1-F0)+(1-cost.parameter)*pi1*F1)



    Look.Up.AUC <- function(xcurve,ycurve,x=0){
      # assumes the curve is monotonic
      result <- NA
      if (all(diff(xcurve) >= 0)){
        ind <- which(xcurve-x>0)[1]
        x1 <- xcurve[ind-1]
        x2 <- xcurve[ind]
        y1 <- ycurve[ind-1]
        y2 <- ycurve[ind]

        if (x2-x1 > 0) {
          pos <- (x2-x)/(x2-x1)
          result <- (1-pos)*y1 + pos*y2         
        } else {result <- y2}
      }
      return(result)
    }

    SensFixed <- matrix(NA,1,length(level))
    SpecFixed <- matrix(NA,1,length(level))
    temp <- array(NA,length(level))
    for (l in 1:length(level)){
      SensFixed[l] <- c(Look.Up.AUC(F0,1-F1,x=level[l]))
      temp[l] <- paste('Sens.Spec',floor(level[l]*100),sep='')
    } 
    SensFixed <- as.data.frame(SensFixed)
    colnames(SensFixed) <- temp

    for (l in 1:length(level)){
      SpecFixed[l] <- Look.Up.AUC(F1,F0,x=1-level[l])
      temp[l] <- paste('Spec.Sens',floor(level[l]*100),sep='')
    }
    SpecFixed <- as.data.frame(SpecFixed)
    colnames(SpecFixed) <- temp

    # restrict to upper convex hull by considering ROC above diagonal only
    chull.points <- chull(1-F0,pmax(1-F1,1-F0))
    G0 <- 1-F0[chull.points]
    G1 <- 1-F1[chull.points] 
    hc <- length(chull.points)
    sG0 <- c(0,G0[c(2:length(G0))] - G0[c(1:(length(G0)-1))])
    sG1 <- c(0,G1[c(2:length(G1))] - G1[c(1:(length(G1)-1))])
    AUCH <- sum(sG0 * (G1 - 0.5 * sG1))


    # get sorted scoring densities
    s.class0 <- sort(s[y==0])
    s.class1 <- sort(s[y==1])


    # Calculate the LHshape1 value
    cost <- c(1:(hc+1))
    b0 <- c(1:hc+1)
    b1 <- c(1:hc+1)
    
    # extract shape
    if (severity.ratio > 0){
      shape1 <- 2
      shape2 <- 1+(shape1-1)*1/severity.ratio
    }
    if (severity.ratio < 0){
      shape1 <- pi1+1
      shape2 <- pi0+1
    }
    cost[1] <- 0
    cost[hc+1] <- 1

    b00 <- beta(shape1,shape2)
    b10 <- beta(1+shape1,shape2)
    b01 <- beta(shape1,1+shape2)


    b0[1] <-
      pbeta(cost[1], shape1=(1+shape1), shape2=shape2)*b10/b00

    b1[1] <-
      pbeta(cost[1], shape1=shape1, shape2=(1+shape2))*b01/b00

    b0[hc+1] <-
      pbeta(cost[hc+1], shape1=(1+shape1), shape2=shape2)*b10/b00

    b1[hc+1] <-
      pbeta(cost[hc+1], shape1=shape1, shape2=(1+shape2))*b01/b00

    ### NB: can become massively faster
    for (i in 2:hc){
      cost[i] <- pi1*(G1[i]-G1[i-1]) / 
      (pi0*(G0[i]-G0[i-1]) + pi1*(G1[i]-G1[i-1]))

      b0[i] <-
        pbeta(cost[i], shape1=(1+shape1), shape2=shape2)*b10/b00

      b1[i] <-
        pbeta(cost[i], shape1=shape1, shape2=(1+shape2))*b01/b00
    }

    LHshape1 <- 0
    for (i in 1:hc){
      LHshape1 <- LHshape1 + pi0*(1-G0[i])*(b0[(i+1)]-b0[i]) + pi1*G1[i]*(b1[(i+1)]-b1[i])
    }

    B0 <- 
      pbeta(pi1, shape1=(1+shape1), shape2=shape2)*b10/b00

    B1 <-
      pbeta(1, shape1=shape1, shape2=(1+shape2))*b01/b00 -
      pbeta(pi1, shape1=shape1, shape2=(1+shape2))*b01/b00

    H <- 1 - LHshape1/(pi0*B0 + pi1*B1)

    data <- list(F0=F0, F1=F1, G0=G0, G1=G1, cost=cost,
                 pi1=pi1, pi0=pi0, n0=n0, n1=n1, n=n, hc=hc,
                 s.class0=s.class0, s.class1=s.class1,
                 severity.ratio=severity.ratio)

    metrics <- data.frame(H=H, Gini=Gini, AUC=AUC, AUCH=AUCH, KS=KS, MER=MER, MWL=MWL)
    metrics <- cbind(metrics,SpecFixed,SensFixed)
    metrics <- cbind(metrics,misclass.metrics,misclass.conf)

    return(list(data=data,metrics=metrics))
  }


  #############################################
  ### SINGLE CLASSIFIER DEFINITION complete ###
  #############################################


  ######################################
  ### PROCESS CLASSIFIERS ONE BY ONE ###
  ######################################

  	
  data <- list()
  for (count in 1:k){
    name.now <- colnames(scores)[count]
    s <- scores[,count]
    threshold.now <- threshold[count]
    output <- HMeasure.single(y=true.class, s=s, classifier.name=name.now,
                              severity.ratio=severity.ratio, 
                              threshold=threshold.now,level=level)

    if (count == 1){
      metrics <- output$metrics
    }
    if (count > 1){metrics <- rbind(metrics,output$metrics)}
    
    # retrieve data for plotting purposes
    data[[count]] <- output$data
        
  }

  # name the rows by classifier
  rownames(metrics) <- colnames(scores)
  # name the data output by classifier
  names(data) <- colnames(data)

  # construct output
  hmeasure <- list(metrics=metrics)
  attr(hmeasure,'data') <- data
  class(hmeasure) <- 'hmeasure'
  return(hmeasure)

} 


### SUMMARY METHOD ###
summary.hmeasure <- function(object, show.all=FALSE,...){
  new.object <- unclass(object)
  new.object$data <- NULL
  
  noof <- dim(new.object$metrics)[2]
  
  if (show.all){
  	print(as.data.frame(new.object$metrics))

  	} else {
  	print(as.data.frame(
  	    new.object$metrics[,1:(noof-13)]
  		))  	

  }
}

### PRINT METHOD ###
summary.hmeasure <- function(object, show.all=FALSE,...){
  new.object <- unclass(object)
  new.object$data <- NULL
  
  noof <- dim(new.object$metrics)[2]
  
  if (show.all){
  	print(as.data.frame(new.object$metrics))

  	} else {
  	print(as.data.frame(
  	    new.object$metrics[,1:(noof-13)]
  		))  	

  }
}

### PLOT METHOD ###

plotROC <- function(results, which=1,bw='nrd0'){

  if (class(results)!='hmeasure'){
    stop('The first argument of plotROC() must be of class \'hmeasure\'')
  }
  data <- attr(results,'data')
  pi1 <- data$pi1
  pi0 <- data$pi0
  no.of <- length(data)
  classifier.names <- rownames(results$metrics)
  
  # ROC CURVES
  plot.types <- c(1:4)
  if (!any(which %in% plot.types)){
    stop('Type of plot not recognised. Use a numeric value between 1 and 4.')
  }
  colorlist <-  c('red','blue','green','magenta','yellow','forestgreen')
    
  if (which == 1){ # ROC

    legend.colors <- array(data=NA, no.of)
    legend.lty <- array(data=NA, dim=no.of)
    legend.names <- array(data=NA, dim=no.of)

    # Plot the diagonal
    plot(c(0,1), c(0,1), type='l', lty=4, col='black',
         xlab='FPR (1-F0)', ylab='TPR (1-F1)',
         main= 'ROC (continuous) and ROCH (dotted)')

    # Define legend entry
    legend.lty[1] <- 4
    legend.colors[1] <- 'black'
    legend.names[1] <- 'trivial'

    for (count in (1:no.of)){
      data.now <- data[[count]]  
     
      # Plot the ROC
      lines(1-data.now$F0,1-data.now$F1, type='l', lty=1, col=colorlist[count])
      
      # Plot the convex hull
      lines(data.now$G0,data.now$G1, type='l', lty=3, col=colorlist[count])
      
      # Define legend entry
      legend.lty[count+1] <- 1
      legend.colors[count+1] <- colorlist[count]
      legend.names[count+1] <- classifier.names[count]
      
    }
    # put in legend
    legend('bottomright', legend=legend.names,
           lty=legend.lty, col=legend.colors)
  }
  
  if (which==2){ # H cost weights

    legend.colors <- array(data=NA, no.of)
    legend.lty <- array(data=NA, dim=no.of)
    legend.names <- array(data=NA, dim=no.of)

    # prior parameters are the same for all classifiers
    severity.ratio <- data[[1]]$severity.ratio
    
    if (is.na(severity.ratio)){
      severity.ratio <- pi1/pi0
    }
    if (severity.ratio > 0){
      shape1 <- 2
      shape2 <- 1+(shape1-1)*1/severity.ratio
    }
    if (severity.ratio < 0){
      shape1 <- pi1+1
      shape2 <- pi0+1
    }  
    cost.parameter <- severity.ratio*((1+severity.ratio)^(-1))
    
    b <- c(1:100)/100
    y <- dbeta(b,shape1,shape2)
    plot(b, y, type='l', xlab= 'Normalised Cost = SR/(1+SR)', ylab='w(c) ',
         main='H measure w(c)'
         )
    
    lines(rep(cost.parameter,100),c(1:100*max(y))/100,lty=3)
    # put in legend
    legend('topright', lty=c(1,3), col=c('black','black'),
      legend=c('All classifiers', 
      paste('SR = ', format(severity.ratio,digits=3), sep=' ')
      )
    )

  }
  
  if (which==3){ # AUC cost weights
    
    # initialise legend specifications
    legend.colors <- array(data=NA, no.of)
    legend.lty <- array(data=NA, dim=no.of)
    legend.names <- array(data=NA, dim=no.of)

    # iterate over classifiers
    for (count in (1:no.of)){
      temp <- data[[count]]  
      
      aucd <- c((temp$n0*temp$G0 + temp$n1*temp$G1),1)
      aucd2 <- c(1, (temp$n0*temp$G0 + temp$n1*temp$G1))
      aucf <- (aucd-aucd2)/temp$n
  
      if (count == 1){
        plot(temp$cost[2:temp$hc], aucf[2:temp$hc], type='h', 
             xlim=c(0,1), ylim=c(0,1), lwd=2,
             main='AUC w(c)',
             xlab='Cost', ylab='w(c)', col=colorlist[count])
      } else {
        lines(temp$cost[2:temp$hc], aucf[2:temp$hc], type='h', 
             xlim=c(0,1), ylim=c(0,1), lwd=2,
             main='AUC measure weight function of c',
             xlab='Cost', ylab='w(c)', col=colorlist[count])
      } 
    
      # Define legend entry
      legend.lty[count] <- 1
      legend.colors[count] <- colorlist[count]
      legend.names[count] <- classifier.names[count]
    }
    # put in legend
    legend('top', legend=legend.names,
           lty=legend.lty, col=legend.colors)
  }
  
  if (which == 4){ # scoring densities
  
    legend.colors <- array(data=NA, no.of)
    legend.lty <- array(data=NA, dim=no.of)
    legend.names <- array(data=NA, dim=no.of)
    
    
    
    # First compute densities
    ymax <- 0
    to.plot0 <- list()
    to.plot1 <- list()
    for (count in (1:no.of)){
      data.now <- data[[count]]  
      s.min <- min(data.now$s.class0, data.now$s.class1)
      s.max <- max(data.now$s.class0, data.now$s.class1)
            
      to.plot0[[count]] <- density(data.now$s.class0,bw=bw)
      to.plot1[[count]] <- density(data.now$s.class1,bw=bw)
      
      # compute maximum density value for plotting purposes
      ymax <- max(to.plot0[[count]]$y, to.plot1[[count]]$y, ymax)
    }

    # Now plot
    for (count in (1:no.of)){
      
      if (count == 1){
        plot(to.plot0[[count]], xlim=c(s.min,s.max), col=colorlist[count], lty=4,
             	main='Smoothed score distributions \n (class 0: dash-dotted, class 1: dashed) ', 
		xlab= 'Score ', ylim=c(0,ymax))
      } else {
        lines(to.plot0[[count]], xlim=c(s.min,s.max), col=colorlist[count], lty=4, ylim=c(0,ymax))
      }
      lines(to.plot1[[count]], lty=2, col=colorlist[count], ylim=c(0,ymax))
      
      
      # Define legend entry
      legend.lty[count] <- 1
      legend.colors[count] <- colorlist[count]
      legend.names[count] <- classifier.names[count]
        
    }
    # put in legend
    legend('top', legend=legend.names,
           lty=legend.lty, col=legend.colors)

  }
  
}


relabel <- function(labels){


  if (length(levels(as.factor(labels)))==1){
    stop('Only one class is present in the dataset. Need both classes to be represented.')
  }
  if (length(levels(as.factor(labels)))>2){
    stop('More than two classes present, but code can only handle binary classification.')
  }
  


   labels <- as.factor(as.character(labels))
   input.labels <- levels(labels)
   

   cond.temp <- (	identical(input.labels,c('case','non-case')) |
   					identical(input.labels,c('Case','Non-case')) |
   					identical(input.labels,c('case','noncase'))  |
   					identical(input.labels,c('Case','Non-case')) )

   if (cond.temp) {
      levels(labels) <- c('1', '0')
      message('Class labels have been switched from (',
        paste(input.labels[1],input.labels[2], sep=','), ') to (', 
        paste('1', '0', sep=','), ')')   
        labels <- as.factor(labels)
	    labels <- 2-as.numeric(labels) # turn into numeric array of 0s and 1s
   } else {
      levels(labels) <- c('0', '1')
	  if (!(identical(input.labels,c('0', '1')))){
	      message('Class labels have been switched from (',
    	    paste(input.labels[1],input.labels[2], sep=','), ') to (', 
        	paste('0', '1', sep=','), ')')
        	}
      labels <- as.factor(labels)
	  labels <- as.numeric(labels)-1 # turn into numeric array of 0s and 1s

    } 

    return(labels)	
}


misclassCounts <- function(predicted.class,true.class){
 
	true.class <- as.array(true.class)
	predicted.class <- as.array(predicted.class)
	
	# make sure the same convention is employed for both true and predicted
#    check <- relabel(c(true.class,predicted.class))
#	l <- length(check)
#	true.class <- check[1:(l/2)]
#	predicted.class <- check[(l/2+1):l]
	

    TP <- sum(predicted.class == 1 & true.class == 1)
    FP <- sum(predicted.class == 1 & true.class == 0)
    TN <- sum(predicted.class == 0 & true.class == 0)
    FN <- sum(predicted.class == 0 & true.class == 1)

    conf.matrix <- data.frame(pred.1=c(TP,FP),pred.0=c(FN,TN))
    row.names(conf.matrix) <- c('actual.1','actual.0')

    ER <- (FP + FN)/(TP+FP+TN+FN)

    Sens <- TP/(TP+FN)
    Spec <- TN/(TN+FP)

    Precision <- TP/(TP+FP)
    Recall <- Sens 

    TPR <- Recall
    FPR <- 1-Spec

    F <- 2/(1/Precision+1/Sens)
    Youden <- Sens + Spec -1
    metrics <- data.frame(ER=ER,
  		Sens=Sens,Spec=Spec,Precision=Precision,
		Recall=Recall, TPR=TPR, FPR=FPR, F=F, Youden=Youden)
    return(list(conf.matrix=conf.matrix,metrics=metrics))
}

