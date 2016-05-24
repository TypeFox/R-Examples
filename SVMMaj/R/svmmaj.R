svmmaj <- function(X,y,lambda=1,...)  UseMethod('svmmaj')

#==========================================================================
#==========================================================================
#==========================================================================
print.svmmaj <- function(x,...)
    with(x,cat('Model:\n',
        '   method                            ',method$type,'\n',
        '   number of iterations              ',iteration,'\n',
        '   loss value                        ',loss,'\n',
        '   number of support vectors         ',nSV,'\n'
    ))       

#==========================================================================
#==========================================================================
#==========================================================================


summary.svmmaj <- function(object,...)
  structure(object,class='summary.svmmaj')
  
print.summary.svmmaj <- function(x,...) {
    object <- x
    cat('Call:\n   ')
    print(object$call)
   
   splineb <- 'no'
    if(!with(object,splineKnots==0 && splineDegree==1)) splineb <- with(object,
          paste('knots ='  ,splineKnots , ',degree =', splineDegree))
         
    typek  <- attr(object$method$kernel,'class')[1]
    if(typek == 'vanillakernel')
      typek <- 'linear'
       
    cat('\nSettings:\n',
        '   lambda                            ',object$lambda,'\n',
        '   hinge error                       ',attr(object$hinge,'hinge'),'\n',
        '   spline basis                      ',splineb,'\n',
        '   type of kernel                    ',typek,'\n')
              
    if(typek!='linear'){
      cat('    parameters of kernel               ') 
      kernel.param <- attr(object$method$kernel,'kpar')     
      sapply(1:length(kernel.param),function(i) 
          cat(names(kernel.param)[i],'=', kernel.param[[i]],' ')
      )}
    cat('\nData:\n',
        '   class labels                      ',object$classes,'\n',
        '   rank of X                         ',object$method$x,'\n',
        '   number of predictor variables     ',NCOL(object$data)-1,'\n',
        '   number of objects                 ',NROW(object$data),'\n',
        '   omitted objects                   ',length(object$na.output),'\n\n'
    )
    cat('Model:\n',
        '   method                            ',object$method$type,'\n',
        '   number of iterations              ',object$iteration,'\n',
        '   loss value                        ',object$loss,'\n',
        '   number of support vectors         ',object$nSV,'\n'
    )   
    prints <- with(object,classification(q,na.omit(data)[,1],classes=classes,weights=weights.obs))
     
    cat('\nConfusion matrix:\n')
    print(prints$matrix)  
    cat('\nClassification Measures:\n')
    print(prints$overall,digits=3)      
    cat('\n')                          
    print(prints$measures,digits=3)
}
#==========================================================================
#==========================================================================
#==========================================================================

predict.svmmaj <- function(object,X.new,y=NULL,show.plot=FALSE,...)
  {
    k     <- NCOL(object$data) -1
    if(!is.data.frame(X.new)) {
      if(!is.matrix(X.new)) X.new  <- matrix(X.new,ncol=k)
      X.new <- data.frame(X=X.new)
    }
    #TRANSFORM DATA X, IF NEEDED    
	Xnew <- mapply(predict.transDat,X.new,attrib=object$propData,SIMPLIFY=FALSE)
	    if(length(object$weights.var)==k)
		    Xnew <- mapply(`*`,Xnew,object$weights.var)
		    
	Xnew <- data.matrix(data.frame(Xnew))
    
    #PREDICT THE VALUES OF THE PREDICTED VALUE INCLUDING THE INTERCEPT (Q-TILDE)  
    q <- structure(qhat.theta(object$method,object$theta,Xnew),
          class='q.svmmaj', y=y, classes=object$classes)
    attr(q,"yhat") <- factor(q>0 , levels=c(FALSE,TRUE), labels=object$classes)

    if(show.plot & !is.null(y))
        plot.svmmaj(list(q=q,y=y,classes=object$classes))
    else if(show.plot & is.null(y)) 
        plot(density(q),y=NULL,
            main='Density of predicted values',xlab='q',ylab='density') 
    return(q)  
  }

#==========================================================================
#==========================================================================
#==========================================================================
print.q.svmmaj <- function(x,...) {
q <- x
    classnames   <- attr(q,'classes')
    classcount  <- matrix(c( sum(q<0) , sum(q>=0) ),nrow=1,
        dimnames= list(list('frequency'),as.list(classnames)))

    cat('Prediction frequencies:\n')
    print(classcount)
    
    if(!is.null(attr(q,'y'))){
       prints <- classification(q,attr(q,'y'),classes=classnames)
  
      cat('\nConfusion matrix:\n')
      print(prints$matrix)
      cat('\nClassification Measures:\n')
      print(prints$overall,digits=3)
      cat('\n')                                    
      print(prints$measures,digits=3)
    }
}

#==========================================================================
#==========================================================================
#==========================================================================
classification <- function(q,y,classes=c('-1','1'),weights=NULL) {
  #Compute confusion matrix
  tab  <- t(cbind(y==classes[[1]],y==classes[[2]])) %*% cbind(q<=0,q>0)

  tabT        <- array(0,dim=c(3,3))
  tabT[-3,-3] <- tab
  tabT[3,]    <- colSums(tabT)
  tabT[,3]    <- rowSums(tabT)
  dimnames(tabT) <- list('Observed (y)'   =c(classes,'Total'),
                         'Predicted(yhat)'=c(classes,'Total'))

  #Compute the overall rates
  hr <- sum(diag(tab)/sum(tab))
  if(is.null(weights)) {
    overall <- matrix(c(hr,1-hr),ncol=1)
    dimnames(overall) <- list( list('   hit rate                          ',
                                    '   misclassification rate            '),
                               list(''))
  } else {
    hr.w <- (q<=0) == (y==classes[[1]])
    hr.w <- sum(hr.w*weights)/sum(weights)
    overall <- matrix(c(hr,hr.w,1-hr,1-hr.w),ncol=1)
    dimnames(overall) <- list( list('   hit rate                          ',
                                    '   weighted hit rate                 ',
                                    '   misclassification rate            ',
                                    '   weighted missclassification rate  '),
                               list(''))
  }
  #Compute classification measures per class
  positive <- diag(tab)/rowSums(tab)
  measures  <- cbind('       TP'=positive,
                     '       FP'=1-positive,
                     'Precision'=diag(tab)/colSums(tab))
  rownames(measures) <- paste('   ',classes)

  return(list(matrix=tabT,overall=overall,measures=measures))
}

#==========================================================================
#==========================================================================
#==========================================================================
plot.svmmaj <- function(x,titletext=NULL,...)
{
  object <- x
  classnames <- object$classes
  if(all(sort(unique(object$y))!=c(-1,1)))
      object$y <- sign((object$y==classnames[[2]])-.5)
  dens1 <- density(object$q[object$y>0],bw=bw.nrd0(object$q)) ;#dens1$y <- dens1$y / sum(object$y==1)
  dens2 <- density(object$q[object$y<0],bw=bw.nrd0(object$q));#dens2$y <- dens2$y / sum(object$y==-1)
  
  dens1$y <- dens1$y  * length(object$q[object$y>0])
  dens2$y <- dens2$y  * length(object$q[object$y<0])
  
  q <- c(dens1$x,dens2$x)
  n <- c(dens1$y,dens2$y)*1.3
  sub <- paste('Correctly predicted objects: ',
                round(mean(sign(object$q)==object$y)*100,2),'%',sep='')
  
  plot(q,n,'n',main=c('Density of predicted values',sub),ylab='density',
  			xlab=paste('q\n (bandwidth = ',round(bw.nrd0(object$q),4),')'))
  title(titletext)
  lines(dens1,col='blue',lty=1)
  lines(dens2,col='red',lty=1)
  lines( rep(0,2) , c(0,par()$usr[4]),col='black',lty=3)
  classlabels <- c( paste('Class (q<0)',classnames[1]), paste('Class (q>=0)',classnames[2]))
  legend("top",classlabels,horiz=TRUE,lty=c(1,1),col=c('red','blue'),bg = 'white')
}
#==========================================================================
#==========================================================================
#==========================================================================

plotWeights <- function(object, plotdim = c(3,3)){
	nplot <- prod(plotdim)
	
	
	if(!object$method$linear) stop('nonlinear kernel is used')
	Xold  			<- data.matrix(object$data[-1])
	splineLength 	<- sapply(object$propData,`[[`,'dim')[2,]
	beta			<- object$beta
	k				<- NCOL(object$data)-1
	
	#GENERATE GRIDPOINTS        
	grid <- matrix(rep( seq(0,1,length.out=100),k),ncol=k)
	grid <- t(t(grid)*as.vector(diff(apply(Xold,2,range))) +
	                 apply(Xold,2,min))
	grid <- data.frame(grid)

	X <- mapply(predict.transDat,grid,attrib=object$propData,SIMPLIFY=FALSE)
    
	l <-  sapply(object$propData,`[[`,'dim')[2,]  #dim...
	at <- unlist(mapply(rep,1:k,l))
	bspl <- lapply(1:k,function(i) beta[at==i])
	
	
	Splinesplus <- Splinesmin <- Splinesum <- NULL
  for(j in 1:ceiling(k/ nplot)){
		#DETERMINE PLOT LAYOUT
		mini <- nplot*(j-1)+1
		maxi <- min(k,j*nplot)
		noi  <- maxi - mini + 1
		plotlayout <- rep(0,nplot)
		plotlayout[1:noi] <- 1:noi
		plotlayout <- rbind(matrix(plotlayout,byrow=TRUE,nrow=plotdim[1],ncol=plotdim[2]),noi+1)
		layout(plotlayout,heights = c(rep(1,plotdim[1]),.2))
		
		#PLOT PER PAGE
		par(ask=j>1,mar=c(2.1,4.1,3.1,2.1),mgp=c(2, 1, 0))
		Splinesplus <- mapply(`%*%`,X,lapply(bspl,pmax,0))
		Splinesmin  <- mapply(`%*%`,X,lapply(bspl,pmin,0))
		Splinessum  <- mapply(`%*%`,X,bspl)
		ylim        <- range(c(Splinesplus,Splinesmin,Splinessum))
		
		if(is.matrix(Splinessum)){
		  Splinesplus <- data.frame(Splinesplus)
		  Splinesmin  <- data.frame(Splinesmin)
		  Splinessum  <- data.frame(Splinessum)
    }
		for(i in mini:maxi){
			if(!is.null(object$propData[[i]]$values)){
			#FACTOR
				barplot(bspl[[i]],names.arg=object$propData[[i]]$values,
				ylim=ylim, main=colnames(Xold)[i],xlab='',ylab='treatment')
			} else if(!is.null(object$propData[[i]]$type)){
			#LOGICAL
				barplot(c(0,bspl[[i]]),names.arg=c('FALSE','TRUE'),
				ylim=ylim, main=colnames(Xold)[i],xlab='',ylab='treatment')
			} else {
			#NUMERIC / SPLINES
				plot( grid[,i] , Splinessum[[i]]  , col='black',type='l',ylab='',
					ylim=ylim,
					main=colnames(Xold)[i],xlab='')
                         
				if(!is.null(object$propData[[i]]$splineDegree)){
					title(ylab='I-Spline')
					lines( grid[,i] , Splinesplus[[i]], col='blue' ,type='l',lty=3)
					lines( grid[,i] , Splinesmin[[i]] , col='red'  ,type='l',lty=2)
				} else title(ylab = 'x*b')
			}
		}
		par(mar=c(0,0,0,0))
		plot(0,0,type='n',axes=FALSE,bty='n',pty='m',xlab='',ylab='')
		legend('center',c('positive weights','negative weights','Total'),
			horiz=TRUE,
			lty=c(3,2,1),
			col=c('blue','red','black'),
			bg = 'white',
			bty= 'o',
			xpd=T)
		}
	par(mfrow=c(1,1),mgp=c(3,1,0),mar=c(4,5,5,2)+.1)
}
                   

