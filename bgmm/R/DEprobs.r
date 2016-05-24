
DEprobs <- function(model, verbose=FALSE){
  k=model$k
  if (!k%in%c(2,3))
    stop("DEprobs: Wrong number of fitted model components. Only support for 2 or 3-component models is provided")

  ## Get the right component to correspond to the differentially expressed genes -- which looks like the differential one according to the posterior probabilities
  ## For models of size 2 -- should be the one with a broader range (encompassing the other, or the one with higher maximum deviation from 0 (we assume the data are centered around 0 and are both positive and negative)

  if (!is.null(model$knowns)){
    preds=predict(model, model$X, knowns=model$knowns, B=model$B)
	}
else
  preds=predict(model, model$X)

  if (k==2)
    resp=check2res(model$X, model$knowns, preds, verbose)
  else
    resp=check3res(model$X, model$knowns, preds, model$mu, verbose)
#browser()
  return(resp)
}


check2res <- function(X, knowns, preds,verbose){
	d=X
	class=preds$class.X

	if (!is.null(knowns)){
		d=rbind(d,knowns)
		class=c(class, preds$class.knowns)
		}
    #identify the differential component
  d1 <- d[class==1]
  d2 <- d[class==2]
  r1 <- range(d1)
  r2 <- range(d2)
  m1=max(abs(r1))
  m2=max(abs(r2))
  diff=2
  if (r1[1]<r2[1]){ ##1 sticks out on the left
    if(r1[2]>r2[2]){ ##1 sticks out on the right
      diff=1 
    }else{
      if (m1>m2){ ##1 sticks out more from 0 than the 2
        diff=1
        warning(paste("DEprobs: accoding to the data, the differentially expressed component cannot be identified as the one which is broader than the unchanged on both sides. Choosing the component with a larger dev from 0"))
      }
    }
  }## diff either encompasses or has higher dev from 0
  
  if (verbose)
    print(paste("The differential component number:",diff))
    
  p.X <- preds$tij.X[,diff]
  downs.X <- (X<0)
  p.X[downs.X] <- p.X[downs.X]*(-1)
  
  p.knowns=NULL
  if(!is.null(knowns)){
 	p.knowns <- preds$tij.knowns[,diff]
  	downs.knowns <- (knowns<0)
  	p.knowns[downs.knowns] <- p.knowns[downs.knowns]*(-1)
	}
  return(list(diff.p.X=p.X, diff.p.knowns=p.knowns, diff.c=c(diff)))
}


check3res <- function(X, knowns, preds, mu, verbose){
  ##order by means
  ord=order(mu)
 if (verbose){
    print(paste("the downregulated component number:", ord[1]))
    print(paste("the upregulated component number:", ord[3]))
  }

  tij.X=preds$tij.X[,ord]
  p.X <- -1*tij.X[,1] ##probability of down-regulation
  ups.X <- (X>0) 
  p.X[ups.X] <- tij.X[ups.X,3]

  
  p.knowns=NULL
  if(!is.null(knowns)){
	tij.knowns=preds$tij.knowns[,ord]
  	p.knowns <- -1*tij.knowns[,1] ##probability of down-regulation
  	ups.knowns <- (knowns>0) 
  	p.knowns[ups.knowns] <- tij.knowns[ups.knowns,3]
  }
  return(list(diff.p.X=p.X, diff.p.knowns=p.knowns, diff.c=c(ord[1],ord[3])))
}
