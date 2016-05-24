#' Greedy algorithm for prototype selection
#'
#' Selects prototypes for each class in a greedy manner as described in 'Bien and
#' Tibshirani (2011) Prototype Selection for Interpretable Classification. Annals
#' of Applied Statistics. 5(4). 2403-2424.'
#'
#' @param x n by p matrix of training features (optional, see \code{dxz}).
#' @param y n-vector of labels of the training data.
#' @param z set of potential prototypes (optional, see \code{dxz}).
#' @param dxz instead of \code{x} and \code{z}, you can give \code{dxz}, the matrix of
#' pairwise dissimilarities between \code{x} and \code{z}, with ij-th element giving the
#' dissimilarity between training point x_i and prototype-candidate z_j.
#' @param eps size of covering balls.
#' @param lambda cost of adding a prototype.
#' @details
#' It's more efficient to compute dxz just once on your own rather than have protoclass
#' repeatedly compute the pairwise distances on each call.
#' @return
#' An object of class "protoclass," which has the following elements:
#' \itemize{
#'  \item{alpha: }{Matrix of dimensions nrow(z)-by-nclass.  alpha[j,k] indicates whether jth potential prototype has been chosen as a prototype for class k.}
#'  \item{classes: }{Names of classes}
#'  \item{proto.order: }{The sequence of prototypes that were selected.}
#'  \item{ncovered: }{nproto-by-nclass matrix with ncovered[j,k] giving the number of class k training points covered by the jth prototype's ball.}
#' \item{coverlist: }{n-by-nclass matrix with row i giving number of each type of prototype covering point i.}
#' \item{uncovered: }{Indicates whether a training point is not covered by a prototype of its own class.}
#' \item{wrongcover: }{Number of prototypes from other classes covering each training point.}
#' \item{nproto: }{nclass-vector giving the number of prototypes in each class.}
#' }
#' @export
#' @seealso \code{\link{predict.protoclass}}
#' @examples
#' # generate some data:
#' set.seed(1)
#' n <- 200
#' p <- 2
#' x <- matrix(rnorm(n * p), n, p)
#' y <- rep(c("A","B"), each=n/2)
#' x[y=="A", ] <- x[y=="A", ] + 3
#' itr <- sample(n, n/2)
#' xtr <- x[itr, ] # train
#' ytr <- y[itr]
#' xte <- x[-itr, ] # test
#' yte <- y[-itr]
#' 
#' # take prototype candidates identical to training points:
#' z <- xtr
#' dxz <- dist2(xtr, z)
#' # run protoclass:
#' prot <- protoclass(dxz=dxz, y=ytr, eps=2, lambda=1/n)
#' \dontrun{
#' plot(prot,xtr,y=1+(ytr=="A"))
#' }
#' # get predictions on test data:
#' pred1 <- predict(prot, xte, z=xtr)
#' # get predictions on test data using pairwise distances:
#' pred2 <- predictwithd.protoclass(prot, dist2(xte, z))
protoclass <- function(x, y, z, dxz, eps, lambda=1/n)
  {
    if(missing(y)) stop("Missing y")
    if(missing(eps)) stop("Missing epsilon")
    if(missing(dxz))
      {
        if(missing(x)||missing(z)) stop("Must input either dxz or x and z")
        if(nrow(x)==ncol(x) & all(diag(x)==0))
          stop("May be mixing up dxz and x: Use dxz= to be safe")
        x = scale(x)
        cent = attr(x,"scaled:center")
        scale = attr(x,"scaled:scale")
        z = scale(z,center=cent,scale=scale)
        dxz = dist2(x,z)
        print("Using Standardized Euclidean metric.")
      }
    if(!is.matrix(dxz)) stop("Dissimilarites must be a matrix")
    n = nrow(dxz)
    if(length(y)!=n) stop("Dimensions of d and y don't match")
    prot <- greedy(dxz,y,eps,lambda=lambda)
    class(prot) <- "protoclass"
    prot$call <- match.call()
    prot
  }

#' Prints an object of class protoclass
#'
#' @param x object of class protoclass.
#' @param ... ignore this argument.
#' @S3method print protoclass
#' @method print protoclass
print.protoclass <- function(x, ...) {
  if (!is.null(x$call))
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Number of classes:", ncol(x$alpha), "\n")
  cat("Ball size        :", x$eps, "\n")
  cat("-----\n")
  mat <- matrix(x$nproto,byrow=TRUE,nrow=1)
  dimnames(mat) <- list(c("Number of Prototypes:"),x$classes)
  print(mat, quote=FALSE)
  cat("\n")
  invisible(x)
}

#' Compute distances between two sets of points
#'
#' Given (n1)-by-p and (n2)-by-p matrices, computes (n1)-by-(n2) matrix of Euclidean distances
#'
#' @param x a (n1)-by-p matrix
#' @param xx a (n2)-by-p matrix
#' @export
dist2 <- function(x,xx)
  {
    if(nrow(x) < nrow(xx))
      {
        return(t(dist2(xx,x)))
      }
    
    d = apply(xx,1,function(xxrow) rowSums(sweep(x,2,xxrow)^2))
    sqrt(d)
  }

#' Make predictions based on prototypes
#'
#' Performs nearest-prototype classification on newx.
#'
#' @param object output of protoclass function
#' @param newx features from some new data that we want to make predictions on.
#' @param z potential prototype positions z (often taken to be x, the training data)
#' @param ... ignore this argument.
#' @seealso \code{\link{protoclass}} \code{\link{predictwithd.protoclass}}
#' @S3method predict protoclass
#' @method predict protoclass
predict.protoclass <- function(object, newx, z, ...){
  alpha <- object$alpha
  n=nrow(alpha)
  nclass=ncol(alpha)
  classes = object$classes
  
  require(class)

  if(all(colSums(alpha>0) == 0)){pred=rep(NA,nrow(newx));return(pred)}
  
  prototypes=NULL
  protoclass=NULL
  for(k in 1:nclass)
    {
      oo = which(alpha[,k]>0)
      prototypes = rbind(prototypes,z[oo,])
      protoclass = c(protoclass,rep(classes[k],length(oo)))
    }

  pred=knn1(prototypes,newx,protoclass)
  
  return(pred)
}


#' Make predictions based on prototypes
#'
#' Performs nearest-prototype classification like predict.protoclass, except only requires
#' the ntest-by-m matrix of distances between test features and z rather than requiring
#' these matrices themselves. (Using this saves repeated computation of the pairwise distances...)
#'
#' @param object output of protoclass function
#' @param d.test.z ntest-by-m matrix of distances between test features and z
#' @seealso \code{\link{predict.protoclass}}
#' @export
predictwithd.protoclass <- function(object,d.test.z)
  {
    # get the indices of the prototypes
    whichProto = which(rowSums(object$alpha)>0)
    if(length(whichProto)==0) return (rep(NA,nrow(d.test.z)))
    d.Proto = as.matrix(d.test.z[,whichProto])
        
    # for every test point, find the closest prototype:
    closestProto = apply(d.Proto,1,which.min)

    alphaClosest = object$alpha[whichProto[closestProto],]
    # get max and randomly tie-break:
    maxrand <- function(alp)
      {
        w = which(alp == max(alp))
        if(length(w)>1)
          w = sample(w,1)
        
        w
      }
    as.factor(object$classes)[apply(alphaClosest,1,maxrand)]
  }

greedy <- function(dxz, y, eps, lambda=0)
  {
    classes = unique(y)
    nproto = ncol(dxz)
    K = length(classes)
    n = nrow(dxz)
    if(length(y)!=n) stop("y must be of length n.")

    C = NULL # n by K matrix indicating class of each point
    for(cl in classes)
      C = cbind(C,y==cl)
    B = dxz < eps # each col indicates the points covered by this prototype
    
    covered = rep(0,n)# keeps track of points that have been *correctly* covered
    # W_ik = 1{i is of class k and uncovered} - 1{i not of class k}
    W = 2*C-1
    scores = t(B)%*%W
    
    protos=NULL
    ncovered=NULL
    alpha = matrix(0,nproto,K)
    i=0
    while(TRUE)
      {
        i=i+1
        iimax = which.max(scores)
        # convert via:  element (i,j) of a matrix coded as i+(j-1)*nrows
        kmax = ceiling(iimax/nproto)
        pmax = iimax - nproto*(kmax-1) # adding prototype pmax to class kmax's list is best
        
        if(scores[pmax,kmax] > lambda)
          {
            # adding this prototype increases objective by the most, so accept it:
            protos = c(protos,pmax)
            alpha[pmax,kmax] = 1

            ## update scores:

            # identify points that are no longer uncovered:
            justcovered = B[,pmax] & C[,kmax] & !covered
            covered = covered + justcovered
            # determine which (potential) prototypes' scores were affected:
            affectedProtos = which(colSums(matrix(B[justcovered,],ncol=nproto))>0)
            ncovered = rbind(ncovered,
              c(sum(justcovered),sum(B[,pmax] & !C[,kmax])))

            # for each affected prototype, determine # of points
            # that it can cover that no longer need to be covered
            scoreReduction = sapply(affectedProtos,
              function(p) sum(B[,p]&justcovered))
            scores[affectedProtos,kmax] = scores[affectedProtos,kmax] - scoreReduction
          }
        else
          {
            # none of the remaining prototype-candidates increases objective by enough
            break
          }
      }

    aa = list(alpha=alpha,classes=classes,eps=eps,lambda=lambda,
      proto.order=protos,ncovered=ncovered)
    an = analyzeSolution(aa,dxz,y)
    # add items of an to aa
    for(i in seq(length(an)))
      aa[[names(an)[i]]] = an[[i]]
    
    aa
  }

analyzeSolution <- function(aa,dxz,y)
  {
    if(is.null(aa$lambda))
      lambda=0
    else
      lambda=aa$lambda
    
    # j-th col indicates the points covered by prototype z_j.
    B = dxz < aa$eps
    n = nrow(dxz)
    
    alpha = aa$alpha > 0
    coverlist = B %*% alpha # row i gives how many of each type of prototype covers point i
    uncovered = rep(0,n)
    wrongcover = rep(0,n)
    wrongcover.max = rep(0,n)
    for(k in seq(length(aa$classes)))
      {
        wrongcover[y==aa$classes[k]] = rowSums(as.matrix(coverlist[y==aa$classes[k],-k]))
        wrongcover.max[y==aa$classes[k]] = apply(as.matrix(coverlist[y==aa$classes[k],-k]),1,max)
        uncovered[y==aa$classes[k]] = coverlist[y==aa$classes[k],k] == 0
      }
    # in rowe5c, slack incurred if either of these are non-zero
    #psi = pmax(uncovered,wrongcover.max)

    list(coverlist=coverlist,
         #psi=psi,
         uncovered=uncovered,
         wrongcover=wrongcover,
         nproto = colSums(alpha))#,
         #obj.rowe5c = sum(psi)+lambda*sum(alpha),
         #obj.rowe6 = sum(uncovered)+sum(wrongcover)+lambda*sum(alpha))
  }

#' Plots prototypes with data
#'
#' This function is only for two-dimensional data.
#' @param x output of protoclass function
#' @param xtr n by p matrix of training features.
#' @param y n-vector of labels of the training data.
#' @param z set of potential prototypes (only required if protoclass was not run with z=x).
#' @param circleFLAG indicates whether to draw circles around each prototype.
#' @param boundFLAG indicates whether to draw classification boundary.  This can be computationally expensive.
#' @param res resolution of the classification boundary.  Computation increases with res.  Try 100.
#' @param main optional main title for plot.
#' @param ... other arguments to pass to plot.
#' @seealso \code{\link{protoclass}}
#' @S3method plot protoclass
#' @method plot protoclass
plot.protoclass <- function(x,xtr,y,z=NULL,circleFLAG=TRUE, boundFLAG=FALSE,res=NULL,main=NULL,...)
  {
    aa <- x
    x <- xtr
    eps <- aa$eps
    # make sure arguments are correct
    if(is.null(z))
      {
        # make sure z not needed
        if(is.null(aa$generalz) || !aa$generalz)
          z=x
        else
          stop("z needed")
      }
    else
      {
        # make sure that z is needed (or if not, at least =x)
        if(is.null(aa$generalz) || !aa$generalz)
          if(any(dim(x)!= dim(z))|| any(x!=z)) stop("aa used z=x, yet different z given.")
      }
    
    if(nrow(aa$alpha)!=nrow(z)) stop("Wrong z for this aa.")
    if(ncol(x)!=2) stop("Only plots 2 dimensional data.")
    if(ncol(z)!=2) stop("z must have 2 columns.")
    if(is.null(main))
      main = paste("eps =",as.character(round(eps,2)))
    
    # start the plotting

    mins = apply(cbind(x,z),2,min)
    maxs = apply(cbind(x,z),2,max)
    width = max(maxs-mins)
    
    plot(x,
         xlim = c(mins[1],mins[1]+width),
         ylim = c(mins[2],mins[2]+width),
         col=y+1,main=main,axes=F,xlab="",ylab="",...)
    box()

    if(all(rowSums(aa$alpha)==0)) return() # no prototypes case
    
    # draw the prototypes (with tint indicating "strength" of alpha)
    nclass = ncol(aa$alp)
    for(k in 1:nclass)
      {
        sol = aa$alpha[,k]>0
        if(sum(sol)==0)
          next
        
        # size = 1.2/ num prototypes so far at this spot
        size = 1.2 / rowSums(aa$alpha[sol,1:k,drop=F]>0)

        color = rgb2hsv(col2rgb(k+1))
        s = color[2]*aa$alpha[sol,k]
        for(i in seq(length(s)))
          if(s[i]>1)
            s[i]=1
        points(matrix(z[sol,],nrow=sum(sol)),
               col=hsv(h=color[1],s=s,v=color[3]),
               pch=20,
               cex=size)
        points(matrix(z[sol,],nrow=sum(sol)),
               col=1,
               pch=1,
               cex=size)
      }

    if(circleFLAG)
      {
        for(i in which(rowSums(aa$alp)>0))
          {
            # i is the index of a centroid for some class
            addcircle(z[i,1],z[i,2],eps)
          }
      }
    if(boundFLAG)
      {
        isproto = rowSums(aa$alpha)>0
        protox = z[isproto,]
        # which.max below used to pick the element = 1,tie breaking in some way
        protoy = aa$classes[apply(aa$alpha[isproto,],1,which.max)]
        predictPVM <- function(xx) protoy[which.min(rowSums(sweep(protox,2,xx)^2))]
        if(is.null(res)) res = 40
        addclassifierboundary(NULL,predictPVM,lty="solid",res=res)
      }
  }


addcircle <- function(x,y,r)
  {
    tt=0:50/50*2*pi
    lines(x+r*cos(tt),y+r*sin(tt))
  }

# classify is a function that takes a query point and returns a label
addclassifierboundary <- function(x, classify, res=40,...)
  {
    x1 = seq(par("usr")[1],par("usr")[2],length = res)
    x2 = seq(par("usr")[3],par("usr")[4],length = res)
    xgrid = cbind(rep(x1,each=length(x2)),x2)

    predictions = apply(xgrid,1,classify)
    # draw boundary
    for(class in unique(predictions))
      {
        ind = (predictions == class)
        ind = matrix(ind,ncol=res,byrow=TRUE)
        contour(x1, x2, ind, levels = .5, drawlabels=FALSE, add=TRUE,...)
      }
  }
