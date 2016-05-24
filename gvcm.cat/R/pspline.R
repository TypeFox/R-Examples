sp <-
function(x, knots=20, n="L2"){

  a.spline <- function (k) 
      { 
        if (k>1)
        {
          if (k > 2)
          {
            m1 <- cbind(diag(k-2),0,0)
            m2 <- cbind(0,diag(-2,k-2),0)
            m3 <- cbind(0,0,diag(k-2))
            mat <- t(m1 + m2 + m3)
          } else { mat <- matrix(ncol=0,nrow=2) }
        } else {
        mat <- matrix(ncol=0,nrow=1)
        }
        
        return (t(mat))
      }

    if (!is.numeric(x))
      stop("x in 'pspline(x)' must be numeric.")
    
    # definitions
    mini <- min(x)
    maxi <- max(x)
    ki <- knots - 4
    dis <- (maxi - mini)/(knots - 7) # 5, by = ((to - from)/(length.out - 1))
    
    # B-spline design
    knoten <- seq(from=mini-3*dis, to=maxi+3*dis, length.out=knots)
    design <- splineDesign(knots=knoten, x, ord = 4, outer.ok = TRUE) # knots-3 columns, ki

    # transformation
    D <- a.spline(ki)
    PsiPen <- t(D)%*%solve(D%*%t(D))
    PsiUnpen <- matrix(c(rep(1, ki), knoten[1:ki]), ncol=2, byrow=FALSE)
    
    design <- (cbind(design%*%PsiUnpen, design%*%PsiPen))[,-1] # knots-4 columns, global intercept!! = ki-1

    # names
    colnames(design) <- paste(".", 0:(ki-2), sep="") # knots- order - 1 (intercept!!) columns, including 1 free slope parameter
    colnames(design)[1] <- ".slope"
    
    return(design)

}


vspline <- # p-spline with categorical effect modifier
function(x, u, knots=20){

  a.spline <- function (k) 
      { 
        if (k>1)
        {
          if (k > 2)
          {
            m1 <- cbind(diag(k-2),0,0)
            m2 <- cbind(0,diag(-2,k-2),0)
            m3 <- cbind(0,0,diag(k-2))
            mat <- t(m1 + m2 + m3)
          } else { mat <- matrix(ncol=0,nrow=2) }
        } else {
        mat <- matrix(ncol=0,nrow=1)
        }
        
        return (t(mat))
      }

    if (!is.numeric(x))
      stop("x in 'mspline(x, u)' must be numeric.")     
    if (is.factor(x))
       stop("varying coefficient not well defined. \n")
    if (!is.factor(u))
       stop("effect modifier must be nominal or ordinal. \n")

    # options
    old.contr <- getOption("contrasts")
    options(contrasts = c("contr.treatment", "contr.treatment"))
    options(na.action=na.pass)  

    # dummies
    int <- rep(1, times = length(u))
    dummies <- as.matrix(model.matrix(~ u)[,-1])
    dummies <- cbind( int-rowSums(dummies) , dummies)

    # definitions    
    mini <- min(x)
    maxi <- max(x)
    ki <- knots - 4
    dis <- (maxi - mini)/(knots - 7) # by = ((to - from)/(length.out - 1))
    knoten <- seq(from=mini-3*dis, to=maxi+3*dis, length.out=knots)

    # transformation matrices
    D <- a.spline(ki)
    PsiPen <- t(D)%*%solve(D%*%t(D))
    PsiUnpen <- matrix(c(rep(1, ki), knoten[1:ki]), ncol=2, byrow=FALSE)

    # design
    design <- matrix(nrow=length(x), ncol=0)
   
    for (i in 1:ncol(dummies)) {
      
      design.klein <- splineDesign(knots=knoten, x[which(dummies[,i]!=0)], ord = 4, outer.ok = TRUE) # B-spline design      
      design.klein <- (cbind(design.klein%*%PsiUnpen, design.klein%*%PsiPen))[,-1] # transformation 
      design.gross <- matrix(0, ncol=ki-1, nrow=length(x))
      design.gross[which(dummies[,i]!=0),] <- design.klein
      colnames(design.gross)    <- paste(".",levels(u)[i],".", 0:(ki-2), sep="")
      colnames(design.gross)[1] <- paste(".", levels(u)[i], ".slope", sep="")
      design <- cbind(design, design.gross)

    }
       
    options(contrasts = old.contr)
    return(design)

}



