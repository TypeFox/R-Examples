
pencopula <- function(data,d=3,D=d,q=1,base="B-spline",max.iter=20,plot.bsp=FALSE,
                    lambda=NULL,pen.order=2,adapt.grid=FALSE,add=TRUE,
                    alpha=0,symmetric=TRUE,data.frame=parent.frame()) {
  penden.env <- new.env()
  assign("frame",data.frame,penden.env)
  assign("Y",data,penden.env)
  
  if(is.matrix(data)|is.data.frame(data)) den.form <-  pendenForm(penden.env)
  else stop(" 'data' is not a matrix or a data frame. Create 'data' column by column.")

  if(is.null(d)) d <- 3 # Anzahl Halbierungen des Intervals [0,1]
  if(is.null(D)) D <- d # max. Hierachiestufe

  if(base=="B-spline" & D<d) d <- D
  if(base=="Bernstein") {
    D <- 2*d
    pen.order <- 1
  }
  assign("D",D,penden.env)
  assign("d",d,penden.env)
  
  assign("add",add,penden.env)
  assign("pen.order",pen.order,penden.env)
  assign("alpha",alpha,penden.env)
  assign("symmetric",symmetric,penden.env)
  assign("no",FALSE,penden.env)
  assign("max.iter",max.iter,penden.env)
  assign("base",base,penden.env)
  if(base=="Bernstein") adapt.grid <- FALSE
  assign("adapt.grid",adapt.grid,penden.env)
  assign("plot.bsp",plot.bsp,penden.env)
  p <- get("p",penden.env)  #p Dimension der Kovariablen
  
  if(base=="B-spline") {
    dd <- (2**d)+1 # Anzahl Knoten
    ddb <- dd + (q-1) # Anzahl Basisfunktionen
    if(is.null(q)) q <- 1 # Grad des B spline
    assign("q",q,penden.env)
  }
  if(base=="Bernstein") {
    dd <- d
    ddb <- dd+1# Anzahl Basisfunktionen
    DD <- ddb^p
    assign("DD",DD,penden.env)
  }
  assign("dd",dd,penden.env)
  assign("ddb",ddb,penden.env)
  
  assign("D",D,penden.env) #max. Hierarchiestufe

  if(base=="B-spline") dimension <- c(rep(0,q+1),rep(1:d,2**(0:(d-1))))
  if(base=="Bernstein") dimension <- seq(0,D)

  if(base=="B-spline") {
    if(is.null(lambda)) lambda <- rep(10000,p)
    if(!is.null(lambda)) if(length(lambda)<p | length(lambda)>p) stop("length of lambda is wrong")
  }
  if(base=="Bernstein") lambda <- rep(0,p)
  assign("lambda",lambda,penden.env)
  assign("dimension",dimension,penden.env)
  
    
# Dimension gibt die Hierarchiestufe an, aus der der hierarchische B-spline berechnent wird
  
  ##################

  if(base=="B-spline") {  
  #D maximale Hierarchiestufe
    DIMENSION <- dimension
    Index.basis <- matrix(1:ddb)
    index.sparse <- DIMENSION <= D
    Index.basis.D <- matrix(Index.basis[index.sparse,])
    DIMENSION <- DIMENSION[index.sparse]
    
    for ( j in 2:p)
      {
                                        #print(j)
        DIMENSION.j <-  kronecker(matrix(1,ddb,1),DIMENSION) + kronecker( dimension, matrix(1, length(DIMENSION),1))
        Index.basis.plus.1 <- matrix(NA, dim(Index.basis.D)[1] * ddb , j)
        Index.basis.plus.1[,j] <- kronecker(matrix(1:ddb), matrix(1,dim(Index.basis.D)[1],1))
        Index.basis.plus.1[, 1:(j-1)] <-  kronecker(matrix(1, ddb,1),Index.basis.D)
        index.sparse <- DIMENSION.j <= D
        Index.basis.D <- Index.basis.plus.1[index.sparse,]
        DIMENSION <- DIMENSION.j[index.sparse]
      }
    DD <- dim(Index.basis.D)[1] # Dimension of sparse grid basis
    assign("DD",DD,penden.env) # DD Anzahl Koeffizienten
  }
  if(base=="Bernstein") {
    assign("knots",seq(0,1,length=ddb),penden.env)
  
    Index.basis.D <- matrix(NA,DD,p)
    help.val <- rep(seq(1,ddb),ddb)
    help.val2 <- sort(help.val)
    help <- cbind(help.val,help.val2)
    if(p>2) {
      for(j in 3:p) {
        help <- kronecker(matrix(1,ddb,1),help)
        help <- cbind(help,sort(help[,(j-1)]))
      }
    }
    Index.basis.D <- help
  }
  assign("Index.basis.D",Index.basis.D,penden.env)
  
  
  ###################

  # Matrix zur Erstellung der marginalen Spline Koeffizienten
  # 
  j <- 1
  T.marg <- array(NA, dim=c(ddb,DD,p))
  
  for ( j in 1:p)
    {
      for ( l in 1:ddb)
        {
          T.marg[l,,j] <- (Index.basis.D[,j] == l)+0
        }
    }
  
  assign("T.marg",T.marg,penden.env)

  ####################

  #knots
  knots.start(penden.env)
 
  #####################

  tilde.Psi.d <-  array(NA, dim=c(get("n",penden.env),ddb,p))
 
  for (j in 1:p)
    {
      tilde.Psi.d[,,j] <-  hierarch.bs(get("Y",penden.env)[,j], d = d, plot.bsp = plot.bsp,typ=3,penden.env,int=FALSE)$B.tilde
    }

  assign("tilde.Psi.d",tilde.Psi.d,penden.env)

  assign("tilde.PSI.d.D",tilde.Psi.d[,Index.basis.D[,1],1],penden.env)

  for (j in 2:p)
    {
      assign("tilde.PSI.d.D",get("tilde.PSI.d.D",penden.env) * get("tilde.Psi.d",penden.env)[,Index.basis.D[,j],j],penden.env)
    }

  #####################

  #startwerte und gitter berechnen

  start.valgrid(penden.env)

  #####################################

  # Allgemein: Die Restriktionsmatrizen sind als Array aufgefasst:
  
  A <- array(NA, dim=c(get("ddb",penden.env),DD,p))

  for ( j in 1:p)
    {
      A[,,j] <- get("tilde.Psi.knots.d",penden.env) %*% T.marg[,,j]
    }

  assign("A.Restrict",A,penden.env)
  
  # Nun muss gelten A[,,j] %*% c = 1 fuer alle j

  #############################

  penalty.matrix(penden.env=penden.env)
  
  #############################

  liste <- matrix(0,1,3+DD+p)

  lam <- coef <- c()
  for(i in 1:p) lam[i] <- paste("lambda.",i,sep="")
  for(j in 1:DD) coef[j] <- paste("b.",j,sep="")
 
  colnames(liste) <- c("pen.log.like","log.like","marg.log.like",lam,coef)
   
  #print("fitted values at the beginning")
  help.str <- paste("d=",get("d",penden.env),"D=",get("D",penden.env),"lambda=",get("lambda",penden.env)[1],sep="")
  assign("help.str",help.str,penden.env)
  assign("liste",liste,penden.env)
 
  calculate(penden.env)
  
  class(penden.env) <- "pencopula"
  return(penden.env)
}
