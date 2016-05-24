randEEG <-
function(n.classes=2,n.rec=10,n.channels = 20,n.signals=250,ar="default", ma="default", order="default", vars = c(1,2)) {
	
n.elec = n.channels

if(ar=="default")
	ar = mat.or.vec(n.classes,1)
if(ma=="default")
	ma = mat.or.vec(n.classes,1)
if(order=="default")
	order = mat.or.vec(n.classes,3)
  
#--------------------------------------------------------------
#INI: transforming vector into matrices to avoid errors  
  if (n.classes==1){
    if (is.null(nrow(order))) order<-t(as.matrix(order))
    if (is.null(nrow(ar))) ar<-t(as.matrix(ar))
    if (is.null(nrow(ma))) ma<-t(as.matrix(ma))
  } else {
    if (is.null(nrow(order))) order<-t(as.matrix(order))
    if (is.null(nrow(ar))) ar<-as.matrix(ar)
    if (is.null(nrow(ma))) ma<-as.matrix(ma)  
  }
#END: transforming vector into matrices to avoid errors 
#--------------------------------------------------------------

  
#--------------------------------------------------------------
#INI: input control  
  for(i in 1:length(vars))
  {
    if(vars[i]<=0) stop("Parameter 'vars': All entries must be greater than 0.")}
    if(n.classes<=0) stop("Parameter 'n.classes': This number must be greater than 0.")
    if(n.rec<=0) stop("Parameter 'n.rec': This number must be greater than 0.")
    if(n.elec<=0) stop("Parameter 'n.channels': This number must be greater than 0.")
    if(n.signals<=0) stop("Parameter 'n.signals': This number must be greater than 0.")
    
    if (nrow(order)!=n.classes & nrow(order)!=(n.classes*n.elec)) stop(
      "Parameter 'order': \nThe number of rows in 'order' must be equal to 'n.classes' or 'n.classes*n.channels'")
    if (nrow(order)!=nrow(ar)) stop("Parameter 'ar': \nMatrix 'ar' must have the same number of rows of 'order'.")
    if (nrow(order)!=nrow(ma)) stop("Parameter 'ma': \nMatrix 'ma' must have the same number of rows of 'order'.")
    if (nrow(order)!=length(vars)) stop("Parameter 'vars': \nThe length of 'vars' must be equal to the number of classes.")

    m<-max(order[,1])
    if (ncol(ar)>m & m>0) {ar<-ar[,1:m] ; warning("Matrices 'order' and 'ar' do not conform. Some parameters were deleted in 'ar' matrix to proceed the simulation.")}
    repeat{ if (ncol(ar)<m) ar<-cbind(ar,numeric(nrow(ar))) else break;}
    m<-max(order[,3])
    if (ncol(ma)>m & m>0) {ma<-ma[,1:m]; warning("Matrices 'order' and 'ma' do not conform. Some parameters were deleted in 'ma' matrix to proceed the simulation.")}
    repeat{ if (ncol(ma)<m) ma<-cbind(ma,numeric(nrow(ma))) else break;}
    
    if (length(vars)!=n.classes & length(var)!=(n.classes*n.elec)) stop("Parameter 'vars': \nlength of 'vars' must be equal to the number of classes or equal to 'n.classes'*'n.channels'.")
#END: input control
#--------------------------------------------------------------
 
  
    data <- mat.or.vec(n.signals*n.classes*n.rec,n.elec)
    reps <- numeric(n.signals*n.classes*n.rec)
    class <- numeric(n.signals*n.classes*n.rec)
    for (i in 1:n.rec) 
    {
      reps[((i-1)*n.signals+1):(i*n.signals)]<-rep(i,n.signals)
    }
    for (i in 1:n.classes) 
    {
      class[((i-1)*n.signals*n.rec+1):(i*n.signals*n.rec)] <- rep(i,n.signals*n.rec) ;
                         if (i>1) reps[((i-1)*n.signals*n.rec+1):(i*n.signals*n.rec)]<-reps[1:(n.signals*n.rec)]
    }
  
    if (nrow(order)==n.classes) 
    {
      for (i in 1:n.elec)
      {
        for (j in 1:n.classes)
        {
          for (k in 1:n.rec)
          {
            ord <- order[j,]
            if (ord[1]==0) 
            {
              AR<-c()
            } else 
            {
              AR<-ar[j,1:ord[1]]
            }
            if (ord[3]==0) 
            {
              MA<-c()
            } else 
            {
              MA<-ma[j,1:ord[3]]
            }
            data[((j-1)*n.signals*n.rec+(k-1)*n.signals+1):((j-1)*n.signals*n.rec+k*n.signals),i] <-  tail(arima.sim(n = n.signals, list(ar=AR,ma=MA,order=ord),sd = sqrt(vars[j])),n.signals)
          } 
        }
      }
    } else 
    {
      for (i in 1:n.elec)
      {
        for (j in 1:n.classes)
        {
          for (k in 1:n.rec)
          {
            ind <- (j-1)*n.elec+i
            ord <- order[ind,]
            if (ord[1]==0) {AR<-c()} else {AR<-ar[ind,1:ord[1]]}
            if (ord[3]==0) {MA<-0} else {MA<-ma[ind,1:ord[3]]}
            data[((j-1)*n.signals*n.rec+(k-1)*n.signals+1):((j-1)*n.signals*n.rec+k*n.signals),i] <-  arima.sim(n = n.signals, list(ar=AR,ma=MA,order=ord),sd = sqrt(vars[ind]))
          } 
        }
      }
    
    }
  
  
  result <- list(
    data = data , classes.Id = class , rec.Id = reps , n.classes=n.classes, n.rec = n.rec, n.channels=n.elec,n.signals=n.signals,vars=vars
    )
  class(result) <- "RandEEG"
  return(result) 
}


print.RandEEG <-
  function(x, ...){
    cat("Class: RandEEG \n \n")
    cat("Description: Simulation of EEG data \n")
    cat("Number of classes:",x$n.classes,"\n")
    cat("Number of recordings of each class:",x$n.rec,"\n")
    cat("Number of channels:",x$n.channels,"\n")
    cat("Number of signals in each recording:",x$n.signals,"\n")
    cat("Variance for each class:",x$vars,"\n")
  }


summary.RandEEG <-
  function(object, ...){
    x<- object
    cat("Class: RandEEG \n \n")
    cat("Description: Simulation of EEG data \n")
    cat("Number of classes:",x$n.classes,"\n")
    cat("Number of recordings of each class:",x$n.rec,"\n")
    cat("Number of channels:",x$n.channels,"\n")
    cat("Number of signals in each recording:",x$n.signals,"\n")
    cat("Variance for each class:",x$vars,"\n")  
  }

