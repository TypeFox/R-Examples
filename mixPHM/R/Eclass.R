`Eclass` <-
function(x, old, K, method, Sdist,p, cutpoint)
{

shape <- matrix(NA, K, p)	           # K x p Matrix with Shape-Parameter
scale <- matrix(NA, K, p)            # K x p Matrix with Scale-Parameter

#-------------sanity checks during EM-iteration-----------------
if (any(as.vector(table(old))==1)) {   #if a frequency equals 1
  outlier <- (1:length(old))[which(old==which(table(old)==1))]
  cat("Cluster contains only one observation! Subject",outlier,"may be an outlier!\n") 
  stop("Cannot proceed with estimation!")
}
if (length(unique(old))!=K) {         #if a cluster doesn't contain any element
  stop("Cluster contains 0 elements! Re-run with less components!")
}

for (j in 1:K) {
  y <- as.matrix(x[old==j,]) 
  ttab <- apply(y,2,table,exclude=0)         #table of dwell-times (list)
  lvec <- sapply(ttab,length)                #vector with different dwell-times for each page
  ind0 <- which(lvec<=1)                     #column index for those with less than 2 values
  rep.el <- sort(unique(as.vector(y)))[2:3]  #elements for 0-replacement (2 smallest except 0)
  if (length(ind0) >= 1) {
       for (i in ind0) y[,i][which(y[,i]==0)][1:2] <- rep.el
       warning("Complete 0 survival times in cluster occured. Two of them are replaced by minimum survival times in order to proceed with estimation!")
  }
  x[old==j,] <- y
}
#-------------end sanity checks----------------------

priorl <- by(x,old,function(y) {                           #list of prior probabilities
                    y <- as.matrix(y)
                    nses <- length(y[,1])                  #number of sessions in group
                    apply(y,2,function(z){ lz <- length(z[z>0])/nses})  
                    })
prior <- matrix(unlist(priorl),ncol=p,byrow=TRUE)                             #matrix of prior probabilities

#------------- separate ----------------------
if (method=="separate") {  
   parlist <- tapply(1:dim(x)[1],old, function(ind) {
                         y <- as.matrix(x[ind,])
                         apply(y,2,function(z) {
                            censvec <- rep(1, length(z))   
                            censvec[z > cutpoint] <- 0     #vector for censored data (set to 0)
                            wphm <- survreg(Surv(z[z>0], censvec[z>0])~1,dist=Sdist)        #wphm for each page within group
                            shapep <- 1/wphm$scale
                            scalep <- exp(wphm$coefficients[1])
                            list(scalep,shapep)
                            }) })
                            
   shsclist <- tapply(unlist(parlist),rep(1:2,length(unlist(parlist))/2),function(y){
                   matrix(y,nrow=K,byrow=TRUE)})                             #reorganizing parlist
   shape <- shsclist[[2]]                                                      #shape matrix K,p
   scale <- shsclist[[1]]                                                      #scale matrix K,p
   anzpar <- 2*K*p
}

#---------------------- group contrast ----------------------
if (method=="main.g") {
    for (i in 1:p) {
       datreg  <- as.vector(x[,i])			#VD-vektor i-te Seite
       datreg  <- datreg[x[,i] > 0]
       censvec <- rep(1, length(datreg))   
       censvec[datreg > cutpoint] <- 0     #vector for censored data (set to 0)                 
       xold    <- old[x[,i] > 0]				#Gruppenvektor i-te Seite

       wphm <- survreg(Surv(datreg, censvec)~factor(xold),dist=Sdist)
       scalebase <- as.vector(wphm$coefficients[1])	#scale parameter group 1 (reference group)
       scalevec1 <- as.vector(exp(wphm$coefficients[2:K]+scalebase)) #scale parameter of the remaining groups
       scale [,i] <- c(exp(scalebase),scalevec1)
       shape [,i] <- 1/wphm$scale			#shape gruppen konstant
   }
anzpar <- K*p+p
}

#------------- page constasts -----------------

if (method=="main.p") {
    for (j in 1:K)  {
       datregmat <- as.matrix(x[old == j,])
       nsess <- dim(datregmat)[1]				#sessionanzahl in gruppe j
       pagevek <- rep(1:p,rep(nsess,p))				#Seitenvektor sessions in gruppe j
       datreg <- as.vector(datregmat)
       xold <- pagevek[datreg > 0]				#VD > 0
       datreg <- datreg[datreg > 0]
       censvec <- rep(1, length(datreg))   
       censvec[datreg > cutpoint] <- 0     #vector for censored data (set to 0)                 
           
       wphm <- survreg(Surv(datreg, censvec)~factor(xold),dist=Sdist)  		#xold bezieht sich auf seiten
       scalebase <- as.vector(wphm$coefficients[1])
       scalevec1 <- as.vector(exp(wphm$coefficients[2:p]+scalebase))
       scale[j,] <- c(exp(scalebase),scalevec1)
       shape[j,] <- 1/wphm$scale				#shape bleibt seiten konstant
     }
anzpar <- K*p+K
}

#------------ page*group interaction  ----------------

if (method=="int.gp") {
    datreg <- as.vector(x)
    nsess <- dim(x)[1]
    pagevek <- rep(1:p,rep(nsess,p))
    oldall <- rep(old,p)
    xoldg <- oldall[datreg > 0]				#Gruppencontrast
    xoldp <- pagevek[datreg > 0]			#Seitencontrast
    datreg <- datreg[datreg > 0]
    censvec <- rep(1, length(datreg))   
    censvec[datreg > cutpoint] <- 0     #vector for censored data (set to 0)                 
     
    wphm <- survreg(Surv(datreg, censvec)~factor(xoldg)*factor(xoldp),dist=Sdist)
    scalebase <- as.vector(exp(wphm$coefficients[1]))
    scaleg <- exp(c(0,wphm$coefficient[2:K]))		#group contrast
    scalep <- exp(c(0,wphm$coefficient[(K+1):(K+p-1)])) #page contrast
    scaleimat <- matrix(exp(wphm$coefficient[(K+p):(K*p)]),(K-1),(p-1)) #interaction effects
    scaleimat <- rbind(rep(1,p),cbind(rep(1,K-1),scaleimat))
    scaletemp <- outer(scaleg,scalep)*scalebase
    scale <- scaletemp*scaleimat
    shape <- matrix((1/wphm$scale),K,p)
    anzpar <- K*p+1
}

#------------------ page + group main effects ----------------

if (method=="main.gp") {
    datreg <- as.vector(x)
    nsess <- dim(x)[1]
    pagevek <- rep(1:p,rep(nsess,p))
    oldall <- rep(old,p)
    xoldg <- oldall[datreg > 0]				#Gruppencontrast
    xoldp <- pagevek[datreg > 0]			#Seitencontrast
    datreg <- datreg[datreg > 0]
    censvec <- rep(1, length(datreg))   
    censvec[datreg > cutpoint] <- 0     #vector for censored data (set to 0)                 
    
    wphm <- survreg(Surv(datreg, censvec)~factor(xoldg)+factor(xoldp),dist=Sdist)
    scalebase <- as.vector(exp(wphm$coefficients[1]))
    scaleg <- exp(c(0,wphm$coefficient[2:K]))		#group contrast
    scalep <- exp(c(0,wphm$coefficient[(K+1):(K+p-1)]))	#page contrast
    scale <- outer(scaleg,scalep)*scalebase
    shape <- matrix((1/wphm$scale),K,p)
    anzpar <- K+p
}

list (scale = scale, shape = shape, prior = prior, anzpar = anzpar)
}

#returns matrices with shape and scale parameters as well as prior matrix
