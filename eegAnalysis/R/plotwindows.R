plotwindows <-
function(data, classes.Id , rec.Id , which.classes = "ALL", which.rec="ALL",
 which.channels ="ALL", win=10, stat="sum", power = 2, abs=FALSE, log=FALSE,complete = FALSE, mintomax=FALSE) {

  which.elec = which.channels
  #par(mfrow=c(1,1))
    
  if (!(stat%in%c("sum","mean","var","sd","max","prod","median","min","geometric","harmonic"))) {
   stop("Parameter stat: \n'stat' must be one of 'sum','mean','var','sd','max','prod','median','min','geometric','harmonic'.")
  }

  if (nrow(data)!=length(rec.Id) ||  length(rec.Id)!=length(classes.Id)) stop("The number of rows in 'data' must be equal to the length of the vectors rec.Id and classes.Id.")
    
  #Caso "ALL" redefinindo os vetores
  if (which.elec[1]=="ALL") {
    if (is.null(ncol(data))) {which.elec<-1} else which.elec <- 1:ncol(data)
  }
  if (which.classes[1]=="ALL") which.classes <- unique(classes.Id)
  
  #Numerando as classes.
  cl<-numeric(length(classes.Id))
  lwc <- length(which.classes)
  for (i in 1:lwc){
    cl[which(classes.Id==which.classes[i])]=i
  }
  which.classes<-1:lwc
  classes.Id<-cl
  
  if (which.rec[1] == "ALL") which.rec<-lapply(1:lwc,function (g) unique(rec.Id[which(classes.Id==which.classes[g])]))
  
  if (is.null(ncol(data))){
   L<-length(data[which(rec.Id%in%which.rec[[1]][1] & classes.Id == which.classes[1])])
  }else{
   L<-length(data[which(rec.Id%in%which.rec[[1]][1] & classes.Id == which.classes[1]),1])
  }
  
  #Definindo o banco de dados que deve ser utilizado.
  v<-c()
  for (g in 1:lwc){
    v<-c(v,which(rec.Id%in%which.rec[[g]] & classes.Id == which.classes[g]))
  }
  if (is.null(ncol(data))){
    dados<-data[v] } else{
    dados<-data[v,which.elec]
  }
  rec.Id<-rec.Id[v]
  classes.Id<-classes.Id[v]
  

  
  Rest<-L%%win
  Tot<- length(rec.Id)/L
  
  
  if (complete){
  
    if (abs) dados<-abs(dados)
  
    dados<-dados^power
  
    if (min(dados)<0 & log) 
    {
      dados<-abs(dados)
      warning("Parameter log: \nSome negative values were found. Using abs=T instead.")
    }
    if (log) dados<-log(dados+1)
    
    Dados <- numeric(length(which.elec)*(L-win+1)*Tot)
    el <- numeric(length(which.elec)*(L-win+1)*Tot)
    Reps <- numeric(length(which.elec)*(L-win+1)*Tot)
    Classes <- numeric(length(which.elec)*(L-win+1)*Tot)
    cont <- 0
    for (E in 1:length(which.elec)){
      for (CL in which.classes) {
        for (R in which.rec[[CL]]){
          cont <- cont+1
          if (length(which.elec)==1) prov <- dados[which(classes.Id==CL & rec.Id==R)] else prov <- dados[which(classes.Id==CL & rec.Id==R),E]
          for (i in 1:(L-win+1)){
            el[(cont-1)*(L-win+1)+i] <- E
            Reps[(cont-1)*(L-win+1)+i] <- R
            Classes[(cont-1)*(L-win+1)+i] <- CL
            
            if (stat=="sum") Dados[(cont-1)*(L-win+1)+i] <- sum(prov[i:(i+win-1)])
            if (stat=="mean") Dados[(cont-1)*(L-win+1)+i] <- mean(prov[i:(i+win-1)])
            if (stat=="var") Dados[(cont-1)*(L-win+1)+i] <- var(prov[i:(i+win-1)])
            if (stat=="sd") Dados[(cont-1)*(L-win+1)+i] <- sd(prov[i:(i+win-1)])
            if (stat=="median") Dados[(cont-1)*(L-win+1)+i] <- median(prov[i:(i+win-1)])
            if (stat=="max") Dados[(cont-1)*(L-win+1)+i] <- max(prov[i:(i+win-1)])
            if (stat=="min") Dados[(cont-1)*(L-win+1)+i] <- min(prov[i:(i+win-1)])
            if (stat=="prod") Dados[(cont-1)*(L-win+1)+i] <- prod(prov[i:(i+win-1)])
            if (stat=="geometric") Dados[(cont-1)*(L-win+1)+i] <- (prod(prov[i:(i+win-1)]))^(1/(win))
            if (stat=="harmonic") Dados[(cont-1)*(L-win+1)+i] <- 1/mean(1/prov[i:(i+win-1)])
          }
        }
      }
    }

    dados <- Dados
    rec.Id <- Reps
    classes.Id <- Classes

    LG <- L-win+1

  } else {
    
    if (Rest!=0) {  
      cei<-ceiling(Rest/2)
      fl<-floor(Rest/2)
      discart <- numeric(Tot*Rest)
      for (i in 1:Tot) {
      discart[((i-1)*Rest+1):((i-1)*Rest+cei)] <- c(((i-1)*L+1):((i-1)*L+cei))
      discart[((i-1)*Rest+1+cei):(i*Rest)] <- c((i*L-fl+1):(i*L))
    }
    
    if (length(which.elec)==1) dados<-dados[-discart] else dados<-dados[-discart,]
      rec.Id<-rec.Id[-discart]
      classes.Id<-classes.Id[-discart]
    }
    
    COLS<-length(which.elec)
    ROWS<-length(rec.Id)
    
    if (abs) dados<-abs(dados)
    dados<-dados^power
    if (min(dados)<0 & log) 
    {
      dados<-abs(dados)
      warning("Parameter log: \nSome negative values were found. Using abs=T instead.")
    }
    if (log) dados<-log(dados+1)
    
    rec.Id<-matrix(rec.Id,COLS*ROWS/win,win,byrow=TRUE)[,1]
    classes.Id<-matrix(classes.Id,COLS*ROWS/win,win,byrow=TRUE)[,1]
    
    if (stat=="sum") {dados<-rowSums(matrix(dados,COLS*ROWS/win,win,byrow=TRUE))}
    if (stat=="mean")  {dados<-rowMeans(matrix(dados,COLS*ROWS/win,win,byrow=TRUE))}
    if (stat=="min")  {dados<-rowMins(matrix(dados,COLS*ROWS/win,win,byrow=TRUE))}
    if (stat=="max")  {dados<-rowMaxs(matrix(dados,COLS*ROWS/win,win,byrow=TRUE))}
    if (stat=="var")  {dados<-t(colVars(t(matrix(dados,COLS*ROWS/win,win,byrow=TRUE))))}
    if (stat=="sd")  {dados<-sqrt(t(colVars(t(matrix(dados,COLS*ROWS/win,win,byrow=TRUE)))))}
    if (stat=="median")  {
      m<-matrix(dados,COLS*ROWS/win,win,byrow=TRUE)
      dados<-sapply(1:nrow(m),function (g) median(m[g,]))
    }
    if (stat=="prod")  {
      m<-matrix(dados,COLS*ROWS/win,win,byrow=TRUE)
      dados<-sapply(1:nrow(m),function (g) prod(m[g,]))
    }
    if (stat=="geometric")  {
      m<-matrix(dados,COLS*ROWS/win,win,byrow=TRUE)
      dados<-sapply(1:nrow(m),function (g) (prod(m[g,]))^(1/win))
    }
    if (stat=="harmonic")  {
      m<-matrix(dados,COLS*ROWS/win,win,byrow=TRUE)
      dados<-sapply(1:nrow(m),function (g) 1/mean(1/m[g,]) )
    }

    
    el <- numeric(length(rec.Id))
    NN <- length(rec.Id)/length(which.elec)
    for (i in 1:(length(which.elec))) {
       el[((i-1)*NN+1):(i*NN)]<-rep(i,NN)
    }

    
    LG <- ROWS/win/Tot

  }

  
  if (mintomax) {
    
    cont<-0
    for (E in 1:length(which.elec)){
      min <- min(dados[which(el==E)])
      max <- max(dados[which(el==E)])
      plot(1:LG,1:LG,type='l',col="white",ylim=c(min,max),xlab="",ylab="",
      main=paste("Channel",which.elec[E],sep=" "))
      for (CL in which.classes) {
       for (R in which.rec[[CL]]){
          lines(1:LG,sort(dados[which(el==E & classes.Id==CL & rec.Id==R)]),col=CL)
       }
      }
      cont<-cont+1
      if (which.elec[cont]==which.elec[length(which.elec)]) break;
      if (length(which.elec)>1){ cat("\n","Enter s to stop. Enter any other key to continue.","\n") # prompt
              key<-readline()}
      if (key=="s") break;  
    }

  } else{

    cont<-0
    for (E in 1:length(which.elec)){
      min <- min(dados[which(el==E)])
      max <- max(dados[which(el==E)])
      plot(1:LG,1:LG,type='l',col="white",ylim=c(min,max),xlab="",ylab="",
      main=paste("Channel",which.elec[E],sep=" "))
      for (CL in which.classes) {
        for (R in which.rec[[CL]]){
          lines(1:LG,dados[which(el==E & classes.Id==CL & rec.Id==R)],col=CL)
        }
      }
      cont<-cont+1
      if (which.elec[cont]==which.elec[length(which.elec)]) break;
      if (length(which.elec)>1){ 
          cat("\n","Enter s to stop. Enter any other key to continue.","\n") # prompt
          key<-readline()
      }
          if (key=="s") break;  
    }

  }#mintomax
  
  L0<-(length(classes.Id)/length(which.elec))
  classes.Id<-classes.Id[1:L0]
  rec.Id<-rec.Id[1:L0]
  
  cont<-0
  for(CL in which.classes) {
    for(R in 1:length(unique(rec.Id[which(classes.Id==CL)]))) {
      cont<-cont+1
      rec.Id[((cont-1)*LG+1):(cont*LG)]<-rep(R,LG)
    }
  }

  prov <- mat.or.vec(L0,length(which.elec))
  
  if (length(which.elec)!=1){
  
    for (E in 1:length(which.elec)){
      prov[,E]<-dados[((E-1)*L0+1):(E*L0)]
    }
    dados<-prov
  }
  
  result <- list(data=dados, classes.Id=classes.Id, rec.Id=rec.Id, which.classes=which.classes,
  which.rec=which.rec, which.channels =which.elec, stat=stat, power=power, 
  abs=abs,log=log, win=win,mintomax=mintomax, length=LG)
  class(result) <- "plotwindows"
  return(result)
}




print.plotwindows <-
  function(x,...){
    cat("Features by windows using the statistic ", x$stat,".\n",sep="")
    cat("Window size: ",x$win,".\n",sep="")
    if (x$abs){
      if (x$log) {cat("Transformation in data: log(abs(data^",x$power,")+1).\n",sep="")} else {cat("Transformation in data: abs(data^",x$power,").\n",sep="")}
    } else {
      if (x$log) {cat("Transformation in data: log(data^",x$power,"+1).\n",sep="")} else {cat("Transformation in data: data^",x$power,".\n",sep="")}
    }
  }


summary.plotwindows <-
  function(object,...){
    x<-object
    cat("Features by windows using the statistic ", x$stat,".\n",sep="")
    cat("Window size: ",x$win,".\n",sep="")
    if (x$abs){
      if (x$log) {cat("Transformation in data: log(abs(data^",x$power,")+1).\n",sep="")} else {cat("Transformation in data: abs(data^",x$power,").\n",sep="")}
    } else {
      if (x$log) {cat("Transformation in data: log(data^",x$power,"+1).\n",sep="")} else {cat("Transformation in data: data^",x$power,".\n",sep="")}
    }
  }


plot.plotwindows <-
  function(x,...){
    
    which.elec<-x$which.channels
    which.classes<-x$which.classes
    which.rec<-x$which.rec
    mintomax<-x$mintomax
    LG<-x$length
    
    for (i in 1:length(which.classes)){
      which.rec[[i]]<-1:(length(which.rec[[i]]))
    }
    
    if(is.null(ncol(x$data))){ 
      
      dados<-x$data
      rec.Id<-x$rec.Id
      classes.Id<-x$classes.Id
      el<-rep(1,length(dados))
      
    }else{
      NR<-nrow(x$data)
      L0<-nrow(x$data)*ncol(x$data)
      dados<-numeric(L0)
      classes.Id<-numeric(L0)
      rec.Id<-numeric(L0)
      el<-numeric(L0)
      
      for (g in 1:ncol(x$data)) {
        dados[((g-1)*NR+1):(g*NR)]<-x$data[,g] 
        classes.Id[((g-1)*NR+1):(g*NR)]<-x$classes.Id
        rec.Id[((g-1)*NR+1):(g*NR)]<-x$rec.Id
        el[((g-1)*NR+1):(g*NR)]<- rep(g,NR)
      } #for
      
    }
    
    
    if (mintomax) {
      
      cont<-0
      for (E in 1:length(which.elec)){
        min <- min(dados[which(el==E)])
        max <- max(dados[which(el==E)])
        plot(1:LG,1:LG,type='l',col="white",ylim=c(min,max),xlab="",ylab="",
             main=paste("Channel",which.elec[E],sep=" "))
        for (CL in which.classes) {
          for (R in which.rec[[CL]]){
            lines(1:LG,sort(dados[which(el==E & classes.Id==CL & rec.Id==R)]),col=CL)
          }
        }
        cont<-cont+1
        if (which.elec[cont]==which.elec[length(which.elec)]) break;
        if (length(which.elec)>1){ cat("\n","Enter s to stop. Enter any other key to continue.","\n") # prompt
                                   key<-readline()}
        if (key=="s") break;  
      }
      
    } else{
      
      cont<-0
      for (E in 1:length(which.elec)){
        min <- min(dados[which(el==E)])
        max <- max(dados[which(el==E)])
        plot(1:LG,1:LG,type='l',col="white",ylim=c(min,max),xlab="",ylab="",
             main=paste("Channel",which.elec[E],sep=" "))
        for (CL in which.classes) {
          for (R in which.rec[[CL]]){
            lines(1:LG,dados[which(el==E & classes.Id==CL & rec.Id==R)],col=CL)
          }
        }
        cont<-cont+1
        if (which.elec[cont]==which.elec[length(which.elec)]) break;
        if (length(which.elec)>1){ cat("\n","Enter s to stop. Enter any other key to continue.","\n") # prompt
                                   key<-readline()}
        if (key=="s") break;  
      }
      
    }#mintomax
  }
