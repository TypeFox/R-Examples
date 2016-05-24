################################################
#### Copied bbplot (from uroot package)
################################################


bbplot <- function(wts, colour=c("SlateBlue","SeaGreen","red","magenta"))
{
  s <- frequency(wts)
  if(s==4){
    nl <- 4; np <- 1; fr <- 2
    snames <- c("Qrt1", "Qrt2", "Qrt3", "Qrt4")
  }
  else if(s==12){
    nl <- 3; np <- fr <- 4
    snames <- month.abb
  }
  else{
    if(identical(c(s/2 - as.integer(s/2)), 0)){
      nl <- 4; np <- s/nl
    }
    else{
      np <- floor((s+1)/4); nl <- c(rep(4, np-1), s/np+1)
    }
  }

  Mseas <- Msts(wts)

  aux <- Mseas[nrow(Mseas),]
  labaux <- c(na.omit(aux), Mseas[(nrow(Mseas)-1),which(is.na(aux))])

  xlim <- c(start(wts)[1], end(wts)[1]+1.5)
  ylim <- range(wts, na.rm=TRUE)
  aux1 <- seq(1,s,nl); aux2 <- seq(nl,s,nl)

  if(s==4)
    opar <- par(mar=c(3,3.5,2,2), las=1)
  if(s==12)
    opar <- par(mfrow=c(fr/2, fr/2), mar=c(3,3.5,2,2), las=1)
  for(i in 1:np){
    ref <- aux1[i]:aux2[i]
    ts.plot(ts(Mseas[,ref], start=start(wts)[1]), xlab="", ylab="",
      xlim=xlim, ylim=ylim , col=colour)

    #text(xlim[2]-0.2, as.matrix(Mseas[,ref])[nrow(Mseas),], month.abb[ref])
    text(xlim[2]-0.2, labaux[ref], snames[ref], col=colour)
  }
  par(opar)
}

################################################
#### Auxiliary function Msts (from uroot package)
################################################


Msts <- function(wts)  ##~ hacer que trate a las columnas como ts, ts.union.
{
  s <- frequency(wts)
  seas.data <- split(wts, cycle(wts))
  Mseas <- matrix(NA, nrow=ceiling(length(wts)/s), ncol=s)
  ref1 <- c(rep(2, start(wts)[2]-1), rep(1, s-start(wts)[2]+1))
  ref2 <- c(rep(0, end(wts)[2]), rep(1, s-end(wts)[2]))
  for(i in 1:s)
    Mseas[ref1[i]:(nrow(Mseas)-ref2[i]),i] <- seas.data[[i]]  # unlist(seas.data[[i]])

  ynames <- as.integer(time(wts))
  if(s==4)
    snames <- c("Qtr1", "Qtr2", "Qtr3", "Qtr4")
  else if(s==12)
    snames <- month.abb
  else
    snames <- paste("Seas", 1:s, sep="")

  Mseas <- matrix(Mseas, nrow=nrow(Mseas), ncol=ncol(Mseas),
           dimnames=list(as.character(ynames[1]:ynames[length(ynames)]), snames))

#  for(i in 1:ncol(Mseas))
#    Mseas[,i] <- ts(Mseas[,i], frequency=1, start=start(wts))
#  as.matrix(Mseas)
  Mseas
}

################################################
#### Deprecated bbplot (old code)
################################################

##~ Ver hacer método hasta obtener seas.data.
##~ ver usar cbind.ts.
# bbplot <- function(wts)
# {
#   colour <- c("SlateBlue","SeaGreen","red","magenta")
#   seas.data <- split(wts, cycle(wts))
# 
#   t0 <- c(rep(1, start(wts)[2]-1), rep(0, frequency(wts)-(start(wts)[2]-1)))
#   for(i in 1:frequency(wts))
#     assign(paste("season", i, sep=""), ts(seas.data[[i]], frequency=1, start=start(wts)[1]+t0[i]))
# 
#   if(frequency(wts) == 4){
#     seas.labels <- c("Q1","Q2","Q3","Q4")
#     seas.data <- list(season1, season2, season3, season4)
#     opar <- par(las=1)
#     ts.plot(season1, season2, season3, season4, lty=c(1,2,3,4), xlab="", ylab="",
#             xlim=c(start(wts)[1], end(wts)[1]+2.5)) ##~ col=colour
#     par(opar)
#     for(i in 1:4)
#       text(end(seas.data[[i]])[1]+2, seas.data[[i]][length(seas.data[[i]])], seas.labels[i])
#       ##~ , col = colour[i])
#   }
# 
#   if(frequency(wts) == 12){
#     seas.data <- list(season1, season2, season3, season4, season5, season6,
#                       season7, season8, season9, season10, season11, season12)
#     xlim <- c(start(wts)[1], end(wts)[1]+1.5); ylim <- c(min(wts), max(wts))
# 
#     opar <- par(mfrow=c(2,2), las=1)
#     ts.plot(season1, season2, season3, lty=c(1,2,3), xlab="", ylab="", xlim=xlim, ylim=ylim) ##~ col=colour)
#     for(i in 1:3)
#       text(end(seas.data[[i]])[1]+1, seas.data[[i]][length(seas.data[[i]])], month.abb[i]) ##~ , col =  colour[i])
#     ts.plot(season4, season5, season6, lty=c(1,2,3), xlab="", ylab="", xlim=xlim, ylim=ylim) ##~ col=colour)
#     for(i in 4:6)
#       text(end(seas.data[[i]])[1]+1, seas.data[[i]][length(seas.data[[i]])], month.abb[i]) ##~ , col = colour[i])
#     ts.plot(season7, season8, season9, lty=c(1,2,3), xlab="", ylab="", xlim=xlim, ylim=ylim) ##~ col=colour)
#     for(i in 7:9)
#       text(end(seas.data[[i]])[1]+1, seas.data[[i]][length(seas.data[[i]])], month.abb[i]) ##~ , col = colour[i])
#     ts.plot(season10, season11, season12, lty=c(1,2,3), xlab="", ylab="", xlim=xlim, ylim=ylim) ##~ col=colour)
#     for(i in 10:12)
#       text(end(seas.data[[i]])[1]+1, seas.data[[i]][length(seas.data[[i]])], month.abb[i]) ##~ , col = colour[i])
#     par(opar)
#   }
# }
