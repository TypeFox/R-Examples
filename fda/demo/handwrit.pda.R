#  load the data

load('../data/handwrit.rda')
  
#  array with coordinates and vector of times

fdaarray = handwrit
fdatime  <- seq(0, 2.3, len=1401)

#  basis for coordinates

fdarange <- c(0, 2.3)
breaks = seq(0,2.3,length.out=501)
norder = 6
fdabasis = create.bspline.basis(fdarange,norder=norder,breaks=breaks)

#  parameter object for coordinates

fdaPar = fdPar(fdabasis,int2Lfd(4),1e-8)

#  coordinate functions and a list tontaining them

Xfd = smooth.basis(fdatime, fdaarray[,,1], fdaPar)$fd
Yfd = smooth.basis(fdatime, fdaarray[,,2], fdaPar)$fd

xfdlist = list(Xfd, Yfd)

#  basis and parameter object for weight functions

fdabasis2 = create.bspline.basis(fdarange,norder=norder,nbasis=51)
pdaPar = fdPar(fdabasis2,1,1e-8)

#  set up a first order equation

pdaParlist = list(pdaPar)

bwtlist = list( list(pdaParlist,pdaParlist), list(pdaParlist,pdaParlist) )

#  do the pda

pdaList = pda.fd(xfdlist, bwtlist)

#  extract the weight functions

bwtlistout = pdaList$bwtlist
bwtfd11    = bwtlistout[[1]][[1]][[1]]$fd
bwtfd21    = bwtlistout[[2]][[1]][[1]]$fd
bwtfd12    = bwtlistout[[1]][[2]][[1]]$fd
bwtfd22    = bwtlistout[[2]][[2]][[1]]$fd

#  plot the weight functions

par(mfrow=c(2,2))
plot(bwtfd11)
title("Weight for variable 1 in equation 1")
plot(bwtfd21)
title("Weight for variable 2 in equation 1")
plot(bwtfd12)
title("Weight for variable 1 in equation 2")
plot(bwtfd22)
title("Weight for variable 2 in equation 2")

#  plot residual functions

reslist = pdaList$resfdlist
par(mfrow=c(2,1))
plot(reslist[[1]])
title("Residual function for variable 1")
plot(reslist[[2]])
title("Residual function for variable 2")

#  set up a second order equation

pdaParlist = list(pdaPar, pdaPar)

bwtlist = list( list(pdaParlist,pdaParlist), list(pdaParlist,pdaParlist) )

#  do the second order pda

pdaList = pda.fd(xfdlist, bwtlist)

#  extract and plot the weight coefficients for X and Y

bwtlistout = pdaList$bwtlist
bwtfd11    = bwtlistout[[1]][[1]][[1]]$fd
bwtfd21    = bwtlistout[[2]][[1]][[1]]$fd
bwtfd12    = bwtlistout[[1]][[2]][[1]]$fd
bwtfd22    = bwtlistout[[2]][[2]][[1]]$fd
par(mfrow=c(2,2))
plot(bwtfd11)
title("Weight for variable 1 in equation 1")
plot(bwtfd21)
title("Weight for variable 2 in equation 1")
plot(bwtfd12)
title("Weight for variable 1 in equation 2")
plot(bwtfd22)
title("Weight for variable 2 in equation 2")

#  extract and plot the weight coefficients for DX and DY

bwtlistout = pdaList$bwtlist
bwtfd11    = bwtlistout[[1]][[1]][[2]]$fd
bwtfd21    = bwtlistout[[2]][[1]][[2]]$fd
bwtfd12    = bwtlistout[[1]][[2]][[2]]$fd
bwtfd22    = bwtlistout[[2]][[2]][[2]]$fd
par(mfrow=c(2,2))
plot(bwtfd11)
title("Weight for variable 1 in equation 1")
plot(bwtfd21)
title("Weight for variable 2 in equation 1")
plot(bwtfd12)
title("Weight for variable 1 in equation 2")
plot(bwtfd22)
title("Weight for variable 2 in equation 2")

#  display residual functions

reslist = pdaList$resfdlist
par(mfrow=c(2,1))
plot(reslist[[1]])
title("Residual function for variable 1")
plot(reslist[[2]])
title("Residual function for variable 2")
