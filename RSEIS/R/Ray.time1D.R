`Ray.time1D` <-
function(indelta, inhpz, instaz, inlay , ztop ,  vel)
  {
    ###  One Dimensional travel times with low velocity layers
    ##### tt1d  travel.time1D calculate the travel time, dtdr, dtdz and takeoff angle for the given vel, dist(km), hpz and stz
    
    ######  input:
    ######  indelta = distance from source to reciever
    ######  inhpz   = depth of source
    ######   instaz  = elevation of reciever
    ######    inlay  =number of layers
    ######    ztop  = tops of the vel model
    ######    vel   =  velocity in the layers

    ##  output:  dtdr=derivative with respect to r   dtdz=derivative with respect to z
    ###  angle=angle of takeoff measured from nadir (down vecotor)  tt=travel time

    dtdr=0;
    dtdz=0; angle=0; outt=0;
    nnod = inlay*4
    znod = rep(0,length=nnod)
    rnod = rep(0,length=nnod)

    slness = 1/vel
    
    TTout = .C("CALL_DTTray",PACKAGE = "RSEIS",
      as.double(indelta),
      as.double(inhpz),as.double(instaz),
      as.integer(inlay), as.double(ztop) ,  as.double(slness),
      as.double(dtdr), as.double(dtdz), as.double(angle),  as.double(outt),
      as.integer(nnod) , as.double(znod), as.double(rnod) )

    dtdr=as.numeric(unlist(TTout[7]));
    dtdz=as.numeric(unlist(TTout[8]));
    angle=as.numeric(unlist(TTout[9]));
    tt=as.numeric(unlist(TTout[10]));

    nnod = as.numeric(unlist(TTout[11]));
    znod = as.numeric(unlist(TTout[12]));
    rnod = as.numeric(unlist(TTout[13]));

    angle = 180*atan2(rnod[2]-rnod[1], znod[2]-znod[1] )/pi

    return(list( dtdr=dtdr, dtdz=dtdz, angle=angle, tt=tt, nnod=nnod, znod=znod, rnod=rnod)) 

  }

