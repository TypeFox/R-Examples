`travel.time1D` <-
function(indelta, inhpz, instaz, inlay , ztop ,  vel)
  {
    ###  One Dimensional travel times with low velocity layers
    ##### tt1d  travel.time1D calculate the travel time, dtdr, dtdz
    #######   and takeoff angle for the given vel, dist(km), hpz and stz
    
    ######  input:
    ######  indelta = distance from source to reciever
    ######  inhpz   = depth of source
    ######   instaz  = elevation of reciever
    ######    inlay  =number of layers
    ######    ztop  = tops of the vel model
    ######    vel   =   slness  in the layers

    ##  output:  dtdr=derivative with respect to r   dtdz=derivative with respect to z
    ###  angle=angle of takeoff measured from nadir (down vecotor)  tt=travel time

    dtdr=0;
    dtdz=0; angle=0; outt=0;
    

    if(length(ztop)<inlay) { print("ERROR: BAD velocity model"); return(NULL) }
    if(length(vel)<inlay) { print("ERROR: BAD velocity model"); return(NULL) }
    if( any(indelta <0) )  { print("ERROR: BAD Distance in travel.time1D "); return(NULL) }

    
    if(any(!is.numeric(ztop))){ print("ERROR: BAD velocity model"); return(NULL) }
    if(any(!is.numeric(vel))){ print("ERROR: BAD velocity model"); return(NULL) }
    
    TTout = .C("CALL_DTT1",PACKAGE = "RSEIS",
      as.double(indelta), as.double(inhpz),as.double(instaz), as.integer(inlay), as.double(ztop) ,  as.double(vel),
      as.double(dtdr), as.double(dtdz), as.double(angle),  as.double(outt) )

    dtdr=as.numeric(unlist(TTout[7]));
    dtdz=as.numeric(unlist(TTout[8]));
    angle=as.numeric(unlist(TTout[9]));
    tt=as.numeric(unlist(TTout[10]));

    return(list( dtdr=dtdr, dtdz=dtdz, angle=angle, tt=tt)) 

  }

