read.magn<-function(data)
{ 
  read.fwf(data, skip=3,fill=T,widths=c(5,7,1,6,1,3,3,3,5,5,9,3,6,5,rep(6,13),5),
           col.names=c("IMOcode","long","EW","lat","NS","day","month","year","start",
                        "stop","sollong","Shw","lmg","m6","m5","m4","m3","m2","m1",
                        "zero","p1","p2","p3","p4","p5","p6","p7","N"))
 }
