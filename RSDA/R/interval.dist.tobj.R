interval.dist.tobj <-
function(sym.obj.x,sym.obj.y,distance=c('hausdorff','centers','interscal'),p=2) {
    distance<-match.arg(distance)
    idn1<-all(sym.obj.x$var.types=='$I')
    idn2<-all(sym.obj.y$var.types=='$I')
    if((idn1==FALSE)||(idn2==FALSE))
      stop("The two variables have to be interval type")         
    if(distance=='hausdorff') {
      pos<-1
      sum<-0
      for(j in 1:(sym.obj.x$M)) {
        a<-abs(sym.obj.x$obj.data.vector[pos]-sym.obj.y$obj.data.vector[pos])
        b<-abs(sym.obj.x$obj.data.vector[pos+1]-sym.obj.y$obj.data.vector[pos+1])
        m<-max(a,b)
        pos<-pos+2
        sum<-sum+m^p
      }
      return(sum^(1/p))
    }
    if(distance=='centers') {
        pos<-1
        sum<-0
        for(j in 1:(sym.obj.x$M)) {
          a<-(sym.obj.x$obj.data.vector[pos]+sym.obj.x$obj.data.vector[pos+1])/2
          b<-(sym.obj.y$obj.data.vector[pos]+sym.obj.y$obj.data.vector[pos+1])/2
          m<-abs(a-b)
          pos<-pos+2
          sum<-sum+m^p
        }
        return(sum^(1/p))
    }    
    if(distance=='interscal') { 
      #  Daneaux and Masson Distance
      pos<-1
      suma<-0
      sumb<-0
      for(j in 1:(sym.obj.x$M)) {
        a<-(sym.obj.x$obj.data.vector[pos+1]-sym.obj.x$obj.data.vector[pos]+
              sym.obj.y$obj.data.vector[pos+1]-sym.obj.y$obj.data.vector[pos]-
              2*abs((sym.obj.x$obj.data.vector[pos+1]+sym.obj.x$obj.data.vector[pos])/2-
                      (sym.obj.y$obj.data.vector[pos+1]+sym.obj.y$obj.data.vector[pos])/2)-
              abs(sym.obj.x$obj.data.vector[pos+1]-sym.obj.x$obj.data.vector[pos]+
                    sym.obj.y$obj.data.vector[pos+1]-sym.obj.y$obj.data.vector[pos]-
                    2*abs((sym.obj.x$obj.data.vector[pos+1]+sym.obj.x$obj.data.vector[pos])/2-
                            (sym.obj.y$obj.data.vector[pos+1]+sym.obj.y$obj.data.vector[pos])/2)))^2
        suma<-suma+a                
        b<-(sym.obj.x$obj.data.vector[pos+1]-sym.obj.x$obj.data.vector[pos]+
            sym.obj.y$obj.data.vector[pos+1]-sym.obj.y$obj.data.vector[pos]+
            2*abs((sym.obj.x$obj.data.vector[pos+1]+sym.obj.x$obj.data.vector[pos])/2-
                  (sym.obj.y$obj.data.vector[pos+1]+sym.obj.y$obj.data.vector[pos])/2))^2
        sumb<-sumb+b        
        pos<-pos+2
      }
      return(c((1/4)*sqrt(suma),(1/2)*sqrt(sumb)))
    }  
}
