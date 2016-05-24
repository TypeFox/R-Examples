knots.start <- function(penden.env) {
  if(get("base",penden.env)=="B-spline") {
    assign("knots",knots.transform(get("d",penden.env),get("alpha",penden.env),get("symmetric",penden.env)),penden.env)
    q <- get("q",penden.env)
    
    if(q==1) {
    x.help <- c()
    for(j in 1:(length(get("knots",penden.env))-1)) {
      x.help <- c(x.help,(get("knots",penden.env)[j+1]-get("knots",penden.env)[j])/2+get("knots",penden.env)[j])
    }
  }
    if(q==2) {
      x.help <- c()
      for(j in 1:(length(get("knots",penden.env))-1)) {
        x.help <- c(x.help,(get("knots",penden.env)[j+1]-get("knots",penden.env)[j])/3+get("knots",penden.env)[j],2*(get("knots",penden.env)[j+1]-get("knots",penden.env)[j])/3+get("knots",penden.env)[j])
      }
    }
    assign("knots.help",sort(c(x.help,get("knots",penden.env))),penden.env)
    
    if(q==2) {
      assign("knots.t",seq(0,1,length=get("ddb",penden.env)),penden.env)
      if(get("d",penden.env)==3) assign("k.order",c(1,5,10,6,3,8,2,4,7,9),penden.env)
      if(get("d",penden.env)==4) assign("k.order",c(1,9,18,10,5,14,3,7,12,16,2,4,6,8,11,13,15,17),penden.env)
    }
    
    if(q==1) {
      assign("k.order",knots.order(penden.env),penden.env)
      assign("knots.t",get("knots",penden.env),penden.env)
    }
                                        #Knoten nach hierarchischer Basis anordnen
    
    assign("tilde.Psi.knots.d",hierarch.bs(get("knots.t",penden.env), d = get("d",penden.env),plot.bsp=get("plot.bsp",penden.env),typ=3,penden.env,int=FALSE)$B.tilde,penden.env)
    assign("tilde.Psi.knots.d.help",hierarch.bs(get("knots.help",penden.env), d = (get("d",penden.env)+1),plot.bsp=get("plot.bsp",penden.env),typ=3,penden.env,int=FALSE)$B.tilde,penden.env)
    assign("knots.t",get("knots.t",penden.env)[get("k.order",penden.env)],penden.env)
  }
  if(get("base",penden.env)=="Bernstein") {
    assign("knots.t",get("knots",penden.env),penden.env)
    #assign("k.order",knots.order(penden.env),penden.env)
    assign("tilde.Psi.knots.d",hierarch.bs(get("knots.t",penden.env), d = get("d",penden.env),plot.bsp=get("plot.bsp",penden.env),typ=3,penden.env,int=FALSE)$B.tilde,penden.env)
    #assign("tilde.Psi.knots.d.help",hierarch.bs(get("knots.help",penden.env), d = (get("d",penden.env)+1),plot.bsp=get("plot.bsp",penden.env),typ=3,penden.env,int=FALSE)$B.tilde,penden.env)
    #assign("knots.t",get("knots.t",penden.env)[get("k.order",penden.env)],penden.env)
  }
}
