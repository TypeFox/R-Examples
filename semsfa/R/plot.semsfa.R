plot.semsfa<-function(x,g.type="reg",mod=NULL,...){
  dati<-data.frame(x$data)
  if(g.type=="eff" | g.type=="reg"){
     if(g.type=="eff"){
     if(is.null(mod)){stop("argument mod non specified")}
      if(mod=="hist" | mod=="dens"){
      if(is.null(x$efficiencies)) {stop("efficiency not yet calculated")}
      else{
         if(mod=="hist") hist(x$eff,...)
         if(mod=="dens") plot(density(x$eff),...)
         }
      }else{stop("only histogram or density plot admitted")} 
     }else{
         if(x$sem.method!="loess") {plot(x$reg,...)}else{stop("plot is avalible only for gam or kernel methods")}
      }

  } else {stop("only efficiencies or regression results may be plotted")}
}
