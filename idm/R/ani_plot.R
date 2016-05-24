ani_plot <- function(outmo,nfrs=25,moname="mymovie.gif",labs,att=TRUE,pca=TRUE,contrib,zoom){

  oopt <- animation::ani.options(interval = 0.1)
  if(att==TRUE){
    if(zoom==FALSE){
      xr=range(unlist(lapply(1:outmo$nchunk,function(jj) {outmo[[jj]]$xr })))
      yr=range(unlist(lapply(1:outmo$nchunk,function(jj) {outmo[[jj]]$yr })))
    }else{xyr=active_zoom(outmo,nframes=nfrs,att=T)}
  }else{
    if(zoom==FALSE){
      xr=range(unlist(lapply(1:outmo$nchunk,function(jj) {outmo[[jj]]$uxr })))
      yr=range(unlist(lapply(1:outmo$nchunk,function(jj) {outmo[[jj]]$uyr })))  
    }else{xyr=active_zoom(outmo,nframes=nfrs,att=F)}
  }
  
  FUN2 <- function(nfrs,nchunk,xr,yr,contrib) {
    for (chu in 1:nchunk){
      if(zoom==FALSE){
        lapply(seq(1,nfrs, by = 1), function(i) {
          plot_fun(outmo,chu,i,xr,yr,lab=labs,att=att,pca=pca,contrib=contrib)
          animation::ani.pause()})
      }else{
        lapply(seq(1,nfrs, by = 1), function(i) {
          plot_fun(outmo,chu,i,xyr[[chu]]$xr[i,],xyr[[chu]]$yr[i,],lab=labs,att=att,pca=pca,contrib=contrib)
          animation::ani.pause()})
        
      }
    }
  } 
  
  if (file.exists(moname)) file.remove(moname)
 # saveGIF(FUN2(nfrs,outmo$nchunk,xr,yr,contrib), interval = 0.1,movie.name=moname)
  saveLatex(FUN2(nfrs,outmo$nchunk,xr,yr,contrib), interval = 0.1,movie.name=moname)
}


