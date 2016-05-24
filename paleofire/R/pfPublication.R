#' Get citations for charcoal sites
#'
#' @author O. Blarquez
#' @param x A "pfSiteSel" object
#' @param output Defines the output as a "list" or a "data.frame" (default).
#' @return A list or data frame with citation informations related to charcoal sites.
#' @examples
#' x=pfSiteSel(id_site %in% c(1:4))
#' pfPublication(x,output="list")


pfPublication <- function(x,output="data.frame") {
  
  pub_key<-pub<-NULL
  data(pub,envir = environment())
  data(pub_key,envir = environment())
  

  
  ## data.frame output
  if(output=="data.frame"){ 
    
    a=b=id_pub_list=list()
    for(i in 1:length(x$id_site)){
      id_pub_list[[i]]=pub_key[pub_key[,1]==x$id_site[i],2]
      a[[i]]=rep(x$id_site[i],length(pub_key[pub_key[,1]==x$id_site[i],2]))
      b[[i]]=rep(x$site_name[i],length(pub_key[pub_key[,1]==x$id_site[i],2]))
    }
    
    c=cbind(unlist(a),unlist(b),unlist(id_pub_list))
    for(i in 1:nrow(c))
      c[i,3]=c(as.character(pub[pub[,1]== c[i,3],2]))
    
    colnames(c)=c("id_site","site_name","citation")
    pub_frame=data.frame(c)
    return(pub_frame)
    
  } else {
    ## List output
    pub_list=NULL
    for(i in 1:length(x$id_site)){
      pub_list[[i]]=list(id_site=x$id_site[i],
                         site_name=x$site_name[i],
                         citation=as.character(
                           pub[pub[,1] %in% pub_key[pub_key[,1]==x$id_site[i],2],2]) )
    }
    return(pub_list)
  }
  
}
