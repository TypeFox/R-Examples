setGeneric("equalProj", function(x){standardGeneric("equalProj")})
setMethod("equalProj", 
          signature=c(x="list"),
          definition=function(x){
		  if(all(is.na(unlist(lapply(x, proj4string)))))
		  {
			  return(T)
		  }else{
		  if(length(x)==1){
			  return(T)}else{
				  return(all(unlist(lapply(x[-1], identicalCRS, x[[1]]))))
		  }}
#             if(!isClass(Class="Raster",x)) stop("The list contains objects that are not Rasters.")
 #           tmp <- unlist(strsplit(gsub(pattern="+", replacement="", unlist(strsplit(as.vector(unlist(lapply(x, projection))), " ")), fixed=T), "="))
 #           if(!all(tmp=="NA")){
 #             nom <- as.character(tmp[seq(1,length(tmp),2)])
 #             dat <- data.frame(nom=rep(unique(nom),length(x)), ID=sort(rep(1:length(x), length(unique(nom)))))
 #             vals <- data.frame(val=as.character(tmp[seq(2,length(tmp),2)]), nom=nom)
 #             vals$ID=rep(1:length(x), unlist(lapply(lapply(lapply(as.vector(lapply(x, projection)), strsplit, " "), unlist),length)))
 #             dat <- merge(dat, vals, by=c("nom", "ID"), all=T)
 # 
 #             if(any(is.na(dat$val))){
 #               warning(paste("Projections have different number of arguments (differences in: ",paste(as.character(dat$nom[is.na(dat$val)]), collapse=", "),").", sep=""))}
 #             if(!nrow(unique(dat[complete.cases(dat), c("val", "nom")])) == length(unique(nom))){
 #               return(FALSE)} else {return(TRUE)}
 #             } else {return(TRUE)}
            })
