pvx.super <- function(theta_v, thres=NULL, dat=NULL){
  # func. by joerg-henrik heine jhheine(at)googlemail.com
  # theta_v: ein vector oder zahl; oder ein pers_obj  
  # thres_m: matrix thurstonian thresholds der items ggf sind NAs drin
  # dat : benutz daten um die jeweilige P der gewählten antwort auszugeben
  # funct. needs pvx.matrix in i.pvx.matrix.R
  
  resp <- NULL
  
  if(any( (class(theta_v)=="pers") & (is.null(thres)) ) ){
    ## wenn nur pers_obj übergeben wird -- OK
    resp  <- theta_v$pair$resp
    thres <- (theta_v$pair$threshold)
    namen <- rownames(theta_v$pair$resp)
    theta_v <- (theta_v$pers$WLE)
    names(theta_v) <- namen
    
  }
  if(any( (class(theta_v)=="pers") & (!is.null(thres)) ) ){
    ## wenn pers_obj übergeben wird  mit separatem threshold 
    resp  <- theta_v$pair$resp
    thres <- thres
    theta_v <- (theta_v$pers$WLE)
  }
  
# hier könnte man noch weiter bauen 21.11.2015  
#   if((class(theta_v)=="numeric") & (length(thres)!=0) ){
#     ## wenn  
#     theta_v <- theta_v
#     thres <- thres
#     if(class(dat)!="logical"){resp <- dat}
#   }
  
  thresL <- lapply(1:nrow(thres), function(i) {na.omit(thres[i,])})
  names(thresL) <- rownames(thres)
  # do.call(cbind , lapply(thresL,function(x){t(pvx.matrix(theta_v,x ))}) )

    
#   if(length(resp)==0){
#     suppressWarnings(erg <- data.frame(lapply(thresL,function(x){t(pvx.matrix(theta_v,x ))}), row.names=names(theta_v),check.rows=F) )
#   } currently not supported

  if(length(resp)!=0){
    respL <- lapply(1:ncol(resp), function(i) {(resp[,i])})
    names(respL) <- colnames(resp)
    #erg <- mapply(function(x,y){t(pvx.matrix(theta_v,x ,y+1))}, x=thresL, y=respL)
    erg <- mapply(function(x,y){t(pvx.matrix(theta_v = theta_v,thres = x ,xm_v = (y+1) ))}, x=thresL, y=respL)
  }

  return(erg)
}

