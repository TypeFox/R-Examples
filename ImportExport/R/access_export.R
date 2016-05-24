access_export <- function(file, x, tablename = as.character(1:length(x)), uid = "", pwd = "", ...){
  mycon <- RODBC::odbcConnectAccess(file, uid = uid, pwd = pwd, ...)
  if(class(x)=="data.frame")
    RODBC::sqlSave(mycon, x, tablename=tablename[1], ...)
  else
    for(i in 1:length(x))
      RODBC::sqlSave(mycon,x[[i]],tablename=tablename[i],...)
  close(mycon)
}
