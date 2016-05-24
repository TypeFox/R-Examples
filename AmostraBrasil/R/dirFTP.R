#' List Files in FTP Site
#' @param myURL A string
#' @param codibge A string
#' @return A dataframe with IBGE files for one municipality
#' @export
dirFTP <- function (myURL="", codibge=""){
  ibgecode <- getURL(myURL,ftp.use.epsv=F,dirlistonly=T)
  ibgecode <- paste(strsplit(ibgecode, "\r*\n")[[1]], sep = "")
  listaIBGE <- as.data.frame(ibgecode)
  listaIBGE <- subset(listaIBGE,codibge == substring(listaIBGE$ibgecode,1,7))
  return(listaIBGE)
}
