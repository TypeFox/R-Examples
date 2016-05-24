#' Get a Shapefile From IBGE Website
#'
#' @param codibge A string
#' @param nomemun A string
#' @return Shapefile of Municipality in Brazil
#' @export
getIBGEMunSHP <- function(codibge="",nomemun="")
{
  Encoding(MUNICIPIOS.IBGE$MUNICIPIO)<-"latin1"
  myMUN<-subset(MUNICIPIOS.IBGE,MUNICIPIOS.IBGE$CODIBGE==codibge | MUNICIPIOS.IBGE$MUNICIPIO==nomemun)
  # testar length myMUN -> se maior 1 msg usuario
  codibge<-myMUN$CODIBGE
  nomemun<-myMUN$MUNICIPIO
  uf<-myMUN$UF
  ufCode<-substring(codibge,1,2)
  temp <- tempfile()
  #testar length listasetor
  ##ftp://geoftp.ibge.gov.br/malhas_digitais/municipio_2013/SP/
  MYfilename<-paste("ftp://geoftp.ibge.gov.br/malhas_digitais/municipio_2013/",uf,"/",tolower(myMUN$UF),"_municipios.zip",sep="")
  message("")
  message(iconv("Buscando mapa do munic\u00EDpio (shp file)","UTF-8","latin1"))
  download.file(MYfilename,temp)
  unzip(temp)
  ## dados <- read.table(unz(temp, paste(codibge,"0500.TXT",sep="")),head=F,sep="%",as.is=T)
  unlink(temp)
  UFSHP<-readOGR(dsn=".",layer=paste(ufCode,"MUE250GC_SIR",sep=""),verbose = F)
  MUNSHP <- UFSHP[UFSHP$CD_GEOCMU == codibge, ]
  MUNSHP <- spTransform(MUNSHP, CRS("++proj=longlat +ellps=WGS84"))
  message(paste(iconv("Criando arquivo SHP do munic\u00EDpio de ","UTF-8","latin1"),nomemun,sep=""))
  writeOGR(MUNSHP, dsn = ".", layer = paste(nomemun,"_lmts",sep=""),driver = "ESRI Shapefile", overwrite_layer = TRUE)
  return(MUNSHP)
}
