#' Generates Brazil's IBGE Householding Sample
#'
#' @param codibge A string
#' @param municipio A string
#' @param listasetor A string vector
#' @param N A number
#' @param geocod A boolean
#' @param shape A boolean
#' @param saveall A boolean
#' @param dados A null
#' @return A dataframe with N addresses sampled from IBGE census
#' @examples
#' amostraBrasil(municipio="Pindoba", N=5, geocod=TRUE, shape=TRUE)
#'
#' @import stringr
#' @import rgdal
#' @import RJSONIO
#' @import foreign
#' @import RCurl
#' @import sp
#' @importFrom methods as
#' @importFrom utils download.file memory.limit read.table unzip
#' @importFrom sp over spTransform
#' @export
amostraBrasil <- function(codibge="",municipio="",listasetor=c(),N=0,geocod=F,shape=F,saveall=T,dados=NULL)
{
  stopifnot(require(AmostraBrasil, quietly=TRUE))
  Sys.setlocale("LC_CTYPE","Portuguese_Brazil.1252")
  memory.limit(size = 4095)
  Encoding(MUNICIPIOS.IBGE$MUNICIPIO)<-"latin1"
  myMUN<-subset(MUNICIPIOS.IBGE,MUNICIPIOS.IBGE$CODIBGE==codibge | MUNICIPIOS.IBGE$MUNICIPIO==municipio)
  # testar length myMUN -> se maior 1 msg usuario
  xmun<-1
  if (nrow(myMUN)==0){
    message(iconv("Cidade n\u00E3o encontrada! Verifique o nome ou tente o c\u00F3digo do munic\u00EDpio no IBGE","UTF-8","latin1"))
    return (0)
  }
  if (nrow(myMUN)>1){
    message(paste("Foram encontrados mais de um munic\u00EDpio chamados ",municipio," (exibindo UF e c\u00F3digo do IBGE):",sep=""))
    for (i in 1:nrow(myMUN)){
      message(paste(i,") ",myMUN$UF[i], " - ",myMUN$CODIBGE[i],sep=""))
    }
    xmun<-readline(iconv("Qual voc\u00EA deseja ? (tecle o n\u00FAmero ou 0 para sair e ENTER) ?","UTF-8","latin1"))
    xmun<-as.numeric(unlist(strsplit(xmun, ",")))
  }
  if (xmun==0) {
    message(iconv("Opera\u00E7\u00E3o cancelada!","UTF-8","latin1"))
    return (0)
  }
  if (xmun>nrow(myMUN)) {
    message(iconv("N\u00FAmero inv\u00E1lido! opera\u00E7\u00E3o cancelada!","UTF-8","latin1"))
    return (0)
  }
  codibge<-myMUN$CODIBGE[xmun]
  municipio<-myMUN$MUNICIPIO[xmun]
  uf<-myMUN$UF[xmun]


  if (class(dados)!="data.frame") {
  message(paste("Aguarde! Conectando IBGE... Buscando dados de ",municipio,"/",uf," (",codibge,")" ,sep=""))
  message("")
  myURL<-paste("ftp://ftp.ibge.gov.br/Censos/Censo_Demografico_2010/Cadastro_Nacional_de_Enderecos_Fins_Estatisticos/",uf,"/",sep="")
  mylist<-dirFTP(myURL=myURL,codibge=codibge)
  options(stringAsFactors=FALSE)
  n<-0
  narqs <-nrow(mylist)
  #narqs<-3
  i<-1
  while (i < narqs+1) {
    temp <- paste(getwd(),"/zipIBGE.ZIP",sep="")             ##tempfile()
    MYfilename<-paste(myURL,mylist$ibgecode[i],sep="")
    message(paste("\nAcessando arquivo ",mylist$ibgecode[i]," (",i,"/",narqs,")",sep=""))
    options(timeout=60)
    sair=9
    while (sair!=0){
      if (try(download.file(MYfilename,temp,quiet=T))!=0){
        message(iconv("\nErro na conex\u00E3o com IBGE... Tentando novamente... ","UTF-8","latin1"))
        sair<-readline("Sair ? (tecle 0 para sair e ENTER)")
        sair<-as.numeric(unlist(strsplit(sair, ",")))
      } else {
        sair<-0
      }
    }
    #message("abrindo zip")

    dadost <- try(read.table(unz(temp, paste(substr(mylist$ibgecode[i],1,11),".TXT",sep="")),header=F,sep="%",as.is=T,stringsAsFactors = F))
    unlink(temp)
    rm(temp)
    if (class(dadost)=="try-error"){
      message(iconv("Erro na descompacta\u00E7\u00E3o do arquivo. Tentando novamente...","UTF-8","latin1"))
    } else {
#      message("Acrescentando registros do arquivo ",mylist$ibgecode[i],sep="")
      Encoding(dadost$V1) <- "latin1"
      dadost<-data.frame(setor=paste(gsub(" ","0",substr(dadost$V1,8,9)),gsub(" ","0",substr(dadost$V1,10,11)),gsub(" ","0",substr(dadost$V1,12,15)),sep=""),
                        urbrur=substr(dadost$V1,16,16),
                        espend=substr(dadost$V1,472,473),
#                        desend=substr(dadost$V1,474,500),
                        endIBGE=paste(str_trim(substr(dadost$V1,17,36))," ",ifelse(str_trim(substr(dadost$V1,37,66))=="","",paste(str_trim(substr(dadost$V1,37,66))," ",sep="")),str_trim(substr(dadost$V1,67,126)),", ",str_trim(substr(dadost$V1,127,134))," ",ifelse(codibge=="3550308" & dadost$V1>"",paste("0",substr(dadost$V1,551,554),sep=""),substr(dadost$V1,551,555)),"-",ifelse(codibge=="3550308" & dadost$V1>"",substr(dadost$V1,555,557),substr(dadost$V1,556,558)), sep="")
      )
      dadost<-subset(dadost,dadost$espend=="01" | dadost$espend=="02")
#      dadost<-subset(dadost,dadost$espend=="06" & (grep("BAR ",dadost$desend)>0 | grep("BUTECO ",dadost$desend)>0))
#      dadost<-subset(dadost,dadost$espend=="06" & (substr(dadost$desend,1,4)=="BAR " | substr(dadost$desend,1,7)=="BUTECO "))
      nrdadost<-nrow(dadost)
      n<-n+nrdadost
      message("Acrescentados ",nrdadost,iconv(" domic\u00EDlios (total: ","UTF-8","latin1"), n,")",sep="")
      if (i==1){
        dados<-dadost
        }
      else {
      dados<-rbind(dados,dadost)
#      message(memory.size())
#      message(object.size(dados)/1024)
#      message(nrow(dados))
      }
      i<-i+1
    }
  }

  message("\nConfigurando o arquivo...",municipio)
  dados$endIBGE<-paste(dados$endIBGE,", ",municipio," - ",uf,", Brasil",sep="")

  if (saveall==T){
    write.dbf(dados,file=paste(municipio,"_domicilios.dbf",sep=""))
  }
  message("\nGerados ",n, " registros com sucesso!",sep="")




  if (length(listasetor)>0) {
    dados<-subset(dados, dados$setor %in% listasetor)
  }

  #testar amostra se amostra >0
  if (N>0){
    csample<-sample(seq(1:nrow(dados)),N,replace=F)
    dados<-dados[csample,]
    n<-N
  }


  }

  #testar geocod<>F
  if (geocod==TRUE){
    message(iconv("\nAguarde... geocodificando endere\u00E7os via Google Maps.","UTF-8","latin1"))
    #message(dados$endIBGE)
    n<-nrow(dados)
    dados<-setLatLong(dados,dados$endIBGE)
    write.dbf(dados,file=paste(municipio,"_geo.dbf",sep=""))
    dados<-subset(dados,dados$Status=="OK")
    message(paste("\nForam descartados ",n-nrow(dados), iconv(" registros n\u00E3o geocodificados!","UTF-8","latin1"),sep=""))
    if (nrow(dados)==0) {
      message("Nenhum registro geocodificado! Verifique o arquivo ",municipio,"_geo.dbf",sep="")
    }
    if (shape==TRUE & nrow(dados)>0) {
      lmts<-getIBGEMunSHP(codibge = codibge)
      LLcoor<-data.frame(Lng=dados$Lng,Lat=dados$Lat)
      coordinates(LLcoor)=~Lng+Lat
      proj4string(LLcoor)=CRS("++proj=longlat +ellps=WGS84") # set it to dec Long Lat
      dados<-SpatialPointsDataFrame(LLcoor,data.frame(SC=dados$setor,ENDIBGE=dados$endIBGE,UR=dados$urbrur,TIPO=dados$espend))
      dados<-setDentroFora(dados,lmts)
      writeOGR(dados, dsn="." ,layer=paste(municipio,"_pts",sep=""),driver="ESRI Shapefile", overwrite_layer = TRUE)
      message(paste("\nForam criados o arquivo ",municipio, "_pts.shp e seus componentes em ",getwd() ,sep=""))
    }
  }
  return(dados)
}
