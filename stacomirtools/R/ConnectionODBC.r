# Nom fichier :        ConnectionODBC (classe)
# Organisme :          IAV
# Auteur :             Cedric Briand
# Contact :            cedric.briand"at"eptb-vilaine.fr
# Date de creation :   06/02/2007 10:58:37
# Etat :               OK
# Description          Classe de connexion a la base de donnee
#**********************************************************************

#validation function
validite_ODBC=function(object)
{
	rep1= class(object@baseODBC[1])=="Character"
	rep2=class(object@baseODBC[2])=="Character"
	rep3=class(object@baseODBC[3])=="ANY"
	rep4=length(object@baseODBC)==3
	return(ifelse(rep1 & rep2 & rep3 & rep4,TRUE,c(1:4)[!c(rep1, rep2, rep3, rep4)]))
}

#' @title ConnectionODBC class 
#' @note Mother class for connection, opens the connection but does not shut it
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
#' @slot baseODBC="vector" (of length 3, character)
#' @slot silent="logical"
#' @slot etat="ANY" # can be -1 or string
#' @slot connection="ANY" # could be both string or S3
#' @slot sql="character"
#' @slot query="data.frame"
#' @return connectionODBC an S4 object of class connectionODBC
#' @examples 
#' objet=new("ConnectionODBC")
#' objet@baseODBC=c("myODBCconnection","myusername","mypassword")
#' objet@silent=FALSE
#' objet<-connect(objet)
#' odbcClose(objet@connection)
setClass(Class="ConnectionODBC",
		representation= representation(baseODBC="vector",silent="logical",etat="ANY",connection="ANY"),
		prototype = list(silent=TRUE),
		validity=validite_ODBC)
setGeneric("connect",def=function(objet,...) standardGeneric("connect"))

#' connect method for ConnectionODBC class
#' @returnType ConnectionODBC S4 object
#' @return a connection with slot filled
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
#' @examples objet=new("ConnectionODBC")
#' objet@baseODBC=baseODBC
#' connect(objet)
setMethod("connect",signature=signature("ConnectionODBC"),definition=function(objet) {     
			if (length(objet@baseODBC)!=3)  {
				if (exists("baseODBC",envir=.GlobalEnv)){ 
					objet@baseODBC<-get("baseODBC",envir=.GlobalEnv) 
				} else {
					if(exists("envir_stacomi")){# the program is called within stacomiR
						funout(get("msg",envir_stacomi)$ConnectionODBC.1,arret=TRUE)
					} else	  {
						stop("you need to define a vector baseODBC with the ODBC link, user and password")
					}
				}
			}
			e=expression(channel <-odbcConnect(objet@baseODBC[1],
							uid = objet@baseODBC[2],
							pwd = objet@baseODBC[3],
							case = "tolower",
							believeNRows = FALSE))
			if (!exists("odbcConnect")) {
				if(exists("envir_stacomi")){
					funout(get("msg",envir_stacomi)$ConnectionODBC.2,arret=TRUE)
				} else	  {
					stop("the RODBC library is necessary, please load the package")
				}
			}
			if (!objet@silent) {
				if(exists("envir_stacomi")){
					print(paste(get("msg",envir_stacomi)$ConnectionODBC.3,objet@baseODBC[1]))
				} else {
					print(paste("connection trial, warning this class should only be used for test: ",objet@baseODBC[1]))
				}
			}	
			# sends the result of a trycatch connection in the
			#l'object (current connection), e.g. a character vector
			connection_error<-if(exists("envir_stacomi")){
						error=paste(get("msg",envir_stacomi)$ConnectionODBC.4,objet@baseODBC[1])
					} else {
						error="impossible connection"
					}
			currentConnection<-tryCatch(eval(e), error=connection_error) 
			if (class(currentConnection)=="RODBC") {
				if (!objet@silent){
					if(exists("envir_stacomi")){
						print(get("msg",envir_stacomi)$ConnectionODBC.5)
					} else {
						print("connection successful")
					}
				} 
				objet@connection=currentConnection  # an objet S3 RODBC
				if(exists("envir_stacomi")){
					state<-get("msg",envir_stacomi)$ConnectionODBC.6
				} else {
					state<-"connected"
				}
				objet@etat=state
			} else {
				funout(currentConnection)
				objet@etat=currentConnection # reporting error
			}
			return(objet)
		})
