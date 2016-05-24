
############################################################

setClass(
		Class="menuitem", 
		representation(
			nombre= "character",
			previo = "function",
			ayuda = "function",
			opciones = "character",
			acciones = "list",
			fi= "logical"
		)
)

##### Get
setMethod(
	f="[",
	signature=c("menuitem","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			nombre		= return(x@nombre),
			previo		= return(x@previo),
			ayuda 		= return(x@ayuda),
			opciones 	= return(x@opciones),
			acciones 	= return(x@acciones),
			fi 			= return(x@fi)
		)
	}
)

###### Set
setMethod(
	f="[<-",
	signature=c("menuitem","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			nombre		= x@nombre <- value,
			previo		= x@previo <- value,
			ayuda 		= x@ayuda <- value,
			opciones 	= x@opciones <- value,
			acciones 	= x@acciones <- value,
			fi 			= x@fi <- value
		)
		return(x)
	}
)

############################################################

setClass(
		Class="menu", 
		representation(
			items = "list"
		)
)

##### Get
setMethod(
	f="[",
	signature=c("menu","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			items	= return(x@items)
		)
	}
)

###### Set
setMethod(
	f="[<-",
	signature=c("menu","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			items = x@items <- value
		)
		return(x)
	}
)

############################################################

setClass(
		Class="report", 
		representation(
			report 		= "logical",
			graph 		= "numeric",
			trueVal 	= "logical",
			actRep 		= "function",
			desGraph 	= "function",
			desRep 		= "function",
			commentt	= "logical",
			files		= "logical",
			contRep 	= "list"
		)
)

##### Get
setMethod(
	f="[",
	signature=c("report","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			report		= return(x@report),
			graph		= return(x@graph),
			trueVal		= return(x@trueVal),
			actRep		= return(x@actRep),
			desGraph	= return(x@desGraph),
			desRep		= return(x@desRep),
			commentt	= return(x@commentt),
			files		= return(x@files),
			contRep		= return(x@contRep)
		)
	}
)

###### Set
setMethod(
	f="[<-",
	signature=c("report","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			report		= x@report <- value,
			graph		= x@graph <- value,
			trueVal		= x@trueVal <- value,
			actRep		= x@actRep <- value,
			desGraph	= x@desGraph <- value,
			desRep		= x@desRep <- value,
			commentt	= x@commentt <- value,
			files		= x@files <- value
		)
		return(x)
	}
)

############################################################

setClass(
		Class="cami", 
		representation(
			cami = "character"
		)
)

##### Get
setMethod(
	f="[",
	signature=c("cami","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			cami	= return(x@cami)
		)
	}
)

###### Set
setMethod(
	f="[<-",
	signature=c("cami","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			cami	= x@cami <- value,
			
		)
		return(x)
	}
)

############################################################

setClass(
		Class="lineal", 
		representation(
			crit = "numeric",
			ls = "logical",
			atip = "ANY"
		)
)

##### Get
setMethod(
	f="[",
	signature=c("lineal","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			crit		= return(x@crit),
			ls		= return(x@ls),
			atip 		= return(x@atip)
		)
	}
)

###### Set
setMethod(
	f="[<-",
	signature=c("lineal","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			crit		= x@crit <- value,
			ls		= x@ls <- value,
			atip 		= x@atip <- value
		)
		return(x)
	}
)

############################################################

setClass(
		Class="serie",
		representation(
			nom = "character",
			serie = "ts",
			orig = "numeric",
			sact = "logical",
			trans = "numeric",
			est = "numeric",
			reg = "numeric",
			stac = "logical",
			lin = "lineal"
		)
)

##### Get
setMethod(
	f="[",
	signature=c("serie","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			nom		= return(x@nom),
			serie	= return(x@serie),
			orig 	= return(x@orig),
			sact 	= return(x@sact),
			trans 	= return(x@trans),
			est 	= return(x@est),
			reg 	= return(x@reg),
			stac 	= return(x@stac),
			lin 	= return(x@lin)
			
		)
	}
)

###### Set
setMethod(
	f="[<-",
	signature=c("serie","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			nom		= x@nom <- value,
			serie	= x@serie <- value,
			orig 	= x@orig <- value,
			sact 	= x@sact <- value,
			trans 	= x@trans <- value,
			est 	= x@est <- value,
			reg 	= x@reg <- value,
			stac 	= x@stac <- value,
			lin 	= x@lin <- value
		)
		return(x)
	}
)

############################################################

setClass(
		Class="modelo", 
		representation(
			modelo = "ANY",
			mact = "logical",
			ser = "numeric",
			int = "logical",
			valid = "logical",
			est = "logical",
			eqm = "numeric",
			best = "logical"
		)
)

##### Get
setMethod(
	f="[",
	signature=c("modelo","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			modelo	= return(x@modelo),
			mact	= return(x@mact),
			ser 	= return(x@ser),
			int 	= return(x@int),
			valid 	= return(x@valid),
			est 	= return(x@est),
			eqm 	= return(x@eqm),
			best 	= return(x@best)
		)
	}
)

###### Set
setMethod(
	f="[<-",
	signature=c("modelo","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			modelo	= x@modelo <- value,
			mact	= x@mact <- value,
			ser 	= x@ser <- value,
			int 	= x@int <- value,
			valid 	= x@valid <- value,
			est 	= x@est <- value,
			eqm 	= x@eqm <- value,
			best 	= x@best <- value
		)
		return(x)
	}
)

############################################################

setClass(
		Class="datos",
		representation(
			lserie = "list",
			sident = "numeric",
			lmodelo = "list",
			mident = "numeric",
			modif = "logical"
			)
)

##### Get
setMethod(
	f="[",
	signature=c("datos","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			lserie	= return(x@lserie),
			sident	= return(x@sident),
			lmodelo	= return(x@lmodelo),
			mident 	= return(x@mident),
			modif 	= return(x@modif)
		)
	}
)

###### Set
setMethod(
	f="[<-",
	signature=c("datos","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			lserie	= x@lserie <- value,
			sident	= x@sident <- value,
			lmodelo	= x@lmodelo <- value,
			mident 	= x@mident <- value,
			modif 	= x@modif <- value
		)
		return(x)
	}
)
	
############################################################	

setClass(
		Class="sessio.ts",
        representation(
			menu = "menu",
			cami = "cami",
			datos = "datos",
			student = "logical",
			report = "report"
		)
)

##### Get
setMethod(
	f="[",
	signature=c("sessio.ts","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			menu		= return(x@menu),
			cami		= return(x@cami),
			datos 		= return(x@datos),
			student 	= return(x@student),
			report		= return(x@report)
		)
	}
)

###### Set
setMethod(
	f="[<-",
	signature=c("sessio.ts","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			menu		= x@menu <- value,
			cami		= x@cami <- value,
			datos 		= x@ayuda <- value,
			student 	= x@student <- value,
			report 		= x@report <- value
		)
		return(x)
	}
)
	
############################################################
	
setClass(
		Class="tex",
        representation(
			menu = "matrix",
			opciones = "matrix",
			previo = "matrix",
			ayuda = "matrix",
			ini  = "matrix",
			gesser = "matrix",
			trans = "matrix",
			iden = "matrix",
			gesmod = "matrix",
			estim = "matrix",
			valid = "matrix",
			cap = "matrix",
			atip = "matrix",
			prev = "matrix",
			salir = "matrix",
			plot = "matrix",
			drawser = "matrix",
			drawmod = "matrix",
			writelistser = "matrix",
			writemod = "matrix",
			writelistmod = "matrix",
			writecoef = "matrix",
			drawarma = "matrix",
			atipics = "matrix",
			drawatip = "matrix",
			previsiones = "matrix",
			enterComment = "matrix",
			latexCntrl = "matrix"
		)
)

###### Get
setMethod("[",
	signature=c("tex","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			menu		= return(x@menu),
			opciones	= return(x@opciones),
			previo		= return(x@previo),
			ayuda		= return(x@ayuda),
			ini			= return(x@ini),
			gesser		= return(x@gesser),
			trans		= return(x@trans),
			iden		= return(x@iden),
			gesmod		= return(x@gesmod),
			estim		= return(x@estim),
			valid		= return(x@valid),
			cap			= return(x@cap),
			atip		= return(x@atip),
			prev		= return(x@prev),
			salir		= return(x@salir),
			plot		= return(x@plot),
			drawser		= return(x@drawser),
			drawmod		= return(x@drawmod),
			writelistser	= return(x@writelistser),
			writemod	= return(x@writemod),
			writelistmod= return(x@writelistmod),
			writecoef	= return(x@writecoef),
			drawarma	= return(x@drawarma),
			atipics		= return(x@atipics),
			drawatip	= return(x@drawatip),
			previsiones = return(x@previsiones),
			enterComment = return(x@enterComment),
			latexCntrl = return(x@latexCntrl)
		)
	}
)

###### Set
setMethod("[<-",
	signature=c("tex","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			menu		=x@menu <- value,
			opciones	= x@opciones <- value,
			previo		= x@previo <- value,
			ayuda		= x@ayuda <- value,
			ini			= x@ini <- value,
			gesser		= x@gesser <- value,
			trans		= x@trans <- value,
			iden		= x@iden <- value,
			gesmod		= x@gesmod <- value,
			estim		= x@estim <- value,
			valid		= x@valid <- value,
			cap			= x@cap <- value,
			atip		= x@atip <- value,
			prev		= x@prev <- value,
			salir		= x@salir <- value,
			plot		= x@plot <- value,
			drawser		= x@drawser <- value,
			drawmod		= x@drawmod <- value,
			writelistser	= x@writelistser<- value,
			writemod	= x@writemod <- value,
			writelistmod= x@writelistmod <- value,
			writecoef	= x@writecoef <- value,
			drawarma	= x@drawarma <- value,
			atipics		= x@atipics <- value,
			drawatip	= x@drawatip <- value,
			previsiones = x@previsiones <- value,
			enterComment = x@enterComment <- value,
			latexCntrl = x@latexCntrl <- value
		)
		return(x)
	}
)
	
setGeneric("makeMenu", function(x, opciones, ...) standardGeneric("makeMenu"))
setGeneric("initializerRep", function(report, contRep, ...) standardGeneric("initializerRep"))
setGeneric("initializer", function(object, student, report, contRep, ...) standardGeneric("initializer"))
setGeneric("ender", function(session, ...) standardGeneric("ender"))
setGeneric("activate", function(x, ...) standardGeneric("activate"))
setGeneric("addcami", function(x, ...) standardGeneric("addcami"))
setGeneric("getcami", function(x, ...) standardGeneric("getcami"))
setGeneric("nextGraphic", function(x, name, twoplot, pos, contRep, ...) standardGeneric("nextGraphic"))
setGeneric("getserie", function(x, ...) standardGeneric("getserie"))
setGeneric("addserie", function(x, ...) standardGeneric("addserie"))
setGeneric("getsident", function(x, ...) standardGeneric("getsident"))
setGeneric("addsident", function(x, ...) standardGeneric("addsident"))
setGeneric("getmodelo", function(x, ...) standardGeneric("getmodelo"))
setGeneric("addmodelo", function(x, ...) standardGeneric("addmodelo"))
setGeneric("getmident", function(x, ...) standardGeneric("getmident"))
setGeneric("addmident", function(x, ...) standardGeneric("addmident"))
setGeneric("getmodif", function(x, ...) standardGeneric("getmodif"))
setGeneric("addmodif", function(x, ...) standardGeneric("addmodif"))
setGeneric("resumen", function(x, ...) standardGeneric("resumen"))
setGeneric("texinitializer", function(x, ...) standardGeneric("texinitializer"))



#######Funciones de realizacion del informe#######
setMethod("nextGraphic",signature(x = "numeric",name="character",twoplot="ANY",pos="ANY", contRep="list"),
          function(x,name,twoplot,pos,contRep,...)
		{
	if(missing(twoplot)) twoplot=F
    x <- x+1
    aux <-paste("Grafh",x,".jpg", sep="")
    jpeg(filename=aux)
	if(twoplot){
		if(pos==1){
			cat("\n\\begin{figure}[H]")
			cat("\n\\centering")
			cat("\n\\captionsetup{listformat=empty}")
			
			cat(paste("\n\\subfloat[","*",name,"]{\\includegraphics[width=",contRep$twograph,"\\linewidth]{",aux,"}}",sep=""))
		}else{
			cat(paste("\n\\subfloat[","*",name,"]{\\includegraphics[width=",contRep$twograph,"\\linewidth]{",aux,"}}",sep=""))			
			cat("\n\\end{figure}\n")
		}
	}else{
		cat("\n\\begin{figure}[H]")
		cat("\n\\centering")
		cat(paste("\n\\includegraphics[width=",contRep$maingraph,"\\linewidth]{",aux,"}\n",sep=""))
		cat(paste("\n\\caption{",name,"}",sep=""))
		cat("\n\\end{figure}\n")	
	}
	return(x)

  }
)
###########

setMethod("addcami",signature(x = "cami"),
          function(x, id=NULL,...)
		{
		if(is.null(id)){
			length(x@cami)=length(x@cami)-1
		}else{
			x@cami=c(x@cami,id)
		}
		x
	}
)


setMethod("getcami",signature(x = "cami"),
          function(x,...)
		{
		a=x@cami[length(x@cami)]
		a
	}
)

setMethod("getserie", signature(x = "datos"),
          function(x, id=NULL,...)
		{
		if(is.null(id)){
			x@lserie[[x@sident]]
		}else{
			x@lserie[[id]]
		}
	}
)

setMethod("addserie", signature(x = "datos"),
          function(x, s,id=NULL,...)
		{
		if(is.null(id)){
			x@lserie[[length(x@lserie)+1]]=s
			x=addsident(x,length(x@lserie))
		}else{
			x@lserie[[id]]=s
		}
		x
	}
)



setMethod("getsident", signature(x = "datos"),
          function(x,...)
		{
		x@sident
	}
)

setMethod("addsident", signature(x = "datos"),
          function(x, a,...)
		{
		x@sident=a
		x
	}
)


setMethod("getmodelo", signature(x = "datos"),
          function(x, id=NULL,...)
		{
		if(is.null(id)){
			x@lmodelo[[x@mident]]
		}else{
			x@lmodelo[[id]]
		}
	}
)

setMethod("addmodelo", signature(x = "datos"),
          function(x, m, id=NULL,...)
		{
		if(is.null(id)){
			if (x@mident == 0) {
				x@lmodelo[[1]]=m
			} else {
				x@lmodelo[[length(x@lmodelo)+1]]=m
			}
			x=addmident(x,length(x@lmodelo))
		}else{
			x@lmodelo[[id]]=m
		}
		x
	}
)


setMethod("getmident", signature(x = "datos"),
          function(x,...)
		{
		x@mident
	}
)

setMethod("addmident", signature(x = "datos"),
          function(x, a,...)
		{
		x@mident=a
		x
	}
)


setMethod("getmodif", signature(x = "datos"),
          function(x,...)
		{
		x@modif
	}
)

setMethod("addmodif", signature(x = "datos"),
          function(x, a,...)
		{
		x@modif=a
		x
	}
)


setMethod("resumen", signature(x = "datos"),
          function(x, a,...)
		{
		for(n in 1:length(x@lserie)){
			s=getserie(x,n)
			s@sact=T
			x=addserie(x,s,n)

		}
		if(getmident(x) != 0){
			for(n in 1:length(x@lmodelo)){
				s=getmodelo(x,n)
				s@mact=T
				x=addmodelo(x,s,n)
			}
		}
		x
	}
)








slotcreator=function(tex,pos,nlines,name){
	aux=matrix(NA,nrow=max(nlines),ncol=length(pos),dimnames=list(NULL,name))
	for(i in 1:length(pos)){
		aux[,i]=c(tex[(pos[i]+1):(pos[i]+nlines[i])],rep(NA,max(nlines)-nlines[i]))
	}
	return(aux)
}


findtext=function(keynames,tex){
	m=matrix(c(keynames,rep(NA,length(keynames)*2)),nrow=length(keynames),ncol=3)
	j=1
	for(i in 1:(length(keynames)-1)){
		find=F
		while(!find){
			if(m[i,1]==tex[j]) find=T
			j=j+1
		}
		m[i,2]=j-1
	}
	m[length(keynames),2]=length(tex)
	aux1=as.numeric(m[,2][-1])
	aux2=as.numeric(m[,2][-length(m[,2])])
	a=aux1-aux2-1
	a[length(a)+1]=0
	m[,3]=a
	return(m)
}


setMethod(	"texinitializer",
			signature="tex",
			function(x, ...){
				# for(i in 1:length(.libPaths())){
					# fil=try(file(paste(.libPaths()[i],"/TSTutorial/","english",".txt",sep="")),silent=T)
					# af=suppressWarnings(try(readLines(fil),silent=T))
					# if(!is(af,"try-error")) tex=af
					# close(fil)
				# }
				tex=englishLanguage(NULL)
				keynames=c("--menu",
							paste("--opciones.",1:10,sep=""),
							paste("--previo.",1:3,sep=""),
							paste("--ayuda.",1:10,sep=""),
							paste("--ini.",1,sep=""),
							paste("--gesser.",1:2,sep=""),
							paste("--trans.",4:6,sep=""),
							paste("--iden.",2:3,sep=""),
							paste("--gesmod.",1:2,sep=""),
							paste("--estim.",2:4,sep=""),
							paste("--valid.",2:5,sep=""),
							paste("--cap.",2:4,sep=""),
							paste("--atip.",2,sep=""),
							paste("--prev.",2,sep=""),
							"--salir",
							"--plot",
							"--drawser",
							"--drawmod",
							"--writelistser",
							"--writemod",
							"--writelistmod",
							"--writecoef",
							"--drawarma",
							"--atipics",
							"--drawatip",
							"--previsiones",
							"--enterComment",
							"--latexCntrl",
							"--end"
				)
				m=findtext(keynames,tex)
				slotpos=c(1,10,3,10,1,2,3,2,2,3,4,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
				slotname=c("menu","opciones", "previo", "ayuda", "ini", "gesser", "trans", "iden", "gesmod", "estim", "valid", "cap", "atip",
							"prev", "salir","plot", "drawser", "drawmod", "writelistser", "writemod", "writelistmod", "writecoef", "drawarma",
							"atipics", "drawatip", "previsiones", "enterComment", "latexCntrl"
				)
				listname=list(
					c(1),c("ini","trans","gesser","iden","estim","gesmod","valid","cap","atip","prev"),c("previotit1","previotit2","previo"),
					c("ini","gesser","trans","iden","gesmod","estim","valid","cap","atip","prev"),c(1), c(1:2), c(4:6), c(2:3), c(1:2),c(2:4),c(2:5), c(2:4), c(2),
					c(2), c(1), c(1), c(1), c(1), c(1), c(1), c(1), c(1), c(1), c(1), c(1), c(1), c(1), c(1)
				)
				aux=1
				pos=as.numeric(m[,2][aux:slotpos[1]])
				nlines=as.numeric(m[,3][aux:slotpos[1]])
				name=listname[[1]]
				x[slotname[1]]<-slotcreator(tex,pos,nlines,name)
				aux=slotpos[1]+1
				for(i in 2:length(slotpos)){
					pos=as.numeric(m[,2][aux:(aux+slotpos[i]-1)])
					nlines=as.numeric(m[,3][aux:(aux+slotpos[i]-1)])
					name=listname[[i]]
					x[slotname[i]]<-slotcreator(tex,pos,nlines,name)
					aux=aux+slotpos[i]
				}
				return(x)
			}
)







setMethod(	"makeMenu", 
			signature(x = "menu", opciones = "list"),
			function(x, opciones, ...){
			co=1
			for (nom in names(opciones)){
				k=length(opciones[[nom]])
				x@items[[co]]=new(Class="menuitem",
					nombre=nom,
					previo = get(paste(nom,"previo",sep=".")),
					ayuda = get(paste(nom,"ayuda",sep=".")),
					opciones=opciones[[nom]],
					acciones=apply(matrix(paste(nom,1:k,sep=".")),1,function(el) get(el,pos=environment(TSTutorial))),
					fi=FALSE)
				co=co+1
				}
			names(x@items)=names(opciones)
			return(x)
		}
)

latexCntrl=function(...){
	latex=F
	for(i in 1:length(.libPaths())){
			fil=try(file(paste(.libPaths()[i],"/TSTutorial/test/LatexTest.tex",sep="")),silent=T)
			af=suppressWarnings(try(readLines(fil),silent=T))
			if(!is(af,"try-error")) text=af
			close(fil)
	}
	suppressWarnings(try(close(fil),silent=T))
	aux="LatexTest.tex"
	sink(aux)
	for(i in 1:length(text)){
		cat(paste(text[i],"\n",sep=""))
	}
	sink()
	suppressWarnings(sink())					
	# b=suppressWarnings(try(ifelse(test=as.logical(Sys.info()["sysname"] == "Windows"),
		# shell("pdfLatex LatexTest.tex",intern=T),
		# system("pdflatex LatexTest.tex",intern=T)),silent=T))
				
	b=suppressWarnings(try(ifelse(test=as.logical(Sys.info()["sysname"] == "Windows"),
		system("pdfLatex LatexTest.tex",intern=T),
		system("pdflatex LatexTest.tex",intern=T)),silent=T))
	div=apply(as.matrix(rep(1:nchar(b[1]))),1,function(el) substr(b[1],el,el))
	if(paste(div[1],div[2],div[3],div[4],sep="")=="This") latex=T
	suppressWarnings(try(unlink(paste("LatexTest",c(".pdf",".aux",".log",".tex"),sep="")),silent=T))
	if(!latex){
		prettyprint(tex["latexCntrl"][1,"1"])
		cat("\n")
		return(F)
	}else{
		return(T)
	}
}

reportCntrl=function(report=T, comment=F, files=F){
	stopifnot(is.logical(report),is.logical(comment),is.logical(files))
	if(report){
			report=latexCntrl()
	}
	if(!report){
		comment=F
		files=F
	}
	list(
		report = as.logical(report),
		comment = as.logical(comment),
		files = as.logical(files)
	)
}

contRepCntrl=function(fil=F, tex=NULL, name=NULL, maingraph=0.7, twograph=0.47){
	stopifnot(is.logical(fil),is.numeric(maingraph),is.numeric(twograph))
	if((fil)&(is.null(tex))) stop("tex needs to have a file .tex conection")
	list(
		fil = fil,
		tex = tex,
		name = name,
		maingraph = maingraph,
		twograph = twograph
	)
}

nameser=function(name){
	print(name)
	if(substr(name, start=11, stop=17)=="series"){
		cont=1
		nam=""
		while((substr(name, start=17+cont, stop=17+cont)!=",")|(substr(name, start=11+cont, stop=11+cont)!=")")){
			print(substr(name, start=17+cont, stop=17+cont))
			nam=paste(nam,substr(name, start=17+cont, stop=17+cont),sep="")
			cont=cont+1
			print(cont)
		}
	}else{
		cont=1
		nam=""
		while((substr(name, start=11+cont, stop=11+cont)!=",")|(substr(name, start=11+cont, stop=11+cont)!=")")){
			print(substr(name, start=11+cont, stop=11+cont))
			nam=paste(nam,substr(name, start=11+cont, stop=11+cont),sep="")
			cont=cont+1
			print(cont)
		}		
	}
	return(nam)
}

setMethod(	"initializerRep", 
			signature(report = "list", contRep = "list"),
			function(report, contRep, ...){
					
				actRep=function(...){
					sink(paste(contRep$name,".tex",sep=""),append=T)	
					return(invisible())
				}	
				desGraph=function(...){
					try(dev.off(),silent=T)
				}
				desRep=function(...){
					sink()
					suppressWarnings(sink())
					return(invisible())
				}
				aux=new(Class="report",report=report[["report"]],graph=0,trueVal=report[["report"]],actRep=actRep,desGraph=desGraph,
					desRep=desRep,commentt=report[["comment"]],files=report[["files"]],contRep=contRep)
				return(aux)
			}
)

setMethod(	"initializer",
			signature(object = "ts", student="logical",  report="list", contRep="list"),
			function(object, student, report, contRep, ...){
				
				menus=new(Class="menu",items=list(NULL))
				menus=makeMenu(menus,opciones)
				creport=initializerRep(report,contRep)
				cami=new(Class="cami",cami="ini")
				lineal=new(Class="lineal",crit=0,ls=F,atip=NULL)
				serie=new(Class="serie",nom=contRep$name ,serie=(object),orig=1,sact=T,trans=1,est=0,reg=0,stac=F,lin=lineal)
				datos=new(Class="datos",lserie=list(serie),sident=1,lmodelo=list(NULL),mident=0,modif=T)

				session=new(Class="sessio.ts",
					menu=menus,
					cami=cami,
					datos=datos,
					student = student,
					report=creport
				)
				if(report[["report"]]){
					if(!contRep$fil){
						sink(paste(contRep$name,".tex",sep=""))	
						cat("\\documentclass{article}\n")
						cat("\n\\usepackage[colorlinks=true,urlcolor=blue]{hyperref}\n")
						cat("\n\\usepackage{color}\n")
						cat("\n\\usepackage[cp1252]{inputenc}\n")
						cat("\n\\usepackage{amscd}\n")
						cat("\n\\usepackage{graphicx}\n")
						cat("\n\\usepackage{float}\n")
						cat("\n\\usepackage{subfig}\n")				
						cat("\n\\addtolength{\\oddsidemargin}{-.875in}\n")
						cat("\n\\addtolength{\\evensidemargin}{-.875in}\n")
						cat("\n\\addtolength{\\textwidth}{1.75in}\n")
						cat("\n\\addtolength{\\topmargin}{-.875in}\n")
						cat("\n\\addtolength{\\textheight}{1.75in}\n")
						cat("\n\\makeindex\n")
						cat("\n\\begin{document}\n")
						cat(paste("\n\\title{",tex["previsiones"][1,"1"],": ",contRep$name,"}\n",sep=""))
						user=Sys.info()["user"]
						cat(paste("\n\\author{",user,"}\n",sep=""))
						cat("\n\\maketitle\n")
						cat("\n\\tableofcontents\n")
						cat("\n\\newpage\n")
						creport["desRep"]()
					}else{
						tex=suppressWarnings(readLines(contRep$tex,n=-1))
						sink(paste(contRep$name,".tex",sep=""))
						cat(tex,sep="\n")
						creport["desRep"]()
					}
				}
				return(session)
		}
)

setMethod(	"ender",
			signature(session = "sessio.ts"),
			function(session, ...){
				session["report"]["actRep"]()
				cat("\n\\section{ References }\n")
				cat("\nThe report is made with the library 'TSTutorial' from the software R.\n")
				cat("\n\\end{document}\n")
				session["report"]["desRep"]()

				# b=try(ifelse(test=as.logical(Sys.info()["sysname"] == "Windows"),
				# shell(cmd=paste("pdfLatex ",session["report"]@contRep$name,".pdf",sep=""),intern=T),
				# system(paste("pdflatex ",session["report"]@contRep$name,".tex",sep=""),intern=T)),silent=T)	
				#Se vuelve a realizar ya que con solo una ejecucion, el indice, en la gran mayoria de ocasiones, no aparece.

				# b=try(ifelse(test=as.logical(Sys.info()["sysname"] == "Windows"),
				# shell(cmd=paste("pdfLatex",session["report"]@contRep$name,sep=" "),intern=T),
				# system(paste("pdflatex",session["report"]@contRep$name,sep=" "),intern=T)),silent=T)	
			

				# b=try(ifelse(test=as.logical(Sys.info()["sysname"] == "Windows"),
				# shell.exec("pdfLatex ",paste(session["report"]@contRep$name,".pdf",sep="")),
				# system(paste("xpdf ",session["report"]@contRep$name,".pdf&",sep=""))),silent=T)

				b=try(ifelse(test=as.logical(Sys.info()["sysname"] == "Windows"),
				system(paste("pdflatex ",session["report"]@contRep$name,".tex",sep="")),
				system(paste("xpdf ",session["report"]@contRep$name,".pdf&",sep=""))),silent=T)

				b=try(ifelse(test=as.logical(Sys.info()["sysname"] == "Windows"),
				system(paste("pdflatex ",session["report"]@contRep$name,".tex",sep="")),
				system(paste("xpdf ",session["report"]@contRep$name,".pdf&",sep=""))),silent=T)

				
				if(!session["report"]["files"]){
						LatexFiles=c(".tex",".aux",".idx",".log",".out",".toc")
						deleteList=c(paste("Grafh",1:session["report"]["graph"],".jpg",sep=""),paste(session["report"]@contRep$name,LatexFiles,sep="")
						)
						unlink(deleteList)
				}
				prettyprint(tex["previsiones"][2,"1"])				
				return(invisible())
		}
)

setMethod(	"activate", 
			signature(x = "sessio.ts"),
			function(x, ...){
				men=x@menu@items[[getcami(x@cami)]]
				n=length(men@opciones)
				valores=men@previo(men=list(cami=x@cami,datos=x@datos,student=x@student,report=x@report))
				x@datos=valores$datos
				res=0
				res=menu(c(men@opciones,tex["menu"][1,"1"],tex["menu"][2,"1"]),graphics=F)
				while((res<=0)|(res>(n+2))){
					res=menu(c(men@opciones,tex["menu"][1,"1"],tex["menu"][2,"1"]),graphics=F)
				}
				if (res<=n)	valores=men@acciones[[res]](men=list(cami=x@cami,datos=x@datos,student=x@student,report=x@report))
				if (res==n+1) valores=men@ayuda(men=list(cami=x@cami,datos=x@datos,student=x@student,report=x@report))
				if (res==n+2) valores=salir(men=list(cami=x@cami,datos=x@datos,student=x@student,report=x@report))
				
				x@cami=valores$cami
				x@datos=valores$datos
				x@report=valores$report
				return(x)
			}
)



tex=new(Class="tex",opciones=matrix(NA),previo=matrix(NA),ayuda=matrix(NA),ini=matrix(NA),gesser=matrix(NA),trans=matrix(NA),iden=matrix(NA),gesmod=matrix(NA),
			estim=matrix(NA),valid=matrix(NA),cap=matrix(NA),atip=matrix(NA),prev=matrix(NA),salir=matrix(NA),plot=matrix(NA),drawser=matrix(NA),drawmod=matrix(NA),
			writelistser=matrix(NA),writemod=matrix(NA),writelistmod=matrix(NA),writecoef=matrix(NA), drawarma=matrix(NA),atipics=matrix(NA),drawatip=matrix(NA),
			enterComment=matrix(NA),latexCntrl=matrix(NA)
)

englishLanguage=function(a){
	return(c(
	"--menu"
	,"Help"
	,"Exit"
	,"--opciones.1"
	,"Exploratory analysis"
	,"Transformations"
	,"Tutorial"
	,"--opciones.2"
	,"Series manager"
	,"Brief explanation of a Stationary Time Series"
	,"Exploratory analysis"
	,"Non-constant variance"
	,"Seasonal component"
	,"Non-constant mean"
	,"Model detection"
	,"Back to previous menu"
	,"--opciones.3"
	,"Choose series"
	,"Delete series"
	,"Back to main menu"
	,"--opciones.4"
	,"Series manager"
	,"ACF/PACF"
	,"Insert model"
	,"Model estimation"
	,"Back to previous menu"
	,"--opciones.5"
	,"Model manager"
	,"Estimation"
	,"Fix part to zero"
	,"Modify model"
	,"Models validity"
	,"Back to previous menu"
	,"--opciones.6"
	,"Choose model"
	,"Delete model"
	,"Back to main menu"
	,"--opciones.7"
	,"Model manager"
	,"Residual analysis"
	,"Compare theoretical-sample ACF/PACF"
	,"AR(infinite)/MA(infinite)"
	,"Delete model"
	,"Prediction capacity"
	,"Back to previous menu"
	,"--opciones.8"
	,"Model manager"
	,"Stability"
	,"Prediction capacity"
	,"Choose best model"
	,"Atypical treatment"
	,"Long-term predictions"
	,"Back to previous menu"
	,"--opciones.9"
	,"Model manager"
	,"Atypical treatment"
	,"Back to previous menu"
	,"--opciones.10"
	,"Model manager"
	,"Long-term predictions"
	,"Back to previous menu"
	,"--previo.1"
	,"PREDICTIONS"
	,"TRANSFORMATIONS"
	,"MODEL DETECTION"
	,"MODEL ESTIMATION"
	,"MODEL VALIDITY"
	,"PREDICTION CAPACITY"
	,"ATYPICAL TREATMENT"
	,"LONG-TERM PREDICTIONS"
	,"PREVIOUS"
	,"--previo.2"
	,"Predictions"
	,"Transformations"
	,"Model Detection"
	,"Model Estimation"
	,"Model Validity"
	,"Prediction Capacity"
	,"Atypical Treatment"
	,"Long-term Predictions"
	,"Previous"
	,"--previo.3"
	,"The objective of the program is to get predictions of the series"
	,"which starts in"
	,"and with a period"
	,". Firstly, it is recommended to make an exploratory analysis to get a previous idea of the series behaviors and properties. Then, the next step is 'Transformations'."
	,"To can obtain models to get good predictions, the series has to be stationary (see the vignette 'Stationary' or select the option 2). For this, if it necessary, you have to make certain transformations. Firstly it has to check the uniform of the variance, if it is not uniform, in most cases, it is solved applying logarithm."
	,"Once we have a stationary series we can identify its possible models. To can do it, you have to look its ACF/PACF plot and, then, introduce the found models."
	,"Once we have some models, we have to estimate them and, in case it has no significant coefficients, depends on the case, you have to fix these coefficients or modify the models until they have significant all their coefficients."
	,"When we have significant models, the last step is the validation where we may rule out the invalid models."
	,"With the obtained valid models you can obtain predictions but it can be some differences between the quality (Reality approximation) in these predictions. Also, these models could not be steady (The same model fitted without the latest observations has coefficients very different). For this, it is recommended analyze the stability and the prediction capacity of the different valid models to then choose the best model and proceed to the 'Atypical Analysis' or to get 'Long-term Predictions'. The advisable number of reservations is one period of the series."
	,"The atypical treatment allows obtaining the series that would exist if the unexpected events, which have caused changes in the series, hadn't happened. In case to use the atypical treatment, with the obtained series, the program readdress to the menu 'Transformations' to restart (The linearized series maintains the transformation applied to obtain constant variance). For the atypical detection be optimal, it is necessary that the active model was the best model created."
	,"Once you have a good model, you can obtain long-term predictions. Remember that the confident intervals are at 95% of confidence."
	,"--ayuda.1"
	,"'Initial' help menu"
	,"1: The 'Exploratory analysis' option makes a data analysis showing numeric and graphical results."
	,"2: The 'Transformations' option readdress to the following menu where you can make transformations to the series for making it a stationary one."
	,"3: The 'Tutorial' option opens a .pdf file which contains a brief tutorial of TSTutorial."
	,"--ayuda.2"
	,"'Series manager' help menu"
	,"1: The 'Choose series' option let choose, from a list which contains all the series, the series which you want to work."
	,"2: The 'Delete series' option let delete one series from the list of series which contains all the created series."
	,"--ayuda.3"
	,"'Transformations' help menu"
	,"1: The 'Series manager' option lets to manage the series choosing with which work or delete."
	,"2: The 'Brief explanation of an Stationary Time Series' option opens a .pdf file which shows with examples how is an stationary series and the different non-stationary series that it could be founded."
	,"3: The 'Exploratory analysis' option makes a data analysis showing numeric and graphical results."
	,"4: The 'Non-constant variance' option shows, through some analysis, whether the series has constant variance and, if it is necessary, you have to indicate what transformation you want to make."
	,"5: The 'Stationality' option shows whether the series can have seasonal component for, then, make it a seasonal differentiation, in case it was necessary."
	,"6: The 'Non-constant mean' option shows whether the series has non-constant mean and, if it is necessary, makes a regular differentiation."
	,"7: The 'Model detection' readdress to the next menu where you can identify and introduce the possible models of the series."
	,"--ayuda.4"
	,"'Model detection' help menu"
	,"1: The 'Series manager' option lets to manage the series choosing with which work or delete."
	,"2: The 'ACF/PACF' option shows the ACF/PACF of the active series to let you to deduce its possible models."
	,"3: The 'Insert model' lets to introduce one possible model of the active series."
	,"4: The 'Model estimation' readdress to the next menu where you can estimate and modify the created models."
	,"--ayuda.5"
	,"'Model manager' help menu"
	,"1: The 'Choose model' option let choose, from a list which contains all the models, the model which you want to work."
	,"2: The 'Delete model' option let delete one model from the list of models which contains all the created models."
	,"--ayuda.6"
	,"'Model estimation' help menu"
	,"1: The 'Model manager' option lets to manage the models choosing with which work or delete."
	,"2: The 'Estimation' option shows the active model estimation; suggest whether it has to work with constant and whether it contains non-significant coefficients."
	,"3: The 'Fix part to zero' option lets to modify one coefficient on each time you select the option. Then, it contrasts the fixed model with the previous one to confirm the action."
	,"4: The 'Modify model' option modify the active model letting to choose what components you want to modify (p, q, P o Q)."
	,"5: The 'Models validity' readdress to the next menu where, through some analysis, you can determine whether the active model is valid to make predictions."
	,"--ayuda.7"
	,"'Model validity' help menu"
	,"1: The 'Model manager' option lets to manage the models choosing with which work or delete."
	,"2: The 'Residual analysis' option shows through some plots whether the analyzed model has good residuals and, in this case, if the model is valid."
	,"3: The 'Compare theoretical-sample ACF/PACF' option compares the ACF/PACF's active model with the theoretical to see if their behaviors are similar."
	,"4: The 'AR(infinite)/MA(infinite)' option shows a comparative of the AR(infinite)/MA(infinite) from all models to can decide the existence of equivalent models."
	,"5: The 'Delete model' option let delete the active model."
	,"6: The 'Prediction capacity' option readdress to the next menu where you can calculate the stability and the prediction capacity of the models to decide which is the best model."
	,"--ayuda.8"
	,"'Prediction capacity' help menu"
	,"1: The 'Model manager' option lets to manage the models choosing with which work or delete."
	,"2: The 'Stability' option calculates if the active model is stable."
	,"3: The 'Prediction capacity' option calculates the model prediction capacity calculating its Mean Squared Error."
	,"4: The 'Choose best model' option assign the chosen model as the best model."
	,"5: The 'Atypical treatment' option readdress to the next menu where you can apply the atypical treatment to linearize the series."
	,"6: The 'Long-term predictions' option readdress to the next menu in which you can obtain long-term predictions."
	,"--ayuda.9"
	,"'Atypical treatment' help menu"
	,"1: The 'Model manager' option lets to manage the models choosing with which work or delete."
	,"2: The 'Atypical treatment' option indicates what criterions can be the best to use and you can linearize the series with the desired criterion."
	,"--ayuda.10"
	,"'Long-term predictions' help menu"
	,"1: The 'Model manager' option lets to manage the models choosing with which work or delete."
	,"2: The 'Long-term predictions' makes long-term predictions."
	,"--ini.1"
	,"It shows the series numeric descriptive; the series plot, its histogram, the qqnorm plot and its ACF; and the descompose plot which shows the series without its seasonal component and trend, the series without its seasonal component, its trend and its remember plot (reordered ACF)."
	,"EXPLORATORY ANALYSIS OF THE SERIES:"
	,"Exploratory analysis of the series:"
	,"General numeric description:"
	,"The mean of the series is"
	,"and its standard deviation is"
	,"Its median value is"
	,", its minimum value is"
	,"and its maximum value is"
	,"--gesser.1"
	,"Series list"
	,"Variance"
	,"Indicate the list number which contains the series with which you want to work."
	,"The number that you have chosen it doesn't exist in the list."
	,"--gesser.2"
	,"Indicate the list number which contains the series with which you want to delete."
	,"The number that you have chosen it doesn't exist in the list."
	,"(Remember that you can't delete the original one)"
	,"The deleted series was the active series, for this reason you have to choose the series with which you want to work."
	,"You have to have at least 2 series to can choose one to delete."
	,"--trans.4"
	,"It shows a plot which contains the boxplots for every period, if you see that the initials are enough different from the final ones, it's very possible that the series doesn't have constant variance. Moreover, it shows a plot with the means vs. the Standard Deviations that if it is a flat line it means that the variance is constant. Finally, it shows a plot with the Box-Cox transformation that indicates the better value of lambda. In case to transform the series, the program will make automatically an exploratory analysis."
	,"Non-constant variance series:"
	,"SUGGESTION: It is suggested to transform the data."
	,"SUGGESTION: It seems no necessary transform the data."
	,"Write 0 if you want to make a logarithm transformation or the lambda value of the transformation. If you don't want to make a transformation, write 1."
	,"The transformation done is logarithmic."
	,"It hasn't made a transformation."
	,"It has made a transformation with lambda"
	,"Do you consider that the series is stationary?"
	,"The series"
	,"it is considered stationary."
	,"it isn't considered stationary."
	,"The active series has been transformed yet. Please, select a not transformed series."
	,"--trans.5"
	,"As well as helping with the 'Exploratory analysis' monthplot, it also creates a new plot which contain one period with all the observations separated by position of the period that it take up. With it, you can see if it contains different heights is the key of the possible existence of stationary component. In case to want to differentiate the series, the seasonal order is the period of the series."
	,"Series seasonal:"
	,"Do you want to differentiate?"
	,"What is the seasonal order of the series?"
	,"It is made a seasonal differentiation with order"
	,"to the series."
	,"Do you consider that the series is stationary?"
	,"It is decided not make a seasonal differentiation to the series."
	,"The series has been seasonal differenced before. For R specifications, R can't work with more than one seasonal differentiation."
	,"You have written a negative number or it is not integer."
	,"--trans.6"
	,"You can be helped using the monthplot from 'Exploratory analysis'. If the trend doesn't seems like a horizontal straight line is possible that you have to make a regular differentiation."
	,"Non-constant mean series:"
	,"Do you want to differentiate?"
	,"It is made a regular differentiation to the series."
	,"It is decided not to make a regular differentiation to the series."
	,"Do you consider that the series is stationary?"
	,"--iden.2"
	,"It generates a plot which contains the ACF and PACF of the series."
	,"ACF/PACF of the series:"
	,"--iden.3"
	,"Write p"
	,"Write q"
	,"Write P"
	,"Write Q"
	,"Insertion the model of the series:"
	,"The model that you have introduced doesn't achieve the stationary properties, try with another model."
	,"It is created the model"
	,"You have written a negative number or it is not an integer."
	,"--gesmod.1"
	,"Model List"
	,"Indicate the list number which contains the model with which you want to work."
	,"The number that you have chosen isn't in the list."
	,"You have to create at least one model to choose."
	,"--gesmod.2"
	,"Indicate the list number which contains the model with which you want to work."
	,"The number that you have chosen isn't in the list."
	,"The delete model was the active model, for this reason you have to choose the model with which you want to work."
	,"You need at least 1 model to can choose one to delete."
	,"--estim.2"
	,"Model estimation of"
	,"The model constant isn't significant and this is the model estimation without constant:"
	,"SUGGESTION: The model AIC without constant is less than with constant, for this reason, the model is better."
	,"SUGGESTION: The model AIC without constant is bigger than with constant, for this reason, the model is worse."
	,"Do you want to work without constant?"
	,"It is chosen to work without constant."
	,"It is chosen to work with constant."
	,"SUGGESTION: The"
	,"model has at least one final coefficient not significant. In case you want to fix it, is better to modify the model subtracting one of the components which has this coefficient."
	,"model has at least one final and middle coefficients not significant. In case you want to fix one of the final coefficients, is better to modify the model subtracting one of the components which has this coefficient. On the contrary, you may fix the intermediate non-significant coefficients."
	,"model contains middle coefficients not significant. It is recommended to fix it."
	,"You have to have at least 1 model to can estimate it."
	,"new"
	,"--estim.3"
	,"Sometimes the model contains non-significant coefficients (|Coef/s.e.|<1.96) that it has to fix to zero to make the model better."
	,"Really do you want to fix one model coefficient?"
	,"Fix part of the model"
	,"Active model"
	,"SUGGESTION: It is recommended fixing the coefficient less significant (the closest to zero)."
	,"Indicate the position number that occupies the coefficient you want fix (left->right)."
	,"You have written an incorrect position."
	,"The chosen coefficient is already fixed to zero."
	,"It is decided to fix the coefficient that occupies the position"
	,"Modify model"
	,"The modified model AIC is less than the active, for this reason the model has improved."
	,"The modified model AIC is bigger than the active, for this reason the model has deteriorated."
	,"Do you want to work with the modified model?"
	,"It has chosen to work with the modified model."
	,"It has chosen not to work with the modified model."
	,"SUGGESTION: The"
	,"model has at least one final coefficient non-significant. In case you want to fix it, is better to modify the model subtracting one to the component which has this coefficient."
	,"model has at least one final and middle coefficients not significant. In case you want to fix one of the final coefficients, is better to modify the model subtracting one of the components which has this coefficient. On the contrary, you may fix the intermediate non-significant coefficients."
	,"model contains middle coefficients not significant. It is recommended to fix it."
	,"You cannot fix more coefficients because the model only contains one."
	,"You need at least 1 model to can estimate it."
	,"new"
	,"--estim.4"
	,"Do you want to modify the model?"
	,"Modify the model"
	,"What parameter do you want to modify?"
	,"What parameter do you want to modify?"
	,"It is decided to modify the parameter"
	,"The actual value of p is"
	,". Write its new value."
	,"The new value of p is"
	,"The actual value of q is"
	,"The new value of q is"
	,"The actual value of P is"
	,"The new value of P is"
	,"The actual value of Q is"
	,"The new value of Q is"
	,"The model that you have introduced doesn't achieve the stationary properties, try with another model."
	,"You need at least 1 created model to can modify it."
	,"You have written a negative number or it is not integer."
	,"--valid.2"
	,"Generates 6 plots: the model residuals, the data approximation to the Normal straight line, the residual histogram, the residual ACF/PACF, the squared residual ACF/PACF and the Tsdiag"
	,"SUGGESTION: The analysis indicates that the model is valid."
	,"SUGGESTION: The analysis indicates that the model is not valid."
	,"Do you consider that the model is valid?"
	,"The model"
	,"is considered valid."
	,"is not considered valid."
	,"You need at least 1 created model to can analyze the residuals."
	,"Residual analysis of"
	,"--valid.3"
	,"Comparison of the theoretical ACF/PACF with the sample"
	,"Generates 2 plots, the first contains the ACF and the second the PACF."
	,"You need at least 1 created model to can compare the theoretical ACF/PACF with the sample."
	,"--valid.4"
	,"Shows in the console the invertible coefficients and creates a plot that contains the model expressions like AR and MA infinite."
	,"Model AR/MA infinites"
	,"Invertibility coefficients"
	,"SUGGESTION: The model is invertible."
	,"SUGGESTION: The model is not invertible, for this reason, it is advisable to work with another model because its predictions will not be very accurate."
	,"You need at least 1 created model to can see the AR(infinite)/MA(infinite)."
	,"Your model is an AR model. By definition, your model is invertible."
	,"--valid.5"
	,"Do you want to delete the active model?"
	,"--cap.2"
	,"How many observations do you want to reserve?"
	,"You cannot reserve more observations than the length of the series or less than one and a non-integer number."
	,"Model stability of"
	,"It has been chosen reserve"
	,"observations."
	,"A stable model is when the same model without the reserved observations maintains similar values of the original model coefficients. In  case to be a quite a lot of different it can be an indication of the atypical presence in the last part of the series."
	,"SUGGESTION: It seems that the model is a quite a lot of unstable because many coefficients have changed considerably in its values. You should choose another model to analyze its prediction capacity because its MSE will be very high. Another option is to make an 'Atypical treatment' and then draw conclusions."
	,"SUGGESTION: It seems that the model is a little bit unstable, it could be for the existence of atypical values in the latest observations from the series."
	,"The linearize criterion may be smaller."
	,"SUGGESTION: All the tests indicate that the model is stable."
	,"Do you consider that the model is stable?"
	,"It is considered that the model is stable"
	,"It is not considered that the model is stable"
	,"You need at least one created model to can study its stability."
	,"--cap.3"
	,"Generates a list shown in the console and a plot that contains the obtained reserved predictions with their confident intervals and the real data of the observations."
	,"Model predict capacity of"
	,"How many observations do you want to reserve?"
	,"You cannot reserve more observations that the length of the series or less than one and a non-integer number."
	,"It has been chosen reserve"
	,"observations."
	,"Obtained predictions with the original data"
	,"You have"
	,"You need at least one created model to can check its prediction capacity."
	,"--cap.4"
	,"To select the best model, it's useful use the plot which shows the model collection to can compare it."
	,"Best model selection"
	,"Do you want to continue with the best model selection?"
	,"SUGGESTION: It's useful analyze the prediction capacity of all the models to have an optimal comparison between it."
	,"SUGGESTION: The model with best AIC is"
	,"Write the position taken up in the model list that you want to fix as best model."
	,"You have written an incorrect position."
	,"It have been chosen as best model"
	,"Do you want that the best model was the active model?"
	,"SUGGESTION: The model with best AIC and MSE is"
	,"SUGGESTION: The model with best MSE is"
	,"The unique created model that you have is already selected as best model."
	,"You have only one created model. Do you want to choose it as best model?"
	,"You need at least one created model to can choose it as best model."
	,"--atip.2"
	,"Creates two plots to make a comparison making different treatments: one have the atypical list with LS and the other plot without LS. You have to compare between the plots and between the results in the plots."
	,"This procedure could take some time...."
	,"Do you think that you have to make an atypical treatment to the series?"
	,"Atypical treatment of the model"
	,"Write the treatment criterion:"
	,"Do you want that the treatment works with LS?"
	,"The treatment criterion is"
	,"with"
	,"without"
	,"LS."
	,"Detected atypical list"
	,"With the chosen criterion, it hasn't found any atypical."
	,"The linearized series variance is"
	,"You cannot make the 'Atypical Treatment' with a model already linearized."
	,"You have to be at least one created model to can apply it the 'Atypical Treatment'."
	,"The criterion has to be bigger than 0.1."
	,"--prev.2"
	,"Generates a list shown in the console and a plot that contains the obtained predictions with their confident intervals and the previous series data."
	,"Lon-term predictions of the model"
	,"How many predictions you want to obtain?"
	,"You need at least one created model to can obtain its long-term predictions."
	,"Long-term predictions"
	,"The number has to be a positive integer."
	,"--salir"
	,"Are you sure that you want to exit?"
	,"Session Data Summary"
	,"--plot"
	,"(Plot, Histogram, QQnorm, Acf)"
	,"(Plot descompose)"
	,"(Boxplot)"
	,"(Monthplot)"
	,"(ACF/PACF plot)"
	,"(Residual plot)"
	,"(Residual QQnorm plot)"
	,"(Residual histogram)"
	,"(Residual ACF/PACF plot)"
	,"(Squared residual ACF/PACF plot)"
	,"(Tsdiag plot)"
	,"(ACF sample vs. theoretical plot)"
	,"(PACF sample vs. theoretical plot)"
	,"(Capacity of Prevision plot)"
	,"(Series minus linearized series plot)"
	,"(Superimpose the series with the linearized series plot)"
	,"(Long-term predictions plot)"
	,"(List of series plot)"
	,"(List of models plot)"
	,"(AR/MA infinite)"
	,"(Plot Mean vs. Standard deviation)"
	,"(Plot Box-Cox transformation)"
	,"--drawser"
	,"SERIES COLLECTION"
	,"Series"
	,"Variance"
	,"Stat."
	,"Lineal"
	,"no"
	,"and"
	,"The active series is"
	,"--drawmod"
	,"MODEL COLLECTION"
	,"Arima model"
	,"Int."
	,"AIC"
	,"Valid"
	,"Sta."
	,"MSE"
	,"Lineal"
	,"no"
	,"Crit"
	,"The active model is"
	,"The best model is"
	,"MODEL COLLECTION"
	,"It doesn't exist any created model."
	,"--writelistser"
	,"(Stationary)"
	,"(Linearized)"
	,"--writemod"
	,"with constant."
	,"--writelistmod"
	,": Arima model"
	,"(Validate)"
	,"(Linearized)"
	,"(Best model)"
	,"--writecoef"
	,"Model coefficients of:"
	,"Significant?"
	,"no"
	,"nearly"
	,"yes"
	,"--drawarma"
	,"with constant."
	,"--atipics"
	,"SUGGESTION: It is not necessary the atypical treatment"
	,"SUGGESTION: With LS it is recommended to work with"
	,"SUGGESTION: Without LS it is recommended to work with"
	,"the criterions"
	,"the criterion"
	,"--drawatip"
	,"Atypical analysis with LS"
	,"Atypical analysis without LS"
	,"--previsiones"
	,"FITTING TEMPORAL SERIES"
	,"You have your report already in a .pdf file in your working directory."
	,"--enterComment"
	,"Do you want to write a comment?"
	,"Following, write your comment."
	,"COMMENT"
	,"--latexCntrl"
	,"You have chosen that TSTutorial make a report, but you don't have installed any .tex compiler. For this reason, you cannot get it. Automatically, the function starts to work but without making the report."
	,"--end"
	))
}
tex=texinitializer(tex)

