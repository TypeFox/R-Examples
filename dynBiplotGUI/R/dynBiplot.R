dynBiplot <-


function (lang="es")
{
	#library(tcltk)
	#library(tcltk2)
	#library(tkrplot)
	#
	# 	b.x 	- Archivo de entrada
	# 	b.x2 	- Matriz de 2 vias
	# 	b.x3 	- Cubo de 3 vias
	#	b.f* 	- Archivos de formatos
	#	bt.*	- Archivos/variables temporales
	
	#	Variables
	#
	bt.leer  <- tclVar(0)		# tipo de archivo a cargar
	bt.le  <- 1					# niveles de situaciones
	tb  <- tclVar("1")			# tipo de biplot HJ
	vf  <- tclVar()				# para formato
	i3v <- tclVar("0")			# indicador de datos de 3 vias
	ifo <- tclVar("0")			# indicador de carga de formatos
	if1 <- tclVar("0")			# indicador de campos de formatos
	icar <- tclVar("0")			# indicador de datos cargados
	ifd <- tclVar("0")			# indicador de ya mostrado formato datos
	it1 <- tclVar("0")			# indicador de mostrar titulo
	it2 <- tclVar("0")			# indicador de mostrar subtitulo
	iagr <- 1					# indicador de boton de formato
	iyax <- 0					# indicador de ya rotado x
	iyay <- 0					# indicador de ya rotado y
	isd <- tclVar("0")			# indicador de ya mostrado selector de datos
	imos <- tclVar("1")			# indicador de mostrar en dibujo
	ieje <- tclVar("1")			# indicador mostrar ejes.
	ipca <- tclVar("0")			# indicador de PCA
	icen <- tclVar("1")			# indicador de centrado
	iesc <- tclVar("1")			# indicador de escalado
	iref <- tclVar("1")			# indicador estandarizar por #dimension.
	ibg  <- tclVar("0")			# indicador de biplot global
	ietx <- tclVar("1")			# indicador de etiquetas para x
	iety <- tclVar("1")			# indicador de etiquetas para y
	ietzr <- tclVar("1")		# indicador de etiquetas para z filas
	ietzc <- tclVar("1")		# indicador de etiquetas para z columnas
	itrr <- tclVar("1")			# indicador de mostar trayectoria filas
	itrc <- tclVar("1")			# indicador de mostar trayectoria columas
	ittp <- tclVar(0)			# indicador de existe ventana grafica
	igx  <- tclVar(0)			# indicador rotar x en el grafico
	igy  <- tclVar(0)			# indicador rotar y en el grafico
	ipv  <- tclVar(1)			# indicador de p-valor para trayectorias
	ivs  <- tclVar(0)			# indicador de concatenacion etiquetas trayectorias
	pval <- tclVar(0.05)		# margen para p-valor
	ex1 <- tclVar()				# variable rotar -x
	ex2 <- tclVar()				# variable rotar +x
	ey1 <- tclVar()				# variable rotar -y
	ey2 <- tclVar()				# variable rotar +y
	veti <- tclVar()			# variable que tiene las etiquetas 
	vsit <- tclVar()			# variable que tiene las situaciones 
	vagr <- tclVar()			# variable de agrupacion para formato 
	vag1 <- tclVar()			# valor de agrupacion seleccionada
	vref <- tclVar()			# variable que tiene la referencia 
	vinx <- tclVar("0")			# variable cantidad de inercia a mostrar para x
	viny <- tclVar("0")			# variable cantidad de inercia a mostrar para y	
	vesc <- tclVar("1")			# variable para reescalado
	vtit <- tclVar("Biplot")	# variable para recoger el titulo
	vsub <- tclVar()			# variable para recoger el subtitulo
	neje <- tclVar("2")			# captura ejes a tratar
	di1 <- tclVar("1")			# captura dim1 a tratar
	di2 <- tclVar("2")			# captura dim2 a tratar
	nejes <- 2					# numero de ejes a tratar
	dim1 <- 1;	dim2 <- 2		# para tratar las dimensiones
	label.ejes <- array()		# etiquetas "Eje" para los factores
	wout1 <- tclVar("1.8")		# ancho ventana grafica de salida
	wout2 <- tclVar("1.8")		# alto ventana grafica de salida
	p4cbx <- tclVar()			# Combobox del panel 4 
	#
	#		Declaracion de variables para evitar NOTE en R CMD check
	#
	bt.hoja <- bt.hc <- bt.hr <- bt.hs <- NULL	
	b.x2 <- b.x3 <- NULL
	bt.cubo.b <- bt.nl <- bt.leve <- NULL
	bt.t <- bt.x <- bt.x2 <- bt.x3 <- bt.fx <- bt.fy <- bt.tf <- bt.fxg <- bt.fyg <- NULL
	bt.mean <- bt.sd <- bt.x2m <- bt.x2sd <- NULL
	bt.c <- bt.e <- bt.xce <- bt.zc <- bt.zr <- bt.r2c <- bt.r2r <- NULL
	bt.Fc <- bt.Prc <-NULL
	bt.svd <- bt.U <- bt.V <- bt.D <- bt.a <- bt.b <- NULL
	bt.varexpl <- bt.res.vp <- bt.res.a <- bt.res.b <- NULL
	bt.res.cr <- bt.res.cc <- bt.res.ty <- bt.res.tx <- bt.res.r2y <- NULL
	bt.res.Fy <- bt.res.Pry <- NULL
	bt.en1 <- bt.micol <- bt.cv1 <- bt.cbb1 <- bt.en2 <- bt.cbb2 <- NULL
	bt.micol2 <- bt.cv2 <- bt.typ1 <- bt.bu3 <- tb3.t2 <- NULL
	bt.limx <- bt.limy <- bt.limx1 <- bt.limy1 <- NULL
	bt.ttp <- bt.ttp1 <- bt.ttp2 <- bt.img <- NULL
	
	#
	#	Seleccion de lenguaje
	#
	lit0 <- read.csv(file.path(path.package("dynBiplotGUI"),"lang",
				"Language.csv",fsep=.Platform$file.sep),
				header=T,as.is=1,sep=";",encoding="latin1")
	if(lang=="es")		lit <- lit0["es"]		# espanol
	else if(lang=="fr") lit <- lit0["fr"]		# frances
	else if(lang=="pt") lit <- lit0["pt"]		# portugues
		else {			lit <- lit0["en"]		# ingles, en los demas casos
			if(lang!="en") print("Language not implemented. We use English.")
			}
	bt.lit <- apply(lit,2,paste)				# forzar comillas
	#
	#	Funciones	<<<<<<<==========================
	#
	#	Listar archivos en R (basado en la funcion lsos()
	#
	lsos <- function(pos=1) {
		napply <- function(nn, fn) sapply(nn, function(x) fn(get(x)))
		nn <- ls(".GlobalEnv")
		clase <- napply(nn, function(x) as.character(class(x))[1])
		modo <- napply(nn, mode)
		tipo <- ifelse(is.na(clase), modo, clase)
		out <- data.frame(nn,tipo,stringsAsFactors =F)
		n2 <- out["tipo"]=="data.frame" | out["tipo"]=="matrix" | out["tipo"]=="array"
		out["nn"][n2]
	}
	#
	#		leer archivos	------------
	#
	leer.archivos <- function()	{
	if (tclvalue(icar)=="1") {		# Datos ya cargados
		print(bt.lit[20,])			# ERROR datos ya cargados
		return()
	}
	if(1==2) b.x <- b.x				# para evitar error en R CMD check
	if (tclvalue(bt.leer)=="1") leer.df()

        else if (tclvalue(bt.leer)=="2") leer.excel()
			else if (tclvalue(bt.leer)=="3") b.x <<- leer.csv()
				else if (tclvalue(bt.leer)=="4") b.x <<- leer.txt()
					else if (tclvalue(bt.leer)=="5") b.x <<- leer.spss()
						else b.x <<- leer.clipboard()
	if (length(b.x)==0) {print(bt.lit[99,])		# -Selecciona tipo de fichero-
						return() }
	tkconfigure(la, text=mens.leer, foreground="black", background="yellow2")	# para mostrar a pie de pagina
	tk2tip(la, bt.lit[21,])			# Fichero cargado
	cubo(b.x)
	tclvalue(icar) <- "1"			# indicador de datos cargados
	}
	#		leer excel
	#
	leer.excel <- function()	{	
		# require(xlsx)
		if(1==2) b.x <- NULL				# para evitar error en R CMD check
		hojaex <- tclvalue(tkgetOpenFile(filetypes = "{{Excel files} {.xls .xlsx}}"))
		if (!length(hojaex))   return()
		tmp <- xlsx::loadWorkbook(hojaex)
		hojas <- names(xlsx::getSheets(tmp))		# lista de hojas
		whoja <- tktoplevel()
		tkwm.title(whoja,bt.lit[22,])					# Selecciona
		tl <- tk2listbox(whoja, height = min(length(hojas),15),
				values=hojas, selectmode = "browse",background = "white")
		tkpack(tk2label(whoja, text = bt.lit[23,]))		# Selecciona tabla
		tkpack(tl)
		tkselection.set(tl, 0)
		OnOK <- function()	{
			bt.hoja <- hojas[as.numeric(tkcurselection(tl)) + 1]
			b.x <<- xlsx::read.xlsx(hojaex,sheetName=bt.hoja,encoding="UTF-8")
			mens.leer <<- bt.hoja
			tkdestroy(whoja)
		}

		OK.but <- tk2button(whoja,text=bt.lit[25,],command=OnOK)	# OK
		tkpack(OK.but)
		tkfocus(whoja) 
		tkwait.window(whoja)
		}

	#		leer csv, sep=";", dec="."
	#
	leer.csv <- function()	{	
		mens.leer <<- tclvalue(tkgetOpenFile(filetypes = "{{CSV files} {.csv}}"))
		if (mens.leer=="") return()
		read.csv(mens.leer, header=T, sep=";", dec=".")
		}

	#		leer dataframe
	# 
	leer.df <- function()	{
		if(1==2) b.x <- b.x				# para evitar error en R CMD check
		if (length(ls(.GlobalEnv))==0) 
							{print(bt.lit[111,])		# Error en tipo de archivo
							b.x <<- NULL
							return() }	# no hay ficheros
		h0 <- lsos()
		aa <- h0 %in% ls(.GlobalEnv,pattern=".f")	# para no mostrar formatos
		h0 <- h0[!aa]
		if (length(h0)==0) {print(bt.lit[111,])		# Error en tipo de archivo
							b.x <<- NULL
							return() }	# no hay ficheros
		t0 <- tktoplevel()
		tkwm.title(t0,bt.lit[22,])		# Selecciona
		tl <- tk2listbox(t0, height = min(length(h0),15),
				values=h0, selectmode = "browse",background = "white")
		tkgrid(tk2label(t0,text=bt.lit[24,]))	# Selecciona objeto
		tkgrid(tl)
		tkselection.set(tl,0)
		OnOK <- function()	{
			bt.h1 <- h0[as.numeric(tkcurselection(tl)) + 1]
			b.x <<- get(bt.h1)
			mens.leer <<- bt.h1

			tkdestroy(t0)

		}
		OK.but <- tk2button(t0,text=bt.lit[25,],command=OnOK)	# OK



		tkgrid(OK.but) 
		tkfocus(t0) 
		tkwait.window(t0)
		}









	#		leer txt
	#
	leer.txt <- function()	{
		mens.leer <<- tclvalue(tkgetOpenFile(filetypes = "{{Text files} {.txt}}"))
		if (mens.leer=="") return()
		read.table(mens.leer, header=T)
		}
	#		leer SPSS
	#
	leer.spss <- function()	{
		# library(foreign)
		mens.leer <<- tclvalue(tkgetOpenFile(filetypes = "{{SPSS Files} {.sav}}"))
		if (mens.leer=="") return()
		foreign::read.spss(mens.leer, use.value.labels=F,to.data.frame=T)
		}
	#		leer clipboard
	#
	leer.clipboard <- function()	{
		mens.leer <<- bt.lit[26,]		# Portapapeles
		read.table("clipboard",header=T)
		}
	#
	#	Generacion del cubo de 3 vias
	# 		fichero de entrada b.x
	cubo <- function(x) {
		# variable de etiquetas y referencia:
		#	panel de datos - para generar el cubo
		fr12 <- tk2frame(fr.d2,relief="raised", borderwidth=2,padding="2")
		fr12.2 <- tk2frame(fr12)
		fr12.3 <- tk2frame(fr12)


		vari <- colnames(x)
		tclvalue(veti) <- colnames(x)[[1]]		# valor inicial a mostrar




		tclvalue(vsit) <- colnames(x)[[2]]		# valor inicial a mostrar
		bt.cubo.b <<- tk2button(fr12,text=bt.lit[27,], command=function() cubo.gen(x))	# Generar matrices
		tkpack(tk2label(fr12.2, text=bt.lit[28,],background="lightyellow",width=12),	# Etiquetas
			tk2combobox(fr12.2,values=vari,textvariable=veti,width=15),side="left")
		tk2tip(fr12.2,bt.lit[29,])				# Variable que tiene las etiquetas
		if (tclvalue(i3v)=="1") {
			tkpack(tk2label(fr12.3,text=bt.lit[30,],background="lightyellow",width=12),	# Situaciones
				tk2combobox(fr12.3,values=vari,textvariable=vsit,width=15),side="left")
				tk2tip(fr12.3, bt.lit[31,])	}	# Variable que tiene las situaciones
		# los valores seleccionados se recogen al pulsar el boton
		tkpack(fr12.2, fr12.3, side="top")
		tkpack(fr12, bt.cubo.b)	
	}
	
	cubo.gen <- function(x) {
		ifelse (tclvalue(i3v)=="1", {b.x3 <<- cubo3(x)
									b.x2 <<- b.x3[,,1]}, 	# provisional para etiquetas
									b.x2 <<- cubo2(x))
		tkconfigure(bt.cubo.b,state="disable")
		formato()
		llena.lbr()					# Llamada a la funcion seleccion de variables
		VeOpciones ()
		tk2notetab.select(nb,bt.lit[5,])	# Variables
		}
		
	cubo2 <- function(x) {
		ve <- tclvalue(veti)
		rownames(x) <- x[,ve]
		x[colnames(x)!=ve]					# elimino la columna de etiquetas
		}
		
	cubo3 <- function(x) {
		ve <- tclvalue(veti)
		vs <- tclvalue(vsit)
		bt.nl <<- length(levels(as.factor(x[[vs]])))
		bt.leve <<- levels(as.factor(x[[vs]]))
		bt.t <- max(bt.leve) 		# Referencia provisional hasta seleccionarla
		# elimino las columnas de etiquetas y situaciones
		bt.x <- x[sapply(x,is.numeric)]	# solo datos numericos
		if(nrow(bt.x)%%bt.nl==0) 	# ckeck bloques completos
			bt.x3 <- array(0,c(nrow(x)/bt.nl,ncol(bt.x),bt.nl))
			else stop(bt.lit[133,])	# ERROR: bloques incompletos
			

		# 	carga la via 3
		for (k in 1:bt.nl) bt.x3[,,k] <- as.matrix(subset(bt.x,x[vs]==bt.leve[k]))
		#
		#	etiqueta el cubo
		colnames(bt.x3) <- colnames(bt.x)
		rownames(bt.x3) <- x[ve][x[vs]==bt.t]	
		dimnames(bt.x3)[[3]] <- bt.leve
		#	Chequeamos trayectorias nulas:
		#	filas
		apply(bt.x3,c(1,2),sum)==0 -> hay0
		if(any(hay0)) {bt.x3[,,1][hay0] <- 0.1
				print(bt.lit[32,])}	# AVISO: una fila tiene todos 0. Se pone 0.1 a una celda
		#	columnas
		apply(bt.x3,c(2,3),sum)==0 -> hay0
		if(any(hay0)) {bt.x3[1,,] <- 0.1
				print(bt.lit[33,])}	# AVISO: una fila tiene todos 0. Se pone 0.1 a una celda
		bt.x3
		}
	#
	#		Formato de individuos, variables y ocasiones
	#
	#			Data Frame de formato: bt.fx, bt.fy
	#				2 dimensiones, x,y. 
	#			Todas con las mismas columnas:
	#				Etiqueta: 		eti
	#				Color eti:		ecol
	#				Simbolo:		pch
	#				Tipo linea:		lty	
	#				Grosor linea:	lwd
	#				Posicion:		pos
	#				color trayect.:	tcol
	#				Se pinta:		type
	#				Inercia del plano: ine		se modifica en run
	#				mostrar variable: inl		se modifica en run

	formato <- function() {
		if (tclvalue(ifd)=="1") return()		# Termina la funcion
		if (1==2) b.fx <- b.fy <- NULL			# para evitar error en R CMD check
		tix <- tiy <- 1							# indicadores de cargar formato
		
		if (tclvalue(ifo)=="1") {				# chequear formatos leidos
			if(!exists("b.fx")) {print(bt.lit[122,])	# ERROR: formato filas
								tclvalue(ifo) <- "0"	# reset para dar formato
								tix <- 0		# indicador de error formato x
								}
			else 
				if (any(rownames(b.x2)!=rownames(b.fx))) 
								{print(bt.lit[122,])	# ERROR: formato filas
								tclvalue(ifo) <- "0"	# reset para dar formato
								tix <- 0		# indicador de error formato x
								}
			if(!exists("b.fy")) {print(bt.lit[123,])	# ERROR: formato columnas
								tclvalue(ifo) <- "0"	# reset para dar formato
								tiy <- 0		# indicador de error formato x
								}
			else 					
				if (any(colnames(b.x2)!=rownames(b.fy))) 
								{print(bt.lit[123,])	# ERROR: formato columnas
								tclvalue(ifo) <- "0"	# reset para dar formato
								tiy <- 0		# indicador de error formato x
								}							
			}
		else {tix <- tiy <- 0}					# indicador para preparar formatos
		#	Inicializacion de los data frame de formatos
		if (tix==0) {							# No cargamos formato desde archivo
			b.fx <<- data.frame(eti=rownames(b.x2),ecol="#0000ff",pch=16,lty=1,lwd=1,
						pos=2,tcol="#0000ff",type=1, stringsAsFactors = F)
			rownames(b.fx) <<- rownames(b.x2)
					}
		if (tiy==0) {							# No cargamos formato desde archivo
			b.fy <<- data.frame(eti=colnames(b.x2),ecol="#000000",pch=18,lty=1,lwd=2,
						pos=2,tcol="#000000",type=1, stringsAsFactors = F)
			rownames(b.fy) <<- colnames(b.x2)
			}

		#	Formato de datos					
		row0.frm <- tk2frame(tb1,relief="sunken",padding="2")
		tkpack(tk2label(row0.frm,text=bt.lit[34,],
					background="lightyellow"),fill="x")	# Editar formatos
		b.x <- b.x									# para evitar error en R CMD check
		row.erow <- function()	fix(b.fx)
		row.ecol <- function()	fix(b.fy)
		row.b <- tk2frame(row0.frm, relief="sunken")		
		row.b1 <- tk2button(row.b,text=bt.lit[10,], command=row.erow)	# Individuos
		row.b2 <- tk2button(row.b,text=bt.lit[5,], command=row.ecol)	# Variables

		tkpack(row.b1,row.b2, side="left")
		tkpack(row0.frm,row.b, side="bottom")
		tclvalue(ifd) <- "1"
	}
	
	#	llenado de listbox
	#		salida: bt.hr, bt.hc, bt.hs = lista de variables seleccionadas
	#
	llena.lbr <- function() {

		# 	Carga de iconos
		picon <- paste(path.package("dynBiplotGUI"),"images",sep="/")
		all.img <- tkimage.create("photo", file=paste(picon,"accept_list.gif",sep="/"))
		none.img <- tkimage.create("photo", file=paste(picon,"delete_list.gif",sep="/"))
		b.nada <- function(x,y)	{tkselection.clear(x,0,length(y)-1)}
		b.todo <- function(x,y)	{tkselection.set(x,0,length(y)-1)}
		
		#		Seleccion de filas
		row.frm <- tk2frame(fr.a1, relief="sunken",padding="2")
		row1.frm <- tk2frame(row.frm, relief="sunken",padding="2")
		row0 <- paste(1:nrow(b.x2), " - ", rownames(b.x2))
		lbr <- tk2listbox(row1.frm, values=row0, height = 16, 
				selectmode = "extended", scroll = "y",autoscroll = "y")
		tkpack(tk2label(row.frm,text=bt.lit[37,],background="lightcyan"	# Selecciona filas
						,tip=bt.lit[38,]),fill="x")	# Marca las filas que quieras
		
		row.b <- tk2frame(row.frm, relief="sunken")		
		row.b1 <- tk2button(row.b, tip=bt.lit[39,], text = "None", 	# Quita todas
					image=none.img, command = function() b.nada(lbr,row0))
		row.b2 <- tk2button(row.b, text = "All", tip=bt.lit[40,],	# Marca todas
					image=all.img, command = function() b.todo(lbr,row0))
		
		tkpack(lbr,fill="y", side="left")
		tkpack(row.b1,row.b2, side="right")
		tkpack(row1.frm,row.b,row.frm)
		
		#		Seleccion de columnas
		col.frm <- tk2frame(fr.a2, relief="sunken",padding="2")		
		col1.frm <- tk2frame(col.frm, relief="sunken",padding="2")
		col0 <- paste(1:ncol(b.x2), " - ", colnames(b.x2))
		lbc <- tk2listbox(col1.frm, values=col0, height = 16, 
				selectmode = "extended", scroll = "y",autoscroll = "y")
		tkpack(tk2label(col.frm,text=bt.lit[41,],background="lightcyan"	# Selecciona columnas
					,tip=bt.lit[42,]),fill="x")	# Marca las variables
		
		col.b <- tk2frame(col.frm, relief="sunken")		
		col.b1 <- tk2button(col.b,tip=bt.lit[39,],text = "None",	# Quita todas
							image=none.img, command = function() b.nada(lbc,col0))
		col.b2 <- tk2button(col.b,tip=bt.lit[40,],text = "All",		# Marca todas
							image=all.img, command = function() b.todo(lbc,col0))
		
		tkpack(lbc,fill="y", side="left")
		tkpack(col.b1,col.b2, side="right")
		tkpack(col1.frm,col.b,col.frm)
							
		#		Seleccion de situaciones
		if (tclvalue(i3v)=="1") {
			sit.frm <- tk2frame(fr.a3, relief="sunken",padding="2")		
			sit1.frm <- tk2frame(sit.frm, relief="sunken",padding="2")
			sit0 <- paste(1:dim(b.x3)[3], " - ", dimnames(b.x3)[[3]])
			lbs <- tk2listbox(sit1.frm, values=sit0, height = 16, 
					selectmode = "extended", scroll = "y",autoscroll = "y")
			tkpack(tk2label(sit.frm,text=bt.lit[137,],background="lightcyan"	# Selecciona situaciones
						,tip=bt.lit[138,]),fill="x")	# Marca las situaciones
			
			sit.b <- tk2frame(sit.frm, relief="sunken")		
			sit.b1 <- tk2button(sit.b,tip=bt.lit[39,],text = "None",	# Quita todas
								image=none.img, command = function() b.nada(lbs,sit0))
			sit.b2 <- tk2button(sit.b,tip=bt.lit[40,],text = "All",		# Marca todas
								image=all.img, command = function() b.todo(lbs,sit0))
			
			tkpack(lbs,fill="y", side="left")
			tkpack(sit.b1,sit.b2, side="right")
			tkpack(sit1.frm,sit.b,sit.frm)
		}
		
		col.sel <- function()	{bt.hc <<- as.numeric(tkcurselection(lbc)) + 1
								bt.hr <<- as.numeric(tkcurselection(lbr)) + 1
			# Validacion de elementos marcados:
			if (length(bt.hr)==0) {print(bt.lit[43,])	# ERROR: filas no seleccionadas
									return()}
			if (length(bt.hc)==0) {print(bt.lit[44,])	# ERROR: columnas no seleccionadas
									return()}
			if (tclvalue(i3v)=="1") {bt.hs <<- as.numeric(tkcurselection(lbs)) + 1
				if (length(bt.hs)<=1) {print(bt.lit[136,])	# ERROR: situaciones no seleccionadas
									return()}
									}
			if (length(bt.hr)<length(bt.hc)) {print(bt.lit[121,])	# ERROR: filas < columnas
									return()}
			tkconfigure(tb4.cb0,values= c(2:length(bt.hc)))
			tkconfigure(tb4.cb1,values= c(1:(length(bt.hc)-1)))
			tkconfigure(tb4.cb2,values= c(2:length(bt.hc)))
			
			tkconfigure(run.but,state="normal")
			tk2notetab.select(nb,bt.lit[6,])		# Analisis
			if (tclvalue(isd)!="1") {panel4.1()		# opciones del analisis
				if (tclvalue(i3v)=="1") panel4.3v()	# opciones para 3 vias
						tclvalue(isd) <- "1"}		# ya mostrado panel
			else if (tclvalue(i3v)=="1") {			# reconfigura lista de referencia
				tmp1 <- dimnames(b.x3)[[3]][bt.hs]
				tclvalue(vref) <- max(tmp1)
				tkconfigure(p4cbx,values=tmp1)
				}
			###	la captura de los valores se hace en el boton 'Run Biplot' ###
		}
													
		bok <- tk2button(tb3, tip=bt.lit[25,],text = "OK",			# OK
							 command = col.sel)
		tkpack(row.frm,col.frm, side="left", fill="y",expand=TRUE)
		tkpack(bok)
	}
		
	#		============================
	#		Analisis Biplot	<<<<<<<<------------
	#		datos de entrada: b.x2
	
	run.biplot <- function()	{	
	#	captura ejes y plano
		b.fx <- b.fx							# para evitar error en R CMD check
		b.fy <- b.fy							# para evitar error en R CMD check
		nejes <<- as.numeric(tclvalue(neje))
		dim1 <<- as.numeric(tclvalue(di1))
		dim2 <<- as.numeric(tclvalue(di2))
		if (tclvalue(i3v)=="1") {bt.t <<- tclvalue(vref)	# 3 vias
						b.x2 <<- b.x3[,,bt.t]
						bt.x3 <<- b.x3[bt.hr,bt.hc,bt.hs]	
						nl <- dim(bt.x3)[[3]]	# numero niveles seleccionado
			}
		bt.x2  <<- b.x2[bt.hr,bt.hc]	# Datos seleccionados para analizar
		bt.fxg <<- b.fx[bt.hr,]			# formatos de datos seleccionados
		bt.fyg <<- b.fy[bt.hc,]
		
		#	Calculo de medias y sd para estandarizacion
		bt.x2m <<- apply(bt.x2,2,mean)	# media de matriz de referencia
		bt.x2sd <<- apply(bt.x2,2,sd)	# sd de matriz de referencia
		if (tclvalue(i3v)=="1") {
			if (tclvalue(iref)=="2") {	# iref 1,2,3 = 1,2,3 vias
				bt.x2m <<- apply(bt.x3,2,mean)	# media del filete
				# bt.x2sd <<- apply(bt.x3,2,sd)		# sd global	DA ERROR
				bt.x2sd <<- bt.x2m		# prepara estructura
				for (i in 1:ncol(bt.x3)) bt.x2sd[i] <<- sd(as.vector(bt.x3[,i,]))
			}	else {	
					if (tclvalue(iref)=="3") {	# media del cubo
					bt.x2m[] <<- mean(bt.x3)
					bt.x2sd[] <<- sd(bt.x3)
			}	}
		}

			# si biplot con todos los datos:
		if (tclvalue(ibg)=="1")	{
			bt.x2 <<- bt.x3[,,1]
			for(i in 2:nl) bt.x2 <<- rbind(bt.x2,bt.x3[,,i])
			temp <- NULL
			for(i in 1:nl) {temp1 <- rep(dimnames(bt.x3)[[3]][i],length(bt.hr))
								temp <- c(temp,temp1)}
			rownames(bt.x2) <<- paste(bt.fxg$eti,temp,sep="")
			for(i in 2:nl) bt.fxg <<- rbind(bt.fxg,b.fx[bt.hr,])
			bt.fxg$eti <<- rownames(bt.x2)
		}	
		label.ejes <<- paste(bt.lit[45,],1:ncol(bt.x2),sep = "")	# Eje
		bt.c <<- as.numeric(tclvalue(icen)) 
		bt.e <<- as.numeric(tclvalue(iesc))
		if (bt.c != 1) bt.x2m[] <<- 0
		if (bt.e != 1) bt.x2sd[] <<- 1

		bt.xce <<- scale(bt.x2,bt.x2m,bt.x2sd)
		bt.svd <<- svd(bt.xce) 
		bt.U <<- bt.svd$u
		bt.V <<- bt.svd$v
		bt.D <<- diag(bt.svd$d)
		#
		if 	(tclvalue(tb)=="2") biplot.gh() 
			else if (tclvalue(tb)=="3") biplot.jk() 
			else biplot.hj()
		#
		dimnames(bt.a) <<- list(bt.fxg$eti,label.ejes[1:ncol(bt.x2)])
		dimnames(bt.b) <<- list(bt.fyg$eti,label.ejes[1:ncol(bt.x2)])

		# Trayectorias:
		if (tclvalue(i3v)=="1") {
			bt.zc <<-array(,c(nl,length(bt.hc),ncol(bt.x2)))	#inicializamos matriz de trayectorias
			bt.zr <<-array(,c(nl,length(bt.hc),nrow(bt.x2)))	#inicializamos matriz de trayectorias
			bt.r2c <<-array(,c(ncol(bt.x2),nl))			#inicializamos matriz de R2
			bt.Fc <<- bt.Prc <<- bt.r2c 				# matriz anova F, p-value
			bt.r2r <<-array(,c(nrow(bt.x2),nl))			#inicializamos matriz de R2
			# solo si no global
			if (tclvalue(ibg)!="1") {
				for (i in 1:length(bt.hc)) tray(i)		# Trayectorias de variables				
				for (i in 1:length(bt.hr)) trax(i)		# Trayectoria de individuos		
				leve <- dimnames(bt.x3)[[3]]			# niveles seleccionados
				dimnames(bt.zc) <<- list(leve,label.ejes[1:ncol(bt.x2)],colnames(bt.x2))
				dimnames(bt.zr) <<- list(leve,label.ejes[1:ncol(bt.x2)],rownames(bt.x2))
				dimnames(bt.r2c) <<- list(colnames(bt.x2),leve)
				dimnames(bt.Fc) <<- dimnames(bt.Prc) <<- dimnames(bt.r2c)
				dimnames(bt.r2r) <<- list(rownames(bt.x2),leve)
			}
		}		
		# Funcion de resultados numericos
		Resultados ()							# crea la matriz con los resultados
		tkconfigure(res.but,state="normal")		# Habilita boton de resultados
		# Funcion de dibujar
		plotBiplot ()
	}
	#
	#		analisis HJ
	biplot.hj <- function() {
		bt.a <<- bt.U %*% bt.D
		bt.b <<- bt.V %*% bt.D
	}
	# 
	#		analisis GH
	biplot.gh <- function() {					# Revisar todo si dejamos PCA
		if (tclvalue(ipca)=="1") {
			bt.a <<- bt.U * sqrt(nrow(bt.xce) - 1)	# PCA
			bt.b <<- bt.V %*% bt.D / sqrt(nrow(bt.xce) - 1)}
		else {
			bt.a <<- bt.U
			bt.b <<- bt.V %*% bt.D
	}	}
	#
	#		analisis JK
	biplot.jk <- function() {
		bt.a <<- bt.U %*% bt.D
		bt.b <<- bt.V
	}
	#
	#		Trayectorias
	#
	#	Para variables:
	#	Entrada: Numero de variable a dibujar
	#	
	tray <- function(v) {
		lmp <- function (xx) {						# funcion para p-value
			f <- summary(xx)$fstatistic
			p <- pf(f[1],f[2],f[3],lower.tail=F)
			attributes(p) <- NULL
			return(p)
		}
		hj <- tclvalue(tb)
		va <- colnames(bt.x2)[v]
		x <- bt.x3[,va,]
		A <- bt.a
		tmp1 <- bt.x2m[va]
		tmp2 <- bt.x2sd[va]		
		x2 <- scale(x,rep(tmp1,ncol(x)),rep(tmp2,ncol(x)))
		Z <- t(solve(t(A) %*% A) %*% t(A) %*% x2)
		if (hj=="1") Z <- Z %*% bt.D				# reescalado de HJ
		bt.zc[,,v] <<- Z
		reg <- apply(x2,2,function(x) lm(x~A[,c(dim1,dim2)]))
		bt.r2c[v,] <<- sapply(reg, function(x) summary.lm(x)$r.squared)		# todos los R2 de la variable
		bt.Fc[v,] <<- sapply(reg, function(x) summary.lm(x)$fstatistic)[1,]		# estadistico F
		bt.Prc[v,] <<- sapply(reg, function(x) lmp(x))			# p-valor de F
	}
	#	Para individuos
	#	Entrada: Numero de individuo a dibujar
	#
	trax <- function(v)	{
		hj <- tclvalue(tb)
		va <- rownames(bt.x2)[v]
		x <- t(bt.x3[va,,])
		x2 <- matrix(,dim(bt.x3)[[3]],ncol(bt.x3))
		B <- bt.b	
		x2 <- scale((x),bt.x2m,bt.x2sd)				# siempre con la referencia
			
		Z <- t(solve(t(B) %*% B) %*% t(B) %*% t(x2))	# traspuesta de x2
		if (hj=="1") Z <- Z%*%bt.D						# reescalado de HJ
		bt.zr[,,v] <<- Z	
	}
	#	======================================================
	    # Construye la matriz de resultados
    #  
    Resultados <- function()
    { 
	#	Se guardan los datos integros. Los redondeos a la hora de mostrarlos.
		vartotal <- sum(bt.svd$d^2)
        bt.varexpl <<- (bt.svd$d^2/vartotal) *100
        sc.a <- rowSums(bt.a^2)
        CRFE.a <- round(((bt.a^2) * 1000)/sc.a,0)
        sc.b <- rowSums(bt.b^2)
        CRFE.b <- round(((bt.b^2) * 1000)/ sc.b,0)
		
		bt.res.vp <<- cbind(bt.svd$d, bt.varexpl, cumsum(bt.varexpl))	# poner etiqueta
		bt.res.a  <<- bt.a[,1:nejes]		# coordenadas de filas
		bt.res.b  <<- bt.b[,1:nejes]		# coordenadas de columnas
		bt.res.cr <<- CRFE.a[,1:nejes]		# contribuciones de filas
		bt.res.cc <<- CRFE.b[,1:nejes]		# contribuciones de columnas
		if (tclvalue(i3v)=="1" & tclvalue(ibg)!="1") {
			bt.res.ty <<- round(bt.zc[,1:nejes,],3)
			bt.res.tx <<- round(bt.zr[,1:nejes,],3)
			bt.res.r2y <<- round(bt.r2c,4)
			bt.res.Fy <<- round(bt.Fc,3)
			bt.res.Pry <<- round(bt.Prc,4)
		}
	}

	#	======================================================
	    # Funcion fichero de resultados
    #  
	# 	Mostrar los resultados
	#
    ShowRes <- function()
      {                        
        cat(bt.lit[2,], "\n", file = "Results.txt")				#cabecera: Biplot Dinamico
		if(tclvalue(icen)==1) cat (bt.lit[46,],file="Results.txt", append=TRUE)	# Centrado
			else cat(bt.lit[47,], file="Results.txt", append=TRUE)		# No centrado
		if (tclvalue(iesc)==1) cat(bt.lit[48,],"\n", file="Results.txt", append=TRUE)	# Escalado
			else cat(bt.lit[49,],"\n", file="Results.txt", append=TRUE)	# No escalado.
		if 			(tclvalue(tb)=="2") cat("GH Biplot","\n", file="Results.txt", append=TRUE) 
			else if (tclvalue(tb)=="3") cat("JK Biplot","\n", file="Results.txt", append=TRUE) 
					else cat("HJ Biplot","\n", file="Results.txt", append=TRUE)
			
		if (tclvalue(i3v)=="1") cat(bt.lit[50,],bt.t,"\n","\n",file="Results.txt",append=TRUE)	# Referencia
		# varianza	
		cat(bt.lit[51,],"\n", file="Results.txt", append=TRUE) #valores propios
		cab1 <- c(bt.lit[52,], bt.lit[53,], bt.lit[54,])		# vp, var, acum
        write.table(cbind(round(bt.svd$d,3),round(bt.varexpl,3),
			round(cumsum(bt.varexpl),3)),
			file="Results.txt",quote=FALSE,sep="\t",dec=",", 
			append=TRUE, col.names=cab1)
        # coordenadas
		cat("\n",bt.lit[55,],"\n", file="Results.txt",append=TRUE)  #coordenadas filas  
        write.table(round(bt.a[,1:nejes],3),file="Results.txt",quote=FALSE,
			sep="\t",dec=",", append=TRUE)
        cat("\n",bt.lit[56,],"\n",file="Results.txt",append=TRUE)   #coordenadas columnas  
        write.table(round(bt.b[,1:nejes],3),file="Results.txt",quote=FALSE,
			sep="\t",dec=",",append=TRUE)        
        cat("\n",file="Results.txt",append=TRUE)
        # contribuciones
		cat(bt.lit[57,],"\n",file="Results.txt",append=TRUE)	# contribuc. relativas
        cat(bt.lit[58,],"\n",file="Results.txt",append=TRUE)	# contribuciones filas
        write.table(bt.res.cr,file="Results.txt",quote=FALSE,sep="\t",dec=",",append=TRUE)
        cat(bt.lit[59,],"\n",file="Results.txt",append=TRUE)	#contribuciones columnas
        write.table(bt.res.cc,file="Results.txt",quote=FALSE,sep="\t",dec=",",append=TRUE)
		# Trayectorias
		if (tclvalue(i3v)=="1" & tclvalue(ibg)!="1") {
			cat("\n",file="Results.txt",append=TRUE)
			cat(bt.lit[60,],"\n", file="Results.txt", append=TRUE)	# Coordenadas de las TRAYECTORIAS
			#	variables
			for (i in 1:dim(bt.res.ty)[3]) {
			cat("\n",bt.lit[61,],colnames(bt.x2)[i],"\n", file="Results.txt",append=TRUE) 	# Variable
			write.table(bt.res.ty[,,i],file="Results.txt",quote=FALSE,sep="\t",
			dec=",", append=TRUE, col.names=TRUE,row.names=TRUE)
			}
			#	individuos
			for (i in 1:dim(bt.res.tx)[3]) {
			cat("\n",bt.lit[62,],rownames(bt.x2)[i],"\n",file="Results.txt",append=TRUE) 	# Individuo
			write.table(bt.res.tx[,,i],file="Results.txt",quote=FALSE,sep="\t",
			dec=",", append=TRUE, col.names=TRUE,row.names=TRUE)
			}
			cat("\n",paste(bt.lit[124,],dim1,"-",dim2,sep=""),sep="",
				"\n",file="Results.txt",append=TRUE)    	# R2 - columnas, plano
			write.table(bt.res.r2y,file="Results.txt",quote=FALSE,sep="\t",dec=",",append=TRUE)	
			cat("\n",paste(bt.lit[127,],dim1,"-",dim2,sep=""),sep="",
				"\n",file="Results.txt",append=TRUE)    	# F anova
			write.table(bt.res.Fy,file="Results.txt",quote=FALSE,sep="\t",dec=",",append=TRUE)	
			cat("\n",paste(bt.lit[128,],dim1,"-",dim2,sep=""),sep="",
				"\n",file="Results.txt",append=TRUE)    	# p-valor de F
			write.table(bt.res.Pry,file="Results.txt",
				quote=FALSE,sep="\t",dec=",",append=TRUE)
		}
		file.show("Results.txt", title=bt.lit[63,])	# Biplot Dinamico - Resultados
	}
    #
    # 	Funcion de opciones de grafico
    #
	VeOpciones <- function()
	{
	if (1==2) b.fx <- b.fy <- NULL		# para evitar error en R CMD check
	#	Titulos
	frf1 <- tk2frame(fr.f2,borderwidth=1)
	frf2 <- tk2frame(fr.f2,borderwidth=1)
	frf3 <- tk2frame(fr.f2,borderwidth=1)
	tkpack(tk2checkbutton(frf1,variable=it1,tip=bt.lit[110,]),		# Mostrar
		tk2label(frf1,text=bt.lit[64,], width="11"), 	# Titulo
		tk2entry(frf1, width="40", textvariable=vtit), side="left")
	tkpack(tk2checkbutton(frf2,variable=it2,tip=bt.lit[110,]),		# Mostrar
		tk2label(frf2,text=bt.lit[65,], width="11"), 	# Subtitulo
		tk2entry(frf2, width="40", textvariable=vsub), side="left")
	tkpack(tk2label(frf3,text=paste(bt.lit[135,],", H: "),	# Escala de la ventana
		tip="Horizontal (1.4 < h < 2.5)"),
		tk2entry(frf3,textvariable=wout1,width="3"), 
		tk2label(frf3,text=" V: ",tip="Vertical (1.4 < v < 2.5)"),
		tk2entry(frf3,textvariable=wout2,width="3"), 
		side="left")
	tkpack(frf1, side="top")
	# if (tclvalue(i3v)==0) tkpack(frf2)
	tkpack(frf2,frf3)

	#	Para formato de datos
		# Solapa individuos tb1
	ff21 <- tk2frame(tb21, relief="sunken",padding="2")
	tkpack(ff21,fill="x", expand=0)
	tkpack(ttklabel(ff21, text=bt.lit[69,],	# Selecciona filas
			background="peachpuff"),fill="x")
			
	ff21.lb2 <- tk2listbox(ff21, values=b.fx$eti, height = 10, 
			tip=bt.lit[66,],			# Puedes seleccionar mas de 1
			selectmode = "extended",
			scroll = "y",autoscroll = "y")
	tkpack(ff21.lb2,fill="x")
	
	sel.boton <- function () { 
		tclvalue(vf) <- "1"				# filas
		tmplb <- as.numeric(tkcurselection(ff21.lb2)) + 1
		bt.fx <<- b.fx
		bt.tf <<- b.fx[tmplb,]
		Formatos(bt.tf)	
	}
	tb21.but <- tk2button(ff21,text=bt.lit[67,],command=sel.boton)	# Seleccionar
	tk2tip(tb21.but, bt.lit[68,])		# Marca fila y selecciona
	tkpack(tb21.but)
	Formatos1()
	
			# Solapa variables tb2
	ff22 <- tk2frame(tb22, relief="sunken",padding="2")
	tkpack(ff22,fill="x", expand=0)
	tkpack(ttklabel(ff22, text=bt.lit[69,],	# Selecciona filas
			background="peachpuff"),fill="x")
	
	ff22.lb2 <- tk2listbox(ff22, values=b.fy$eti, height = 10, 
			tip=bt.lit[66,],			# Puedes seleccionar mas de 1
			selectmode = "extended",
			scroll = "y",autoscroll = "y")
	tkpack(ff22.lb2,fill="x")
	
	sel2.boton <- function () { 
		# b.fy <- NULL					# para evitar error en R CMD check
		tclvalue(vf) <- "2"				# columnas
		tmplb <- as.numeric(tkcurselection(ff22.lb2)) + 1
		
		bt.fx <<- b.fy
		bt.tf <<- b.fy[tmplb,]
		Formatos(bt.tf)	
	}
	tb22.but <- tk2button(ff22, text=bt.lit[67,],command=sel2.boton)	# Selecciona
	tk2tip(tb22.but, bt.lit[68,])		# Marca fila y selecciona
	tkpack(tb22.but)
	
	#	Opciones grafico del panel 4
	#	etiquetas
	fr.o1 <- tk2frame(fr.s5,padding="2",relief="sunken")
	tkpack(tk2label(fr.o1,text=bt.lit[70,],background="honeydew"),	# Etiquetas para filas
		tk2checkbutton(fr.o1, variable=ietx), 
		tk2label(fr.o1,text=bt.lit[71,],background="honeydew"),	# para columnas
		tk2checkbutton(fr.o1, variable=iety), 
		side="left")
	fr.o2 <- tk2frame(fr.s5,padding="2",relief="flat")
	fr.o3 <- tk2frame(fr.s5,padding="2",relief="flat")	
	# 	inercia 
	la1 <- tk2label(fr.o2,text=bt.lit[72,],width=20,background="honeydew")	# Inercia filas
	sc1 <- tk2scale(fr.o2,tip=bt.lit[73,],from = 0, to = 1000,	# filas
			variable=vinx, length=200)
	e41 <- tk2entry(fr.o2, textvariable=vinx, width=4)
	tkpack(la1,sc1,e41, side="left", fill="x")
	la2 <- tk2label(fr.o3,text=bt.lit[74,], width=20,background="honeydew")	# Inercia columnas
	sc2 <- tk2scale(fr.o3,tip=bt.lit[75,],from = 0, to = 1000,	# columnas
			variable=viny, length=200)
	e42 <- tk2entry(fr.o3, textvariable=viny, width=4)
	tkpack(la2,sc2,e42, side="left", fill="x")
	tkpack(fr.o1, fr.o2, fr.o3, side="top", fill="x")
	}

	# opciones de 3 vias
	panel4.3v <- function() {
		fr27 <- tk2frame(fr.s4, relief="raised",padding="2")
		tmp1 <- dimnames(b.x3)[[3]][bt.hs]			# etiquetas de situaciones
		tclvalue(vref) <- max(tmp1)
		p4cbx <<- tk2combobox(fr27,values=tmp1,textvariable=vref,width=8) # situaciones seleccionadas
		tkpack(tk2label(fr27, text=bt.lit[76,],background="palegreen"),	# Referencia
			p4cbx,
			tk2label(fr27, text=bt.lit[77,],background="palegreen"),	# Biplot global
			tk2checkbutton(fr27, variable=ibg),
			side="left")
		fr27.2 <- tk2frame(fr.s4, padding="2")
		tkpack(tk2label(fr27.2, text=bt.lit[78,], width=15,	# Trayectoria - filas
							background="palegreen"),
			tk2checkbutton(fr27.2, variable=itrr), 
			tk2label(fr27.2,text=bt.lit[79,],background="palegreen"),	# columnas
			tk2checkbutton(fr27.2, variable=itrc), 
			tk2label(fr27.2,text=bt.lit[130,],background="palegreen"),	# p-valor
			tk2checkbutton(fr27.2, variable=ipv), 
			tk2entry(fr27.2, textvariable=pval, width=4),
			side="left", fill="x")
		fr27.3 <- tk2frame(fr.s4, padding="2")
		tkpack(tk2label(fr27.3, text=bt.lit[80,], width=15,	# Etiquetas - filas
							background="palegreen"),
			tk2checkbutton(fr27.3, variable=ietzr), 
			tk2label(fr27.3,text=bt.lit[79,],background="palegreen"),	# columnas
			tk2checkbutton(fr27.3, variable=ietzc), 
			tk2label(fr27.3,text=bt.lit[131,],background="palegreen",	# Concatenar
					tip=bt.lit[132,]),	# Etiquetas = nombre de variable + situacion
			tk2checkbutton(fr27.3, variable=ivs), 
			side="left", fill="x")
		tkpack(fr27, fr27.2, fr27.3, side="top")
		}
	#
	#	Funcion formato de datos
	#
	Formatos <- function (x) {
		# Entrada: x = b.fx[i:j,]	# todos los datos de formato
		# Esta definido:
		#	b.fx <<- data.frame(eti=rownames(b.x2),ecol="#0000ff",pch=16,lwd=1,pos=2,
		#			tcol="#0000ff",type=1, stringsAsFactors = F)
	
		eti1 <- tclVar(x$eti[1])
		if (nrow(x)!=1) tkconfigure(bt.en1,state="disable")
			else {	tkconfigure(bt.en1,textvariable=eti1)
				tkconfigure(bt.en1,state="enable")}
		bt.micol <<- x$ecol[1]
		tkconfigure(bt.cv1,bg=x$ecol[1])
		pch1 <- tclVar(x$pch[1])
		tkconfigure(bt.cbb1,textvariable=pch1)
		lwd1 <- tclVar(x$lwd[1])
		tkconfigure(bt.en2,textvariable=lwd1)
		pos1 <- tclVar(x$pos[1])
		tkconfigure(bt.cbb2,textvariable=pos1)
		bt.micol2 <<- x$tcol[1]
		tkconfigure(bt.cv2,bg=x$tcol[1])
		tclvalue(imos) <- x$type[1]
		tkconfigure(bt.typ1,variable=imos)
		tkconfigure(bt.bu3,state="enable")
	}
	
	camb <- function () { 
		x <- bt.tf
		b.fx <- b.fx					# para evitar errores en R CMD check
		b.fy <- b.fy
		if (nrow(x)==1) x$eti <- tclvalue(tkget(bt.en1))
		x$ecol <- bt.micol
		x$pch <- as.numeric(tkget(bt.cbb1))
		x$lwd <- tclvalue(tkget(bt.en2))
		x$pos <- tclvalue(tkget(bt.cbb2))
		x$type <- tclvalue(imos)
		if (tclvalue(i3v)=="1") x$tcol <- bt.micol2
		rx <- rownames(x)
		tbf <- t(bt.fx)					# formatos de fila o columna
		tbf[,rx] <- t(x)
		tmp <- as.data.frame(t(tbf), stringsAsFactors = F)
		if (tclvalue(vf)=="1") b.fx <<- tmp
			else b.fy <<- tmp
		tkconfigure(bt.bu3,state="disable")
		}
		
		# Dibuja los campos la primera vez
		#	Para formato de datos
	Formatos1 <- function ()
	{
			# Solapa individuos tb1
		ff20 <- tk2frame(fr.f12, relief="sunken",padding="2")
		tkpack(ff20,fill="x", expand=0)
		ff21.2 <- tk2frame(ff20, relief="sunken",padding="2")
		ff21.3 <- tk2frame(ff20, relief="sunken",padding="2")
		tkpack(ff21.3,ff21.2,side="top",fill="x")
		bt.lf1 <-tk2labelframe(ff21.2,text=bt.lit[81,])	# Cambia formato
		bt.lf2 <-tk2labelframe(ff21.3,text=bt.lit[82,])	# Multiple
		tkpack(bt.lf1,bt.lf2, side="left",fill="x")

		fagr <- function() {
			b.x <- b.x					# para evitar error en R CMD check
			b.fx <- b.fx				# para evitar error en R CMD check
			if(iagr==1) {				# para alternar el uso del boton
				tkconfigure(eagr,values=sort(unique(b.x[[tclvalue(vagr)]])))
				tkconfigure(bagr,text="sel>")
				iagr <<- 2
			} else {
				tmplb <- rownames(b.x2)[b.x2[,tclvalue(vagr)]==tclvalue(vag1)]
				tkconfigure(bagr,text="<sel")
				iagr <<- 1
				bt.fx <<- b.fx
				bt.tf <<- b.fx[tmplb,]
				tclvalue(vf) <- "1"
				Formatos(bt.tf)
		}	}
		bagr <- tk2button(bt.lf2,text="<sel",command=fagr,width=4)
		eagr <- tk2combobox(bt.lf2,textvariable=vag1,width=7)
		tkpack(tk2label(bt.lf2,text=bt.lit[61,],width=7,background="peachpuff"),	# Variable
				tk2combobox(bt.lf2,values=colnames(b.x2),textvariable=vagr,width=5),
				bagr,eagr,side="left",fill="x")
		tk2tip(bt.lf2, bt.lit[83,])			# Variable de agrupacion
		
			# Frames para los campos 
		fr1 <- tk2frame(bt.lf1)
		fr2 <- tk2frame(bt.lf1)
		fr3 <- tk2frame(bt.lf1)
		fr4 <- tk2frame(bt.lf1)
		fr5 <- tk2frame(bt.lf1)
		fr6 <- tk2frame(bt.lf1)
		fr7 <- tk2frame(bt.lf1)
		tkpack(fr1,fr2,fr3,fr4,fr5,fr6,fr7, side="top",fill="x")
		
		eti1 <- tclVar("")
		bt.en1 <<- tk2entry(fr1,textvariable=eti1,width="20")
		tkpack(tk2label(fr1,text=bt.lit[84,],width = 10,background="papayawhip"),	# Etiqueta
				bt.en1, side="left",fill="x")
		
		colb <- function () { 
			bt.micol <<- tclvalue(tcl("tk_chooseColor",initialcolor="#0000ff"))
			tkconfigure(bt.cv1,bg=bt.micol)
		}
		
		bt.micol <<- "#0000ff"
		bt.cv1 <<- tk2canvas(fr2,bg=bt.micol,width="20",height="20",relief="raised")
		cv1b <- tk2button(fr2, text="+", width = 1,command = colb)
		tkpack(tk2label(fr2,text=bt.lit[85,],width = 10,background="papayawhip"), 	# Color
				bt.cv1, cv1b, side="left",fill="both")
		pch1 <- tclVar(16)
		bt.cbb1 <<- tk2combobox(fr3,values=0:25,width="3", textvariable=pch1)
		tkpack(tk2label(fr3,text=bt.lit[86,],width = 10,background="papayawhip"),	# Forma
				bt.cbb1, side="left")
		
		lwd1 <- tclVar(1)
		bt.en2 <<- tk2entry(fr4,textvariable=lwd1,width="4")
		tkpack(tk2label(fr4,text=bt.lit[87,],width = 10,background="papayawhip"), bt.en2, side="left")	# Peso
		
		pos1 <- tclVar(2)		
		bt.cbb2 <<- tk2combobox(fr5,values=1:4,width="2", textvariable=pos1)
		tkpack(tk2label(fr5,text=bt.lit[88,], tip="1-sur,2-o,3-n,4-e",		# Posicion
				width = 10,background="papayawhip"), bt.cbb2, side="left")
		
		bt.typ1 <<- tk2checkbutton(fr7, variable=imos)
		tkpack(tk2label(fr7,text=bt.lit[89,],width = 10,background="papayawhip"),	# Mostrar
				bt.typ1, side="left")
				
		bt.bu3 <<- tk2button(bt.lf1, text=bt.lit[90,],width=8,	# Cambia
				command = camb )
		tkconfigure(bt.bu3,state="disable")
		tkpack(bt.bu3)	
		
		#	Si datos de 3 vias, mostramos color de las trayectorias
		col2b <- function () { 
			bt.micol2 <<- tclvalue(tcl("tk_chooseColor",initialcolor="#0000ff"))
			tkconfigure(bt.cv2,bg=bt.micol2)
		}
		col2 <- "#0000ff"
		bt.micol2 <<- "#0000ff"
		bt.cv2 <<- tk2canvas(fr6,bg=col2,width="20",height="20",relief="raised")
		cv2b <- tk2button(fr6, text="+", width=1,command = col2b)
		if (tclvalue(i3v)=="1") {
			tkpack(tk2label(fr6,text=bt.lit[91,],width = 10,background="papayawhip"),	# Trayectoria
					bt.cv2, cv2b, side="left",fill="both")
		}
	}
	#
	#	Funcion de ayuda
	#
	ayuda <- function(topic) {			# por ejemplo: topic <- "panelAnalysis"
		if(lang=="es") topic <- paste(topic,"_",lang,sep="")
			else topic <- paste(topic,"_en",sep="")
		print(help(topic))
	}

	#
	#	Generarcion de paneles
	#
	panel1 <- function()
	{
		fr.d0 <-tk2labelframe(tb1,text=bt.lit[92,],relief="raised",	# Opcion
					borderwidth=2,padding="2")			# para 3 vias
		fr.d1 <-tk2labelframe(tb1,text=bt.lit[93,])		# leer datos
		tb1.t0 <- tk2frame(tb1)
		tb1.t1 <- tk2label(tb1.t0,text=bt.lit[94,],background="yellow",width=65)	# Carga de datos
		topic <- "panelData"
		tb1.t2 <- tk2button(tb1.t0,text="?",width=4,
					command=function() ayuda(topic),tip=bt.lit[103,])	# Ayuda
		tkpack(tb1.t1, side="left", fill="x")
		tkpack(tb1.t2, side="right")
		tkpack(tb1.t0,fr.d0, fr.d1, fr.d2, side="top", fill="x")

		fr11 <- tk2frame(fr.d0)	
		tkpack(tk2label(fr11,text=bt.lit[95,],background="lightyellow"), 	# Son datos de 3 vias
				tk2checkbutton(fr11, variable=i3v), side="left",fill="both")
		tk2tip(fr11, bt.lit[96,])	# Marca antes de buscar el fichero
		fr10 <- tk2frame(fr.d0) 
		tkpack(tk2label(fr10,text=bt.lit[97,],background="lightyellow"),	# Carga formatos desde R
				tk2checkbutton(fr10, variable=ifo), side="left")
		tk2tip(fr10, bt.lit[98,])	# Archivo filas: b.fx \nArchivo columnas: b.fy
		#
		fr1 <- tk2frame(fr.d1,relief="sunken",borderwidth=2,padding="2")
		tkpack(tk2label(fr1,text=bt.lit[99,],background="lightyellow"))	# Selecciona tipo de fichero
		# if (.Platform$OS.type=="windows") tmp="enable"	
			# else tmp="disable"
		fr1a <- tk2frame(fr1)
		fr1b <- tk2frame(fr1)
		tkpack(fr1a,fr1b,side="top")
		tkpack(
				tk2radiobutton(fr1a, command=leer.archivos, text="R",
							 value=1, variable=bt.leer), 
				tk2radiobutton(fr1a, command=leer.archivos, text="Excel",
							 # value=1, variable=bt.leer,state=tmp), 
							 value=2, variable=bt.leer), 
				tk2radiobutton(fr1b, command=leer.archivos, text="txt",
							 value=4, variable=bt.leer),			 
				tk2radiobutton(fr1b, command=leer.archivos, text="CSV",
							 value=3, variable=bt.leer),	
				tk2radiobutton(fr1a, command=leer.archivos, text="SPSS",
							 value=5, variable=bt.leer),			 
				tk2radiobutton(fr1b, command=leer.archivos, text=bt.lit[100,],	# portapapeles
							 value=6, variable=bt.leer), 
				side="left")	
		#
		tkpack(fr11, fr10, fr1, side="top")
	}
	
	panel2 <- function()
	{
		sep1 <- tk2separator(tb2)
		tb2.t0 <- tk2frame(tb2)
		tb2.t1 <- tk2label(tb2.t0,text=bt.lit[101,],background="salmon",width=65)	# Formato de datos
		topic <- "panelFormat"
		tb2.t2 <- tk2button(tb2.t0,text="?",width=4,
						command=function() ayuda(topic),tip=bt.lit[103,])	# Ayuda
		tkpack(tb2.t1, side="left", fill="x")
		tkpack(tb2.t2, side="right")
		tkpack(tb2.t0, fr.f2, sep1, fr.f1, side="top", fill="x")
		tkpack(fr.f11, fr.f12, side="left")
		tkpack(nb2,fill="both")
	}
	
	panel3 <- function()
	{
		tb3.t0 <- tk2frame(tb3)
		tb3.t1 <- tk2label(tb3.t0,text=bt.lit[102,],	# Seleccion de filas y columnas
						background="lightblue",width=65)
		topic <- "panelVariables"
		tb3.t2 <<- tk2button(tb3.t0,text="?",width=4, 
						command=function() ayuda(topic),tip=bt.lit[103,])	# Ayuda
		tkpack(tb3.t1, side="left", fill="x")
		tkpack(tb3.t2, side="right")
		tkpack(tb3.t0, fill="x")
		tkpack(fr.a0,side="top")
		tkpack(fr.a1, fr.a2, fr.a3, side="left")
	}
	
	panel4 <- function()
	{
		tb4.t0 <- tk2frame(tb4)
		tb4.t1 <- tk2label(tb4.t0,text=bt.lit[104,],	# Opciones de analisis
					background="lightgreen",width=65)
		topic <- "panelAnalysis"			
		tb4.t2 <- tk2button(tb4.t0,text="?",width=4,
						command=function() ayuda(topic),tip=bt.lit[103,])	# Ayuda
		tkpack(tb4.t1, side="left", fill="x")
		tkpack(tb4.t2, side="right")
		tkpack(tb4.t0, fill="x")
	}
	
	panel4.1 <- function()
	{
		#	panel 4 - para seleccion de salida
		#
		tkpack(fr.s1, fr.s2, fr.s3, fr.s5, fr.s4, side="top", fill="x")
		fr21 <- tk2frame(fr.s1,relief="groove", borderwidth=2,padding="2")
		tkpack(tk2label(fr21,text=bt.lit[105,],background="honeydew"),	# Centrado
				tk2checkbutton(fr21, variable=icen), side="left")
		tkpack(tk2label(fr21,text=bt.lit[106,],background="honeydew"),	# Escalado
				tk2checkbutton(fr21, variable=iesc), side="left")
			
		if (tclvalue(i3v)!=0)
			tkpack(tk2label(fr21,text=bt.lit[107,],background="palegreen"),	# Sobre Ref
					tk2radiobutton(fr21, text="1", value=1, variable=iref), 
					tk2radiobutton(fr21, text="2", value=2, variable=iref), 
					tk2radiobutton(fr21, text="3", value=3, variable=iref), side="left")

		fr2 <- tk2frame(fr.s2,relief="sunken",borderwidth=2,padding="2") 
		tkpack(tk2radiobutton(fr2, text="HJ", value=1, variable=tb),
			tk2radiobutton(fr2, text="GH", value=2, variable=tb),
			tk2radiobutton(fr2, text="JK", value=3, variable=tb), side="left",fill="x")
		
		#	la captura de los valores se hace en el boton 'Biplot'.
		#
		tkpack(fr21)
		tkpack(fr2)
		
		# 	Panel 4 - para la seleccion de ejes
		#	se dan valores iniciales hasta obtener los reales
			#	Seleccion de ejes y plano
		tkpack(tk2label(fr25,text=bt.lit[108,],background="honeydew"), tb4.cb0, side="left")	# Ejes
		tkpack(tk2label(fr25,text=bt.lit[109,],background="honeydew"), tb4.cb1,tb4.cb2, side="left")	# Plano
		tkpack(tk2label(fr25,text=bt.lit[110,],background="honeydew"), tk2checkbutton(fr25, variable=ieje), side="left")	# Mostrar
		tkpack(fr25)
	}
	#
	#	Funciones de grafico
	#
	#	elementos para grafico
	#		x = datos;	fx = formato de datos
	#
	Flechas	<- function (x,fx){
		arrows(0,0, x[,1], x[,2],
			col=fx$ecol,length=0.1,angle=30,lty=as.numeric(fx$lty),lwd=fx$lwd)
	}
	Texto <- function (x,fx) {
		text(x[,1], x[,2], labels=fx$eti, col=fx$ecol,
			cex=0.7, pos=fx$pos)
	}
	Puntos <- function (x,fx) {
		points(x[,1], x[,2], col=fx$ecol,
			cex=.7, pch=as.numeric(fx$pch))
	}
	Trayect <- function (x,fx) {
		lines(x[,1], x[,2], col=fx$tcol,lty="dotted",
			cex=.5, type="o")
	}
	Textot <- function (x,fx,e) {
		text(x[,1], x[,2], labels=e, col=fx$tcol,cex=0.5, pos=3)
	}

	#	dibujo biplot
	#
    
    plotBiplot1 <- function(screen = TRUE) 
    {     
		if (tclvalue(it1)==1) wintitle <- tclvalue(vtit)    
			else {wintitle <- NULL
				tmp <- par("mar") 
				tmp[3] <- 1
				par(mar=tmp)
				}
		col.title <- "black"
		background <- "white"    
		params <- par(bg="white")
		col.col <- "black"
		col.row <- "blue"

		if (tclvalue(it2)==1)     
			if (tclvalue(i3v)=="1" & tclvalue(ibg)!="1") subtitulo <- paste("Ref.: ",bt.t)
				else	subtitulo <- tclvalue(vsub)
			else {subtitulo <- NULL
				tmp <- par("mar") 
				tmp[1] <- 4
				par(mar=tmp)
				}
		plot(bt.a, main = wintitle, type = "n",
			col.main = col.title,  family="sans",
			font.main=4,		
			sub=subtitulo, col.sub=col.title, cex.sub=0.8, cex.lab=.8,
			xlab=paste(label.ejes[dim1],": ",round(bt.varexpl[dim1],2),"%",sep=" "),
		    ylab=paste(label.ejes[dim2],": ",round(bt.varexpl[dim2],2),"%",sep=" "),
			xlim=bt.limx1 * 1.1,		# incremento 10%
			ylim=bt.limy1 * 1.1,
			cex.axis=0.8
			)

		if (tclvalue(ieje)== "1") abline (h=0,v=0,lty="dotted")     
	  
		# inercia del plano:	### chequear que al menos hay un elemento ###
		bt.fxg$ine <<- rowSums(bt.res.cr[,c(dim1,dim2)])
		bt.fyg$ine <<- rowSums(bt.res.cc[,c(dim1,dim2)])
		vix <- as.numeric(tclvalue(vinx))
		viy <- as.numeric(tclvalue(viny))
		if (any(bt.fxg$ine > vix)== TRUE) bt.fxg$inl <<- bt.fxg$ine > vix
			else {print(bt.lit[125,])		# La inercia escogida excluye todas las filas
					return()}
		if (any(bt.fyg$ine > viy)== TRUE) bt.fyg$inl <<- bt.fyg$ine > viy
			else {print(bt.lit[126,])		# La inercia escogida excluye todas las columnas
					return()}
		bt.fxg$inl[bt.fxg$type==0] <- F		# si opcion no pintar
		bt.fyg$inl[bt.fyg$type==0] <- F
	  
		# individuos
		tmp <- bt.a[,c(dim1,dim2)]			# Datos
		tmp <- t(t(tmp)[,bt.fxg$inl])
		tmpf <- as.data.frame(t(t(bt.fxg)[,bt.fxg$inl]),stringsAsFactors=F)	# Formato
		if (tclvalue(igx)==1) tmp[,1] <- -1*tmp[,1]		# rota eje x
		if (tclvalue(igy)==1) tmp[,2] <- -1*tmp[,2]		# rota eje y
		tmp <- as.numeric(tclvalue(vesc))* tmp			# reescala x
		Puntos(tmp,tmpf)
		if (tclvalue(ietx)=="1" & tclvalue(ibg)!="1") 
			Texto(tmp,tmpf)	
	  
		# variables
		tmp <- bt.b[,c(dim1,dim2)]			# Datos
		tmp <- t(t(tmp)[,bt.fyg$inl])
		tmpf <- as.data.frame(t(t(bt.fyg)[,bt.fyg$inl]),stringsAsFactors=F)	# Formato
		if (tclvalue(igx)==1) tmp[,1] <- -1*tmp[,1]		# rota eje x
		if (tclvalue(igy)==1) tmp[,2] <- -1*tmp[,2]		# rota eje y
		Flechas(tmp,tmpf)
		if (tclvalue(iety)=="1") 
			Texto(tmp,tmpf)

		#
		#	Trayectorias
		if (tclvalue(i3v)=="1" & tclvalue(ibg)!="1")  {
		#	de variables
			if (tclvalue(itrc)=="1")
			for (i in (1:dim(bt.res.ty)[3])[bt.fyg$inl]) {
				if(tclvalue(ipv)=="1") 		# indicador de p-valor
					{tpval <- as.numeric(tclvalue(pval))
					tmp <- bt.res.ty[,c(dim1,dim2),i][bt.Prc[i,]<tpval,]	# p-valor significativo
					if(length(tmp)==0) next
					}
				else	tmp <- bt.res.ty[,c(dim1,dim2),i]
				tmpf <- bt.fyg[i,]
				if (tclvalue(igx)==1) tmp[,1] <- -1*tmp[,1]		# rota eje x
				if (tclvalue(igy)==1) tmp[,2] <- -1*tmp[,2]		# rota eje y
				Trayect(tmp,tmpf)
				if (tclvalue(ietzc)=="1") {
					if(length(tmp) == 0) return()
					else {
						tmpe <- bt.fyg$eti[i]
						if (tclvalue(ivs)==1) tmpe <- paste(tmpe,rownames(tmp),sep="")	# concatenar etiquetas
							else tmpe <- rownames(tmp)
						Textot(tmp,tmpf,tmpe) }
			}}
		#	de individuos
			if (tclvalue(itrr)=="1")
			for (i in (1:dim(bt.res.tx)[3])[bt.fxg$inl]) {
				tmp <- bt.res.tx[,c(dim1,dim2),i]
				tmpf <- bt.fxg[i,]
				if (tclvalue(igx)==1) tmp[,1] <- -1*tmp[,1]		# rota eje x
				if (tclvalue(igy)==1) tmp[,2] <- -1*tmp[,2]		# rota eje y
				tmp <- as.numeric(tclvalue(vesc))* tmp			# reescala x
				Trayect(tmp,tmpf)
				if (tclvalue(ietzr)=="1") {
					tmpe <- bt.fxg$eti[i]
					if (tclvalue(ivs)==1) tmpe <- paste(tmpe,rownames(tmp),sep="")	# concatenar etiquetas
						else tmpe <- rownames(tmp)
					Textot(tmp,tmpf,tmpe)			
			}}
		}

	#	de individuos si biplot global
	if (tclvalue(ibg)=="1")
		if (tclvalue(itrr)=="1") {
			tmp <- length(bt.hr)
			tmp3 <- length(bt.hs)
			for (i in (1:tmp)[bt.fxg$inl[1:tmp]]) {	
				tma <- matrix(,tmp3,2)
				for (j in 1:tmp3) {
					tma[j,] <- bt.res.a[i+(j-1)*tmp,c(dim1,dim2)]}
				if (tclvalue(igx)==1) tma[,1] <- -1*tma[,1]		# rota eje x
				if (tclvalue(igy)==1) tma[,2] <- -1*tma[,2]		# rota eje y
				tma <- as.numeric(tclvalue(vesc))* tma			# reescala trayectoria
				Trayect(tma,bt.fxg)
				# lines(tma, lty="dotted", type="l", col=bt.fxg$tcol[i])
				if (tclvalue(ietzr)=="1")
					text(tma[,1],tma[,2], labels=bt.fxg$eti[i],col=bt.fxg$tcol[i],cex=.6, pos=3)
				if (tclvalue(ietx)=="1")			# solo sacamos 1 etiqueta
					text(tma[1,1],tma[1,2], labels=bt.fxg$eti[i],col=bt.fxg$tcol[i],
						cex=.6, pos=bt.fxg$pos[i])
	}	}	}
	
	#	Funcion para Menu: Guardar	
	#
	GuardaArchivo <- function() 
	{
		FileName <- tclvalue(tkgetSaveFile(filetypes = "
				{{Imagen png} {.png}} 
				{{Imagen jpeg} {.jpg .jpeg}}
				{{Imagen svg} {.svg}} 
				{{Imagen wmf} {.wmf}} 
				{{Imagen pdf} {.pdf}}
				{{Imagen eps} {.eps .ps}}
				{{All files} *}
				"	))
		# por defecto, png
		if (!grepl("\\.",FileName)) 
			FileName <- paste(FileName, "png", sep = ".")
		nn <-  tolower(unlist(strsplit(FileName,"\\."))[[2]])
		
		if (nn=="pdf") pdf(FileName, width = 7, height = 7)
		else if (nn=="eps" | nn=="ps")
			postscript(file = FileName, width = 7, height = 7, horizontal = FALSE,
			onefile = FALSE, paper = "special", 
			family = "URWHelvetica",fonts=c("sans","serif"))
		else if (nn=="jpg" | nn=="jpeg")
			jpeg(FileName, width = 7, height = 7, units = "in", 
				restoreConsole = FALSE, res = 96, quality = 100)
		else if (nn=="svg")
			svg(FileName, width = 7, height = 7)
		else if (nn=="wmf")
			{plotBiplot1 (screen = FALSE)
			savePlot(FileName, type="wmf")	}
		else if (nn=="png")
			png(FileName, width = 7, height = 7, units = "in", 
					restoreConsole = FALSE, res = 96)
		else print(bt.lit[111,])		# ERROR en tipo de archivo

		plotBiplot1 (screen = FALSE)
        dev.off()
	}
	
	#	Menu para el grafico
	#
	fMenu <- function()
	{
		topmenu <- tk2menu(bt.ttp)
		tkconfigure(bt.ttp,menu=topmenu)
		filemenu <- tk2menu(topmenu,tearoff=FALSE)
		viewmenu <- tk2menu(topmenu,tearoff=FALSE)
		tkadd(filemenu,"command",label=bt.lit[112,],	# Copiar imagen
				command=function() tkrreplot(bt.img))
		tkadd(filemenu,"command",label=bt.lit[113,],command=function()	# Guardar Imagen
				GuardaArchivo() )
		tkadd(filemenu,"separator")
		tkadd(filemenu,"command",label=bt.lit[114,],command=function() tkdestroy(bt.ttp))	# Salir
		tkadd(topmenu,"cascade",label=bt.lit[115,],menu=filemenu)	# Archivo
    	# Control de la ventana, variable ittp
	    tkbind(bt.ttp,"<Destroy>", function() tclvalue(ittp)<-0) # ponerlo tambien en quit
	}
	
	#	Zoom para grafico
	#
	fZoom <- function()
	{
		bt.ex1 <- tk2entry(bt.ttp2,textvariable=ex1,width="5")
		bt.ex2 <- tk2entry(bt.ttp2,textvariable=ex2,width="5")
		bt.ey1 <- tk2entry(bt.ttp2,textvariable=ey1,width="5")
		bt.ey2 <- tk2entry(bt.ttp2,textvariable=ey2,width="5")
		
		fxy <- function(){				# boton cambio limites
			bt.limx1 <<- as.numeric(c(tclvalue(ex1),tclvalue(ex2)))
			bt.limy1 <<- as.numeric(c(tclvalue(ey1),tclvalue(ey2)))
			if (tclvalue(igx)==1) {
				if(iyax==0){bt.limx1 <<- sort(-1*bt.limx1)	# reajusta eje x
					iyax <<- 1
					tclvalue(ex1) <-min(bt.limx1)
					tclvalue(ex2) <-max(bt.limx1)}
				}
			else if(iyax==1){bt.limx1 <<- sort(-1*bt.limx1)	# reajusta eje x
					iyax <<- 0
					tclvalue(ex1) <-min(bt.limx1)
					tclvalue(ex2) <-max(bt.limx1)
					}
			if (tclvalue(igy)==1) {
				if(iyay==0) {bt.limy1 <<- sort(-1*bt.limy1)	# reajusta eje y
					iyay <<- 1
					tclvalue(ey1) <-min(bt.limy1)
					tclvalue(ey2) <-max(bt.limy1)}}
			else if(iyay==1) {bt.limy1 <<- sort(-1*bt.limy1)	# reajusta eje y
					iyay <<- 0
					tclvalue(ey1) <-min(bt.limy1)
					tclvalue(ey2) <-max(bt.limy1)}
		tkrreplot(bt.img)
		}
		bxy <- tk2button(bt.ttp2,text=bt.lit[116,],command=fxy)		# Refrescar
		bt.ttp3 <- tk2frame(bt.ttp2,relief="raised", borderwidth=2,padding="2")
		a<-tk2checkbutton(bt.ttp3,variable=igx,tip=bt.lit[117,])	# invierte x
		b<-tk2label(bt.ttp3, text="-x   -y ")
		c<-tk2checkbutton(bt.ttp3,variable=igy,tip=bt.lit[118,])	# invierte y
		tkpack(a,b,c,side="left")
		if(tclvalue(tb)!="3") {				# Si HJ y GH
			bt.ttp4 <- tk2scale(bt.ttp2,tip=bt.lit[119,],from=1,to=10,	# Escala x
					variable=vesc)
			bt.ttp5 <- tk2entry(bt.ttp2, textvariable=vesc, width=3) }
		else {	bt.ttp4 <- NULL
				bt.ttp5 <- NULL
			}
        tkpack(tk2label(bt.ttp2, text=bt.lit[120,]),bt.ex1,bt.ex2,	# Limites... x
				tk2label(bt.ttp2, text=" y:"),bt.ey1,bt.ey2,
				tk2label(bt.ttp2, text=" "),bxy,
				tk2separator(bt.ttp2,orientation="vertical"),bt.ttp3,
				bt.ttp4,bt.ttp5, side="left")
	}
	
	plotBiplot <- function(screen = TRUE) 
		{   
   		if (tclvalue(i3v)=="1" & tclvalue(ibg)!="1")  {		
			bt.limx <<- range(0,range(bt.res.ty[,dim1,]),range(bt.res.tx[,dim1,]))
			bt.limy <<- range(0,range(bt.res.ty[,dim2,]),range(bt.res.tx[,dim2,])) 
			}	else {	
			bt.limx <<- range(0,rbind(bt.a,bt.b)[,dim1])
			bt.limy <<- range(0,rbind(bt.a,bt.b)[,dim2])
			}
    	bt.limx1 <<- bt.limx
    	bt.limy1 <<- bt.limy
    	tclvalue(ex1) <-min(bt.limx)	# preparo para zoom
		tclvalue(ex2) <-max(bt.limx)
		tclvalue(ey1) <-min(bt.limy)
		tclvalue(ey2) <-max(bt.limy)
		tclvalue(igx) <- 0				# Comenzamos siempre sin rotar
		tclvalue(igy) <- 0
		tclvalue(vesc) <- 1				# reescalado
		iyax <<- iyay <<- 0				# ya rotado
        if (tclvalue(ittp)==0) 			# dibuja grafico
		{
			bt.ttp <<- tktoplevel()                
			tkwm.title(bt.ttp,bt.lit[2,]) 	# "Biplot Dinamico"
			bt.ttp1 <<- tk2frame(bt.ttp)	# para figura
			bt.ttp2 <<- tk2frame(bt.ttp,relief="raised",borderwidth=2,padding="2")	# pie para coordenadas
			tkpack(bt.ttp1,bt.ttp2, side="top",fill="x")

			fMenu()
			fZoom()						# dibuja campos para zoom
			bt.img <<- tkrplot(bt.ttp1, fun = plotBiplot1,
				hscale = as.numeric(tclvalue(wout1)),
				vscale = as.numeric(tclvalue(wout2)))
        	tclvalue(ittp)<-1
        	tkpack(bt.img, expand = "TRUE", fill = "both")
        } else {
        	tkrreplot(bt.img) 
        	}
	}	
	#
	#	Fin funciones	==========================
	#
	if (!is.ttk()) stop(bt.lit[1,])			# "Tcl/Tk >= 8.5 is required"

	#======		Ventana principal	==========
	#
	tt <- tktoplevel()
	Sys.sleep(0.1)							# para bug PR#15150
	# tkwm.title(tt,bt.lit[2,]) 				# "Biplot Dinamico"
	tktitle(tt)<-(bt.lit[2,]) 				# "Biplot Dinamico"
	fontTextLabel <- tkfont.create(family="times",size=12)
	
	nb <- tk2notebook(tt,tabs=c(bt.lit[3,],bt.lit[4,],bt.lit[5,],bt.lit[6,]))
	tkpack(nb,fill="both")
	
	#	Panel 1
	tb1 <- tk2notetab(nb, bt.lit[3,])
	fr.d2 <- tk2labelframe(tb1,text=bt.lit[7,])	# para variables
	tkconfigure(fr.d2, relief="groove")
	panel1()

	#	Panel 2
	tb2 <- tk2notetab(nb, bt.lit[4,])
	fr.f1 <-tk2labelframe(tb2,text=bt.lit[8,])	# Formato de datos:
	fr.f11 <- tk2frame(fr.f1)
	fr.f12 <- tk2frame(fr.f1)
	tkpack(fr.f11, fr.f12, side="left")
	fr.f2 <-tk2labelframe(tb2,text=bt.lit[9,])	# Titulos:
	nb2 <- tk2notebook(fr.f11,tabs=c(bt.lit[10,],bt.lit[5,]))	# individuos
	tb21 <- tk2notetab(nb2, bt.lit[10,])		# Individuos
	tb22 <- tk2notetab(nb2, bt.lit[5,])			# Variables
	panel2()
	
	#	Panel 3
	tb3 <- tk2notetab(nb, bt.lit[5,])
	fr.a0 <-tk2frame(tb3)
	fr.a1 <-tk2frame(fr.a0)
	fr.a2 <-tk2frame(fr.a0)
	fr.a3 <-tk2frame(fr.a0)
	panel3()
	
	#	panel 4
	tb4 <- tk2notetab(nb, bt.lit[6,])
	fr.s1 <-tk2labelframe(tb4,text=bt.lit[11,])	# Estandarizacion
	fr.s2 <-tk2labelframe(tb4,text=bt.lit[12,])	# Analisis Biplot
	fr.s3 <-tk2labelframe(tb4,text=bt.lit[13,])	# Ejes
	fr.s4 <-tk2labelframe(tb4,text=bt.lit[14,],relief="raised")	# Trayectorias
	fr.s5 <-tk2labelframe(tb4,text=bt.lit[15,],relief="sunken")	# Opciones de grafico
		#	Seleccion de ejes y plano
	fr25 <- tk2frame(fr.s3,relief="groove", borderwidth=2,padding="2")
	tb4.cb0 <- tk2combobox(fr25, values= 2,  width="2", textvariable=neje)
	tb4.cb1 <- tk2combobox(fr25, values= 1, width="2", textvariable=di1)
	tb4.cb2 <- tk2combobox(fr25, values= 2,  width="2", textvariable=di2)
	panel4()

	#	abajo
	down.frm <- tk2frame(tt,padding="2")	
	tkpack(down.frm, side="bottom", fill="x")	
	#
	
	#
	##	 Botones de la ventana general:
	frame3 <- tk2frame(down.frm,relief="sunken",borderwidth=2,padding="2")
	mens.leer <- bt.lit[93,]			# Leer datos
	la <- tk2label(frame3,text=mens.leer, foreground="red")
	tkpack(la)							# para mostrar a pie de pagina
	tkpack(frame3, side="left") 

	frame31 <- tk2frame(down.frm, padding="2")
    run.but <- tk2button(frame31,text=bt.lit[16,],command=run.biplot,	# Run Biplot
				state="disable")
    res.but <- tk2button(frame31,text=bt.lit[17,],command=ShowRes,	# Resultados,
				state="disable")
	tkpack(frame31, run.but, res.but, side="left")
	
	frame32 <- tk2frame(down.frm,padding="2")
    q.but <- tk2button(frame32,text=bt.lit[18,], 	# Quit
				command=function() {tkdestroy(tt)	
				rm(list=ls(.GlobalEnv,pattern="bt"),pos=1) } )
	tk2tip(q.but, bt.lit[19,])			# Cierra y borra
    tkpack(frame32, q.but)
	tkpack(frame32, frame31, side="right") 
	
	# tkfocus(tt) 
	a<-1
}
