#################################################################################################
##
## ripa: R Image Processing and Analysis
##
##
## Copyright (c) 2014 Talita Perciano
## For complete license terms see file LICENSE
##################################################################################################

##################################################################################################################
############################### Initialization of environment ripaEnv ############################################
##################################################################################################################

ripaEnv <- new.env()

##################################################################################################################

##################################################################################################################
####################################### Aviris functions #########################################################
##################################################################################################################

		       ##############
		       #            #
         		 #            #
	     #########    ####    ############
	     #       #    #   #   #          #
	################  ####    #   2006   #
	#              #  #  #    #          #
	#    AVIRIS    #  #   #   ############
	#              # .#   #   #
	################          #
	             #            #
	             ##############

########################################################################
###                                                                  ###
###       Marcelo G. Almiron		        Adrian E. Muract         ###
###                                                                  ###
###   <almiron.marcelo@gmail.com>     <amuract@dc.exa.unrc.edu.ar>   ###
###                                                                  ###
########################################################################






	#########################
######### Class aviris_band ############################################################
	#########################



setClass("aviris_band",representation(scene="character",band="numeric",type="character",
	numberOfLines="numeric",samples="numeric",data="matrix",min="numeric",
	max="numeric",mean="numeric",sd="numeric"))



setMethod("initialize",
	"aviris_band",
	function(.Object,scene=character(),band=numeric(),type=character(),
	numberOfLines="numeric()",samples=numeric(),data=matrix()){
		.Object@scene=scene
		.Object@band=band
		.Object@type=type
		.Object@numberOfLines=numberOfLines
		.Object@samples=samples
		.Object@data=data
		.Object@min=min(data)
		.Object@max=max(data)
		.Object@mean=mean(data)
		.Object@sd=sd(as.vector(data))
		.Object
	})



print_information.aviris_band <- function(Object){
	cat("  Scene: ",Object@scene,"\n")
	cat("  Band : ",Object@band,"\n")
	cat("  Number of lines: ",Object@numberOfLines,"\n")
	cat("  Number of samples: ",Object@samples,"\n")		
	cat("  Minimum: ",Object@min,"\n")
	cat("  Maximum: ",Object@max,"\n")
	cat("  Mean: ",Object@mean,"\n")
	cat("  Standard desviation:",Object@sd,"\n")
}



# Interfaz para la Carga de una banda #
#######################################

lband <- function(scene,b){
	path <- paste(scene@path,scene@name,sep="/")
	Object <- new("aviris_band",scene@name,b,scene@type,scene@numberOfLines,
		scene@samples,loadBand(path,X=b,C=scene@samples,F=scene@numberOfLines,
		B=scene@bands))
	Object@data <- clineal(Object@data,0,1)
	return(Object)
}

lbandsample <- function(scene,b){
	path <- paste(scene@path,scene@name,sep="/")
	Object <- new("aviris_band",scene@name,b,scene@type,scene@numberOfLines,
		scene@samples,loadBandSample(path,X=b,C=scene@samples,F=30,
		B=scene@bands))
	Object@data <- clineal(Object@data,0,1)
	return(Object)
}


# Interfaz para Escribir en una banda #
#######################################

wband <- function(scene,band){
	path <- paste(scene@path,scene@name,sep="/")
	writeBand(path,band@data,band@band,band@samples,band@numberOfLines,scene@bands)
}




# Algor�mo de Carga de una banda #
###################################

loadBand <- function(I,X=5,C=614,F=512,B=224){
    f <- 0
    conexionParaLeer <- file(I,open="r+b")   
    Z <- matrix(nrow=F,ncol=C)
    while (f!=F){
        c <- 0
        while (c!=C){
            seek(conexionParaLeer,where=(((B*(f*C+c))+(X-1))*2),rw="read")
	    Z[f+1,c+1] <- readBin(conexionParaLeer,"integer",size=2,endian="swap")
            c <- c+1
        }
        f <- f+1
    }
    close(conexionParaLeer,rw="read")
    Z
}

loadBandSample <- function(I,X=5,C=614,F=30,B=224){
    f <- 0
    conexionParaLeer <- file(I,open="r+b")   
    Z <- matrix(nrow=F,ncol=30)
    while (f!=F){
        c <- 0
        while (c!=30){
            seek(conexionParaLeer,where=(((B*(f*C+c))+(X-1))*2),rw="read")
	    Z[f+1,c+1] <- readBin(conexionParaLeer,"integer",size=2,endian="swap")
            c <- c+1
        }
	f <- f+1
    }
    close(conexionParaLeer,rw="read")
    Z
}

# Algor�mo de Escritura en una banda #
#######################################

writeBand <- function(I,Z,X=NA,C=614,F=512,B=224){
    f <- 0
    conexionParaEscribir <- file(I,open="r+b")
    while (f!=F){
        c <- 0
        while (c!=C){
            seek(conexionParaEscribir,where=(((B*(f*C+c))+(X-1))*2),rw="write") 
            writeBin(as.integer(Z[f+1,c+1]),conexionParaEscribir,size=2,endian="swap")
            c <- c+1
        }
        f <- f+1
    }
    close(conexionParaEscribir,rw="write")
}




# Interfaz de Contrastes #
##########################

contrast <- function(band,type=c("gauss","lineal"),...){
	if(type=="gauss"){
		c <- new("aviris_band",band@scene,band@band,band@type,
			band@numberOfLines,band@samples,cgauss(band@data))	
	}else{
		c <- new("aviris_band",band@scene,band@band,band@type,
			band@numberOfLines,band@samples,clineal(band@data,...))
	}
}




# Algor�mo de Expansion Lineal de Contraste #
##############################################

clineal <- function(Z,A,B){
	if (max(Z)==0) return(Z)
	minZ <- min(Z)
	rango <- (B - A)
	rangoZ <- (max(Z)-minZ)
	pendiente <- (rango/rangoZ)
	res <- pendiente*Z+(A-(pendiente*minZ))
	
	if (is.nan(res[1])){
		for (i in 1:length(res)) res[i] <- 0
	}
	return(res)
}




# Algor�mo de Expansion Gausseano de Contraste #
#################################################

cgauss <- function(Z){
	s <- sd(as.vector(Z))
	x <- mean(Z)
	dnorm(Z,mean=x,sd=s)
}




# Interfaz de Visualizacion #
#############################

plot_band.aviris_band <- function (R=NULL,G=NULL,B=NULL,type=NULL,x0=1,y0=1,...){
	if (is.null(R) && is.null(G) && is.null(B)){ 
      	stop("Missing objects")
	}else if (class(R)[1]=="aviris_band" && is.null(G) && is.null(B)){ 
      	type <- "grey"
	}else if (class(R)[1]=="aviris_band" && class(G)[1]=="aviris_band" && class(B)[1]=="aviris_band"){
      	type <- "rgb"
	}else{ 
		stop("Incorrect parameters")
	}
	if (type == "rgb"){
		RGB(R,G,B,x0,y0,...)
	}else{
		Grey(R,x0,y0,...)	
	}
}




# Interfaz de Visualizacion en escala de Grises #
#################################################

Grey <- function(band,x0,y0,...){
	w <- 6.39
	h <- 5.33
	sw <- band@samples/96
	sh <- band@numberOfLines/96
	if((w-sw)<(h-sh)){
		c <- w/sw
	}else{
		c <- h/sh
	}
	dev.new(width=sw*c,height=sh*c)
	opar <- par(no.readonly=TRUE)
	opar$usr <- c(x0,x0+band@samples-1,y0+band@numberOfLines-1,y0)
	image(1:band@samples,1:band@numberOfLines,t(band@data[band@numberOfLines:1,])     
		,col=gray(0:255/256),xlab="Samples",ylab="Lines",axes=FALSE,...)
	par(cex=0.7)
	title(paste(paste("Scene",band@scene,sep=": "),paste("Band",band@band,sep=": "),sep=" - "))
	par(opar)
}




# Interfaz de Visualizacion en RGB #
####################################

RGB <- function(red,green,blue,x0,y0,...){
	rc <- clineal(red@data,0,1)
	gc <- clineal(green@data,0,1)
	bc <- clineal(blue@data,0,1)
	colvec <- rgb(round(rc,7),round(gc,7),round(bc,7))
    	colors <- unique(colvec)
    	colmat <- array(match(colvec,colors), dim = dim(red@data)[1:2])
	dev.new(width=6.39,height=5.33)
	opar <- par(no.readonly=TRUE)
	opar$usr <- c(x0,x0+dim(colmat)[2],y0+dim(colmat)[1],y0)
	image(1:dim(colmat)[2],1:dim(colmat)[1],t(colmat[nrow(colmat):1,])     
		,col=colors,xlab="Samples",ylab="Lines",axes=FALSE,...)
	r <- paste("R: Band",red@band,sep=" ")
	g <- paste("G: Band",green@band,sep=" ")
	b <- paste("B: Band",blue@band,sep=" ")
	subt <- paste(r,g,sep=" , ")
	t <- paste(subt,b,sep=" , ")
	par(cex=0.7)
	title(paste(paste("Scene",red@scene,sep=" "),t,sep=" - "))
	par(opar)
}



# ROIs #
########

zoom <- function (R=NULL,G=NULL,B=NULL){
	if (is.null(R) && is.null(G) && is.null(B)){ 
      	stop("Missing objects")
	}else if (class(R)[1]=="aviris_band" && is.null(G) && is.null(B)){ 
      	type <- "grey"
	}else if (class(R)[1]=="aviris_band" && class(G)[1]=="aviris_band" && class(B)[1]=="aviris_band"){
      	type <- "rgb"
	}else{ 
		stop("Incorrect parameters")
	}
	if (type == "rgb"){
		zoomRGB(R,G,B)
	}else{
		zoomGrey(R)	
	}
}



zoomGrey <- function(band){
	pos <- locator(2)
	if(pos$x[1]>pos$x[2]){
		aux <- pos$x[1]
		pos$x[1] <- pos$x[2]
		pos$x[2] <- aux  
	}
	if(pos$y[1]>pos$y[2]){
		auy <- pos$y[1]
		pos$y[1] <- pos$y[2]
		pos$y[2] <- auy  
	}
	pos$x <- round(pos$x)
	pos$y <- round(pos$y)
	rect(pos$x[1],pos$y[1],pos$x[2],pos$y[2],col="red",density=0.1)
	newBand <- band
	newBand@data <- band@data[pos$y[1]:pos$y[2],pos$x[1]:pos$x[2]]
	newBand@samples <- ncol(newBand@data)
	newBand@numberOfLines <- nrow(newBand@data)
	plot(newBand,x0=pos$x[1],y0=pos$y[1])
}

zoomRGB <- function(Red, Green, Blue){
	pos <- locator(2)
	if(pos$x[1]>pos$x[2]){
		aux <- pos$x[1]
		pos$x[1] <- pos$x[2]
		pos$x[2] <- aux  
	}
	if(pos$y[1]>pos$y[2]){
		auy <- pos$y[1]
		pos$y[1] <- pos$y[2]
		pos$y[2] <- auy  
	}
	pos$x <- round(pos$x)
	pos$y <- round(pos$y)
	rect(pos$x[1],pos$y[1],pos$x[2],pos$y[2],col="red",density=0.1)
	newBandR <- Red
	newBandG <- Green
	newBandB <- Blue
	newBandR@data <- Red@data[pos$y[1]:pos$y[2],pos$x[1]:pos$x[2]]
	newBandG@data <- Green@data[pos$y[1]:pos$y[2],pos$x[1]:pos$x[2]]
	newBandB@data <- Blue@data[pos$y[1]:pos$y[2],pos$x[1]:pos$x[2]]
	newBandR@samples <- ncol(newBandR@data)
	newBandG@samples <- ncol(newBandG@data)
	newBandB@samples <- ncol(newBandB@data)
	newBandR@numberOfLines <- nrow(newBandR@data)
	newBandG@numberOfLines <- nrow(newBandG@data)
	newBandB@numberOfLines <- nrow(newBandB@data)
	plot(newBandR, newBandG, newBandB,x0=pos$x[1],y0=pos$y[1])
}





	##########################
######### Class aviris_scene ############################################################
	##########################



setClass("aviris_scene",representation(name="character",numberOfLines="numeric",
	samples="numeric",bands="numeric",imageName="character",type="character",path="character"))



setMethod("initialize",
	"aviris_scene",
	function(.Object,name=character(),numberOfLines=numeric(),samples=numeric(),
		bands=numeric(),imageName=character(),type=character(),path=character()){
		.Object@name=name
		.Object@numberOfLines=numberOfLines
		.Object@samples=samples
		.Object@bands=bands
		.Object@imageName=imageName
		.Object@type=type
		.Object@path=path
		.Object
	})


	

print_information.aviris_scene <- function(Object){
	cat("  Name: ",Object@name,"\n")
	cat("  Number of lines: ",Object@numberOfLines,"\n")
	cat("  Number of samples :",Object@samples,"\n")		
	cat("  Number of bands: ",Object@bands,"\n")
}




# Interfaz para cargar una escena de la imagen #
################################################

lscene <- function(image,n){
	if (n == image@numberOfScenes){
		numberOfLines <- image@linesInLastScene
	}else if (n < image@numberOfScenes && n > 0){
		numberOfLines <- 512
	}else{
		stop("La escena solicitada no existe")
	}
	sc <- paste("_sc0",as.character(n),sep="")
	if (image@type=="reflectance"){
		name <- sub(".a",paste(sc,".a.rfl",sep=""), image@name)
	}else{
		name <- sub(".c",paste(sc,".c.img",sep=""), image@name)
	}
	Object <- new("aviris_scene",name,numberOfLines,614,224,image@name,image@type,image@path)
}



	##########################
######### Class aviris_image ############################################################
	##########################



setClass("aviris_image",representation(name="character",numberOfScenes="numeric",
	linesInLastScene="numeric",type="character",path="character"))



setMethod("initialize",
	"aviris_image",
	function(.Object,name=character(),numberOfScenes=numeric(),linesInLastScene=numeric(),
	type=character(),path=character()){
		.Object@name=name
		.Object@numberOfScenes=numberOfScenes
		.Object@linesInLastScene=linesInLastScene
		.Object@type=type
		.Object@path=path
		.Object
	})




print_information.aviris_image <- function(Object){
	cat("  Name: ",Object@name,"\n")
	cat("  Number of scenes: ",Object@numberOfScenes,"\n")
	cat("  Lines in last scene: ",Object@linesInLastScene,"\n")		
	cat("  Image type: ",Object@type,"\n")
}




# Interfaz para cargar el archivo cabecera de una imagen #
##########################################################

limage <- function(H,type){
	if (type!="reflectance" & type!="radiance"){ 
		stop("Debe especificar el tipo de imagen: reflectance o radiance")
	}
	con <- file(paste(H,".log",sep=""),"r+b")
	txt <- readLines(con)
	close(con,rw="read")
	lpath <- unlist(strsplit(H, "\\/"))
	nameImg <- lpath[length(lpath)] 
	path <- sub(lpath[length(lpath)],"",H)
	if (length(lpath)==1) path <- getwd()
	e <- txt[5]
	u <- txt[6]
	e <- as.numeric(sub("number of scenes","",e))#,fixed=TRUE))  
	u <- as.numeric(sub("lines in last scene","",u))#,fixed=TRUE))
	avirisImage <- new("aviris_image",nameImg,as.numeric(e),as.numeric(u),type,path)     	
	}


	###############################
######### Class aviris_training ############################################################
	###############################

	
setClass("aviris_training",representation(category="character",color="character"
			,scene="aviris_scene",bands="list",posX="numeric",posY="numeric"))


setMethod("initialize",
	"aviris_training",
	function(.Object, category=character(),color=character(),scene=aviris_scene,
			band=list(),posX=vector(),posY=vector()){
		.Object@category <- category
		.Object@color <- color
		.Object@scene <- scene
#### tratamiento por si las bandas pertenecen a diferentes escenas ##################
		i <- 1									#
		bool<-TRUE								#
		while (i!=(length(band))){						#
			bool <- all.equal(band[[i]]@scene, band[[i+1]]@scene)		#
####		cat(band[[i]]@scene," = ", band[[i+1]]@scene," es ", bool,"\n")		#
			i <- i+1							#
		}									#
		if (!bool) {								#
			cat("OJO QUE LAS BANDAS PERTENECEN A DIFERENTES ESCENAS \n")	#
		}else {									#
			cat("TODO BIEN LOCO \n")					#
		}									#
###	la decision que tome es solo advertirle al usuario que las bandas pertenecen a 	#
###	distintas escenas. 								#
#####################################################################################
		.Object@bands <- band
		.Object
	}
)

print_information.aviris_training <- function(Object){
	cat("  Category: ",Object@category,"\n")
	cat("  Color: ",Object@color,"\n")
	cat("  Scene: ",Object@scene@name,"\n")
	i <- 0
	cat("BANDAS INVOLUCRADAS \n")
	while (i!=(length(Object@bands))){
		i <- i+1
		cat("  Band : ",Object@bands[[i]]@band,"\n")
	}
}



#################################################################################
# falta hacer analisis por casos porque falta cubrir algunos casos de parametros 
# incorrectos. por ejemplo "takeSamples(clase,5, 56)->clase"
#################################################################################
takeSamples <- function(t, n=NULL, Sample=NULL, Line=NULL){
	if (is.null(Sample) && is.null(Line) && is.null(n)){ 
		p <- locator(1,type="p",pch=20, col=t@color)
		t@posY <- c(t@posY, round(p$y))
		t@posX <- c(t@posX, round(p$x))
	}else if (is.null(Sample) && is.null(Line) && !is.null(n)){ 
		p <- locator(n,type="p",pch=20, col=t@color)
		t@posY <- c(t@posY, round(p$y))
		t@posX <- c(t@posX, round(p$x))
	}else if ((is.null(Sample) && !is.null(Line)) || ((is.null(n) && !is.null(Sample) && is.null(Line))) ){ 
      	stop("Incorrect parameters")
	}else{
		t@posY <- c(t@posY, Line)
		t@posX <- c(t@posX, Sample)
		points(Sample, Line, type="p",pch=20, col="green")
	} 
	cat("  X value: ",t@posX,"\n")
	cat("  Y value: ",t@posY,"\n")
	t
}

#############################################################################################


Zprofile <- function(scene, X=NULL, Y=NULL){
	if (is.null(X) && is.null(Y)){ 
		p <- locator(1)
		Y <- round(p$y)
		X <- round(p$x)
	}else if ((is.null(X) && !is.null(Y)) || ((!is.null(X) && is.null(Y))) ){ 
      	stop("Incorrect parameters")
	} 
	I <- paste(scene@path,scene@name,sep="/")
	conexionParaLeer <- file(as.character(I),open="r+b")   
	seek(conexionParaLeer,where=(224*((Y-1)*614+(X-1))*2),rw="read")
	Z <- vector()
	i <- 1
	while (i!=225) {
		Z <- c(Z,readBin(conexionParaLeer,"integer",size=2,endian="swap"))
		i <- i+1
	}
	close(conexionParaLeer,rw="read")
	dev.new()
	plot(Z, type="l", col="red")	
    Z
}



##################################################################################################################

##################################################################################################################
####################################### Functions ################################################################
##################################################################################################################

# Function to read LAN images
read.lan <- function(arquivo){
	assign('LANfile',arquivo,envir=ripaEnv)
	res <- .C("getImgData",file=as.character(arquivo),numbands=as.integer(0),row=as.integer(0),col=as.integer(0),PACKAGE="ripa")
	row <- res$row
	col <- res$col
	nbands <- res$numbands
	imgMatrix <- array(0,dim=c(row,col,nbands))
	n <- row*col
	res <- .C("readLAN",file = as.character(arquivo),out = as.integer(rep(0,n*nbands)),PACKAGE="ripa")
	mat <- matrix(res$out,ncol=nbands,byrow=T)
	for (i in 1:nbands){
		band <- matrix(mat[,i]/255,nrow=row,byrow=T)
		imgMatrix[,,i] <- band
	}
	return(imgMatrix)
}
#

# Function ro write LAN images
write.lan <- function(arquivo,img){
	mat <- matrix(0,ncol=dim(img)[3],nrow=(nrow(img)*ncol(img)))
	
	for (i in 1:dim(img)[3]){
		mat[,i] = t(img[,,i]*255)
	}
	
	res <- .C("writeLAN",file1 = as.character(get('LANfile',envir=ripaEnv)),file2 = as.character(arquivo), as.integer(as.vector(t(mat))),PACKAGE="ripa")
}
#

# Function to apply the median filter to an image
#medianImg <- function(img,mask){
#	if (is.na(dim(img)[3])) img = array(img,dim=c(nrow(img),ncol(img),1))
#	else img = array(img,dim=dim(img))
#	n <- nrow(img)*ncol(img)
#	nrow <- as.integer(nrow(img))
#	ncol <- as.integer(ncol(img))
#	mat <- img
#	for (i in 1:dim(img)[3]){
#		image <- as.double(as.vector(t(matrix(img[,,i],nrow=nrow(img)))))
#		out <- .C("median",image,nrow, ncol,mask,outImg=as.double(rep(0.0,n)),PACKAGE="ripa")
#		mat[,,i] <- matrix(out$outImg,ncol=ncol(img),byrow=T)
#	}
#	return(mat)
#}
#

# Function to apply the median filter to an image
medianImg <- function(img,mask){
	if (is.na(dim(img)[3])) img = array(img,dim=c(nrow(img),ncol(img),1))
	else img = array(img,dim=dim(img))
	n <- nrow(img)*ncol(img)
	nrow <- as.integer(nrow(img))
	ncol <- as.integer(ncol(img))
	mat <- img

	if (tclvalue(get('useParallel',envir=ripaEnv))==1) {	
		matrix = foreach(i=1:dim(img)[3],.combine=cbind,.packages=c('ripa')) %dopar% {
			image <- as.double(as.vector(t(matrix(img[,,i],nrow=nrow(img)))))
			out <- .C("median",image,nrow, ncol,mask,outImg=as.double(rep(0.0,n)),PACKAGE="ripa")
			mat[,,i] <- matrix(out$outImg,ncol=ncol(img),byrow=T)
		}
		mat <- array(matrix,dim=dim(img))
	}
	else {

		for (i in 1:dim(img)[3]){
			image <- as.double(as.vector(t(matrix(img[,,i],nrow=nrow(img)))))
			out <- .C("median",image,nrow, ncol,mask,outImg=as.double(rep(0.0,n)),PACKAGE="ripa")
			mat[,,i] <- matrix(out$outImg,ncol=ncol(img),byrow=T)
		}
	}
	return(mat)
}
#

# Linear contrast stretch function
stretchImg <- function(img){
	if (is.na(dim(img)[3])) img <- array(img,dim=c(nrow(img),ncol(img),1))
	else img <- array(img,dim=dim(img))
	
	res_final <- array(0,dim=dim(img))
	n <- nrow(img)*ncol(img)
	a <- as.double(0)
	b <- as.double(1)

	if (tclvalue(get('useParallel',envir=ripaEnv))==1) {	
		matrix = foreach(i=1:dim(img)[3],.combine=cbind,.packages=c('ripa')) %dopar% {
			x <- quantile(img[,,i],seq(0,1,by=0.05))
			c <- as.double(x[[2]])
			d <- as.double(x[[20]])
			res <- .C("stretch",img[,,i],a,b,c,d,as.integer(nrow(img)),as.integer(ncol(img)),out=as.double(rep(0,n)),PACKAGE="ripa")
			mat <- matrix(res$out,ncol=ncol(img))
			res_final[,,i] <- imagematrix(mat)
		}
		res_final <- array(matrix,dim=dim(img))
	}
	else {

		for (i in 1:dim(img)[3]){
			x <- quantile(img[,,i],seq(0,1,by=0.05))
			c <- as.double(x[[2]])
			d <- as.double(x[[20]])
			res <- .C("stretch",img[,,i],a,b,c,d,as.integer(nrow(img)),as.integer(ncol(img)),out=as.double(rep(0,n)),PACKAGE="ripa")
			mat <- matrix(res$out,ncol=ncol(img))
			res_final[,,i] <- imagematrix(mat)
		}
	}
	return(res_final)
}

contBriImg <- function(img,cont,bri){

	if (is.na(dim(img)[3])) img <- array(img,dim=c(nrow(img),ncol(img),1))
	else img <- array(img,dim=dim(img))
	
	if (is.na(cont)) cont <- 1
	if (is.na(bri)) bri <- 0
	
	img <- cont*img+bri
	nrow <- as.integer(nrow(img))
	ncol <- as.integer(ncol(img))
		
	index <- which(img<0)
	img[index] <- 0
	index <- which(img>1)
	img[index] <- 1
	
	return(img)
}

# Function to apply the new brightness and contrast
#contBriImg <- function(img,cont,bri){
#
#	if (is.na(dim(img)[3])) img <- array(img,dim=c(nrow(img),ncol(img),1))
#	else img <- array(img,dim=dim(img))
#	
#	if (is.na(cont)) cont <- 1
#	if (is.na(bri)) bri <- 0
#	res_final<-img
#	for (i in 1:dim(img)[3]){
#		res <-as.matrix(img[,,i])
#		
#		res <- cont*res+bri
#		
#		nrow <- as.integer(nrow(res))
#		ncol <- as.integer(ncol(res))
#		
#		#out <- .C("normalize",image<-as.vector(as.double(res)),nrow,ncol,as.double(bri),PACKAGE="RIPA")
#		
#		#res <- matrix(image,nrow=nrow)
#
#		index <- which(res<0)
#		res[index] <- 0
#		index <- which(res>1)
#		res[index] <- 1
#		
#		res_final[,,i]<-res
#	}
#	return(res_final)
#}

#contBriImg <- function(img,cont,bri){
#
#	if (is.na(dim(img)[3])) img <- array(img,dim=c(nrow(img),ncol(img),1))
#	else img <- array(img,dim=dim(img))
#	
#	if (is.na(cont)) cont <- 1
#	if (is.na(bri)) bri <- 0
#	res_final<-img
#
#	if (tclvalue(get('useParallel',envir=ripaEnv))==1) {	
#		matrix = foreach(i=1:dim(img)[3],.combine=cbind) %dopar% {
#
#			res <- as.matrix(img[,,i])
#		
#			res <- cont*res+bri
#		
#			nrow <- as.integer(nrow(res))
#			ncol <- as.integer(ncol(res))
#		
#			#out <- .C("normalize",image<-as.vector(as.double(res)),nrow,ncol,as.double(bri),PACKAGE="RIPA")
#		
#			#res <- matrix(image,nrow=nrow)
#
#			index <- which(res<0)
#			res[index] <- 0
#			index <- which(res>1)
#			res[index] <- 1
#		
#			res_final[,,i]<-res
#		}
#		res_final <- array(matrix,dim=dim(img))
#	}
#	else {
#		
#		for (i in 1:dim(img)[3]){
#			res <-as.matrix(img[,,i])
#		
#			res <- cont*res+bri
#		
#			nrow <- as.integer(nrow(res))
#			ncol <- as.integer(ncol(res))
#		
#			#out <- .C("normalize",image<-as.vector(as.double(res)),nrow,ncol,as.double(bri),PACKAGE="RIPA")
#		
#			#res <- matrix(image,nrow=nrow)
#
#			index <- which(res<0)
#			res[index] <- 0
#			index <- which(res>1)
#			res[index] <- 1
#		
#			res_final[,,i]<-res
#		}
#	}
#	
#	
#	return(res_final)
#}

#Check active tab
checkTab <- function(){
	tn_local <- get('tn',envir=ripaEnv)
	aux <- tclvalue(tcl(tn_local,"select"))
	tab <- unlist(strsplit(aux,'[.]'))
	tab <- tab[length(tab)]
	return(tab)		
}

# Function to read AVIRIS images
read.aviris <- function(fileName,bandsIndexes,bands_local,use_parallel){
	strt<-Sys.time()
	comb <- function(x, ...) {
		lapply(seq_along(x),
		function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
	}	

	#tab <- checkTab()
 	
 	#if (tab=="1") bandsIndexes <- get('AVIRISbands1',envir=ripaEnv)
	#else bandsIndexes <- get('AVIRISbands2',envir=ripaEnv)

	img <- limage(as.character(fileName),"reflectance")
	imgtmpScene2 <- lscene(img,2)
	mat <- array(0,dim=c(512,614,length(bandsIndexes)))
	#assign('bands',vector(length=length(bandsIndexes),mode="list"),envir=ripaEnv)
	#bands_local <- get('bands',envir=ripaEnv)

	tt <- tktoplevel()
	tktitle(tt) <- "Read AVIRIS"
	tkplace(tklabel(tt,text="Reading selected bands..."),x=5,y=1)
	tkconfigure(tt,width=200,height=60)

	if (use_parallel==1) {
		matrix = foreach(i=1:length(bandsIndexes),.combine='comb',.packages=c("ripa")) %dopar% {
			list(lband(imgtmpScene2,bandsIndexes[i]),lband(imgtmpScene2,bandsIndexes[i])@data)

			#bands_local[[i]] <- lband(imgtmpScene2,bandsIndexes[i])
			#mat[,,i] <- clineal(bands_local[[i]]@data,0,1)
								
			#setTkProgressBar(pb, i, label=paste(round(i/length(bandsIndexes)*100, 0),"% done"))
		}
		mat <- array(unlist(matrix[[2]]),dim=dim(mat))
		mat <- clineal(mat,0,1)
		bands_local <- matrix[[1]]
	}
	else {
		#pb <- tkProgressBar(title = "Reading bands...", min = 0, max = length(bandsIndexes), width = 300)
		#setTkProgressBar(pb, 0, label="0% done")
		for (i in 1:length(bandsIndexes)){
			bands_local[[i]] <- lband(imgtmpScene2,bandsIndexes[i])
			#mat[,,i] <- clineal(bands_local[[i]]@data,0,1)
			mat[,,i] <- bands_local[[i]]@data
			#setTkProgressBar(pb, i, label=paste(round(i/length(bandsIndexes)*100, 0),"% done"))
		}
		#close(pb)
	}
	tkdestroy(tt)
	assign('bands',bands_local,envir=ripaEnv)
	print(Sys.time()-strt)
	return(mat)
}

# Function to build a modal dialog
modalDialog <- function(title,question,entryInit,entryWidth=20,returnValOnCancel="ID_CANCEL"){
	dlg <- tktoplevel()
	tkwm.deiconify(dlg)
	tkfocus(dlg)
	tkwm.title(dlg,title)
	tkwm.geometry(dlg,"+350+300")
	textEntryVarTcl <- tclVar(paste(entryInit))
	textEntryWidget <- tkentry(dlg,width=paste(entryWidth),textvariable=textEntryVarTcl)
	tkgrid(tklabel(dlg,text="       "))
	tkgrid(tklabel(dlg,text=question),textEntryWidget)
	tkgrid(tklabel(dlg,text="       "))
	ReturnVal <- returnValOnCancel
	onOK <- function(){
		ReturnVal <<- tclvalue(textEntryVarTcl)
		tkgrab.release(dlg)
		tkdestroy(dlg)
		}
	onCancel <- function(){
		ReturnVal <<- returnValOnCancel
		tkgrab.release(dlg)
		tkdestroy(dlg)
	}
	OK.but <-tkbutton(dlg,text="   OK   ",command=onOK)
	Cancel.but <-tkbutton(dlg,text=" Cancel ",command=onCancel)
	tkgrid(OK.but,Cancel.but)
	tkgrid(tklabel(dlg,text="    "))
	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg)})
	tkbind(textEntryWidget, "<Return>", onOK)
	tkwait.window(dlg)
	return(ReturnVal)
}

##################################################################################################################

##################################################################################################################
################################################## GUI ###########################################################
##################################################################################################################

RIPAgui <- function(){

	##################################################################################################################
	############################# Load shared library with C funcions ################################################
	##################################################################################################################
	dyn.load(system.file("libs/ripa.so",package="ripa"))
	##################################################################################################################

	##################################################################################################################
	############################### Initialization of environment ripaEnv ############################################
	##################################################################################################################

	assign('img1',NULL,envir=ripaEnv)
	assign('img2',NULL,envir=ripaEnv)	
	assign('img3',NULL,envir=ripaEnv)
	assign('visualBands1',c(1,2,3),envir=ripaEnv)
	assign('visualBands2',c(1,2,3),envir=ripaEnv)
	assign('AVIRISbands1',NULL,envir=ripaEnv)
	assign('AVIRISbands2',NULL,envir=ripaEnv)
	assign('imageType1',NULL,envir=ripaEnv)	
	assign('imageType2',NULL,envir=ripaEnv)		
	assign('numbands1',0,envir=ripaEnv)		
	assign('numbands2',0,envir=ripaEnv)		
	assign('regionsList',list(),envir=ripaEnv)			
	assign('regionType',NULL,envir=ripaEnv)			
	assign('regionPointsX',NULL,envir=ripaEnv)	
	assign('regionPointsY',NULL,envir=ripaEnv)
	assign('LANfile',NULL,envir=ripaEnv)	
	assign('tn',NULL,envir=ripaEnv)
	assign('ttregionsbands1',NULL,envir=ripaEnv)
	assign('ttregionsbands2',NULL,envir=ripaEnv)
	assign('regionChoice',NULL,envir=ripaEnv)
	assign('bandsValues',NULL,envir=ripaEnv)
	assign('regionsListAux',NULL,envir=ripaEnv)
	assign('xCoords1',NULL,envir=ripaEnv)
	assign('yCoords1',NULL,envir=ripaEnv)
	assign('bands',NULL,envir=ripaEnv)
	assign('xCoords2',NULL,envir=ripaEnv)
	assign('yCoords2',NULL,envir=ripaEnv)
	assign('ttpoints',NULL,envir=ripaEnv)
	assign('bandSet',NULL,envir=ripaEnv)
	assign('imgtmp',NULL,envir=ripaEnv)
	assign('imgtmp2',NULL,envir=ripaEnv)
	assign('numProcessors',tclVar(detectCores()),envir=ripaEnv)
	assign('useParallel',tclVar(0),envir=ripaEnv)
	assign('cluster',NULL,envir=ripaEnv)
	#registerDoParallel(get('cluster',envir=ripaEnv))

	##################################################################################################################

	##################################################################################################################
	#################################### Required packages ###########################################################
	##################################################################################################################

	if(!is.tclObj(tclRequire("BWidget"))){
		tkmessageBox(title="Error",message="Package BWidget was not found. Please, make sure that this package is installed on your system or add the right path using addTclPath command!",icon="error",type="ok")
		return()
	}
	#if(!is.tclObj(tclRequire("Tktable"))){
	#	tkmessageBox(title="Error",message="Package Tktable was not found. Please, make sure that this package is installed on your system or add the right path using addTclPath command!",icon="error",type="ok")
	#	return()
	#}
	if(!is.tclObj(tclRequire("Img"))){
		tkmessageBox(title="Error",message="Package Img was not found. Please, make sure that this package is installed on your system or add the right path using addTclPath command!",icon="error",type="ok")
		return()
	}

	if(!(require(tkrplot))){
		tkmessageBox(title="Error",message="Package tkrplot was not found. Please, install this package running the command 'install.packages('tkrplot')",icon="error",type="ok")
		return()
	}

	if(!require(e1071)){
		tkmessageBox(title="Error",message="Package e1071 is necessary to complete this action. Please, install this package running the command 'install.packages('e1071')",icon="error",type="ok")
		return()
	}	


	if(!require(rggobi)){
		tkmessageBox(title="Error",message="Package rggobi was not found. Please, install this package running the command 'install.packages('rggobi')",icon="error",type="ok")
		return()
	}

	if(!require(jpeg)){
		tkmessageBox(title="Error",message="Package jpeg is necessary to complete this action. Please, install this package running the command 'install.packages('jpeg')",icon="error",type="ok")
		return()
	}


	if(!require(png)){
		tkmessageBox(title="Error",message="Package png is necessary to complete this action. Please, install this package running the command 'install.packages('png')",icon="error",type="ok")
		return()
	}

	if(!(require(fftw))){
		tkmessageBox(title="Error",message="Package fftw was not found. Please, install this package running the command 'install.packages('fftw')",icon="error",type="ok")
		return()
	}

	if(!(require(foreach))){
		tkmessageBox(title="Error",message="Package foreach was not found. Please, install this package running the command 'install.packages('foreach')",icon="error",type="ok")
		return()
	}

	if(Sys.info()[['sysname']]=='Windows') {	
		if(!(require(doSNOW))){
			tkmessageBox(title="Error",message="Package doSNOW was not found. Please, install this package running the command 'install.packages('doSNOW')",icon="error",type="ok")
			return()
		}
	}
	else {
		if(!(require(doMC))){
			tkmessageBox(title="Error",message="Package doMC was not found. Please, install this package running the command 'install.packages('doMC')",icon="error",type="ok")
			return()
		}
	}

	##################################################################################################################

	#ans <- tkmessageBox(title="Resolution",message="We recomend to use this package with resolution of 1024x768. Do you want to continue?",icon="question",type="yesno")
	#ans <- as.character(ans)
	#if (ans=="no") return()
	
	
	
	##################################################################################################################
	###################################### Main window ###############################################################
	##################################################################################################################
	tt<- tktoplevel()
	tkwm.resizable(tt,0,0)
	tktitle(tt)<-"Image Processing and Analysis in R"

	assign('tn',ttknotebook(tt),envir=ripaEnv)
	tkpack(get('tn',envir=ripaEnv))
	tn_local <- get('tn',envir=ripaEnv)
	tkconfigure(tn_local,width=1015,height=670)
	assign('tn',tn_local,envir=ripaEnv)
	##################################################################################################################
	
	##################################################################################################################
	#################################### Operations Tab ##############################################################
	##################################################################################################################
	
	#################################### Auxiliar Funcitions #########################################################
	

	# Function to build histograms of the region bands
	regionBands <- function(){
		

		tab <- checkTab()
	
		if (tab=="1"){
			tkmessageBox(title="Error",message="Please, use this function only with the second tab!",icon="error",type="ok")
			return()
		}
		
		if (length(get('regionsList',envir=ripaEnv))==0){
			tkmessageBox(title="Error",message="Please, select at least one region in order to use it!",icon="error",type="ok")
			return()
		}
		build1 <- function(){
			assign('ttregionsbands1',tktoplevel(),envir=ripaEnv)
			ttregionsbands1_local <- get('ttregionsbands1',envir=ripaEnv)
			tkwm.geometry(ttregionsbands1_local,"+400+700")
			tkwm.resizable(ttregionsbands1_local,0,0)
			tktitle(ttregionsbands1_local)<-"Regions"
						
			scr1 <- tkscrollbar(ttregionsbands1_local, repeatinterval=5, command=function(...)tkyview(tl1,...))
			tl1<-tklistbox(ttregionsbands1_local,height=6,width=30,selectmode="single",yscrollcommand=function(...)tkset(scr1,...),background="white")
			tkgrid(tklabel(ttregionsbands1_local,text="Choose the region"))
			tkgrid(tl1,scr1)
			tkgrid.configure(scr1,rowspan=6,sticky="nsw")
			assign('ttregionsbands1',ttregionsbands1_local,envir=ripaEnv)
			
			regions <- names(get('regionsList',envir=ripaEnv))
			
			for (i in (1:length(get('regionsList',envir=ripaEnv)))){
				tkinsert(tl1,"end",regions[i])
			}
			tkselection.set(tl1,0)
			
			OnOK1 <- function(){
				assign('regionChoice',as.numeric(tkcurselection(tl1))+1,envir=ripaEnv)
				tkdestroy(get('ttregionsbands1',envir=ripaEnv))
				build2()
			}
			
			OK1.but <-tkbutton(get('ttregionsbands1',envir=ripaEnv),text="   OK   ",command=OnOK1)
			tkgrid(OK1.but)
			tkfocus(get('ttregionsbands1',envir=ripaEnv))
		}
		
		build2 <- function(){
			AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
			assign('ttregionsbands2',tktoplevel(),envir=ripaEnv)

			ttregionsbands2_local <- get('ttregionsbands2',envir=ripaEnv)
			tkwm.geometry(ttregionsbands2_local,"+400+700")
			tkwm.resizable(ttregionsbands2_local,0,0)
			tktitle(ttregionsbands2_local)<-"Bands"
		
			scr2 <- tkscrollbar(ttregionsbands2_local, repeatinterval=5, command=function(...)tkyview(tl2,...))
			tl2<-tklistbox(ttregionsbands2_local,height=6,width=30,selectmode="single",yscrollcommand=function(...)tkset(scr2,...),background="white")
			tkgrid(tklabel(ttregionsbands2_local,text="Choose the band"))
			tkgrid(tl2,scr2)
			tkgrid.configure(scr2,rowspan=6,sticky="nsw")
			assign('ttregionsbands2',ttregionsbands2_local,envir=ripaEnv)
			
			bs <- NULL
			for (i in 1:get('numbands2',envir=ripaEnv)){
				if (is.null(AVIRISbands2_local)) bs <- c(bs,paste("Band ",i,sep=""))
				else bs <- c(bs,paste("Band ",AVIRISbands2_local[i],sep=""))
			}
			
			for (i in (1:get('numbands2',envir=ripaEnv))){
				tkinsert(tl2,"end",bs[i])
			}
			tkselection.set(tl2,0)
			
			OnOK2 <- function(){
				AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
				regionsList_local <- get('regionsList',envir=ripaEnv)
				bandChoice <- as.numeric(tkcurselection(tl2))+1
				tkdestroy(get('ttregionsbands2',envir=ripaEnv))
				
				ttregionhist <- tktoplevel()
				tkwm.geometry(ttregionhist,"+700+0")
				tkwm.resizable(ttregionhist,0,0)
				if (is.null(AVIRISbands2_local)){
					tktitle(ttregionhist)<-paste("Region ",get('regionChoice',envir=ripaEnv)," Band ",bandChoice,sep="")
					histimage <-tkrplot(ttregionhist, function() hist(regionsList_local[[get('regionChoice',envir=ripaEnv)]][[2]][[bandChoice]],main=paste("Band",bandChoice),100,xlab="Values"),vscale=1.04,hscale=0.98)
				}
				else{
					tktitle(ttregionhist)<-paste("Region ",get('regionChoice',envir=ripaEnv)," Band ",AVIRISbands2_local[bandChoice],sep="")
					histimage <-tkrplot(ttregionhist, function() hist(regionsList_local[[get('regionChoice',envir=ripaEnv)]][[2]][[bandChoice]],main=paste("Band",AVIRISbands2_local[bandChoice]),100,xlab="Values"),vscale=1.04,hscale=0.98)
				}
				
				tkpack(histimage)
				build1()
			}
			OK2.but <-tkbutton(get('ttregionsbands2',envir=ripaEnv),text="   OK   ",command=OnOK2)
			tkgrid(OK2.but)
			tkfocus(get('ttregionsbands2',envir=ripaEnv))
		}
		build1()
	}
	#
	
	# Function to build regions bands brushplots
	regionBrush <- function(){
		tab <- checkTab()
		imageType2_local <- get('imageType2',envir=ripaEnv)

		if (tab=="1"){
			tkmessageBox(title="Error",message="Please, use this function only with the second tab!",icon="error",type="ok")
			return()
		}
		
		if (length(get('regionsList',envir=ripaEnv))==0){
			tkmessageBox(title="Error",message="Please, select at least one region in order to use it!",icon="error",type="ok")
			return()
		}
		
		if (imageType2_local=="jpgGrey"){
			tkmessageBox(title="Error",message="Please, use this function with multiple bands images!",icon="error",type="ok")
			return()
		}

		if (imageType2_local=="pngGrey"){
			tkmessageBox(title="Error",message="Please, use this function with multiple bands images!",icon="error",type="ok")
			return()
		}
		
		build1 <- function(){
			assign('ttregionsbands1',tktoplevel(),envir=ripaEnv)
			ttregionsbands1_local <- get('ttregionsbands1',envir=ripaEnv)
			tkwm.geometry(ttregionsbands1_local,"+400+700")
			tkwm.resizable(ttregionsbands1_local,0,0)
			tktitle(ttregionsbands1_local)<-"Regions"
			
			scr1 <- tkscrollbar(ttregionsbands1_local, repeatinterval=5, command=function(...)tkyview(tl1,...))
			tl1<-tklistbox(ttregionsbands1_local,height=6,width=30,selectmode="single",yscrollcommand=function(...)tkset(scr1,...),background="white")
			tkgrid(tklabel(ttregionsbands1_local,text="Choose the region"))
			tkgrid(tl1,scr1)
			tkgrid.configure(scr1,rowspan=6,sticky="nsw")
			assign('ttregionsbands1',ttregionsbands1_local,envir=ripaEnv)
			
			regions <- names(get('regionsList',envir=ripaEnv))
			
			for (i in (1:length(get('regionsList',envir=ripaEnv)))){
				tkinsert(tl1,"end",regions[i])
			}
			tkselection.set(tl1,0)
			
			OnOK1 <- function(){
				assign('regionChoice',as.numeric(tkcurselection(tl1))+1,envir=ripaEnv)
				tkdestroy(get('ttregionsbands1',envir=ripaEnv))
				build2()
			}
			
			OK1.but <-tkbutton(get('ttregionsbands1',envir=ripaEnv),text="   OK   ",command=OnOK1)
			tkgrid(OK1.but)
			tkfocus(get('ttregionsbands1',envir=ripaEnv))
		}
		
		build2 <- function(){
			assign('ttregionsbands2',tktoplevel(),envir=ripaEnv)
			ttregionsbands2_local <- get('ttregionsbands2',envir=ripaEnv)
			tkwm.geometry(ttregionsbands2_local,"+400+700")
			tkwm.resizable(ttregionsbands2_local,0,0)
			tktitle(ttregionsbands2_local)<-"Bands"
			
			scr2 <- tkscrollbar(ttregionsbands2_local, repeatinterval=5, command=function(...)tkyview(tl2,...))
			tl2<-tklistbox(ttregionsbands2_local,height=6,width=30,selectmode="multiple",yscrollcommand=function(...)tkset(scr2,...),background="white")
			tkgrid(tklabel(ttregionsbands2_local,text="Choose the bands"))
			tkgrid(tl2,scr2)
			tkgrid.configure(scr2,rowspan=6,sticky="nsw")
			
			assign('ttregionsbands2',ttregionsbands2_local,envir=ripaEnv)
			AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)

			bs <- NULL
			for (i in 1:get('numbands2',envir=ripaEnv)){
				if (is.null(AVIRISbands2_local)) bs <- c(bs,paste("Band ",i,sep=""))
				else bs <- c(bs,paste("Band ",AVIRISbands2_local[i],sep=""))
			}
			
			for (i in (1:get('numbands2',envir=ripaEnv))){
				tkinsert(tl2,"end",bs[i])
			}
			tkselection.set(tl2,0)
			
			OnOK2 <- function(){
				AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
				regionsList_local <- get('regionsList',envir=ripaEnv)
				bandChoice <- as.vector(as.numeric(tkcurselection(tl2)))+1
				tkdestroy(get('ttregionsbands2',envir=ripaEnv))
				if (length(bandChoice)==1){
					tkmessageBox(title="Error",message="Please, choose at least two bands!",icon="error",type="ok")
					return()
				}
				bs <- NULL
				for (i in 1:length(bandChoice)){
					if (is.null(AVIRISbands2_local)) bs <- c(bs,paste("Band ",bandChoice[i],sep=""))
					else bs <- c(bs,paste("Band ",AVIRISbands2_local[bandChoice[i]],sep=""))
				}
				assign('bandsValues',matrix(nrow=length(regionsList_local[[get('regionChoice',envir=ripaEnv)]][[2]][[1]]),ncol=length(bandChoice)),envir=ripaEnv)
				bandsValues_local <- get('bandsValues',envir=ripaEnv)
				for (i in 1:length(bandChoice)){
					bandsValues_local[,i]<-regionsList_local[[get('regionChoice',envir=ripaEnv)]][[2]][[bandChoice[i]]]
				}
				bandsValues_local <- as.data.frame(bandsValues_local)
				names(bandsValues_local)<-bs
				pairs(bandsValues_local,pch=".")
				assign('bandsValues',bandsValues_local,envir=ripaEnv)
				
				build1()
			}
			OK2.but <-tkbutton(get('ttregionsbands2',envir=ripaEnv),text="   OK   ",command=OnOK2)
			tkgrid(OK2.but)
			tkfocus(get('ttregionsbands2',envir=ripaEnv))
		}
		build1()
	}
	#
	
	# Function to build bands brushplots
	brush <- function(){
		tab <- checkTab()
		imageType1_local <- get('imageType1',envir=ripaEnv)
		imageType2_local <- get('imageType2',envir=ripaEnv)		

		if ((tab=="1" && (imageType1_local=="jpgGrey" || imageType1_local=="pngGrey")) || (tab=="2" && (imageType2_local=="jpgGrey" || imageType2_local=="pngGrey"))){
			tkmessageBox(title="Error",message="Please, use this function with multiple bands images!",icon="error",type="ok")
			return()
		}
			
		build2 <- function(){
			ttbands <- tktoplevel()
			tkwm.geometry(ttbands,"+400+700")
			tkwm.resizable(ttbands,0,0)
			tktitle(ttbands)<-"Bands"
			
			scr <- tkscrollbar(ttbands, repeatinterval=5, command=function(...)tkyview(tl,...))
			tl<-tklistbox(ttbands,height=6,width=30,selectmode="multiple",yscrollcommand=function(...)tkset(scr,...),background="white")
			tkgrid(tklabel(ttbands,text="Choose the bands"))
			tkgrid(tl,scr)
			tkgrid.configure(scr,rowspan=6,sticky="nsw")
			
			bs <- NULL

			tab <- checkTab()
			AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
			AVIRISbands1_local <- get('AVIRISbands1',envir=ripaEnv)
	
			
			if (tab=="1"){
				
				for (i in 1:get('numbands1',envir=ripaEnv)){
					if (is.null(AVIRISbands1_local)) bs <- c(bs,paste("Band ",i,sep=""))
					else bs <- c(bs,paste("Band ",AVIRISbands1_local[i],sep=""))
				}
			}
			else{
				for (i in 1:get('numbands2',envir=ripaEnv)){
					if (is.null(AVIRISbands2_local)) bs <- c(bs,paste("Band ",i,sep=""))
					else bs <- c(bs,paste("Band ",AVIRISbands2_local[i],sep=""))
				}
			}
			if (tab=="1") N=get('numbands1',envir=ripaEnv)
			else N=get('numbands2',envir=ripaEnv)
			for (i in (1:N)){
				tkinsert(tl,"end",bs[i])
			}
			tkselection.set(tl,0)
			
			OnOK <- function(tt){
				AVIRISbands1_local <- get('AVIRISbands1',envir=ripaEnv)
				AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
				bands_local <- get('bands',envir=ripaEnv)
				bandChoice <- as.vector(as.numeric(tkcurselection(tl)))+1
				tkdestroy(tt)
				if (length(bandChoice)==1){
					tkmessageBox(title="Error",message="Please, choose at least two bands!",icon="error",type="ok")
					return()
				}
				bs <- NULL
				if (tab=="1"){
					for (i in 1:length(bandChoice)){
						if (is.null(AVIRISbands1_local)) bs <- c(bs,paste("Band ",bandChoice[i],sep=""))
						else bs <- c(bs,paste("Band ",AVIRISbands1_local[bandChoice[i]],sep=""))
					}
				}
				else{
					for (i in 1:length(bandChoice)){
						if (is.null(AVIRISbands2_local)) bs <- c(bs,paste("Band ",bandChoice[i],sep=""))
						else bs <- c(bs,paste("Band ",AVIRISbands2_local[bandChoice[i]],sep=""))
					}
				}
				if (tab=="1"){
					img1_local <- get('img1',envir=ripaEnv)
					if (is.null(AVIRISbands1_local)) assign('bandsValues',data.frame(a=as.vector(img1_local[,,bandChoice[1]])),envir=ripaEnv)
					else assign('bandsValues',data.frame(a = as.vector(bands_local[[bandChoice[1]]]@data)),envir=ripaEnv)
				}
				else{
					img3_local <- get('img3',envir=ripaEnv)
					if (is.null(AVIRISbands2_local)) assign('bandsValues',data.frame(a=as.vector(img3_local[,,bandChoice[1]])),envir=ripaEnv)
					else assign('bandsValues',data.frame(a = as.vector(bands_local[[bandChoice[1]]]@data)),envir=ripaEnv)
				}
				img1_local <- get('img1',envir=ripaEnv)
				img3_local <- get('img3',envir=ripaEnv)	
				bandsValues_local <- get('bandsValues',envir=ripaEnv)	
				for (i in 2:length(bandChoice)){
					if (tab=="1"){		
						if (is.null(AVIRISbands1_local)) bandsValues_local[,i]<-as.vector(img1_local[,,bandChoice[i]])
						else bandsValues_local[,i] <- as.vector(bands_local[[bandChoice[i]]]@data)
					}
					else{
						if (is.null(AVIRISbands2_local)) bandsValues_local[,i]<-as.vector(img3_local[,,bandChoice[i]])
						else bandsValues_local[,i] <- as.vector(bands_local[[bandChoice[i]]]@data)
					}
				}
				names(bandsValues_local)<-bs
				pairs(bandsValues_local,pch=".")
				assign('bandsValues',bandsValues_local,envir=ripaEnv)
				
			}
			OK.but <-tkbutton(ttbands,text="   OK   ",command=function()OnOK(ttbands))
			tkgrid(OK.but)
			tkfocus(ttbands)
		}
		build2()
	}
	#

	# Function to shows bands histogram and statistics
	imageBands <- function(){
		tab <- checkTab()
	
		if ((tab=="1" && is.null(get('img1',envir=ripaEnv))) || (tab=="2" && is.null(get('img3',envir=ripaEnv))) ){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
		
		ttbands<-tktoplevel()
		tkwm.geometry(ttbands,"+600+700")
		tkwm.resizable(ttbands,0,0)
		tktitle(ttbands)<-"List of bands"
		scr <- tkscrollbar(ttbands, repeatinterval=5, command=function(...)tkyview(tl,...))
		tl<-tklistbox(ttbands,height=6,width=30,selectmode="single",yscrollcommand=function(...)tkset(scr,...),background="white")
		tkgrid(tklabel(ttbands,text="Choose the band"))
		tkgrid(tl,scr)
		tkgrid.configure(scr,rowspan=6,sticky="nsw")
		bands <- NULL

		
		if (tab=="1"){
			AVIRISbands1_local <- get('AVIRISbands1',envir=ripaEnv)
			
			for (i in 1:get('numbands1',envir=ripaEnv)){
				if (is.null(AVIRISbands1_local)) bands <- c(bands,paste("Band ",i,sep=""))
				else bands <- c(bands,paste("Band ",AVIRISbands1_local[i],sep=""))
			}
			for (i in (1:get('numbands1',envir=ripaEnv))){
				tkinsert(tl,"end",bands[i])
			}
		}
		
		if (tab=="2"){
			AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
			for (i in 1:get('numbands2',envir=ripaEnv)){
				if (is.null(AVIRISbands2_local)) bands <- c(bands,paste("Band ",i,sep=""))
				else bands <- c(bands,paste("Band ",AVIRISbands2_local[i],sep=""))
			}
			for (i in (1:get('numbands2',envir=ripaEnv))){
				tkinsert(tl,"end",bands[i])
			}
		}
		
		tkselection.set(tl,0)
		
		OnOK <- function(){
			bandChoice <- as.numeric(tkcurselection(tl))+1
			AVIRISbands1_local <- get('AVIRISbands1',envir=ripaEnv)
			AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
			tab <- checkTab()
		
			ttstatistics<-tktoplevel()
			tkwm.geometry(ttstatistics,"+0+700")
			tkwm.resizable(ttstatistics,0,0)
			if (tab=="1"){
				if (is.null(AVIRISbands1_local)) tktitle(ttstatistics)<-paste("Band ",bandChoice," statistics",sep="")
				else tktitle(ttstatistics)<-paste("Band ",AVIRISbands1_local[bandChoice]," statistics",sep="")
			}
			if (tab=="2"){
				if (is.null(AVIRISbands2_local)) tktitle(ttstatistics)<-paste("Band ",bandChoice," statistics",sep="")
				else tktitle(ttstatistics)<-paste("Band ",AVIRISbands2_local[bandChoice]," statistics",sep="")
			}
			scr <- tkscrollbar(ttstatistics, repeatinterval=5,command=function(...)tkyview(txt,...))
			txt <- tktext(ttstatistics,bg="white",font="courier",width=40,height=10,yscrollcommand=function(...)tkset(scr,...))
			tkgrid(txt,scr)
			tkgrid.configure(scr,sticky="ns")
			
			

			if (tab=="1"){
				ttband <- tktoplevel()
				tkwm.geometry(ttband,"+0+0")
				tkwm.resizable(ttband,0,0)
				if (is.null(AVIRISbands1_local)) tktitle(ttband)<-paste("Band ",bandChoice,sep="")
				else tktitle(ttband)<-paste("Band ",AVIRISbands1_local[bandChoice],sep="")
				img1_local <- get('img1',envir=ripaEnv)
				bandimage <-tkrplot(ttband, function() plot(imagematrix(img1_local[,,bandChoice])),vscale=1.04,hscale=0.98)
				tkpack(bandimage)
				
				nR <- length(img1_local[,,bandChoice])
				minR <- min(img1_local[,,bandChoice])
				maxR <- max(img1_local[,,bandChoice])
				meanR <- mean(img1_local[,,bandChoice])
				medianR <- median(img1_local[,,bandChoice])
				deviationR <- sd(as.vector(img1_local[,,bandChoice]))
				mDeviationR <- mad(img1_local[,,bandChoice])
				kurtosisR <- kurtosis(as.vector(img1_local[,,bandChoice]))
				skewnessR <- skewness(as.vector(img1_local[,,bandChoice]))
				
				ttbandhist <- tktoplevel()
				tkwm.geometry(ttbandhist,"+500+0")
				tkwm.resizable(ttbandhist,0,0)
				if (is.null(AVIRISbands1_local)){
					tktitle(ttbandhist)<-paste("Band ",bandChoice,sep="")
					histimage <-tkrplot(ttbandhist, function() hist(img1_local[,,bandChoice],main=paste("Band",bandChoice),100,xlab="Values"),vscale=1.04,hscale=0.98)
				}
				else{
					tktitle(ttbandhist)<-paste("Band ",AVIRISbands1_local[bandChoice],sep="")
					histimage <-tkrplot(ttbandhist, function() hist(img1_local[,,bandChoice],main=paste("Band",AVIRISbands1_local[bandChoice]),100,xlab="Values"),vscale=1.04,hscale=0.98)
				}
				tkpack(histimage)
			}
			if (tab=="2"){
				img3_local <- get('img3',envir=ripaEnv)
				ttband <- tktoplevel()
				tkwm.geometry(ttband,"+0+0")
				tkwm.resizable(ttband,0,0)
				if (is.null(AVIRISbands2_local)) tktitle(ttband)<-paste("Band ",bandChoice,sep="")
				else tktitle(ttband)<-paste("Band ",AVIRISbands2_local[bandChoice],sep="")
				bandimage <-tkrplot(ttband, function() plot(imagematrix(img3_local[,,bandChoice])),vscale=1.04,hscale=0.98)
				tkpack(bandimage)
				
				nR <- length(img3_local[,,bandChoice])
				minR <- min(img3_local[,,bandChoice])
				maxR <- max(img3_local[,,bandChoice])
				meanR <- mean(img3_local[,,bandChoice])
				medianR <- median(img3_local[,,bandChoice])
				deviationR <- sd(as.vector(img3_local[,,bandChoice]))
				mDeviationR <- mad(img3_local[,,bandChoice])
				kurtosisR <- kurtosis(as.vector(img3_local[,,bandChoice]))
				skewnessR <- skewness(as.vector(img3_local[,,bandChoice]))
				
				ttbandhist <- tktoplevel()
				tkwm.geometry(ttbandhist,"+500+0")
				tkwm.resizable(ttbandhist,0,0)
				if (is.null(AVIRISbands2_local)){
					tktitle(ttbandhist)<-paste("Band ",bandChoice,sep="")
					histimage <-tkrplot(ttbandhist, function() hist(img3_local[,,bandChoice],main=paste("Band",bandChoice),100,xlab="Values"),vscale=1.04,hscale=0.98)
				} else{
					tktitle(ttbandhist)<-paste("Band ",AVIRISbands2_local[bandChoice],sep="")
					histimage <-tkrplot(ttbandhist, function() hist(img3_local[,,bandChoice],main=paste("Band",AVIRISbands2_local[bandChoice]),100,xlab="Values"),vscale=1.04,hscale=0.98)
				}
				tkpack(histimage)
			}

			tkconfigure(txt,state="normal")
			tkinsert(txt,"end",paste("N = ",nR,"\n",sep=""))
			tkinsert(txt,"end",paste("Min = ",minR,"\n",sep=""))
			tkinsert(txt,"end",paste("Max = ",maxR,"\n",sep=""))
			tkinsert(txt,"end",paste("Mean = ",meanR,"\n",sep=""))
			tkinsert(txt,"end",paste("Median = ",medianR,"\n",sep=""))
			tkinsert(txt,"end",paste("Deviation = ",deviationR,"\n",sep=""))
			tkinsert(txt,"end",paste("Median Deviation = ",mDeviationR,"\n",sep=""))
			tkinsert(txt,"end",paste("Kurtosis = ",kurtosisR,"\n",sep=""))
			tkinsert(txt,"end",paste("Skewness = ",skewnessR,"\n\n",sep=""))
			
			tkconfigure(txt,state="disable")
			
		}
		OK.but <-tkbutton(ttbands,text=   " OK ",command=OnOK)
		tkgrid(OK.but)
		tkfocus(ttbands)
	}
	#
	
	# Function to calculate the covariance matrix of a region
	covMatrix <- function(){
		tab <- checkTab()
		imageType2_local <- get('imageType2',envir=ripaEnv)		

		if (tab=="1"){
			tkmessageBox(title="Error",message="Please, use this function only with the second tab!",icon="error",type="ok")
			return()
		}
		
		if (is.null(imageType2_local)){
			tkmessageBox(title="Error",message="Please, open an image in oder to use it!",icon="error",type="ok")
			return()
		}
		
		if (length(get('regionsList',envir=ripaEnv))==0){
			tkmessageBox(title="Error",message="Please, select at least one region in order to use it!",icon="error",type="ok")
			return()
		}
		
		if (imageType2_local=="jpgGrey" || imageType2_local=="pngGrey"){
			tkmessageBox(title="Error",message="Please, use this function only with multiple bands images!",icon="error",type="ok")
			return()
		}
		
		ttregions<-tktoplevel()
		tkwm.geometry(ttregions,"+400+700")
		tkwm.resizable(ttregions,0,0)
		tktitle(ttregions)<-"Regions"
		
		scr <- tkscrollbar(ttregions, repeatinterval=5, command=function(...)tkyview(tl,...))
		tl<-tklistbox(ttregions,height=6,width=30,selectmode="single",yscrollcommand=function(...)tkset(scr,...),background="white")
		tkgrid(tklabel(ttregions,text="Choose the region"))
		tkgrid(tl,scr)
		tkgrid.configure(scr,rowspan=6,sticky="nsw")
		
		regions <- as.vector(names(get('regionsList',envir=ripaEnv)))
		
		for (i in (1:length(get('regionsList',envir=ripaEnv)))){
			tkinsert(tl,"end",regions[i])
		}
		tkselection.set(tl,0)
		
		OnOK <- function(){

			AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
			regionsList_local <- get('regionsList',envir=ripaEnv)
			assign('regionChoice',as.numeric(tkcurselection(tl))+1,envir=ripaEnv)
			tkdestroy(ttregions)
			bands <- NULL
			for (i in 1:get('numbands2',envir=ripaEnv)){
				if (is.null(AVIRISbands2_local)) bands <- c(bands,paste("Band ",i,sep=""))
				else bands <- c(bands,paste("Band ",AVIRISbands2_local[i],sep=""))
			}
			
			Matrix <- matrix(nrow=length(regionsList_local[[get('regionChoice',envir=ripaEnv)]][[2]][[1]]),ncol=get('numbands2',envir=ripaEnv))
			for (i in 1:get('numbands2',envir=ripaEnv)){
				Matrix[,i]<-regionsList_local[[get('regionChoice',envir=ripaEnv)]][[2]][[i]]
			}
			Matrix <- as.data.frame(Matrix)
			names(Matrix)<-bands
			covMat <- cov(Matrix)
			displayInTable <- function(tclarray,title="",height=-1,width=-1,nrow=-1,ncol=-1){
				tttable <- tktoplevel()
				tkwm.title(tttable,title)
				tkwm.geometry(tttable,"+0+70")
				tkwm.resizable(tttable,0,0)
				table <- tkwidget(tttable,"table",rows=nrow,cols=ncol,titlerows=1,titlecols=1,colwidth=20,
					height=height+1,width=width+1,
					xscrollcommand=function(...) tkset(xscr,...),yscrollcommand=function(...) tkset(yscr,...))
				xscr <-tkscrollbar(tttable,orient="horizontal", command=function(...)tkxview(table,...))
				yscr <- tkscrollbar(tttable,command=function(...)tkyview(table,...))
				tkgrid(table,yscr)
				tkgrid.configure(yscr,sticky="nsw")
				tkgrid(xscr,sticky="new")
				tkconfigure(table,variable=tclarray,background="white",selectmode="extended")
				return (table)
			}
			tclArray <- tclArray()
			for (i in 1:get('numbands2',envir=ripaEnv)) tclArray[[0,i]] <- names(Matrix)[i]
			for (i in 1:get('numbands2',envir=ripaEnv)) tclArray[[i,0]] <- names(Matrix)[i]
			for (i in 1:get('numbands2',envir=ripaEnv)){
				for (j in 1:get('numbands2',envir=ripaEnv)){
					tclArray[[i,j]] <- covMat[i,j]
				}
			}
			table <- displayInTable(tclArray,title=paste("Covariance Matrix ","(Region ",get('regionChoice',envir=ripaEnv),")"),nrow=get('numbands2',envir=ripaEnv)+1,ncol=get('numbands2',envir=ripaEnv)+1)
			
			tkconfigure(table, state="disabled")
		}
		
		OK.but <-tkbutton(ttregions,text="   OK   ",command=OnOK)
		tkgrid(OK.but)
		tkfocus(ttregions)
	}
	#
	
	# Function to change the actual visual bands
	changeBands <- function(){
		tab <- checkTab()
		imageType1_local <- get('imageType1',envir=ripaEnv)		
		imageType2_local <- get('imageType2',envir=ripaEnv)		
	
		if ((tab=="1" && is.null(imageType1_local)) || (tab=="2" && is.null(imageType2_local))){
			tkmessageBox(title="Error",message="Please, open an image in oder to use it!",icon="error",type="ok")
			return()
		}
		
		if ((tab=="1" && (imageType1_local=="jpgGrey" || imageType1_local=="pngGrey")) || (tab=="2" && (imageType2_local=="jpgGrey" || imageType2_local=="pngGrey")) || (tab=="1" && get('numbands1',envir=ripaEnv)==3) || (tab=="2" && get('numbands2',envir=ripaEnv)==3)){
			tkmessageBox(title="Error",message="Please, use this function with more than three bands images!",icon="error",type="ok")
			return()
		}
		
		if (tab=="1") assign('visualBands1',NULL,envir=ripaEnv)
		if (tab=="2") assign('visualBands2',NULL,envir=ripaEnv)
		
		build <- function(count){
			ttregions<-tktoplevel()
			tkwm.geometry(ttregions,"+400+300")
			tkwm.resizable(ttregions,0,0)
			tktitle(ttregions)<-paste("Bands ",count,sep="")
			
			scr <- tkscrollbar(ttregions, repeatinterval=5, command=function(...)tkyview(tl,...))
			tl<-tklistbox(ttregions,height=6,width=30,selectmode="single",yscrollcommand=function(...)tkset(scr,...),background="white")
			tkgrid(tklabel(ttregions,text="Choose the band"))
			tkgrid(tl,scr)
			tkgrid.configure(scr,rowspan=6,sticky="nsw")
			
			bands <- NULL
			AVIRISbands1_local <- get('AVIRISbands1',envir=ripaEnv)
			AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
			if (tab=="1"){
				for (i in 1:get('numbands1',envir=ripaEnv)){
					if (is.null(AVIRISbands1_local)) bands <- c(bands,paste("Band ",i,sep=""))
					else bands <- c(bands,paste("Band ",AVIRISbands1_local[i],sep=""))
				}
				for (i in 1:get('numbands1',envir=ripaEnv))
					tkinsert(tl,"end",bands[i])
			}
			if (tab=="2"){
				for (i in 1:get('numbands2',envir=ripaEnv)){
					if (is.null(AVIRISbands2_local)) bands <- c(bands,paste("Band ",i,sep=""))
					else bands <- c(bands,paste("Band ",AVIRISbands2_local[i],sep=""))
				}
				for (i in 1:get('numbands2',envir=ripaEnv))
					tkinsert(tl,"end",bands[i])
				
			}
			
			tkselection.set(tl,0)
			
			OnOK <- function(){
				
				choice <- as.numeric(tkcurselection(tl))+1
				visualBands1_local <- get('visualBands1',envir=ripaEnv)
				visualBands2_local <- get('visualBands2',envir=ripaEnv)
				if (count==2){
					
					if ((tab=="1" && choice==visualBands1_local[1]) || (tab=="2" && choice==visualBands2_local[1])){
						tkmessageBox(title="Error",message="This band was already chosen! Please, select another band",icon="error",type="ok")
						return()
					}
				}
				
				if (count==3){
					if ((tab=="1" && choice==visualBands1_local[1]) || (tab=="2" && choice==visualBands2_local[1]) || (tab=="1" && choice==visualBands1_local[2]) || (tab=="2" && choice==visualBands2_local[2])){
						tkmessageBox(title="Error",message="This band was already chosen! Please, select another band",icon="error",type="ok")
						return()
					}
				}
				if (tab=="1"){
					visualBands1_local <- c(visualBands1_local,choice)
					assign('visualBands1',c(visualBands1_local,choice),envir=ripaEnv)
				}
				if (tab=="2"){ 
					visualBands2_local <- c(visualBands2_local,choice)
					assign('visualBands2',c(visualBands2_local,choice),envir=ripaEnv)
				}
				tkdestroy(ttregions)
				if (count!=3) build(count+1)
				else{
					if (tab=="1"){
						tkdestroy(Frame3)
						Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
						tkconfigure(Frame3,width=485,height=510)
						tkplace(Frame3,x=510,y=30)
						
						Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
						tkplace(Frame4,x=1,y=1)
						
						tkdestroy(Frame1)
						Frame1 <- tkframe(lb1,relief="groove",borderwidth=2)
						tkconfigure(Frame1,width=480,height=510)
						tkplace(Frame1,x=10,y=30)
						
						Frame2 <- tkframe(Frame1,relief="groove",borderwidth=0)
						tkplace(Frame2,x=1,y=1)
						img1_local <- get('img1',envir=ripaEnv)
						auximg <- array(c(img1_local[,,visualBands1_local[1]],img1_local[,,visualBands1_local[2]],img1_local[,,visualBands1_local[3]]),c(nrow(img1_local),ncol(img1_local),3))

						image <-tkrplot(Frame2, function() plot(imagematrix(auximg,type="rgb")),vscale=1.04,hscale=0.98)
						tkpack(image)

						img2_local <- get('img2',envir=ripaEnv)
						auximg <- array(c(img2_local[,,visualBands2_local[1]],img2_local[,,visualBands2_local[2]],img2_local[,,visualBands2_local[3]]),c(nrow(img2_local),ncol(img2_local),3))
						image <-tkrplot(Frame4, function() plot(imagematrix(auximg,type="rgb")),vscale=1.04,hscale=0.98)
						tkpack(image)
					}
					if (tab=="2"){
						img3_local <- get('img3',envir=ripaEnv)
						assign('regionsList',list(),envir=ripaEnv)
						assign('regionsListAux',list(),envir=ripaEnv)
						
						tkconfigure(statisticsTxt,state="normal")
						tkdelete(statisticsTxt,"1.0","end")
						tkconfigure(statisticsTxt,state="disable")
						
						tkdestroy(Frame5)
					
						Frame5 <- tkframe(lb2,relief="groove",borderwidth=2)
						tkconfigure(Frame5,width=480,height=510)
						tkplace(Frame5,x=10,y=30)
				
						Frame6 <- tkframe(Frame5,relief="groove",borderwidth=0)
						tkplace(Frame6,x=1,y=1)
						
						auximg <- array(c(img3_local[,,visualBands2_local[1]],img3_local[,,visualBands2_local[2]],img3_local[,,visualBands2_local[3]]),c(nrow(img3_local),ncol(img3_local),3))
						image <-tkrplot(Frame6, function() plot.imagematrix(imagematrix(auximg,type="rgb",noclipping=TRUE)),vscale=1.04,hscale=0.98)
						tkpack(image)
					}
				}
			}
			
			OK.but <-tkbutton(ttregions,text="   OK   ",command=OnOK)
			tkgrid(OK.but)
			tkfocus(ttregions)
		}
		build(1)
	}
	#
	
	# Function for the brightness and contrast slider
	sliderFunction <- function(){
	
		visualBands1_local <- get('visualBands1',envir=ripaEnv)

		tab <- checkTab()

		if ((is.null(get('img1',envir=ripaEnv)) && tab=="1") || (is.null(get('img3',envir=ripaEnv)) && tab=="2")){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
	
		cont <- (as.numeric(tclvalue(SliderValue2)))
		bri <- as.numeric(tclvalue(SliderValue))
		
		assign('img2',contBriImg(get('img1',envir=ripaEnv),cont,bri),envir=ripaEnv)
		img2_local <- get('img2',envir=ripaEnv)
		
		tkdestroy(Frame3)	
		Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
		tkconfigure(Frame3,width=485,height=510)
		tkplace(Frame3,x=510,y=30)
		Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
		tkplace(Frame4,x=1,y=1)
		
		if (get('numbands1',envir=ripaEnv)>3){
			auximg <- array(c(img2_local[,,visualBands1_local[1]],img2_local[,,visualBands1_local[2]],img2_local[,,visualBands1_local[3]]),c(nrow(img2_local),ncol(img2_local),get('numbands1',envir=ripaEnv)))
			image <-tkrplot(Frame4, function() plot(imagematrix(auximg)),vscale=1.04,hscale=0.98)
		}else{
			image <-tkrplot(Frame4, function() plot(imagematrix(img2_local)),vscale=1.04,hscale=0.98)
		}
		tkpack(image)
	}
	#
	
	# Function to quit the interface
	quit <- function(){
		tkdestroy(tt)
		stopCluster(get('cluster',envir=ripaEnv))
	}
	#
	
	# Function to show the actual pixel value of the mouse
	pixelsValues <- function(){
		plotFunction1 <- function(){

			visualBands1_local <- get('visualBands1',envir=ripaEnv)
			visualBands2_local <- get('visualBands2',envir=ripaEnv)

			params <- par(bg="white")
			
			tab <- checkTab()
			img1_local <- get('img1',envir=ripaEnv)
			img3_local <- get('img3',envir=ripaEnv)
			if (tab=="1" && get('numbands1',envir=ripaEnv)>3){
				auximg <- array(c(img1_local[,,visualBands1_local[1]],img1_local[,,visualBands1_local[2]],img1_local[,,visualBands1_local[3]]),c(nrow(img1_local),ncol(img1_local),get('numbands1',envir=ripaEnv)))
				plot(imagematrix(auximg))
			}
			
			if (tab=="2" && get('numbands2',envir=ripaEnv)>3){
				auximg <- array(c(img3_local[,,visualBands2_local[1]],img3_local[,,visualBands2_local[2]],img3_local[,,visualBands2_local[3]]),c(nrow(img3_local),ncol(img3_local),get('numbands2',envir=ripaEnv)))
				plot(imagematrix(auximg))
			}
			if (tab=="1" && get('numbands1',envir=ripaEnv)<=3)
				plot(imagematrix(img1_local))
			if (tab=="2" && get('numbands2',envir=ripaEnv)<=3)
				plot(imagematrix(img3_local))
			
			parPlotSize1 <<- par("plt")
			usrCoords1   <<- par("usr")
			par(params)
		}
		
		OnLeftClick <- function(x,y){
			img1_local <- get('img1',envir=ripaEnv)
			img3_local <- get('img3',envir=ripaEnv)
			AVIRISbands1_local <- get('AVIRISbands1',envir=ripaEnv)
			AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
			imgtmp_local <- get('imgtmp',envir=ripaEnv)

			xClick <- x
			yClick <- y
			width  <- as.numeric(tclvalue(tkwinfo("reqwidth",imgtmp_local)))
			height <- as.numeric(tclvalue(tkwinfo("reqheight",imgtmp_local)))
			
			xMin <- parPlotSize1[1] * width
			xMax <- parPlotSize1[2] * width
			yMin <- parPlotSize1[3] * height
			yMax <- parPlotSize1[4] * height
			
			rangeX <- usrCoords1[2] - usrCoords1[1]
			rangeY <- usrCoords1[4] - usrCoords1[3]
			
			imgXcoords <- (get('xCoords1',envir=ripaEnv)-usrCoords1[1])*(xMax-xMin)/rangeX + xMin
			imgYcoords <- (get('yCoords1',envir=ripaEnv)-usrCoords1[3])*(yMax-yMin)/rangeY + yMin
			
			xClick <- as.numeric(xClick)+0.5
			yClick <- as.numeric(yClick)+0.5
			yClick <- height - yClick
			
			xPlotCoord <- usrCoords1[1]+(xClick-xMin)*rangeX/(xMax-xMin)
			yPlotCoord <- usrCoords1[3]+(yClick-yMin)*rangeY/(yMax-yMin)
			img1_local <- get('img1',envir=ripaEnv)
			a <- round(xPlotCoord)
			b <- round(nrow(img1_local) - yPlotCoord)
			
			tab <- checkTab()
			
			if (tab=="1"){
				if (a<0 || b<0 || a>ncol(img1_local) || b>nrow(img1_local)) tkmessageBox(title="Error",message="Please, click inside the image!",icon="error",type="ok")
				else{
					aux2 <- NULL
					for (i in 1:get('numbands1',envir=ripaEnv)){
						if (is.null(AVIRISbands1_local)) aux2 <- paste(aux2,"Band ",i,": ",img1_local[b,a,i],"\n",sep="")
						else aux2 <- paste(aux2,"Band ",AVIRISbands1_local[i],": ",img1_local[b,a,i],"\n",sep="")
					}
					tkmessageBox(title="Value",message=aux2,icon="info",type="ok")
				}
			}
			if (tab=="2"){
				if (a<0 || b<0 || a>ncol(img3_local) || b>nrow(img3_local)) tkmessageBox(title="Error",message="Please, click inside the image!",icon="error",type="ok")
				else{
					aux2 <- NULL
					for (i in 1:get('numbands2',envir=ripaEnv)){
						if (is.null(AVIRISbands2_local)) aux2 <- paste(aux2,"Band",i,": ",img3_local[b,a,i],"\n",sep="")
						else aux2 <- paste(aux2,"Band",AVIRISbands2_local[i],": ",img3_local[b,a,i],"\n",sep="")
					}
					tkmessageBox(title="Value",message=aux2,icon="info",type="ok")
				}
			}
		}
		
		tab <- checkTab()
		img1_local <- get('img1',envir=ripaEnv)
		img3_local <- get('img3',envir=ripaEnv)
		if ((tab=="1" && is.null(img1_local)) || (tab=="2" && is.null(img3_local))){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
		
		ttpixels <- tktoplevel()
		tkwm.geometry(ttpixels,"+150+0")
		tkwm.title(ttpixels,"Click on a point to get the pixel value")
		
		if (tab=="1"){
			assign('xCoords1',(1:ncol(img1_local)),envir=ripaEnv)
			assign('yCoords1',(1:nrow(img1_local)),envir=ripaEnv)
		}
		if (tab=="2"){
			assign('xCoords1',(1:ncol(img3_local)),envir=ripaEnv)
			assign('yCoords1',(1:nrow(img3_local)),envir=ripaEnv)
		}
		
		parPlotSize1 <- c()
		usrCoords1 <- c()
	
		assign('imgtmp',tkrplot(ttpixels,fun=plotFunction1,hscale=1.5,vscale=1.5),envir=ripaEnv)
		imgtmp_local <- get('imgtmp',envir=ripaEnv)
		tkgrid(imgtmp_local)
	
		tkbind(imgtmp_local, "<Button-1>",OnLeftClick)
		tkconfigure(imgtmp_local,cursor="hand2")
		assign('imgtmp',imgtmp_local,envir=ripaEnv)
	}
	#
	
	# Function to edit pixels values
	editPixels <- function(){
		tab <- checkTab()
		AVIRISbands1_local <- get('AVIRISbands1',envir=ripaEnv)
	
		if (tab=="2"){
			tkmessageBox(title="Error",message="Please, use this function only with the first tab!",icon="error",type="ok")
			return()
		}
		if (is.null(get('img1',envir=ripaEnv))){
			tkmessageBox(title="Error",message="Please, select an image in order to use it!",icon="error",type="ok")
			return()
		}
		
		ttbands<-tktoplevel()
		tkwm.geometry(ttbands,"+400+700")
		tkwm.resizable(ttbands,0,0)
		tktitle(ttbands)<-"Bands"
		
		scr <- tkscrollbar(ttbands, repeatinterval=5, command=function(...)tkyview(tl,...))
		tl<-tklistbox(ttbands,height=6,width=30,selectmode="single",yscrollcommand=function(...)tkset(scr,...),background="white")
		tkgrid(tklabel(ttbands,text="Choose the band"))
		tkgrid(tl,scr)
		tkgrid.configure(scr,rowspan=6,sticky="nsw")
		
		bands<-NULL
		
		for (i in 1:get('numbands1',envir=ripaEnv)){
			if (is.null(AVIRISbands1_local)) bands <- c(bands,paste("Band ",i,sep=""))
			else bands <- c(bands,paste("Band ",AVIRISbands1_local[i],sep=""))
		}
		
		for (i in (1:get('numbands1',envir=ripaEnv))){
			tkinsert(tl,"end",bands[i])
		}
		tkselection.set(tl,0)
		
		OnOK <- function(tt){
			
			onok <- function(){

				visualBands1_local <- get('visualBands1',envir=ripaEnv)
				img1_local <- get('img1',envir=ripaEnv)
				img2_local <- get('img2',envir=ripaEnv)
				for (i in 1:nrow(img1_local)){
					for (j in 1:ncol(img1_local)){
						img2_local[i,j,bandChoice] <- as.numeric(tclvalue(myarray[[i,j]]))
					}
				}
				
				tkdestroy(Frame3)	
				Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
				tkconfigure(Frame3,width=485,height=510)
				tkplace(Frame3,x=510,y=30)
				Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
				tkplace(Frame4,x=1,y=1)
				
				if (get('numbands1',envir=ripaEnv)>3){
					auximg <- array(c(img2_local[,,visualBands1_local[1]],img2_local[,,visualBands1_local[2]],img2_local[,,visualBands1_local[3]]),c(nrow(img2_local),ncol(img2_local),get('numbands1',envir=ripaEnv)))
					image <-tkrplot(Frame4, function() plot(imagematrix(auximg)),vscale=1.04,hscale=0.98)
				}else{
					image <-tkrplot(Frame4, function() plot(imagematrix(img2_local)),vscale=1.04,hscale=0.98)
				}
				tkpack(image)
			}
			
			bandChoice <- as.numeric(tkcurselection(tl))+1
			tkdestroy(tt)
			
			myarray <- tclArray()
			img1_local <- get('img1',envir=ripaEnv)
			for (i in 1:ncol(img1_local)) myarray[[0,i]] <- i
			for (i in 1:nrow(img1_local)) myarray[[i,0]] <- i
			for (i in 1:nrow(img1_local)){
				for (j in 1:ncol(img1_local)){
					myarray[[i,j]] <- img1_local[i,j,bandChoice]
				}
			}
			
			tttable <- tktoplevel()
			tkwm.title(tttable,paste("Pixels Values (Band ",bandChoice,")",sep=""))
			tkwm.geometry(tttable,"+80+50")
			tkwm.resizable(tttable,0,0)
			table <- tkwidget(tttable,"table",rows=nrow(img1_local)+1,cols=ncol(img1_local)+1,titlerows=1,titlecols=1,colwidth=20,
				height=0,width=0,
				xscrollcommand=function(...) tkset(xscr,...),yscrollcommand=function(...) tkset(yscr,...))
			xscr <-tkscrollbar(tttable,orient="horizontal", command=function(...)tkxview(table,...))
			yscr <- tkscrollbar(tttable,command=function(...)tkyview(table,...))
			tkgrid(table,yscr)
			tkgrid.configure(yscr,sticky="nsw")
			tkgrid(xscr,sticky="new")
			tkconfigure(table,variable=myarray,background="white",selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"",resizeborders="none")
			
			button <-tkbutton(tttable,text="   OK   ",command=onok)
			tkgrid(button)
			
		}
		
		OK.but <-tkbutton(ttbands,text="   OK   ",command=function()OnOK(ttbands))
		tkgrid(OK.but)
		tkfocus(ttbands)

	}
	#
	
	# Function to build dynamic graphics from GGOBI
	dynFunc <- function(func){


		tab <- checkTab()

		AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
		regionsList_local <- get('regionsList',envir=ripaEnv)

		if (tab=="1"){
			tkmessageBox(title="Error",message="Please, use this menu only with the second tab!",icon="error",type="ok")
			return()
		}
		
		if (is.null(get('img3',envir=ripaEnv))){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
		
		if (get('numbands2',envir=ripaEnv)==1){
			tkmessageBox(title="Error",message="Graphics available only for images with more than one band!",icon="error",type="ok")
			return()
		}
		
		if (length(get('regionsList',envir=ripaEnv))==0){
			tkmessageBox(title="Error",message="There are not regions!",icon="error",type="ok")
			return()
		}
		
		names <- vector()
		vars <- vector()
		dataf <- vector()


		if (func!=9){
			for (i in 1:get('numbands2',envir=ripaEnv)){
				auxVector <- vector()
				for (j in 1:length(regionsList_local)){
					if (func==1) auxVector<-c(auxVector,min(regionsList_local[[j]][[2]][[i]]))
					if (func==2) auxVector<-c(auxVector,max(regionsList_local[[j]][[2]][[i]]))
					if (func==3) auxVector<-c(auxVector,mean(regionsList_local[[j]][[2]][[i]]))
					if (func==4) auxVector<-c(auxVector,median(regionsList_local[[j]][[2]][[i]]))
					if (func==5) auxVector<-c(auxVector,sd(regionsList_local[[j]][[2]][[i]]))
					if (func==6) auxVector<-c(auxVector,mad(regionsList_local[[j]][[2]][[i]]))
					if (func==7) auxVector<-c(auxVector,kurtosis(regionsList_local[[j]][[2]][[i]]))
					if (func==8) auxVector<-c(auxVector,skewness(regionsList_local[[j]][[2]][[i]]))
				}
				if (is.null(AVIRISbands2_local)){
					assign(paste("Band",i,sep=""),auxVector)
					names <- c(names,paste("Band",i,sep=""))
				} else{
					assign(paste("Band",AVIRISbands2_local[i],sep=""),auxVector)
					names <- c(names,paste("Band",AVIRISbands2_local[i],sep=""))
				}
				vars <- c(vars,i)
			}
		}else{
			for (i in 1:get('numbands2',envir=ripaEnv)){
				auxVector<-regionsList_local[[1]][[2]][[i]]
				if (is.null(AVIRISbands2_local)){
					assign(paste("Band",i,sep=""),auxVector)
					names <- c(names,paste("Band",i,sep=""))
				} else{
					assign(paste("Band",AVIRISbands2_local[i],sep=""),auxVector)
					names <- c(names,paste("Band",AVIRISbands2_local[i],sep=""))
				}
				vars <- c(vars,i)
			}
		}
		print(auxVector)
		for (i in 1:get('numbands2',envir=ripaEnv)){
			if (i==1){
				if (is.null(AVIRISbands2_local)) dataf <- data.frame(get(paste("Band",i,sep="")))
				else dataf <- data.frame(get(paste("Band",AVIRISbands2_local[i],sep="")))
			}
			else{
				if (is.null(AVIRISbands2_local)) dataf[,i] <- get(paste("Band",i,sep=""))
				else dataf[,i] <- get(paste("Band",AVIRISbands2_local[i],sep=""))
			}
		}
		
		names(dataf) <- names
		g <- ggobi(dataf)		
		d <- display(g[1], "Parallel Coordinates Display")
	}
	#
	
	# Function to show all regions
	showAllRegions <- function(){
		tab <- checkTab()

		visualBands2_local <- get('visualBands2',envir=ripaEnv)
		regionsList_local <- get('regionsList',envir=ripaEnv)

		if (tab=="1"){
			tkmessageBox(title="Error",message="Please, use this menu only with the second tab!",icon="error",type="ok")
			return()
		}
	
		if (length(regionsList_local)==0){
			tkmessageBox(title="Error",message="There are not regions",icon="error",type="ok")
			return()
		}
		for (i in 1:length(regionsList_local)){
			ttregion <- tktoplevel()
			tkwm.resizable(ttregion,0,0)
			tktitle(ttregion)<-paste("Region",i)
			
			if (get('numbands2',envir=ripaEnv)>3){
				auximg <- array(c(regionsList_local[[i]][[3]][,,visualBands2_local[1]],regionsList_local[[i]][[3]][,,visualBands2_local[2]],regionsList_local[[i]][[3]][,,visualBands2_local[3]]),c(nrow(regionsList_local[[i]][[3]]),ncol(regionsList_local[[i]][[3]]),get('numbands2',envir=ripaEnv)))
				cRegion <-tkrplot(ttregion, function() plot(imagematrix(auximg)),vscale=0.8,hscale=0.8)
			}
			if (get('numbands2',envir=ripaEnv)<=3){
				cRegion <-tkrplot(ttregion, function() plot(imagematrix(regionsList_local[[i]][[3]])),vscale=0.8,hscale=0.8)
			}
			tkpack(cRegion)
		}
	}
	#
	
	pca <- function(){
		tab <- checkTab()
		img1_local <- get('img1',envir=ripaEnv)
		if ((tab=="1" && is.null(img1_local)) || (tab=="2" && is.null(get('img3',envir=ripaEnv)))){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
		data <- matrix(nrow=nrow(img1_local[,,1])*ncol(img1_local[,,1]),ncol=get('numbands1',envir=ripaEnv))
		if (tab=="1"){
			for (i in 1:get('numbands1',envir=ripaEnv)){
				data[,i] <- as.vector(t(img1_local[,,i]))
			}
		}
		result <- princomp(data)
		for (i in 1:get('numbands1',envir=ripaEnv)){
			assign(paste("pca",i,sep=""),matrix(result[[6]][,i],nrow=nrow(img1_local[,,1]),ncol=ncol(img1_local[,,1]),byrow=T))
		}
		for (i in 1:get('numbands1',envir=ripaEnv)){
			assign(paste("ttpca",i,sep=""),tktoplevel())
			tkwm.title(get(paste("ttpca",i,sep="")),paste("PCA ",i,sep=""))
			image <-tkrplot(get(paste("ttpca",i,sep="")), function() plot(imagematrix(stretchImg(get(paste("pca",i,sep=""))))),vscale=0.7,hscale=0.7)
			tkpack(image)
		}
	}

	# Function to apply a zoom to an image
	zoomImg <- function(){
		tab <- checkTab()
		img1_local <- get('img1',envir=ripaEnv)
		img3_local <- get('img3',envir=ripaEnv)
		AVIRISbands1_local <- get('AVIRISbands1',envir = ripaEnv)
		AVIRISbands2_local <- get('AVIRISbands2',envir = ripaEnv)
		imageType1_local <- get('imageType1',envir = ripaEnv)		
		imageType2_local <- get('imageType2',envir = ripaEnv)

		bands_local <- get('bands',envir=ripaEnv)
		
		if ((tab=="1" && is.null(img1_local)) || (tab=="2" && is.null(img3_local))){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}

		if (tab=="1"){
			visualBands1_local <- get('visualBands1',envir=ripaEnv)
			if (imageType1_local=="lan"){
				auximg <- array(c(img1_local[,,visualBands1_local[1]],img1_local[,,visualBands1_local[2]],img1_local[,,visualBands1_local[3]]),c(nrow(img1_local),ncol(img1_local),get('numbands1',envir=ripaEnv)))
				plot(imagematrix(auximg))
				zoom_jpg_png_RGB(auximg)
			}

			if (imageType1_local=="jpgRGB" || imageType1_local=="pngRGB"){
				plot(imagematrix(img1_local))
				zoom_jpg_png_RGB(img1_local)
			}
			
			if (imageType1_local=="jpgGrey" || imageType1_local=="pngGrey"){
				plot(imagematrix(img1_local))
				zoom_jpg_png_Grey(img1_local)
			}
				
			if (!is.null(AVIRISbands1_local)){
				plot(bands_local[[1]],bands_local[[2]],bands_local[[3]])
				zoom(bands_local[[1]],bands_local[[2]],bands_local[[3]])
			}
		}

		if (tab=="2"){
			visualBands2_local <- get('visualBands2',envir=ripaEnv)
			if (imageType2_local=="lan"){
				auximg <- array(c(img3_local[,,visualBands2_local[1]],img3_local[,,visualBands2_local[2]],img3_local[,,visualBands2_local[3]]),c(nrow(img3_local),ncol(img3_local),get('numbands2',envir=ripaEnv)))
				plot(imagematrix(auximg))
				zoom_jpg_png_RGB(auximg)
			}

			if (imageType2_local=="jpgRGB" || imageType2_local=="pngRGB"){
				plot(imagematrix(img3_local))
				zoom_jpg_png_RGB(img3_local)
			}
			
			if (imageType2_local=="jpgGrey" || imageType2_local=="pngGrey"){
				plot(imagematrix(img3_local))
				zoom_jpg_png_Grey(img3_local)
			}
				
			if (!is.null(AVIRISbands2_local)){
				plot(bands_local[[1]],bands_local[[2]],bands_local[[3]])
				zoom(bands_local[[1]],bands_local[[2]],bands_local[[3]])
			}
		}

	}

	zoom_jpg_png_RGB <- function(image){
		pos <- locator(2)
		width <- dim(image)[1]
		if(pos$x[1]>pos$x[2]){
			aux <- pos$x[1]
			pos$x[1] <- pos$x[2]
			pos$x[2] <- aux  
		}
		if(pos$y[1]>pos$y[2]){
			auy <- pos$y[1]
			pos$y[1] <- pos$y[2]
			pos$y[2] <- auy  
		}
		
		pos$x <- round(pos$x)
		pos$y <- round(pos$y)
		rect(pos$x[1],pos$y[1],pos$x[2],pos$y[2],col="red",density=0.1)
		plot(imagematrix(image[(width-pos$y[2]):(width-pos$y[1]),pos$x[1]:pos$x[2],]))
	}

	zoom_jpg_png_Grey <- function(image){
		pos <- locator(2)
		width <- dim(image)[1]
		if(pos$x[1]>pos$x[2]){
			aux <- pos$x[1]
			pos$x[1] <- pos$x[2]
			pos$x[2] <- aux  
		}
		if(pos$y[1]>pos$y[2]){
			auy <- pos$y[1]
			pos$y[1] <- pos$y[2]
			pos$y[2] <- auy  
		}
		pos$x <- round(pos$x)
		pos$y <- round(pos$y)
		rect(pos$x[1],pos$y[1],pos$x[2],pos$y[2],col="red",density=0.1)
		plot(imagematrix(image[(width-pos$y[2]):(width-pos$y[1]),pos$x[1]:pos$x[2],]))
	}

	# Function to equalize the image histogram
	equalizeImg <- function(){
		tab <- checkTab()
		img1_local <- get('img1',envir=ripaEnv)
		if (tab=="2"){
			tkmessageBox(title="Error",message="Please, use this menu only with the first tab!",icon="error",type="ok")
			return()
		}
	
		if (is.null(img1_local)){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
		
		tkdestroy(Frame3)
		Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
		tkconfigure(Frame3,width=485,height=510)
		tkplace(Frame3,x=510,y=30)
		Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
		tkplace(Frame4,x=1,y=1)
		
		assign('img2',img1_local,envir=ripaEnv)
		img2_local <- get('img2',envir=ripaEnv)	
		#img2<<-img1_local
		for (i in 1:get('numbands1',envir=ripaEnv)){
			img2_local[,,i] <- equalize(img1_local[,,i])
		}
		image <-tkrplot(Frame4, function() plot(imagematrix(img2_local)),vscale=1.04,hscale=0.98)
		assign('img2',img2_local,envir=ripaEnv)
		tkpack(image)
	}
	#


################################ Comes from package rimage ########################################
##
## rimage: Image Processing Library for R
##
## $Header: /database/repository/rimage/R/rimage.R,v 1.8.2.8 2004/03/17 06:35:18 tomo Exp $
##
## Copyright (c) 2003 Nikon Systems Inc.
## For complete license terms see file LICENSE_rimage
##################################################################################################

	
	##
	## Grey level adjustment
	##

	thresholding <- function(img, mode="fixed", th=0.5) {
	  th.by.discrim <- function(img, L=255) {
	    img.int <- floor(L*img)
	    h <- hist(img.int, breaks=0:(L+1), plot=FALSE)$density
	    lv <- 0:L
	    u.img <- sum(lv * h) / sum(h)
	    s.img <- sum(h * (lv - u.img)^2)
	    Fs <- sapply(1:(L-1), function(k) {
	      w.0 <- sum(h[1:k])
	      w.1 <- sum(h[(k+1):L])
	      u.0 <- sum((1:k) * h[1:k]) / w.0
	      u.1 <- sum(((k+1):L) * h[(k+1):L]) / w.1
	      s.B <- w.0 * (u.0 - u.img)^2 + w.1 * (u.1 - u.img)^2
	      s.B / s.img
	    })
	    lv[rev(order(Fs, na.last = NA))[1]]/L
	  }

	  th <- switch(mode, fixed=th, da=th.by.discrim(img))
	  if (is.null(th)) stop("Either mode or threshold isn't correct.")
	  img[img < th] <- 0
	  img[img >= th] <- 1
	  img
	}

	## image takes up values 0..1    
	equalize <- function(img) {
		img <- (img-min(img))*255 / (max(img)-min(img)) ## normalize it to 0..255
		h <- dim(img)[1]
		w <- dim(img)[2]
		res <- matrix(.C("equalize",
		        as.double(img),
		        as.integer(w),
		        as.integer(h),
		        spec = double(w*h),
		        PACKAGE="ripa"
		        )$spec, nrow=h, ncol=w)
		imagematrix(res / 255)  ## map it to 0..1
	}
	
	clipping <- function(img, low=0, high=1) {
  		img[img < low] <- low
		img[img > high] <- high
  		img
	}

	normalize <- function(img) {
  		(img - min(img))/(max(img) - min(img))
	}

	rgb2grey <- function(img, coefs=c(0.30, 0.59, 0.11)) {
  		if (is.null(dim(img))) stop("image matrix isn't correct.")
  		if (length(dim(img))<3) stop("image matrix isn't rgb image.")
  		imagematrix(coefs[1] * img[,,1] + coefs[2] * img[,,2] + coefs[3] * img[,,3], type="grey")
	}

	##
	## Edge Detection Filters
	##

	sobel.h <- function(img) {
	  w <- dim(img)[2]
	  h <- dim(img)[1]
	  imagematrix(abs(matrix(.C("sobel_h",
		                    as.double(img), as.integer(w), as.integer(h),
		                    eimg = double(w * h),
		                    PACKAGE="ripa")$eimg,
		                 nrow=h, ncol=w)), noclipping=TRUE)
	}

	sobel.v <- function(img) {
	  w <- dim(img)[2]
	  h <- dim(img)[1]
	  imagematrix(abs(matrix(.C("sobel_v",
		                    as.double(img), as.integer(w), as.integer(h),
		                    eimg = double(w * h),
		                    PACKAGE="ripa")$eimg,
		                 nrow=h, ncol=w)), noclipping=TRUE)
	}

	sobel <- function(img) {
	  h.img <- sobel.h(img)
	  v.img <- sobel.v(img)
	  (h.img + v.img)/2
	}

	laplacian <- function(img) {
	  w <- dim(img)[2]
	  h <- dim(img)[1]
	  l.img <- imagematrix(matrix(.C("laplacian",
		                         as.double(img), as.integer(w), as.integer(h),
		                         eimg = double(w * h),
		                         PACKAGE="ripa")$eimg,
		                      nrow=h, ncol=w),
		               noclipping=TRUE)
	}

	##
	## Rank Filters
	##

	meanImg <- function(img) {
	  expand.h <- cbind(img[,1], img, img[,dim(img)[2]])
	  ex.img <- rbind(expand.h[1,], expand.h, expand.h[dim(img)[1],])
	  w <- dim(ex.img)[2]
	  h <- dim(ex.img)[1]
	  f.img <- matrix(.C("meanfilter",
		             as.double(ex.img), as.integer(w), as.integer(h),
		             eimg = double(w * h),
		             PACKAGE="ripa")$eimg,
		          nrow=h, ncol=w)
	  imagematrix(f.img[2:(dim(f.img)[1]-1),2:(dim(f.img)[2]-1)])
	}

	minImg <- function(img) {
	  expand.h <- cbind(img[,1], img, img[,dim(img)[2]])
	  ex.img <- rbind(expand.h[1,], expand.h, expand.h[dim(img)[1],])
	  w <- dim(ex.img)[2]
	  h <- dim(ex.img)[1]
	  f.img <- matrix(.C("minfilter",
		             as.double(ex.img), as.integer(w), as.integer(h),
		             eimg = double(w * h),
		             PACKAGE="ripa")$eimg,
		          nrow=h, ncol=w)
	  imagematrix(f.img[2:(dim(f.img)[1]-1),2:(dim(f.img)[2]-1)])
	}

	maxImg <- function(img) {
	  expand.h <- cbind(img[,1], img, img[,dim(img)[2]])
	  ex.img <- rbind(expand.h[1,], expand.h, expand.h[dim(img)[1],])
	  w <- dim(ex.img)[2]
	  h <- dim(ex.img)[1]
	  f.img <- matrix(.C("maxfilter",
		             as.double(ex.img), as.integer(w), as.integer(h),
		             eimg = double(w * h),
		             PACKAGE="ripa")$eimg,
		          nrow=h, ncol=w)
	  imagematrix(f.img[2:(dim(f.img)[1]-1),2:(dim(f.img)[2]-1)])
	}
	#################################

	# Function to show the chosen region
	region <- function(){


		chooseRegion <- function(){

			img3_local <- get('img3',envir=ripaEnv)
			plotFunction2 <- function(){
				img3_local <- get('img3',envir=ripaEnv)
				visualBands2_local <- get('visualBands2',envir=ripaEnv)
				params <- par(bg="white")
				if (get('numbands2',envir=ripaEnv)>3){
					auximg <- array(c(img3_local[,,visualBands2_local[1]],img3_local[,,visualBands2_local[2]],img3_local[,,visualBands2_local[3]]),c(nrow(img3_local),ncol(img3_local),get('numbands2',envir=ripaEnv)))
					plot(imagematrix(auximg))
				}else{
					plot(imagematrix(img3_local))
				}
				parPlotSize2 <<- par("plt")
				usrCoords2   <<- par("usr")
				par(params)
			}
		
			OnLeftClick2 <- function(x,y){
				img3_local <- get('img3',envir=ripaEnv)
				imgtmp2_local <- get('imgtmp2',envir=ripaEnv)
				xClick <- x
				yClick <- y
				width  <- as.numeric(tclvalue(tkwinfo("reqwidth",imgtmp2_local)))
				height <- as.numeric(tclvalue(tkwinfo("reqheight",imgtmp2_local)))
				
				xMin <- parPlotSize2[1] * width
				xMax <- parPlotSize2[2] * width
				yMin <- parPlotSize2[3] * height
				yMax <- parPlotSize2[4] * height
			
				rangeX <- usrCoords2[2] - usrCoords2[1]
				rangeY <- usrCoords2[4] - usrCoords2[3]
			
				imgXcoords <- (get('xCoords2',envir=ripaEnv)-usrCoords2[1])*(xMax-xMin)/rangeX + xMin
				imgYcoords <- (get('yCoords2',envir=ripaEnv)-usrCoords2[3])*(yMax-yMin)/rangeY + yMin
				
				xClick <- as.numeric(xClick)+0.5
				yClick <- as.numeric(yClick)+0.5
				yClick <- height - yClick
				
				xPlotCoord <- usrCoords2[1]+(xClick-xMin)*rangeX/(xMax-xMin)
				yPlotCoord <- usrCoords2[3]+(yClick-yMin)*rangeY/(yMax-yMin)
				
				a <- round(xPlotCoord)
				b <- round(nrow(img3_local) - yPlotCoord)
				
				if (a<0 || b<0 || a>ncol(img3_local) || b>nrow(img3_local)) tkmessageBox(title="Error",message="Please, click inside the image!",icon="error",type="ok")
				else{
					assign('regionPointsX',c(get('regionPointsX',envir=ripaEnv),a),envir=ripaEnv)
					assign('regionPointsY',c(get('regionPointsY',envir=ripaEnv),b),envir=ripaEnv)
				}
			}
			
			findRegion <- function(){
				img3_local <- get('img3',envir=ripaEnv)
				visualBands2_local <- get('visualBands2',envir=ripaEnv)
				regionValues <- list()
				ttregion <- tktoplevel()
				tkwm.geometry(ttregion,"+0+0")
				tkwm.resizable(ttregion,0,0)
				tktitle(ttregion)<-paste("Region",length(get('regionsList',envir=ripaEnv))+1)
				n = as.integer(length(get('regionPointsX',envir=ripaEnv)))
				#out = .C("grahamMain",n,as.integer(get('regionPointsX',envir=ripaEnv)),as.integer(get('regionPointsY',envir=ripaEnv)),xvectorOut = as.integer(rep(0,n)),yvectorOut = as.integer(rep(0,n)),PACKAGE="ripa")
				aux=NULL
				enter=F
				regionPointsX_local <- get('regionPointsX',envir=ripaEnv)
				regionPointsY_local <- get('regionPointsY',envir=ripaEnv)
				indices = chull(regionPointsX_local,regionPointsY_local)
				x = regionPointsX_local[indices]
				y = regionPointsY_local[indices]
				#x = out$xvectorOut
				#y = out$yvectorOut
				for (i in 1:length(x)){
					if (x[i]==0){
						enter=T
						for (j in 1:(i-1)){
							aux = c(aux,x[j])
						}
					}
					if (enter==T) break
				}
				if (enter==T){
					x = aux
					aux=NULL
					enter=F
				}
				for (i in 1:length(y)){
					if (y[i]==0){
						enter=T
						for (j in 1:(i-1)){
							aux = c(aux,y[j])
						}
					}
					if (enter==T) break
				}
				if (enter==T) y = aux
				
				chosenRegion <- img3_local
				for (i in 1:get('numbands2',envir=ripaEnv)){
					out <- .C("inpolyMain",as.integer(length(x)),as.integer(x),as.integer(y),outregion=as.double(as.vector(t(as.matrix(img3_local[,,i])))),as.integer(nrow(img3_local)),as.integer(ncol(img3_local)),rValues = as.double(rep(0,nrow(img3_local)*ncol(img3_local))),regionLength=as.integer(0),PACKAGE="ripa")
					chosenRegion[,,i] <- matrix(out$outregion,ncol=ncol(img3_local),byrow=T)
					regionValues[[i]] <- out$rValues[1:(out$regionLength-1)]
				}
				
				regionsList_local <- get('regionsList',envir=ripaEnv)

				regionsList_local[[length(regionsList_local)+1]] <- list(Type=get('regionType',envir=ripaEnv),Values=regionValues,Visual=chosenRegion)
					
				names(regionsList_local)[[length(regionsList_local)]] <-paste("Region",length(regionsList_local),sep=" ")
				length(names(regionsList_local)) <- length(regionsList_local)
				assign('regionsList',regionsList_local,envir=ripaEnv)
				
				if (get('numbands2',envir=ripaEnv)>3){
					auximg <- array(c(chosenRegion[,,visualBands2_local[1]],chosenRegion[,,visualBands2_local[2]],chosenRegion[,,visualBands2_local[3]]),c(nrow(chosenRegion),ncol(chosenRegion),get('numbands2',envir=ripaEnv)))
					cRegion <-tkrplot(ttregion, function() plot(imagematrix(auximg,noclipping=TRUE)),vscale=0.8,hscale=0.8)
					tkpack(cRegion)
				}else{
					cRegion <-tkrplot(ttregion, function() plot(imagematrix(chosenRegion)),vscale=0.8,hscale=0.8)
					tkpack(cRegion)
				}
				
				showRegionStatistics(regionValues)
			}
			
			# Auxiliar function to print the statistical summary of the regions
			showRegionStatistics <- function(region){
				AVIRISbands2_local <- get('AVIRISbands2',envir=ripaEnv)
				# Function to print the statistical summary of the regions
				printStatistcs <- function(nR,minR,maxR,meanR,medianR,deviationR,mDeviationR,kurtosisR,skewnessR){
					tkconfigure(statisticsTxt,state="normal")
					tkinsert(statisticsTxt,"end",paste("N = ",nR,"\n",sep=""))
					tkinsert(statisticsTxt,"end",paste("Min = ",minR,"\n",sep=""))
					tkinsert(statisticsTxt,"end",paste("Max = ",maxR,"\n",sep=""))
					tkinsert(statisticsTxt,"end",paste("Mean = ",meanR,"\n",sep=""))
					tkinsert(statisticsTxt,"end",paste("Median = ",medianR,"\n",sep=""))
					tkinsert(statisticsTxt,"end",paste("Deviation = ",deviationR,"\n",sep=""))
					tkinsert(statisticsTxt,"end",paste("Median Deviation = ",mDeviationR,"\n",sep=""))
					tkinsert(statisticsTxt,"end",paste("Kurtosis = ",kurtosisR,"\n",sep=""))
					tkinsert(statisticsTxt,"end",paste("Skewness = ",skewnessR,"\n\n",sep=""))
					tkconfigure(statisticsTxt, state="disabled")
				}
				#
				
				for (i in 1:get('numbands2',envir=ripaEnv)){
					nR <- length(region[[i]])
					minR <- min(region[[i]])
					maxR <- max(region[[i]])
					meanR <- mean(region[[i]])
					medianR <- median(region[[i]])
					deviationR <- sd(as.vector(region[[i]]))
					mDeviationR <- mad(region[[i]])
					kurtosisR <- kurtosis(as.vector(region[[i]]))
					skewnessR <- skewness(as.vector(region[[i]]))
					tkconfigure(statisticsTxt,state="normal")
					if (i==1) tkinsert(statisticsTxt,"end",paste(names(get('regionsList',envir=ripaEnv))[length(get('regionsList',envir=ripaEnv))],"\n"))
					if (is.null(AVIRISbands2_local)) tkinsert(statisticsTxt,"end",paste("Band",i,"\n"))
					else tkinsert(statisticsTxt,"end",paste("Band",AVIRISbands2_local[i],"\n"))
					printStatistcs(nR,minR,maxR,meanR,medianR,deviationR,mDeviationR,kurtosisR,skewnessR)
				}
			}
			#
			
			OnRightClick <- function(){
				tkdestroy(get('ttpoints',envir=ripaEnv))
				findRegion()
			}
			
			rectangleRegion <- function(){
				regionPointsX_local <- get('regionPointsX',envir=ripaEnv)
				regionPointsY_local <- get('regionPointsY',envir=ripaEnv)
				img3_local <- get('img3',envir=ripaEnv)
				visualBands2_local <- get('visualBands2',envir=ripaEnv)
				regionValues <- list()
				tkdestroy(get('ttpoints',envir=ripaEnv))
				
				ttregion <- tktoplevel()
				tkwm.geometry(ttregion,"+0+0")
				tkwm.resizable(ttregion,0,0)
				tktitle(ttregion)<-paste("Region",length(get('regionsList',envir=ripaEnv))+1)
				
				a<-1
				b<-regionPointsY_local[1]
				c<-regionPointsY_local[2]
				d<-nrow(img3_local)
				e<-1
				f<-regionPointsX_local[1]
				g<-regionPointsX_local[2]
				h<-ncol(img3_local)
				
				chosenRegion<-img3_local[-c(a:b,c:d),-c(e:f,g:h),]
				chosenRegion <- array(chosenRegion,dim=c(nrow(chosenRegion),ncol(chosenRegion),dim(img3_local)[3]))
				for (i in 1:get('numbands2',envir=ripaEnv)){
					regionValues[[i]] <- as.vector(chosenRegion[,,i])
				}

				regionsList_local <- get('regionsList',envir=ripaEnv)

				regionsList_local[[length(regionsList_local)+1]] <- list(Type=get('regionType',envir=ripaEnv),Values=regionValues,Visual=chosenRegion)
				
				names(regionsList_local)[[length(regionsList_local)]] <-paste("Region",length(regionsList_local),sep=" ")
				length(names(regionsList_local)) <- length(regionsList_local)
				assign('regionsList',regionsList_local,envir=ripaEnv)
				
				if (get('numbands2',envir=ripaEnv)>3){
					auximg <- array(c(chosenRegion[,,visualBands2_local[1]],chosenRegion[,,visualBands2_local[2]],chosenRegion[,,visualBands2_local[3]]),c(nrow(chosenRegion),ncol(chosenRegion),get('numbands2',envir=ripaEnv)))
					cRegion <-tkrplot(ttregion, function() plot(imagematrix(auximg,noclipping=TRUE)),vscale=0.8,hscale=0.8)
					tkpack(cRegion)
				}else{
					cRegion <-tkrplot(ttregion, function() plot(imagematrix(chosenRegion)),vscale=0.8,hscale=0.8)
					tkpack(cRegion)
				}
				
				showRegionStatistics(regionValues)
				
				
			}
			
			OnLeftClick3 <- function(x,y){
				img3_local <- get('img3',envir=ripaEnv)
				imgtmp2_local <- get('imgtmp2',envir=ripaEnv)
				xClick <- x
				yClick <- y
				width  <- as.numeric(tclvalue(tkwinfo("reqwidth",imgtmp2_local)))
				height <- as.numeric(tclvalue(tkwinfo("reqheight",imgtmp2_local)))
				
				xMin <- parPlotSize2[1] * width
				xMax <- parPlotSize2[2] * width
				yMin <- parPlotSize2[3] * height
				yMax <- parPlotSize2[4] * height
				
				rangeX <- usrCoords2[2] - usrCoords2[1]
				rangeY <- usrCoords2[4] - usrCoords2[3]
				
				imgXcoords <- (get('xCoords2',envir=ripaEnv)-usrCoords2[1])*(xMax-xMin)/rangeX + xMin
				imgYcoords <- (get('yCoords2',envir=ripaEnv)-usrCoords2[3])*(yMax-yMin)/rangeY + yMin
				
				xClick <- as.numeric(xClick)+0.5
				yClick <- as.numeric(yClick)+0.5
				yClick <- height - yClick
				
				xPlotCoord <- usrCoords2[1]+(xClick-xMin)*rangeX/(xMax-xMin)
				yPlotCoord <- usrCoords2[3]+(yClick-yMin)*rangeY/(yMax-yMin)
				
				a <- round(xPlotCoord)
				b <- round(nrow(img3_local) - yPlotCoord)
				
				if (a<0 || b<0 || a>ncol(img3_local) || b>nrow(img3_local)) tkmessageBox(title="Error",message="Please, click inside the image!",icon="error",type="ok")
				else{
					assign('regionPointsX',c(get('regionPointsX',envir=ripaEnv),a),envir=ripaEnv)
					assign('regionPointsY',c(get('regionPointsY',envir=ripaEnv),b),envir=ripaEnv)
				}
				if (length(get('regionPointsX',envir=ripaEnv))==2) rectangleRegion()
			}
			
			circleRegion <- function(){
				img3_local <- get('img3',envir=ripaEnv)
				regionPointsX_local <- get('regionPointsX',envir=ripaEnv)
				regionPointsY_local <- get('regionPointsY',envir=ripaEnv)
				visualBands2_local <- get('visualBands2',envir=ripaEnv)
				tkdestroy(get('ttpoints',envir=ripaEnv))
				
				ttregion <- tktoplevel()
				tkwm.geometry(ttregion,"+0+0")
				tkwm.resizable(ttregion,0,0)
				tktitle(ttregion)<-paste("Region",length(get('regionsList',envir=ripaEnv))+1)
				
				x1<-as.integer(regionPointsX_local[1])
				x2<-as.integer(regionPointsX_local[2])
				y1<-as.integer(regionPointsY_local[1])
				y2<-as.integer(regionPointsY_local[2])
				
				regionValues <- list()
				chosenRegion<-img3_local
				for (i in 1:get('numbands2',envir=ripaEnv)){
					out <- .C("inCircle",x1,x2,y1,y2,outregion=as.double(as.vector(t(as.matrix(img3_local[,,i])))),as.integer(nrow(img3_local)),as.integer(ncol(img3_local)),rValues = as.double(rep(0,nrow(img3_local)*ncol(img3_local))),regionLength=as.integer(0),PACKAGE="ripa")
					regionValues[[i]]<-out$rValues[1:(out$regionLength-1)]
					chosenRegion[,,i]<-matrix(out$outregion,nrow=nrow(img3_local),byrow=T)
				}

				regionsList_local <- get('regionsList',envir=ripaEnv)
				regionsList_local[[length(regionsList_local)+1]] <- list(Type=get('regionType',envir=ripaEnv),Values=regionValues,Visual=chosenRegion)
						
				names(regionsList_local)[[length(regionsList_local)]] <-paste("Region",length(regionsList_local),sep=" ")
				length(names(regionsList_local)) <- length(regionsList_local)
				assign('regionsList',regionsList_local,envir=ripaEnv)
				
				if (get('numbands2',envir=ripaEnv)>3){
					auximg <- array(c(chosenRegion[,,visualBands2_local[1]],chosenRegion[,,visualBands2_local[2]],chosenRegion[,,visualBands2_local[3]]),c(nrow(chosenRegion),ncol(chosenRegion),get('numbands2',envir=ripaEnv)))
					cRegion <-tkrplot(ttregion, function() plot(imagematrix(auximg,noclipping=TRUE)),vscale=0.8,hscale=0.81)
					tkpack(cRegion)
				}else{
					cRegion <-tkrplot(ttregion, function() plot(imagematrix(chosenRegion)),vscale=0.8,hscale=0.8)
					tkpack(cRegion)
				}
				
				showRegionStatistics(regionValues)
				
			}
			
			OnLeftClick4 <- function(x,y){
				img3_local <- get('img3',envir=ripaEnv)
				imgtmp2_local <- get('imgtmp2',envir=ripaEnv)
				xClick <- x
				yClick <- y
				width  <- as.numeric(tclvalue(tkwinfo("reqwidth",imgtmp2_local)))
				height <- as.numeric(tclvalue(tkwinfo("reqheight",imgtmp2_local)))
				
				xMin <- parPlotSize2[1] * width
				xMax <- parPlotSize2[2] * width
				yMin <- parPlotSize2[3] * height
				yMax <- parPlotSize2[4] * height
				
				rangeX <- usrCoords2[2] - usrCoords2[1]
				rangeY <- usrCoords2[4] - usrCoords2[3]
				
				imgXcoords <- (get('xCoords2',envir=ripaEnv)-usrCoords2[1])*(xMax-xMin)/rangeX + xMin
				imgYcoords <- (get('yCoords2',envir=ripaEnv)-usrCoords2[3])*(yMax-yMin)/rangeY + yMin
				
				xClick <- as.numeric(xClick)+0.5
				yClick <- as.numeric(yClick)+0.5
				yClick <- height - yClick
				
				xPlotCoord <- usrCoords2[1]+(xClick-xMin)*rangeX/(xMax-xMin)
				yPlotCoord <- usrCoords2[3]+(yClick-yMin)*rangeY/(yMax-yMin)
				
				a <- round(xPlotCoord)
				b <- round(nrow(img3_local) - yPlotCoord)
				
				if (a<0 || b<0 || a>ncol(img3_local) || b>nrow(img3_local)) tkmessageBox(title="Error",message="Please, click inside the image!",icon="error",type="ok")
				else{
					assign('regionPointsX',c(get('regionPointsX',envir=ripaEnv),a),envir=ripaEnv)
					assign('regionPointsY',c(get('regionPointsY',envir=ripaEnv),b),envir=ripaEnv)
				}
				if (length(get('regionPointsX',envir=ripaEnv))==2) circleRegion()
			}
			
			
			assign('regionPointsX',NULL,envir=ripaEnv)
			assign('regionPointsY',NULL,envir=ripaEnv)
			
			assign('ttpoints',tktoplevel(),envir=ripaEnv)
				
			assign('xCoords2',(1:ncol(img3_local)),envir=ripaEnv)
			assign('yCoords2',(1:nrow(img3_local)),envir=ripaEnv)
			parPlotSize2 <- c()
			usrCoords2 <- c()
			
			assign('imgtmp2',tkrplot(get('ttpoints',envir=ripaEnv),fun=plotFunction2,hscale=1.5,vscale=1.5),envir=ripaEnv)
			imgtmp2_local <- get('imgtmp2',envir=ripaEnv)
			
			tkgrid(imgtmp2_local)
			
			if (get('regionType',envir=ripaEnv) == "Free"){
				tkbind(imgtmp2_local, "<Button-1>",OnLeftClick2)
				tkbind(imgtmp2_local, "<Button-3>",OnRightClick)
				tkwm.title(get('ttpoints',envir=ripaEnv),"Choose the points clicking on the left button of the mouse")
			}
			if (get('regionType',envir=ripaEnv) == "Rectangle"){
				tkbind(imgtmp2_local, "<Button-1>",OnLeftClick3)
				tkwm.title(get('ttpoints',envir=ripaEnv),"Choose the top left and the bottom right points")
			}
			if (get('regionType',envir=ripaEnv) == "Circle"){
				tkbind(imgtmp2_local, "<Button-1>",OnLeftClick4)
				tkwm.title(get('ttpoints',envir=ripaEnv),"Choose the center and the radius of the circle")
			}
			tkconfigure(imgtmp2_local,cursor="hand2")
			assign('imgtmp2',imgtmp2_local,envir=ripaEnv)
			
			
		}
		
		ttkindregions<-tktoplevel()
		tkwm.geometry(ttkindregions,"+400+250")
		tkwm.resizable(ttkindregions,0,0)
		tktitle(ttkindregions)<-"Regions"
		scr <- tkscrollbar(ttkindregions, repeatinterval=5, command=function(...)tkyview(tl,...))
		tl<-tklistbox(ttkindregions,height=6,width=30,selectmode="single",yscrollcommand=function(...)tkset(scr,...),background="white")
		tkgrid(tklabel(ttkindregions,text="Choose the kind of region"))
		tkgrid(tl,scr)
		tkgrid.configure(scr,rowspan=6,sticky="nsw")
		
		tkinsert(tl,"end","Rectangle")
		tkinsert(tl,"end","Circle")
		tkinsert(tl,"end","Free")
		tkselection.set(tl,0)
		
		OnOK <- function(){
			kindChoice <- as.numeric(tkcurselection(tl))+1
			if (kindChoice==1){
				tkdestroy(ttkindregions)
				assign('regionType',"Rectangle",envir=ripaEnv)
				chooseRegion()
			}
			
			if (kindChoice==2){
				tkdestroy(ttkindregions)
				assign('regionType',"Circle",envir=ripaEnv)
				chooseRegion()
			}
			
			if (kindChoice==3){
				tkdestroy(ttkindregions)
				assign('regionType',"Free",envir=ripaEnv)
				chooseRegion()
			}
		}
		
		OK.but <-tkbutton(ttkindregions,text=   " OK ",command=OnOK)
		tkgrid(OK.but)
		tkfocus(ttkindregions)
	}

	# Linear contrast stretch function
	stretch <- function(){
		tab <- checkTab()

		if (tab=="2"){
			tkmessageBox(title="Error",message="Please, use this menu only with the first tab!",icon="error",type="ok")
			return()
		}
		img1_local <- get('img1',envir=ripaEnv)
		if (is.null(img1_local)){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
		
		tkdestroy(Frame3)
		Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
		tkconfigure(Frame3,width=485,height=510)
		tkplace(Frame3,x=510,y=30)
		Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
		tkplace(Frame4,x=1,y=1)
		
		assign('img2',stretchImg(img1_local),envir=ripaEnv)
		img2_local <- get('img2',envir=ripaEnv)
		visualBands1_local <- get('visualBands1',envir=ripaEnv)


		if (get('numbands1',envir=ripaEnv)==1){
			image <-tkrplot(Frame4, function() plot(imagematrix(img2_local)),vscale=1.04,hscale=0.98)
			tkpack(image)
		}else{
			auximg <- array(c(img2_local[,,visualBands1_local[1]],img2_local[,,visualBands1_local[2]],img2_local[,,visualBands1_local[3]]),c(nrow(img2_local),ncol(img2_local),get('numbands1',envir=ripaEnv)))
			image <-tkrplot(Frame4, function() plot(imagematrix(auximg)),vscale=1.04,hscale=0.98)
			tkpack(image)
		}
	}
	
	# Function to choose filters
	filters <- function(opt){
		tab <- checkTab()

		if (tab=="2"){
			tkmessageBox(title="Error",message="Please, use this menu only with the first tab!",icon="error",type="ok")
			return()
		}
		img1_local <- get('img1',envir=ripaEnv)
		if (is.null(img1_local)){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
		
		if (opt==1){
			assign('img2',img1_local,envir=ripaEnv)
			img2_local <- get('img2',envir=ripaEnv)
			for (i in 1:get('numbands1',envir=ripaEnv)){
				img2_local[,,i] <- normalize(sobel(img1_local[,,i]))
			}
			assign('img2',img2_local,envir=ripaEnv)
		}
		else if (opt==2){
			assign('img2',img1_local,envir=ripaEnv)
			img2_local <- get('img2',envir=ripaEnv)		
			for (i in 1:get('numbands1',envir=ripaEnv)){
				img2_local[,,i] <- normalize(laplacian(img1_local[,,i]))
			}
			assign('img2',img2_local,envir=ripaEnv)
		}
		else if (opt==3){
			assign('img2',img1_local,envir=ripaEnv)
			img2_local <- get('img2',envir=ripaEnv)
			for (i in 1:get('numbands1',envir=ripaEnv)){
				img2_local[,,i] <- normalize(meanImg(matrix(img1_local[,,1],nrow=nrow(img1_local))))
			}
			assign('img2',img2_local,envir=ripaEnv)
		}
		else if (opt==4){
			assign('img2',img1_local,envir=ripaEnv)
			img2_local <- get('img2',envir=ripaEnv)
			for (i in 1:get('numbands1',envir=ripaEnv)){
				img2_local[,,i] <- normalize(highpass(img1_local[,,i]))
			}
			assign('img2',img2_local,envir=ripaEnv)
		}
		else if (opt==5){
			assign('img2',img1_local,envir=ripaEnv)
			img2_local <- get('img2',envir=ripaEnv)
			for (i in 1:get('numbands1',envir=ripaEnv)){
				img2_local[,,i] <- normalize(lowpass(img1_local[,,i]))
			}
			assign('img2',img2_local,envir=ripaEnv)
		}
		else if (opt==6){
			assign('img2',img1_local,envir=ripaEnv)
			img2_local <- get('img2',envir=ripaEnv)
			for (i in 1:get('numbands1',envir=ripaEnv)){
				img2_local[,,i] <- normalize(minImg(matrix(img1_local[,,1],nrow=nrow(img1_local))))
			}
			assign('img2',img2_local,envir=ripaEnv)
		}
		else if (opt==7){
			assign('img2',img1_local,envir=ripaEnv)
			img2_local <- get('img2',envir=ripaEnv)
			for (i in 1:get('numbands1',envir=ripaEnv)){
				img2_local[,,i] <- normalize(maxImg(matrix(img1_local[,,1],nrow=nrow(img1_local))))
			}
			assign('img2',img2_local,envir=ripaEnv)
		}
		else if (opt==8){
			ReturnVal <- modalDialog("Mask","Enter the mask length (3, 5, 7...)","")
			if (ReturnVal=="ID_CANCEL")
				return()
			lenMask <- as.integer(ReturnVal)
			assign('img2',img1_local,envir=ripaEnv)			
			assign('img2',medianImg(img1_local,lenMask),envir=ripaEnv)			
		}
		
		tkdestroy(Frame3)
		Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
		tkconfigure(Frame3,width=485,height=510)
		tkplace(Frame3,x=510,y=30)
		Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
		tkplace(Frame4,x=1,y=1)
		
		image <-tkrplot(Frame4, function() plot(imagematrix(get('img2',envir=ripaEnv))),vscale=1.04,hscale=0.98)
		tkpack(image)
	}
	#
	
	# Function to choose segmentation techniques
	segmentation <- function(opt){
		# Function to apply thresholding to an image

		img2_local <- get('img2',envir=ripaEnv)

		thresholdingImg <- function(img,value){
			return(thresholding(imagematrix(img),th=threshValue))
		}
		#
	
		tab <- checkTab()

		if (tab=="2"){
			tkmessageBox(title="Error",message="Please, use this menu only with the first tab!",icon="error",type="ok")
			return()
		}
		img1_local <- get('img1',envir=ripaEnv)
		if (is.null(img1_local)){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
		
		if (opt==1){       #thresholding
			ReturnVal <- modalDialog("Threshold","Enter the threshold value (0.0 - 1.0)","")
			if (ReturnVal=="ID_CANCEL")
				return()
			threshValue <- as.double(ReturnVal)
			assign('img2',img1_local,envir=ripaEnv)
			img2_local <- get('img2',envir=ripaEnv)
			for (i in 1:get('numbands1',envir=ripaEnv)){
				img2_local[,,i] <- thresholdingImg(img1_local[,,i],threshValue)
			}
			assign('img2',imagematrix(img2_local),envir=ripaEnv)
		}
		
		tkdestroy(Frame3)
		Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
		tkconfigure(Frame3,width=485,height=510)
		tkplace(Frame3,x=510,y=30)
		Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
		tkplace(Frame4,x=1,y=1)
		
		image <-tkrplot(Frame4, function() plot(get('img2',envir=ripaEnv)),vscale=1.04,hscale=0.98)
		tkpack(image)
	}
	#
	
	# Function to calculate the image index
	qualityIndex <- function(){
		tab <- checkTab()

		if (tab=="2"){
			tkmessageBox(title="Error",message="Please, use this menu only with the first tab!",icon="error",type="ok")
			return()
		}
		img1_local <- get('img1',envir=ripaEnv)
		if (is.null(img1_local)){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
	
		image_x <- as.double(as.vector(t(as.matrix(get('img1',envir=ripaEnv)))))
		image_y <- as.double(as.vector(t(as.matrix(get('img2',envir=ripaEnv)))))
		ncol <- as.integer(ncol(img1_local))
		nrow <- as.integer(nrow(img1_local))
		ReturnVal <- modalDialog("Window","Enter the window length (common value: 8)","")
		if (ReturnVal=="ID_CANCEL")
			return()
		tamWin <- as.integer(ReturnVal)
		out <- as.double(0.0)
		res <- .C("quality",image_x,image_y,ncol,nrow,tamWin,out,PACKAGE="ripa")
		index <- res[[6]]
		message <- paste("The quality index is:",index)
		tkmessageBox(title="Quality Index Result",message=message)
	}
	#
	
	# Function to save the transformed image
	saveImg <- function(){
		tab <- checkTab()
		imageType1_local <- get('imageType1',envir=ripaEnv)			

		if (tab=="2"){
			tkmessageBox(title="Error",message="Please, use this menu only with the first tab!",icon="error",type="ok")
			return()
		}

		if (tab=="1"){
			if (imageType1_local=="lan"){
				saveLan()
			}

			if (imageType1_local=="jpgGrey" || imageType1_local=="jpgRGB"){
				saveJPG()
			}

			if (imageType1_local=="pngGrey" || imageType1_local=="pngRGB"){
				savePNG()
			}
		}
		
	}
	#
	
	# Function to read AVIRIS images
	openAVIRIS <- function(){
		tab <- checkTab()

		if (tab=="1"){
			assign('AVIRISbands1',NULL,envir=ripaEnv)
			.Tcl(paste(slider, "set",0))
			.Tcl(paste(slider2, "set",1))
	
			tkdestroy(Frame3)
			Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
			tkconfigure(Frame3,width=485,height=510)
			tkplace(Frame3,x=510,y=30)
			
			Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
			tkplace(Frame4,x=1,y=1)
			
			tkdestroy(Frame1)
			Frame1 <- tkframe(lb1,relief="groove",borderwidth=2)
			tkconfigure(Frame1,width=485,height=510)
			tkplace(Frame1,x=10,y=30)
			
			Frame2 <- tkframe(Frame1,relief="groove",borderwidth=0)
			tkplace(Frame2,x=1,y=1)
		}
		
		if (tab=="2"){
			assign('AVIRISbands2',NULL,envir=ripaEnv)
			assign('regionsList',list(),envir=ripaEnv)
			assign('regionsListAux',list(),envir=ripaEnv)
			
			tkconfigure(statisticsTxt,state="normal")
			tkdelete(statisticsTxt,"1.0","end")
			tkconfigure(statisticsTxt,state="disable")
			
			tkdestroy(Frame5)
		
			Frame5 <- tkframe(lb2,relief="groove",borderwidth=2)
			tkconfigure(Frame5,width=480,height=510)
			tkplace(Frame5,x=10,y=30)
	
			Frame6 <- tkframe(Frame5,relief="groove",borderwidth=0)
			tkplace(Frame6,x=1,y=1)	
		}
		
		## Bands samples
		assign('bandsSet',c(rep(FALSE,224)),envir=ripaEnv)
 		fileName <- tkgetOpenFile(filetypes="{{AVIRIS Files} {.a.log}}")
 		if (!nchar(tclvalue(fileName))){
 			tkmessageBox(message="No file was selected!")
 			return()
 		}
 		fileName <- tclvalue(fileName)
		fileName <- strsplit(as.character(fileName),".log")

		ttsamples <- tktoplevel()
		tktitle(ttsamples) <- "Bands Samples"
		sw <- tkwidget(ttsamples,"ScrolledWindow",relief="sunken",borderwidth=2)
		sf <- tkwidget(sw,"ScrollableFrame")
		tcl(sw,"setwidget",sf)
		subfID <- tclvalue(tcl(sf,"getframe"))
		lab <- tcl("label",paste(subfID,".lab",sep=""),text="Select the bands you want:")
		tkpack(sw,fill="both",expand="yes")
		tkgrid(lab,row=1,column=1)
		sampleList <<- list()

		OnClick <- function(W){
			bandsSet_local <- get('bandsSet',envir=ripaEnv)
			if (tab=="1") {
				bands_selection <- get('AVIRISbands1',envir=ripaEnv)
			} else {
				bands_selection <- get('AVIRISbands2',envir=ripaEnv)
			}
			band <- as.integer(strsplit(W,paste(subfID,".lab",sep=""))[[1]][2])
			if (bandsSet_local[band]==TRUE){
				tkmessageBox(title="Band number",message=paste("Band ",band," deselected!",sep=""),type="ok")
				bandsSet_local[band] <- FALSE
				bands_selection = setdiff(bands_selection, band)
			}else{
				tkmessageBox(title="Band number",message=paste("Band ",band," selected!",sep=""),type="ok")
				bandsSet_local[band] <- TRUE
				bands_selection = c(bands_selection,band)
			}
			assign('bandsSet',bandsSet_local,envir=ripaEnv)
			if (tab=="1") {
				assign('AVIRISbands1',bands_selection,envir=ripaEnv)
			} else {
				assign('AVIRISbands2',bands_selection,envir=ripaEnv)
			}
		}
		r <- 2
		imagetmp <- limage(as.character(fileName),"reflectance")
		imagetmpScene2 <- lscene(imagetmp,2)
		pb <- tkProgressBar(title = "Generating samples...", min = 0, max = 224, width = 300)		
		for (i in (1:224))
		{
			## Take samples from band i and put into img
			band_temp <- lbandsample(imagetmpScene2,i)
			img <- band_temp@data
			#if (is.nan(img[1])){
			#	for (j in 1:length(img)) img[j] <- 0
			#}
			jpeg("imgtmp.jpg",width=200,height=200)
			plot(imagematrix(img))
			dev.off()
			image1 <- tclVar()
			tcl("image","create","photo",image1,file="imgtmp.jpg")
			sampleList[[i]] <- tcl("label",paste(subfID,".lab",i,sep=""),image=image1)
			if (i%%5 == 0){
				tkgrid(sampleList[[i-4]],row=r,column=1)
				tkgrid(sampleList[[i-3]],row=r,column=2)
				tkgrid(sampleList[[i-2]],row=r,column=3)
				tkgrid(sampleList[[i-1]],row=r,column=4)
				tkgrid(sampleList[[i]],row=r,column=5)
				r <- r+1
			}
			if (i==224){
				tkgrid(sampleList[[i-3]],row=r,column=1)
				tkgrid(sampleList[[i-2]],row=r,column=2)
				tkgrid(sampleList[[i-1]],row=r,column=3)
				tkgrid(sampleList[[i]],row=r,column=4)
			}
			tkbind(sampleList[[i]],"<FocusIn>",function() tcl(sf,"see",sampleList[[i]]))
			tkbind(sampleList[[i]], "<Button-1>",OnClick)
			setTkProgressBar(pb, i, label=paste(round(i/224*100, 0),"% done"))
		}
		## End of samples part
		close(pb)

		onOK <- function(){
			tkdestroy(ttsamples)
			countbands <- length(get('AVIRISbands1',envir=ripaEnv))
			#bandsSet_local <- get('bandsSet',envir=ripaEnv)
			#for (i in 1:224){
			#	if (bandsSet_local[i]==TRUE){
			#		countbands <- countbands + 1
			#		if (tab=="1") assign('AVIRISbands1',c(get('AVIRISbands1',envir=ripaEnv),i),envir=ripaEnv)	
			#		else assign('AVIRISbands2',c(get('AVIRISbands2',envir=ripaEnv),i),envir=ripaEnv)
			#	}
			#}
			
			if (tab=="1"){
				bandsIndexes <- get('AVIRISbands1',envir=ripaEnv)
				assign('bands',vector(length=length(bandsIndexes),mode="list"),envir=ripaEnv)
				bands_local <- get('bands',envir=ripaEnv)
				use_parallel<-tclvalue(get('useParallel',envir=ripaEnv))
				assign('img1',read.aviris(as.character(fileName),bandsIndexes,bands_local,use_parallel),envir=ripaEnv)
				image <-tkrplot(Frame2, function() plot(imagematrix(get('img1',envir=ripaEnv))),vscale=1.04,hscale=0.98)
 				tkpack(image)
 				image <-tkrplot(Frame4, function() plot(imagematrix(get('img1',envir=ripaEnv))),vscale=1.04,hscale=0.98)
 				tkpack(image)
				assign('img2',get('img1',envir=ripaEnv),envir=ripaEnv)
 				assign('numbands1',countbands,envir=ripaEnv)
 				assign('imageType1',"jpgRGB",envir=ripaEnv)
			}
			if (tab=="2"){
				bandsIndexes <- get('AVIRISbands2',envir=ripaEnv)
				assign('bands',vector(length=length(bandsIndexes),mode="list"),envir=ripaEnv)
				bands_local <- get('bands',envir=ripaEnv)
				use_parallel<-tclvalue(get('useParallel',envir=ripaEnv))
				assign('img3',read.aviris(as.character(fileName),bandsIndexes,bands_local,use_parallel),envir=ripaEnv)
				image <-tkrplot(Frame6, function() plot(imagematrix(get('img3',envir=ripaEnv))),vscale=1.04,hscale=0.98)
 				tkpack(image)
 				assign('numbands2',countbands,envir=ripaEnv)
 				assign('imageType2',"jpgRGB",envir=ripaEnv)
			}
		}
		
		OK.but <- tkbutton(ttsamples,text="   OK   ",command=onOK)
		tkpack(OK.but)

	}
	
	# Function to read LAN images
	openLan <- function(){
		tab <- checkTab()
				
		if (tab=="1"){
			assign('AVIRISbands1',NULL,envir=ripaEnv)
			.Tcl(paste(slider, "set",0))
			.Tcl(paste(slider2, "set",1))
	
			tkdestroy(Frame3)
			Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
			tkconfigure(Frame3,width=480,height=510)
			tkplace(Frame3,x=510,y=30)
			
			Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
			tkplace(Frame4,x=1,y=1)
			
			tkdestroy(Frame1)
			Frame1 <- tkframe(lb1,relief="groove",borderwidth=2)
			tkconfigure(Frame1,width=485,height=510)
			tkplace(Frame1,x=10,y=30)
			
			Frame2 <- tkframe(Frame1,relief="groove",borderwidth=0)
			tkplace(Frame2,x=1,y=1)
		}
		
		if (tab=="2"){
			assign('AVIRISbands2',NULL,envir=ripaEnv)
			assign('regionsList',list(),envir=ripaEnv)
			assign('regionsListAux',list(),envir=ripaEnv)
			
			tkconfigure(statisticsTxt,state="normal")
			tkdelete(statisticsTxt,"1.0","end")
			tkconfigure(statisticsTxt,state="disable")
			
			tkdestroy(Frame5)
		
			Frame5 <- tkframe(lb2,relief="groove",borderwidth=2)
			tkconfigure(Frame5,width=480,height=510)
			tkplace(Frame5,x=10,y=30)
	
			Frame6 <- tkframe(Frame5,relief="groove",borderwidth=0)
			tkplace(Frame6,x=1,y=1)
			
			
		}
		
		arquivo <- tkgetOpenFile(filetypes="{{LAN Files} {.lan}}")
		if (!nchar(tclvalue(arquivo))){
			tkmessageBox(message="No file was selected!")
			return()
		}
		arquivo<-tclvalue(arquivo)
		
		if (tab=="1"){
			assign('img1',read.lan(arquivo),envir=ripaEnv)
		}
		if (tab=="2"){
			assign('img3',read.lan(arquivo),envir=ripaEnv)
		}
		
		if (tab=="1"){
			image <-tkrplot(Frame2, function() plot(imagematrix(get('img1',envir=ripaEnv))),vscale=1.04,hscale=0.98)
			tkpack(image)
			image <-tkrplot(Frame4, function() plot(imagematrix(get('img1',envir=ripaEnv))),vscale=1.04,hscale=0.98)
			tkpack(image)
			assign('img2',get('img1',envir=ripaEnv),envir=ripaEnv)
			assign('numbands1',7,envir=ripaEnv)
			assign('imageType1',"lan",envir=ripaEnv)
		}
		
		if (tab=="2"){
			image <-tkrplot(Frame6, function() plot(imagematrix(get('img3',envir=ripaEnv))),vscale=1.04,hscale=0.98)
			tkpack(image)
			assign('numbands2',7,envir=ripaEnv)
			assign('imageType2',"lan",envir=ripaEnv)
		}
	}
	#
	
	# Function to write LAN images
	saveLan <- function(){
		tab <- checkTab()
		
		arquivo <- tkgetSaveFile(filetypes="{{LAN Files} {.lan}}")
		if (!nchar(tclvalue(arquivo))){
			tkmessageBox(message="No file was selected!")
			return()
		}
		arquivo<-tclvalue(arquivo)
		
		if (tab=="1") write.lan(arquivo,get('img2',envir=ripaEnv))
	}
	#
	
	# Function to write JPEG images
	saveJPG <- function(){

		

		arquivo <- tkgetSaveFile(filetypes="{{JPEG Files} {.jpg .jpeg}}")
		if (!nchar(tclvalue(arquivo))){
			tkmessageBox(message="No file was selected!")
			return()
		}
		arquivo<-tclvalue(arquivo)
		
		if (get('numbands1',envir=ripaEnv)>3){		
			img2_local <- get('img2',envir=ripaEnv)
			visualBands1_local <- get('visualBands1',envir=ripaEnv)
			auximg <- array(c(img2_local[,,visualBands1_local[1]],img2_local[,,visualBands1_local[2]],img2_local[,,visualBands1_local[3]]),c(nrow(img2_local),ncol(img2_local),3))
			writeJPEG(auximg,arquivo,quality=100)
		} else	writeJPEG(get('img2',envir=ripaEnv),arquivo,quality=100)
	}

	savePNG <- function(){
		arquivo <- tkgetSaveFile(filetypes="{{PNG Files} {.png}}")
		if (!nchar(tclvalue(arquivo))){
			tkmessageBox(message="No file was selected!")
			return()
		}
		arquivo<-tclvalue(arquivo)
		writePNG(get('img2',envir=ripaEnv),arquivo)
	}

	# Funciton to check tab and call openJpg function
	openJPG <- function(){

		
		
		# Function to read JPG images Tab 0
		openJpgTabOne <- function(){
			
			

			.Tcl(paste(slider, "set",0))
			.Tcl(paste(slider2, "set",1))
		
			tkdestroy(Frame3)
			
			Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
			tkconfigure(Frame3,width=480,height=510)
			tkplace(Frame3,x=510,y=30)
			
			Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
			tkplace(Frame4,x=1,y=1)
			
			tkdestroy(Frame1)
			
			Frame1 <- tkframe(lb1,relief="groove",borderwidth=2)
			tkconfigure(Frame1,width=485,height=510)
			tkplace(Frame1,x=10,y=30)
		
			Frame2 <- tkframe(Frame1,relief="groove",borderwidth=0)
			tkplace(Frame2,x=1,y=1)
			
			arquivo <- tkgetOpenFile(filetypes="{{JPEG Files} {.jpg .jpeg}}")
			if (!nchar(tclvalue(arquivo))){
				tkmessageBox(title="Error",message="No file was selected!")
				return()
			}
			assign('img1',readJPEG(tclvalue(arquivo)),envir=ripaEnv)
			image <-tkrplot(Frame2, function() plot(imagematrix(get('img1',envir=ripaEnv))),vscale=1.04,hscale=0.98)
			tkpack(image)
			image <-tkrplot(Frame4, function() plot(imagematrix(get('img1',envir=ripaEnv))),vscale=1.04,hscale=0.98)
			tkpack(image)
			
			if (is.na(dim(get('img1',envir=ripaEnv))[3])){
				assign('imageType1',"jpgGrey",envir=ripaEnv)
				assign('numbands1',1,envir=ripaEnv)
				assign('img1',array(get('img1',envir=ripaEnv),dim=c(nrow(get('img1',envir=ripaEnv)),ncol(get('img1',envir=ripaEnv)),1)),envir=ripaEnv)
				assign('img2',get('img1',envir=ripaEnv))
			}
			else{
				assign('imageType1',"jpgRGB",envir=ripaEnv)
				assign('numbands1',3,envir=ripaEnv)
				assign('img1',array(get('img1',envir=ripaEnv),dim=dim(get('img1',envir=ripaEnv))),envir=ripaEnv)
				assign('img2',get('img1',envir=ripaEnv))
			}
		}
		#
		
		# Function to read jpeg image Tab 1
		openJpgTabTwo <- function(){
			
			assign('regionsList',list(),envir=ripaEnv)
			assign('regionsListAux',list(),envir=ripaEnv)
			
			tkconfigure(statisticsTxt,state="normal")
			tkdelete(statisticsTxt,"1.0","end")
			tkconfigure(statisticsTxt,state="disable")
			
			tkdestroy(Frame5)
			
			Frame5 <- tkframe(lb2,relief="groove",borderwidth=2)
			tkconfigure(Frame5,width=480,height=510)
			tkplace(Frame5,x=10,y=30)
		
			Frame6 <- tkframe(Frame5,relief="groove",borderwidth=0)
			tkplace(Frame6,x=1,y=1)
			
			arquivo <- tkgetOpenFile(filetypes="{{JPEG Files} {.jpg .jpeg}}")
			if (!nchar(tclvalue(arquivo))){
				tkmessageBox(message="No file was selected!")
				return()
			}
			
			assign('img3',readJPEG(tclvalue(arquivo)),envir=ripaEnv)
			#img3_local <- get('img3',envir=ripaEnv)
			image <-tkrplot(Frame6, function() plot(imagematrix(get('img3',envir=ripaEnv))),vscale=1.04,hscale=0.98)
			tkpack(image)
			if (is.na(dim(get('img3',envir=ripaEnv))[3])){
				assign('imageType2',"jpgGrey",envir=ripaEnv)
				assign('numbands2',1,envir=ripaEnv)
				assign('img3',array(get('img3',envir=ripaEnv),dim=c(nrow(get('img3',envir=ripaEnv)),ncol(get('img3',envir=ripaEnv)),1)),envir=ripaEnv)
				#assign('img3',array(img3_local,dim=c(nrow(img3_local),ncol(img3_local),1)))
			}
			else{
				assign('imageType2',"jpgRGB",envir=ripaEnv)
				assign('numbands2',3,envir=ripaEnv)
				assign('img3',array(get('img3',envir=ripaEnv),dim=dim(get('img3',envir=ripaEnv))),envir=ripaEnv)
				#assign('img3',array(img3_local,dim=dim(img3_local)),envir=ripaEnv)
			}
		}
		#
		
		tab <- checkTab()
				
		if (tab=="1"){
			openJpgTabOne()
			assign('AVIRISbands1',NULL,envir=ripaEnv)
		}
		if (tab=="2"){
			openJpgTabTwo()
			assign('AVIRISbands2',NULL,envir=ripaEnv)
		}
	}
	#
	
	# Funciton to check tab and call openPng function
	openPNG <- function(){
		

		# Function to read PNG images Tab 0
		openPngTabOne <- function(){
			
			.Tcl(paste(slider, "set",0))
			.Tcl(paste(slider2, "set",1))
		
			tkdestroy(Frame3)
			
			Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
			tkconfigure(Frame3,width=480,height=510)
			tkplace(Frame3,x=510,y=30)
			
			Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
			tkplace(Frame4,x=1,y=1)
			
			tkdestroy(Frame1)
			
			Frame1 <- tkframe(lb1,relief="groove",borderwidth=2)
			tkconfigure(Frame1,width=485,height=510)
			tkplace(Frame1,x=10,y=30)
		
			Frame2 <- tkframe(Frame1,relief="groove",borderwidth=0)
			tkplace(Frame2,x=1,y=1)
			
			arquivo <- tkgetOpenFile(filetypes="{{PNG Files} {.png}}")
			if (!nchar(tclvalue(arquivo))){
				tkmessageBox(title="Error",message="No file was selected!")
				return()
			}
			assign('img1',readPNG(tclvalue(arquivo)),envir=ripaEnv)

			image <-tkrplot(Frame2, function() plot(imagematrix(get('img1',envir=ripaEnv))),vscale=1.04,hscale=0.98)
			tkpack(image)
			image <-tkrplot(Frame4, function() plot(imagematrix(get('img1',envir=ripaEnv))),vscale=1.04,hscale=0.98)
			tkpack(image)
			
			if (is.na(dim(get('img1',envir=ripaEnv))[3])){
				assign('imageType1',"pngGrey",envir=ripaEnv)
				assign('numbands1',1,envir=ripaEnv)
				assign('img1',array(get('img1',envir=ripaEnv),dim=c(nrow(get('img1',envir=ripaEnv)),ncol(get('img1',envir=ripaEnv)),1)),envir=ripaEnv)
				assign('img2',get('img1',envir=ripaEnv),envir=ripaEnv)				
			}
			else{
				assign('imageType1',"pngRGB",envir=ripaEnv)
				assign('numbands1',3,envir=ripaEnv)
				assign('img1',array(get('img1',envir=ripaEnv),dim=dim(get('img1',envir=ripaEnv))),envir=ripaEnv)
				assign('img2',get('img1',envir=ripaEnv),envir=ripaEnv)
			}
		}
		#
		
		# Function to read png image Tab 1
		openPngTabTwo <- function(){
			
			assign('regionsList',list(),envir=ripaEnv)
			assign('regionsListAux',list(),envir=ripaEnv)
			
			tkconfigure(statisticsTxt,state="normal")
			tkdelete(statisticsTxt,"1.0","end")
			tkconfigure(statisticsTxt,state="disable")
			
			tkdestroy(Frame5)
			
			Frame5 <- tkframe(lb2,relief="groove",borderwidth=2)
			tkconfigure(Frame5,width=480,height=510)
			tkplace(Frame5,x=10,y=30)
		
			Frame6 <- tkframe(Frame5,relief="groove",borderwidth=0)
			tkplace(Frame6,x=1,y=1)
			
			arquivo <- tkgetOpenFile(filetypes="{{PNG Files} {.png}}")
			if (!nchar(tclvalue(arquivo))){
				tkmessageBox(message="No file was selected!")
				return()
			}
			
			assign('img3',readPNG(tclvalue(arquivo)),envir=ripaEnv)
			#img3_local <- get('img3',envir=ripaEnv)
			#image <-tkrplot(Frame6, function() plot(imagematrix(img3_local)),vscale=1.04,hscale=0.98)
			image <-tkrplot(Frame6, function() plot(imagematrix(get('img3',envir=ripaEnv))),vscale=1.04,hscale=0.98)
			tkpack(image)
			if (is.na(dim(get('img3',envir=ripaEnv))[3])){
				assign('imageType2',"pngGrey",envir=ripaEnv)
				assign('numbands2',1,envir=ripaEnv)
				assign('img3',array(get('img3',envir=ripaEnv),dim=c(nrow(get('img3',envir=ripaEnv)),ncol(get('img3',envir=ripaEnv)),1)),envir=ripaEnv)
				#assign('img3',array(img3_local,dim=c(nrow(img3_local),ncol(img3_local),1)),envir=ripaEnv)
			}
			else{
				assign('imageType2',"pngRGB",envir=ripaEnv)
				assign('numbands2',3,envir=ripaEnv)
				assign('img3',array(get('img3',envir=ripaEnv),dim=dim(get('img3',envir=ripaEnv))),envir=ripaEnv)
				#assign('img3',array(img3_local,dim=dim(img3_local)),envir=ripaEnv)
			}
		}
		#
		
		tab <- checkTab()
				
		if (tab=="1"){
			openPngTabOne()
			assign('AVIRISbands1',NULL,envir=ripaEnv)
		}
		if (tab=="2"){
			openPngTabTwo()
			assign('AVIRISbands2',NULL,envir=ripaEnv)
		}
	}
	#

	# Function to apply the negative function to an image
	negative <- function(){
		tab <- checkTab()
		
		if (tab=="2"){
			tkmessageBox(title="Error",message="Please, use this menu only with the first tab!",icon="error",type="ok")
			return()
		}
	
		if (is.null(get('img1',envir=ripaEnv))){
			tkmessageBox(title="Error",message="Please, open an image in order to use it!",icon="error",type="ok")
			return()
		}
		
		tkdestroy(Frame3)
		Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
		tkconfigure(Frame3,width=485,height=510)
		tkplace(Frame3,x=510,y=30)
		Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
		tkplace(Frame4,x=1,y=1)
	
		assign('img2',1-get('img1',envir=ripaEnv),envir=ripaEnv)
	

		img2_local <- get('img2',envir=ripaEnv)
		visualBands1_local <- get('visualBands1',envir=ripaEnv)

		if (get('numbands1',envir=ripaEnv)>3){
			auximg <- array(c(img2_local[,,visualBands1_local[1]],img2_local[,,visualBands1_local[2]],img2_local[,,visualBands1_local[3]]),c(nrow(img2_local),ncol(img2_local),get('numbands1',envir=ripaEnv)))
			image <-tkrplot(Frame4, function() plot(imagematrix(auximg)),vscale=1.04,hscale=0.98)
		}else{
			image <-tkrplot(Frame4, function() plot(imagematrix(img2_local)),vscale=1.04,hscale=0.98)
		}
		tkpack(image)
	}
	#

	showAbout <-function(){
		tkmessageBox(title="About",message="R Image Processing and Analysis\nVersion 2.0-2\n\nTalita Perciano - LBNL, USA \nAlejandro C. Frery - IC/UFAL, Brazil",icon="info",type="ok")
	}


	highPOptions <- function(){
		tt_parallel <- tktoplevel()
		tktitle(tt_parallel)<-"Hight Performance Options"
		cb_parallel <- tkcheckbutton(tt_parallel)
		useParallel_local <- get('useParallel',envir=ripaEnv)
		numProcessors_local <- get('numProcessors',envir=ripaEnv)
		entry.processors <-tkentry(tt_parallel,width="20",textvariable=numProcessors_local)
		tkconfigure(cb_parallel,variable=useParallel_local)
		tkgrid(tklabel(tt_parallel,text="Use parallel process suport"),cb_parallel)
		tkgrid(tklabel(tt_parallel,text="Please enter the number of processors"))
		tkgrid(entry.processors)
		OnOK <- function()
		{
			assign('useParallel',useParallel_local,envir=ripaEnv)	
			if ( (as.integer(tclvalue(useParallel_local))==1) ) {
				if (is.na(as.integer(tclvalue(numProcessors_local)))) {
					tkmessageBox(title="Error",message="Please enter a valid number of processors!",icon="error",type="ok")
				}
				else {
					assign('numProcessors',numProcessors_local,envir=ripaEnv)
					if ( as.integer(tclvalue(useParallel_local))==1 ) {
						numproc <- as.integer(tclvalue(numProcessors_local))
						assign('cluster',makeCluster(numproc),envir=ripaEnv)
						#registerDoParallel(get('cluster',envir=ripaEnv))

						if (Sys.info()[['sysname']]=='Windows'){
							registerDoSNOW(get('cluster',envir=ripaEnv))
						}
						else {
							doMC::registerDoMC(cores=numproc)
						}
					}
					tkdestroy(tt_parallel)
				}
			}
			else {
				tkdestroy(tt_parallel)
			}
		}
		OK.but <- tkbutton(tt_parallel,text="OK",command=OnOK)
		tkgrid(OK.but)
		
	}
	
	##################################################################################################################
	

	fontText <- tkfont.create(family="book",size=10)
	fr1 <- tkframe(get('tn',envir=ripaEnv))
	tkadd(get('tn',envir=ripaEnv),fr1,text="Operations")
	lb1 <- ttklabelframe(fr1)
	tkpack(lb1)
	tkconfigure(lb1,labelanchor="n",text="Operations Tab",width=1010,height=665)

	#OperationsTbn <<- tclvalue(tkadd(tn, label="Operations"))
	
	# Left frame
	Frame1 <- tkframe(lb1,relief="groove",borderwidth=2)
	tkconfigure(Frame1,width=485,height=510)
	tkplace(Frame1,x=10,y=30)
	#
	
	# Frame inside left frame
	Frame2 <- tkframe(Frame1,relief="groove",borderwidth=0)
	tkplace(Frame2,x=1,y=1)
	#
	
	# Right frame
	Frame3 <- tkframe(lb1,relief="groove",borderwidth=2)
	tkconfigure(Frame3,width=485,height=510)
	tkplace(Frame3,x=510,y=30)
	#
	
	# Frame inside right frame
	Frame4 <- tkframe(Frame3,relief="groove",borderwidth=0)
	tkplace(Frame4,x=1,y=1)
	#
	
	Label1 <- tklabel(lb1,text="Original image:")
	tkplace(Label1,x=10,y=10)
	
	Label2 <- tklabel(lb1,text="Changed image:")
	tkplace(Label2,x=510,y=10)
	
	SliderValue <- tclVar("0")
	SliderValueLabel <- tklabel(lb1,text=as.character(tclvalue(SliderValue)))
	Label3 <- tklabel(lb1,text="Brightness Value: ")
	tkplace(Label3,x=10,y=545)
	tkplace(SliderValueLabel,x=940,y=545)
	tkconfigure(SliderValueLabel,textvariable=SliderValue)
	slider <- tkscale(lb1, from=-1, to=1,
                   showvalue=F, variable=SliderValue,
                   resolution=0.01, orient="horizontal",length=810)
	tkplace(slider,x=120,y=545)
	
	tkbind(slider,"<ButtonRelease-1>",function() sliderFunction())
	
	SliderValue2 <- tclVar("0")
	SliderValueLabel2 <- tklabel(lb1,text=as.character(tclvalue(SliderValue2)))
	Label4 <- tklabel(lb1,text="Contrast Value: ")
	tkplace(Label4,x=10,y=575)
	tkplace(SliderValueLabel2,x=940,y=575)
	tkconfigure(SliderValueLabel2,textvariable=SliderValue2)
	
	slider2 <-.Tk.subwin(lb1)
	slider2 <- tkscale(lb1, from=0, to=1,
                   showvalue=F, variable=SliderValue2,
                   resolution=0.01, orient="horizontal",length=810)
	tkplace(slider2,x=120,y=575)
	
	tkbind(slider2,"<ButtonRelease-1>",function() sliderFunction())
	
	##################################################################################################################
	
	##################################################################################################################
	####################################### ROI Tab ##################################################################
	##################################################################################################################
	fr2 <- tkframe(get('tn',envir=ripaEnv))
	tkadd(get('tn',envir=ripaEnv),fr2,text="Regions of Interest")
	lb2 <- ttklabelframe(fr2)
	tkpack(lb2)
	tkconfigure(lb2,labelanchor="n",text="Operations Tab",width=1010,height=665)
	
	# Left frame
	Frame5 <- tkframe(lb2,relief="groove",borderwidth=2)
	tkconfigure(Frame5,width=485,height=510)
	tkplace(Frame5,x=10,y=30)
	#
	
	# Frame inside left frame
	Frame6 <- tkframe(Frame5,relief="groove",borderwidth=0)
	tkplace(Frame6,x=1,y=1)
	#
	
	# Right frame
	Frame7 <- tkframe(lb2,relief="groove",borderwidth=2)
	tkconfigure(Frame7,width=485,height=510)
	tkplace(Frame7,x=510,y=30)
	#
	
	# Configure text space inside right frame
	scr <- tkscrollbar(Frame7, repeatinterval=5,command=function(...)tkyview(statisticsTxt,...))
	statisticsTxt <- tktext(Frame7,bg="white",font="courier",yscrollcommand=function(...)tkset(scr,...))
	tkplace(statisticsTxt,x=1,y=1,width=465,height=502)
	tkconfigure(statisticsTxt, state="disabled")
	tkplace(scr,x=468,y=1)
	tkplace.configure(scr,height=502)
	#
	
	Label3 <- tklabel(lb2,text="Original image:")
	tkplace(Label3,x=10,y=10)
	
	Label4 <- tklabel(lb2,text="Actual Region Statistics:")
	tkplace(Label4,x=510,y=10)
	
	Button1 <- tkbutton(lb2,text="Choose Region",command=region)
	tkplace(Button1,x=10,y=560)

	##################################################################################################################
	
	##################################################################################################################
	######################################### Menu ###################################################################
	topMenu<-tkmenu(tt)
	tkconfigure(tt,menu=topMenu)
	
	# File menu
	fileMenu<-tkmenu(topMenu,tearoff=FALSE)
	#
	
	# Operations menu
	operationsMenu<-tkmenu(topMenu,tearoff=FALSE)
	#
	
	# Regions Menu
	regionsMenu<-tkmenu(topMenu,tearoff=FALSE)
	#
	
	# Bands Menu
	bandsMenu <- tkmenu(topMenu,tearoff=FALSE)
	#
	
	# High Performance Menu
	highPMenu <- tkmenu(topMenu, tearoff=FALSE)

	# Help Menu
	helpMenu <- tkmenu(topMenu,tearoff=FALSE)
	#

	# Open menu inside File menu
	openMenu <- tkmenu(topMenu,tearoff=FALSE)
	#
	
	# Filters menu inside Operations Menu
	filtersMenu <- tkmenu(topMenu,tearoff=FALSE)
	#
	
	# Segmentation menu inside Operations Menu
	segmentationMenu <- tkmenu(topMenu,tearoff=FALSE)
	#
	
	# Configure Regions Menu
	tkadd(regionsMenu,"command",label="Show All Regions",command=showAllRegions,underline=0,accelerator="Ctrl+S")
	tkbind(tt,"<Control-s>",function() showAllRegions())
	tkadd(regionsMenu,"command",label="Covariance Matrix",command=covMatrix)
	tkadd(regionsMenu,"command",label="Regions Brushplots",command=regionBrush)
	#
	
	# Configure Bands Menu
	tkadd(bandsMenu,"command",label="Image Bands",command=imageBands,underline=0,accelerator="Ctrl+M")
	tkbind(tt,"<Control-m>",function() imageBands())
	tkadd(bandsMenu,"command",label="Regions Bands",command=regionBands,underline=0,accelerator="Ctrl+G")
	tkbind(tt,"<Control-g>",function() regionBands())
	tkadd(bandsMenu,"command",label="Bands Brushplots",command=brush)
	#
	
	# Configure Dynamic Graphics Menu
	dynGraph<-tkmenu(topMenu,tearoff=FALSE)
	tkadd(regionsMenu,"cascade",label="Dynamic Graphics",menu=dynGraph)
	tkadd(dynGraph,"command",label="Min",command=function() dynFunc(1),underline=1,accelerator="Ctrl+I")
	tkbind(tt,"<Control-i>",function() dynFunc(1))
	tkadd(dynGraph,"command",label="Max",command=function() dynFunc(2),underline=1,accelerator="Ctrl+A")
	tkbind(tt,"<Control-a>",function() dynFunc(2))
	tkadd(dynGraph,"command",label="Mean",command=function() dynFunc(3),underline=1,accelerator="Ctrl+E")
	tkbind(tt,"<Control-e>",function() dynFunc(3))
	tkadd(dynGraph,"command",label="Median",command=function() dynFunc(4),underline=2,accelerator="Ctrl+D")
	tkbind(tt,"<Control-d>",function() dynFunc(4))
	tkadd(dynGraph,"command",label="Standard Deviation",command=function() dynFunc(5),underline=1,accelerator="Ctrl+T")
	tkbind(tt,"<Control-t>",function() dynFunc(5))
	tkadd(dynGraph,"command",label="Median Deviation",command=function() dynFunc(6),underline=9,accelerator="Ctrl+V")
	tkbind(tt,"<Control-v>",function() dynFunc(6))
	tkadd(dynGraph,"command",label="Kurtosis",command=function() dynFunc(7),underline=0,accelerator="Ctrl+K")
	tkbind(tt,"<Control-k>",function() dynFunc(7))
	tkadd(dynGraph,"command",label="Skewness",command=function() dynFunc(8),underline=3,accelerator="Ctrl+W")
	tkbind(tt,"<Control-w>",function() dynFunc(8))
	tkadd(dynGraph,"command",label="All values",command=function() dynFunc(9))
	#
	
	# Configure Graphics Menu
	#graph<<-tkmenu(topMenu,tearoff=FALSE)
	#tkadd(regionsMenu,"cascade",label="Graphics",menu=graph)
	#tkadd(graph,"command",label="Histograms",command=function() histograms(),underline=0,accelerator="Ctrl+H")
	#tkbind(tt,"<Control-h>",function() histograms())
	#
	
	# Configure File menu
	tkadd(fileMenu,"cascade",label="Open Image",menu=openMenu)
	tkadd(fileMenu,"command",label="Change Visual Bands",command=function() changeBands(),underline=0,accelerator="Ctrl+C")
	tkbind(tt,"<Control-c>",function() changeBands())
	tkadd(fileMenu,"command",label="Save Image",command=saveImg,underline=0,accelerator="Ctrl+S")
	tkbind(tt,"<Control-s>",function() saveImg())
	tkadd(fileMenu,"command",label="Quit",command=function() quit(),underline=0,accelerator="Ctrl+Q")
	tkbind(tt,"<Control-q>",function() quit())
	#
	
	# Configure Operations menu
	tkadd(operationsMenu,"command",label="Stretch Image",command=stretch,underline=2,accelerator="Ctrl+R")
	tkbind(tt,"<Control-r>",function() stretch())
	tkadd(operationsMenu,"command",label="Negative", command=negative,underline=0,accelerator="Ctrl+N")
	tkbind(tt,"<Control-n>",function() negative())
	tkadd(operationsMenu,"command",label="Quality Index", command=qualityIndex,underline=8,accelerator="Ctrl+I")
	tkbind(tt,"<Control-i>",function() qualityIndex())
	tkadd(operationsMenu,"command",label="Equalize", command=equalizeImg,underline=0,accelerator="Ctrl+E")
	tkbind(tt,"<Control-e>",function() equalizeImg())
	#tkadd(operationsMenu,"command",label="Histograms",command=imgHistograms,underline=2,accelerator="Ctrl+G")
	#tkbind(tt,"<Control-g>",function() imgHistograms())
	tkadd(operationsMenu,"command",label="Pixels Values",command=pixelsValues,underline=2,accelerator="Ctrl+X")
	tkbind(tt,"<Control-x>",function() pixelsValues())
	tkadd(operationsMenu,"command",label="Edit Pixels Values",command=editPixels)
	tkadd(operationsMenu,"command",label="Zoom",command=zoomImg)
	tkadd(operationsMenu,"command",label="PCA",command=pca)
	#
	
	# To put Filters menu inside Operations menu
	tkadd(operationsMenu,"cascade",label="Filters",menu=filtersMenu)
	#
	
	# To put Segmentation menu inside Operations menu
	tkadd(operationsMenu,"cascade",label="Segmentation",menu=segmentationMenu)
	#
	
	# Image extensions inside Open menu
	tkadd(openMenu,"command",label="JPG",command=openJPG,underline=0,accelerator="Ctrl+J")
	tkbind(tt,"<Control-j>",function() openJPG())
	tkadd(openMenu,"command",label="PNG",command=openPNG,underline=0)
	tkadd(openMenu,"command",label="LAN",command=openLan,underline=0,accelerator="Ctrl+L")
	tkbind(tt,"<Control-l>",function() openLan())
	tkadd(openMenu,"command",label="AVIRIS",command=openAVIRIS)
	#
	
	# Available filters inside Filters menu
	tkadd(filtersMenu,"command",label="Sobel Filter",command=function() filters(1),underline=2,accelerator="Ctrl+B")
	tkbind(tt,"<Control-b>",function() filters(1))
	tkadd(filtersMenu,"command",label="Laplacian Filter",command=function() filters(2),underline=2,accelerator="Ctrl+P")
	tkbind(tt,"<Control-p>",function() filters(2))
	tkadd(filtersMenu,"command",label="Mean Filter",command=function() filters(3),underline=0,accelerator="Ctrl+M")
	tkbind(tt,"<Control-m>",function() filters(3))
	tkadd(filtersMenu,"command",label="High-pass Filter",command=function() filters(4),underline=0,accelerator="Ctrl+H")
	tkbind(tt,"<Control-h>",function() filters(4))
	tkadd(filtersMenu,"command",label="Low-pass Filter",command=function() filters(5),underline=2,accelerator="Ctrl+W")
	tkbind(tt,"<Control-w>",function() filters(5))
	tkadd(filtersMenu,"command",label="Minimum Filter",command=function() filters(6))
	tkadd(filtersMenu,"command",label="Maximum Filter",command=function() filters(7))
	tkadd(filtersMenu,"command",label="Median Filter",command=function() filters(8))
	#
	
	# Available segmentation techniques inside Segmentation Menu
	tkadd(segmentationMenu, "command", label="Thresholding",command=function() segmentation(1))
	#
	
	# Configure High Performance menu
	tkadd(highPMenu,"command",label="Options",command=function() highPOptions(),underline=0)

	# Configure Help menu
	tkadd(helpMenu,"command",label="About",command=function() showAbout(),underline=0)
	#

	# Insert Menus
	tkadd(topMenu,"cascade",label="File",menu=fileMenu,underline=0)
	tkadd(topMenu,"cascade",label="Operations",menu=operationsMenu,underline=0)
	tkadd(topMenu,"cascade",label="Regions",menu=regionsMenu,underline=0)
	tkadd(topMenu,"cascade",label="Bands",menu=bandsMenu,underline=0)
	tkadd(topMenu,"cascade",label="High Performance",menu=highPMenu,underline=0)
	tkadd(topMenu,"cascade",label="Help",menu=helpMenu,underline=0)
	#
	
	tkselect(get('tn',envir=ripaEnv),0)
	tkfocus(tt)
}
##################################################################################################################
