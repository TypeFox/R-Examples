#Classes et methodes pour les type de donnees specifiques a GENLIB
#1-CG
#2-F
#3-GLgen
#4-GLgroup
#5-GLmultiArray4
#6-GLmultiMatrix
#7-GLmultiNumber
#8-GLmultiVector
#9-GLnoone
#10-GLmultiPhi
#######
#OBJET GLgroup
#
#Objet qui represente un ensemble de proband dans different groupe
#Chaque etement de la liste est le nom d'un groupe avec le no du proband
setClass("GLgroup",representation(), contains="list")
#Methode Length correct

setMethod("show","GLgroup",function(object){
		#Affiche les donnees concernant les groupes
		lapply(1:length(object),function(x,obj) 
		{		 	
			cat("GROUP NAME : ",names(obj)[x],"\n")
			cat("PROBANTS : ",obj[[x]],"\n\n")			
		},obj=object)
})

setMethod("[","GLgroup", function(x,...,drop){
	#Subscript qui ne change pas la classe c'est tout
	x<-as(x,"list")
	ret <- x[...,drop=drop]
	class(ret)<-"GLgroup"
	ret
})

#CG
setClass("GLCGMatrixGroupSingle",representation(group="GLgroup",grindex="list"), contains="matrix")  #Une depth
# la matrice est dans le sens [proband, ancestor]

#FONCTION VIRTUELLE ASSOCIE
#if (isGeneric("group")==F)  group <- function(x){x@group}

setGeneric("group",function(x) standardGeneric("group") )
setMethod("group","GLCGMatrixGroupSingle",function(x) x@group )

#setMethod("Dim","GLCGMatrixGroupSingle",function(object, ...){
Dim <- function(object, ...) {
	param <- list(...)
	if (!is.null(param$drop) && param$drop==F){	
		#Dans ce cas, group seulement
		c(length(object@group),dim(object)[2])
	}else{	
		#Dans ce cas, ces les phi moyen...
		c(length(object@group),dim(object)[2])
	}
}#)

setMethod("show","GLCGMatrixGroupSingle",function(object){
	#S'assure de toujours affiche une matrice et non un vecteur...
	xgrindex <- object@grindex
	xgroupe   <- object@group
	x <- unclass(object) 						
	z=GLapplyCG(object,mean)
	dim(z)<-c(length(xgroupe),dim(x)[2])
	dimnames(z)<-list(names(xgroupe),dimnames(x)[[2]])
	prmatrix(z)
})

#EXTRACTION
setMethod("[","GLCGMatrixGroupSingle",function(x,...,drop) GLPrivExtCG(x=x,...,drop=drop) )

#METHODES MATHEMATIQUE POUR CG
#Operateur mathematique (Pas les operateurs car sa marche pas sur l'autre sens 1 + obj)
setMethod("Ops","GLCGMatrixGroupSingle",function(e1,e2=NULL){
	stop("\nYou can not use operations on a GLCGMatrixGroupSingle object\nUse the empty [] operator to convert into matrix. (ex: z[]+1)")
	#stop("\nVous ne pouvez pas faire d'operation sur un GLCGMatrixGroupSingle\nUtilise l'operateur [] vide pour convertir en matrice EX: z[]+1")
})

setReplaceMethod("[","GLCGMatrixGroupSingle",function(x,...,value){
	stop("\nYou can not directly modify a GLCGMatrixGroupSingle object\nYou must create a new object using gen.gc()")
	#stop("\nVous ne pouvez pas faire de modification directement sur un GLCGMatrixGroupSingle\nVous devez recreer un nouvel objet a l'aide de gen.cggroupe()")
})

setMethod("Math","GLCGMatrixGroupSingle",function(x){
	callGeneric(x[])
})

setMethod("Math2","GLCGMatrixGroupSingle",function(x, digits){
	callGeneric(x[],digits)
})

setMethod("Summary","GLCGMatrixGroupSingle",function(x, ..., na.rm = F){
	callGeneric(x[], ..., na.rm)})

#F
#DECLARATION DES NOUVELLES FONCTION GENERIQUE (DUPPLIQUE POUR EVITE ERREUR)
#if (isGeneric("depth")==F) 
#setGeneric("depth",function(x) standardGeneric("depth") )
#if (isGeneric("group")==F)     
#setGeneric("group",function(x){ methods::standardGeneric("group") })


#class d'inpiration : 
setClass("GLmultiVector",representation(depth="integer"), contains="matrix")
setClass("GLmultiFGroup",representation(group="GLgroup",grindex="list"), contains="GLmultiVector") #Plusieur depth
setClass("GLmultiFGroupSingle",representation(group="GLgroup",grindex="list"), contains="matrix")  #Une depth

#Fonctions virtuelles associees (jusqu'ici pas necessaire)
setMethod("group","GLmultiFGroup"      ,function(x) x@group)
setMethod("group","GLmultiFGroupSingle",function(x) x@group)

setMethod("Dim","GLmultiFGroup",function(object, ...){	
	z=dim(object)	
	c(z[length(z)],length(object@group))
})

setMethod("Dim","GLmultiFGroupSingle",function(object, ...){
	c(length(object@group))
})

setMethod("show","GLmultiFGroup",function(object){
	#Affiche les donnees sur le group
	show(GLapplyF(object,mean,named=T))
})

setMethod("show","GLmultiFGroupSingle",function(object){
	#Affiche les donnees sur le group
	show(GLapplyF(object,mean,named=T))
})

#EXTRACTION 
setMethod("[","GLmultiFGroupSingle",function(x,...,drop) GLPrivExtFSINGLE(x=x,...,drop=drop))

#Methode mathematique
setMethod("Ops","GLmultiFGroupSingle",function(e1,e2=NULL){
	stop("\nYou can not use operations on a GLmultiFGroupSingle object\nUse the empty [] operator to convert into matrix. (ex: z[]+1)")
	#stop("\nVous ne pouvez pas faire d'operation sur un GLmultiFGroupSingle\n	Utilise l'operateur [] vide pour convertir en matrice EX: z[]+1")
})

setReplaceMethod("[","GLmultiFGroupSingle",function(x,...,value){
	stop("\nYou can not directly modify a GLmultiFGroupSingle object")
})

setMethod("Math","GLmultiFGroupSingle",function(x){
	callGeneric(x[])
})

setMethod("Math2","GLmultiFGroupSingle",function(x, digits)
{
	callGeneric(x[],digits)
})

setMethod("Summary","GLmultiFGroupSingle",function(x, ..., na.rm = F)
{
	callGeneric(x[], ..., na.rm)
})

#OBJET 2: Plusieurs depth
setMethod("[","GLmultiFGroup",function(x,...,drop) GLPrivExtF(x=x,...,drop=drop))

#Operateur mathematique (Pas les operateurs car sa marche pas sur l'autre sens 1 + obj)
setMethod("Ops","GLmultiFGroup",function(e1,e2=NULL){
	stop("\nYou can not use operations on a GLmultiFGroup object\nUse the empty [] operator to convert into matrix. (ex: z[]+1)")
	#stop("\nVous ne pouvez pas faire d'operation sur un GLmultiFGroup\nUtiliser l'operateur [] vide pour convertir en matrice EX: z[]+1")
})

setReplaceMethod("[","GLmultiFGroup",function(x,...,value){
	stop("\nYou can not directly modify a GLmultiFGroup object")
})

setMethod("Math","GLmultiFGroup",function(x){
	callGeneric(x[])
})

setMethod("Math2","GLmultiFGroup",function(x, digits){
	callGeneric(x[],digits)
})

setMethod("Summary","GLmultiFGroup",function(x, ..., na.rm = F){
	callGeneric(x[], ..., na.rm)
})

#GLgen
#Declaration de fonctions generiques
#if (isGeneric("depth")==F)
#	setGeneric("depth",function(x) standardGeneric("depth") )
	
#if (isGeneric("group")==F)
#	setGeneric("group",function(x) standardGeneric("group") )

setGeneric("depth",function(x) standardGeneric("depth") )

setMethod("depth","ANY",function(x) {NULL})
setMethod("group","ANY",function(x) {NULL})

setClass("GLgen",representation(.Data="integer",Date="character"),
	validity=	function(object){
		#print(paste(length(x@.Data),x@.Data[4],x@.Data[1],x@.Data[2],x@.Data[3],x@.Data[length(x@.Data)]))
		#Pour la version courante soit la 8
	#print( length(object@.Data))#<16 )
	#print( object@.Data[4])#!=8 )
	#print( object@.Data[1])#!=71)
	#print( object@.Data[2])#!=69)
	#print( object@.Data[3])#!=78)
	#print( object@.Data[1:100])#!=99999999 )
	#print( grep(99999999, object@.Data))
	#print( object@.Data[(length(object@.Data)-30):length(object@.Data)])#!=99999999 )

		if (length(object@.Data)<16 || object@.Data[4]!=1
			|| object@.Data[1]!=71 || object@.Data[2]!=69 || object@.Data[3]!=78 
			|| object@.Data[length(object@.Data)]!=99999999 ){
					stop("Invalid GLgen object")
		}
		return (T)
	}
)

setMethod("show","GLgen",function(object){
		#Affiche les donnees concernant une genealogie
		cat(" GENLIB: Genealogical Data Software version ", object@.Data[4],"\n\n", #Ensemble de genealogies version


		" Number of individuals : ", object@.Data[9], "\n",
		" Number of parent-child relations : ", object@.Data[10], "\n")
		if (object@.Data[12]!=-1)
			cat( " Number of men  : ", object@.Data[12], "\n",
				" Number of women : ", (object@.Data[9]-object@.Data[12]),"\n", 
				" Number of subjects : ", length(gen.pro(gen.genealogy(object))), "\n\n")
		cat("\n Genealogical depth : ", object@.Data[11], "\n\n", " Created on : ", object@Date, "\n", sep = "")
})

setMethod("summary","GLgen",function(object,...) {show(object)} )

setMethod("depth","GLgen",function(x) {x@.Data[11]} )

setMethod("length","GLgen",function(x){
		#la taille correspond au nombre d individu
		#stop("Illegal use of length")
		x@.Data[9]
})

#OBJET GLmultiArray4
#
#Objet qui represente un vecteur de matrice de vecteur
#le 1e vecteur correspond au differente depth
#
################################################

setClass("GLmultiArray4",representation(depth="integer"), contains="array")
#dimension Matrice1,Matrice2,TailleVecteur,depth

#FONCTION VIRTUELLE ASSOCIE
setMethod("depth","GLmultiArray4",function(x) x@depth)
#setMethod("length","GLmultiArray4",function(x) length(x@depth))
setMethod("Dim","GLmultiArray4",function(object, ...) {z=dim(object);l=length(z);z[c(l,1:(l-1))]})

setMethod("show","GLmultiArray4",function(object)
{
		depth 	<- object@depth
		object <- as(object,"array")
		ldim	<-	length(dim(object))

		#Creation de la matrice temporaire
		tmpm   <-  matrix("101010101010101010101010101010100110010101010101010",nrow=dim(object)[1],ncol=dim(object)[2])
		dimnames(tmpm) <- dimnames(object)[1:2]
		
		lapply(1:length(depth),function(x,obj,ldim,depth,tmpm) 
		{
			#Pour chaque depth
		 	cat("depth : ",depth[x],"\n")
			if (dim(obj)[3]>1)
			{
				#CREATION D'UNE MATRICE DE CHARACTER
				t1 <- format( c(obj[,,1,x],obj[,,dim(obj)[3],x]) )
				#CREATION DE LA MATRICE TEMPORAIRE
				t2 <- length(t1)/2
				for(a in 1:(length(t1)/2) )
					tmpm[a] <- paste("[",t1[a]," - ",t1[a+t2],"]",sep="")				
				#AFFICHAGE
				prmatrix(tmpm, quote = F) # prmatrix(x, rowlab =, collab =, quote = TRUE, right = FALSE, na.print = NULL, ...)
			}
			else
				show(obj[,,1,x])		
			cat("\n")
		},obj=object,ldim=ldim,depth=depth,tmpm=tmpm)
})
setMethod("summary","GLmultiArray4",function(object,...)
{
		#Affiche les donnees concernant un vecteur de matrice
		cat(" GENLIB: Matrix on many depths","\n",
		 " Number of lines  : ",  dim(object)[1], "\n",
		 " Number of columns: ",  dim(object)[2], "\n",
		 " Vector size      : ",  dim(object)[3], "\n",
		 " depth(s) included : ",sep = "")
		cat(object@depth,sep=",")
		cat("\n")
})

#METHODE D'EXTRACTION
setMethod("[","GLmultiArray4", function(x,...,drop) {
		#extrait pour certaine depth
		#GLOverGroup3 <- function(pro,dim1,dim2,...,named,abs)
		#GLOverGroup4 <- function(pro,dim1,dim2,dim3,...,named,abs)
		#list(pro=pro,dim1=dim1,dim2=dim2,param=p,named=n,abs=ab)	
		l2 <- GLOverGroup4(...)
		pro <- l2$pro
		dim1 <- l2$dim1
		dim2 <- l2$dim2
		dim3 <- l2$dim3

		if(l2$param==0) return(x)			
		
		if(l2$param==1 || l2$param==4)
		{
			#Cas particulier depth
			if (is(pro,"GLnothing"))	#Proband vide
				class(pro)<-"missing"	
			else
				if (is.numeric(pro) && !l2$abs) #Si c'est pas valeur absolue
				{
					#Recherche comme prevu dans la liste de depth
					pos <- match(pro,x@depth)
					if (any(is.na(pos)))
						stop(cat("Some depths were not found:", pro[is.na(pos)],"\n" ))
						#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees  :",pro[is.na(pos)],"\n" ))
					pro<-pos
				}
				#Else dans ce cas garde pro comme il etait						

				if (is(dim1,"GLnothing")) class(dim1)<-"missing"
				if (is(dim2,"GLnothing")) class(dim2)<-"missing"			
				if (is(dim3,"GLnothing")) class(dim3)<-"missing"			

				#Execution du subscript
				depth <- x@depth
				x <- unclass(x)
				if (!l2$named) 
					dimnames(x)<-NULL
				m <- getMethod("[")(x,dim1,dim2,dim3,pro,drop=drop) #getMethod("[","array")
				GLmulti(m,depth[pro],drop=drop)
		}
		else
			stop("You can only use one or four subscripts (exactly like a 4 dimension array)")
			#stop("Vous ne pouvez utiliser qu'un ou quatre subscript (exactement comme un array de dimension 4)")
})

#REMPLACEMENT
setReplaceMethod("[","GLmultiArray4",function(x,...,value)
	{
		#extrait pour certaine depth
		#GLOverGroup4 <- function(pro,dim1,dim2,dim3,...,named,abs)
		#list(pro=pro,dim1=dim1,dim2=dim2,param=p,named=n,abs=ab)	
		l2 <- GLOverGroup4(...)
		pro <- l2$pro
		dim1 <- l2$dim1
		dim2 <- l2$dim2
		dim3 <- l2$dim3

		if (l2$param==0) return(x)			
		
		if (l2$param==1 || l2$param==4)
		{
			#Cas particulier depth
			if (is(pro,"GLnothing"))	#Proband vide
				class(pro)<-"missing"
			else
				if (is.numeric(pro) && !l2$abs) #Si c'est pas valeur absolue
				{
					#Recherche comme prevu dans la liste de depth
					pos <- match(pro,x@depth)
					if (any(is.na(pos)))
						stop(cat("Some depths were not found:",pro[is.na(pos)],"\n"))
						#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees  :",pro[is.na(pos)],"\n" ))
					pro<-pos
				}
				#Else dans ce cas garde pro comme il etait						

				if (is(dim1,"GLnothing")) class(dim1)<-"missing"
				if (is(dim2,"GLnothing")) class(dim2)<-"missing"			
				if (is(dim3,"GLnothing")) class(dim3)<-"missing"			

				#Execution du subscript
				depth <- x@depth
				x <- as(x,"array")
				if (!l2$named) 
					dimnames(x)<-NULL
				x[dim1,dim2,dim3,pro]<-value
				GLmulti(x,depth)
		}
		else
			stop("You can only use one or four subscripts (exactly like a 4 dimension array)")
			#stop("Vous ne pouvez utiliser qu'un ou quatre subscript (exactement comme un array de dimension 4)")
	}
)

#OBJET GLmultiMatrix
#
#Objet qui represente un vecteur de matrice 
#le vecteur correspond au differente depth
#
#1-GLOverGroup3


setClass("GLmultiMatrix",representation(depth="integer"), contains="array")

#Fonction virtuelles associees
setMethod("depth","GLmultiMatrix",function(x) x@depth)
#setMethod("length","GLmultiMatrix",function(x) length(x@depth))
setMethod("Dim","GLmultiMatrix",function(object, ...) {z=dim(object);l=length(z);z[c(l,1:(l-1))]})

setMethod("show","GLmultiMatrix",function(object){
		depth 	<- object@depth
		object <- as(object,"array")
		ldim	<-	length(dim(object))
		lapply(1:length(depth),function(x,obj,ldim,depth){
		 	cat("depth : ",depth[x],"\n")
			if (ldim==4)
				show(obj[,,x])
			else
				if (ldim==3)
					show(obj[,,x])
				else
						if (ldim==2)
							show(obj[,x])
						else
							if (ldim==1)
								show(obj[x])
			cat("\n")
		},obj=object,ldim=ldim,depth=depth)
})

setMethod("summary","GLmultiMatrix",function(object,...){
		#Affiche les donnees concernant un vecteur de matrice
		cat(" GENLIB: Matrix on many depths","\n",
		 " Number of lines  : ",  dim(object)[1], "\n",
		 " Number of columns: ",  dim(object)[2], "\n",
		 " depth(s) included : ",sep = "")
		cat(object@depth,sep=",")
		cat("\n")
})

#METHODE D'EXTRACTION
setMethod("[","GLmultiMatrix",function(x,i,...,drop){
          #extrait pour certaine depth
          #GLOverGroup3 <- function(pro,dim1,dim2,...,named,abs)
          #list(pro=pro,dim1=dim1,dim2=dim2,param=p,named=n,abs=ab)
          l2 <- GLOverGroup3(i,...)
          pro <- l2$pro
          dim1 <- l2$dim1
          dim2 <- l2$dim2
          #print(paste(dim1,dim2,l2$param, length(x), nargs()))
          if (l2$param==0) return(x)
          
          if (l2$param==1 || l2$param==3){
               #Cas particulier depth
               if (is(pro,"GLnothing")) #Proband vide
                    class(pro)<-"missing"
               else
                    if (is.numeric(pro) && !l2$abs){ #Si c'est pas valeur absolue
                         #Recherche comme prevu dans la liste de depth
                         pos <- match(pro,x@depth)
                         if (any(is.na(pos)))
                              stop(paste("Some depths were not found:",pro[is.na(pos)],"\n" ), call.=F)
                              #stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees :",pro[is.na(pos)],"\n" ))
                         pro<-pos
                    }
                    #Else dans ce cas garde pro comme il etait

               if (is(dim1,"GLnothing")) class(dim1)<-"missing"
               if (is(dim2,"GLnothing")) class(dim2)<-"missing"

               #Execution du subscript
               depth <- x@depth
               x <- unclass(x)
               if (!l2$named) 
                    dimnames(x)<-NULL
               #print(paste(dim1,dim2,pro, drop))
               m <- selectMethod("[","array")(x,dim1,dim2,pro,drop=drop)
               GLmulti(m,depth[pro],drop=drop)
          }
          else
               stop("You can only use one or three subscripts (exactly like a 3 dimension array)")
               #stop("Vous ne pouvez utiliser qu'un ou trois subscript (exactement comme un 'array' de dimension 3)")
})

#METHODE DE REMPLACEMENT
setReplaceMethod("[","GLmultiMatrix",function(x,...,value){
		#extrait pour certaine depth
		#GLOverGroup3 <- function(pro,dim1,dim2,...,named,abs)
		#list(pro=pro,dim1=dim1,dim2=dim2,param=p,named=n,abs=ab)	
		l2 <- GLOverGroup3(...)
		pro <- l2$pro
		dim1 <- l2$dim1
		dim2 <- l2$dim2

		if (l2$param==0) return(x)			
		
		if (l2$param==1 || l2$param==3){
			#Cas particulier depth
			if (is(pro,"GLnothing"))	#Proband vide
				class(pro)<-"missing"	
			else
				if (is.numeric(pro) && !l2$abs){ #Si c'est pas valeur absolue
					#Recherche comme prevu dans la liste de depth
					pos <- match(pro,x@depth)
					if (any(is.na(pos)))
						stop(cat("Some depths were not found:",pro[is.na(pos)],"\n" ))
						#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees  :",pro[is.na(pos)],"\n" ))
					pro<-pos
				}
				#Else dans ce cas garde pro comme il etait						

				if (is(dim1,"GLnothing")) class(dim1)<-"missing"
				if (is(dim2,"GLnothing")) class(dim2)<-"missing"			

				#Execution du subscript
				depth <- x@depth
				x <- as(x,"array")
				if (!l2$named) 
					dimnames(x)<-NULL
				x[dim1,dim2,pro]<-value
				GLmulti(x,depth)
		}else
			stop("You can only use one or three subscripts (exactly like a 3 dimension array)")
			#stop("Vous ne pouvez utiliser qu'un ou trois subscript (exactement comme un array de dimension 3)")					
})
################################################
#OBJET GLmultiNumber
#
#Objet qui represente un vecteur de nombre
#le vecteur correspond au differente depth
#
#1-GLOverNumber1
#

setClass("GLmultiNumber",representation(depth="integer",.Names="character"), contains="numeric")

#Fonctions virtuelles associees
setMethod("depth","GLmultiNumber",function(x) x@depth)
#setMethod("length","GLmultiNumber",function(x) length(x@depth))
setMethod("Dim","GLmultiNumber",function(object, ...) {z=dim(object);l=length(z);z[c(l,1:(l-1))]})

setMethod("show","GLmultiNumber",function(object){
	#Affiche les donnees
	lapply(1:length(object@depth),function(x,obj){
	 	cat("depth : ",obj@depth[x],"\n")
			if(length(obj@.Names)>0) cat(obj@.Names[x],"\n")
			cat(as.numeric(obj)[x],"\n")
			cat("\n")
		},obj=object)
})

setMethod("summary","GLmultiNumber",function(object,...){
		#Affiche les donnees concernant un vecteur de nombre
		cat(" GENLIB: Number on many depths","\n",		 
		 " depth included : ",sep = "")
		cat(object@depth,sep=",")
		cat("\n")
})

#METHODE D'EXTRACTION
setMethod("[","GLmultiNumber", function(x,...,drop){
		#extrait pour certaine depth
		l2 <- GLOverNumber1(...)
		pro <- l2$pro
	
		if (l2$param==0) return(x)	
					
		if (l2$param==1){
			#Cas particulier depth
			if (is(pro,"GLnothing"))	#Proband vide
				class(pro)<-"missing"	
			else
				if (is.numeric(pro) && !l2$abs){ #Si c'est pas valeur absolue			
					#Recherche comme prevu dans la liste de depth
					pos <- match(pro,x@depth)
					if (any(is.na(pos)))
						stop(cat("Some depths were not found:",pro[is.na(pos)],"\n" ))
						#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees  :",pro[is.na(pos)],"\n" ))
					pro<-pos
				}
				#Else dans ce cas garde pro comme il etait						

				#Execution du subscript
				depth <- x@depth
				m <- getMethod("[")(x,pro,drop=drop) #getMethod("[","numeric")
				if (length(x@.Names)>0 && l2$named)	names(m)<- x@.Names[pro]				
				GLmulti(m,depth[pro],drop=drop)
		}else
			stop("You can only use one subscript (exactly like a vector)")
			#stop("Vous ne pouvez utiliser qu'un subscript (exactement comme un vecteur)")							
})

#METHODE DE REMPLACMENT
setReplaceMethod("[","GLmultiNumber", function(x,...,value){
		#extrait pour certaine depth
		l2 <- GLOverNumber1(...)
		pro <- l2$pro
	
		if (l2$param==0) return(x)	
					
		if (l2$param==1){
			#Cas particulier depth
			if (is(pro,"GLnothing"))	#Proband vide
				class(pro)<-"missing"	
			else if (is.numeric(pro) && !l2$abs){ #Si c'est pas valeur absolue				
					#Recherche comme prevu dans la liste de depth
					pos <- match(pro,x@depth)
					if (any(is.na(pos)))
						stop(cat("Some depths were not found:",pro[is.na(pos)],"\n" ))
						#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees  :",pro[is.na(pos)],"\n" ))
					pro<-pos
				}
				#Else dans ce cas garde pro comme il etait						

				#Execution du subscript
				depth <- x@depth
				x<-as(x,"numeric")
				x[pro]<-value
				#adapte les nom au nouvelle valeur				
				GLmulti(x,depth)
		}else
			stop("You can only use one subscript (exactly like a vector)")
			#stop("Vous ne pouvez utiliser qu'un subscript (exactement comme un vecteur)")							
})

setMethod("Ops","GLmultiNumber",function(e1,e2=NULL){
	GLmulti(callGeneric(as(e1,"numeric"),e2),e1@depth)
})
################################################
#OBJET GLmultiVector
#
#Objet qui represente un vecteur de matrice 
#le vecteur correspond au differente depth
#
#1-GLOverVector2

setClass("GLmultiVector",representation(depth="integer"), contains="matrix")

#Fonctions virtuelles associees
setMethod("depth","GLmultiVector",function(x) x@depth)
#setMethod("length","GLmultiVector",function(x) length(x@depth))
setMethod("Dim","GLmultiVector",function(object, ...) {z=dim(object);l=length(z);z[c(l,1:(l-1))]})

setMethod("show","GLmultiVector",function(object){
		depth 	<- object@depth
		object <- as(object,"array")
		ldim	<-	length(dim(object))
		lapply(1:length(depth),function(x,obj,ldim,depth){
		 	cat("depth : ",depth[x],"\n")
			if (ldim==3)
				show(obj[,,x])
				else
					if (ldim==2)
						show(obj[,x])
						else
							if (ldim==1)
								show(obj[x])
			cat("\n")
		},obj=object,ldim=ldim,depth=depth)		
})

setMethod("summary","GLmultiVector",function(object,...){
		#Affiche les donnees concernant un vecteur de vecteur
		cat(" GENLIB: Vector on many depths","\n",
		 " Number of elements  : ",  dim(object)[1], "\n",
		 " depth included : ",sep = "")
		cat(object@depth,sep=",")
		cat("\n")
})

#METHODE D'EXTRACTION
setMethod("[","GLmultiVector",function(x,...,drop){
		#extrait pour certaine depth
		#GLOverVector2 <- function(pro,dim1,...,named,abs)
		#list(pro=pro,dim1=dim1,param=p,named=n,abs=ab)	
		l2 <- GLOverVector2(...)
		pro <- l2$pro
		dim1 <- l2$dim1		
		
		if (l2$param==0) return(x)	
					
		if (l2$param==1 || l2$param==2)		{
			#Cas particulier depth
			if (is(pro,"GLnothing"))	#Proband vide
				class(pro)<-"missing"	
			else
				if (is.numeric(pro) && !l2$abs){ #Si c'est pas valeur absolue						
					#Recherche comme prevu dans la liste de depth
					pos <- match(pro,x@depth)
					if (any(is.na(pos)))
						stop(cat("Some depths were not found:",pro[is.na(pos)],"\n" ))
						#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees  :",pro[is.na(pos)],"\n" ))
					pro<-pos
				}
				#Else dans ce cas garde pro comme il etait						

				if (is(dim1,"GLnothing")) class(dim1)<-"missing"

				#Execution du subscript
				depth <- x@depth
				x <- unclass(x)
				if (!l2$named) 
					dimnames(x)<-NULL
				m <- getMethod("[")(x,dim1,pro,drop=drop) #getMethod("[","matrix")
				GLmulti(m,depth[pro],drop=drop)
		}else
			stop("You can only use one or two subscript (exactly like a matrix)")
			#stop("Vous ne pouvez utiliser qu'un ou deux subscript (exactement comme une matrice)")					
})

#METHODE DE REMPLACEMENT
setReplaceMethod("[","GLmultiVector",function(x,...,value){
		#extrait pour certaine depth
		#GLOverVector2 <- function(pro,dim1,...,named,abs)
		#list(pro=pro,dim1=dim1,param=p,named=n,abs=ab)	
		l2 <- GLOverVector2(...)
		pro <- l2$pro
		dim1 <- l2$dim1		
		
		if (l2$param==0) return(x)	
					
		if (l2$param==1 || l2$param==2){
			#Cas particulier depth
			if (is(pro,"GLnothing"))	#Proband vide
				class(pro)<-"missing"	
			else
				if (is.numeric(pro) && !l2$abs){ #Si c'est pas valeur absolue						
					#Recherche comme prevu dans la liste de depth
					pos <- match(pro,x@depth)
					if (any(is.na(pos)))
						stop(cat("Some depths were not found:",pro[is.na(pos)],"\n" ))
						#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees  :",pro[is.na(pos)],"\n" ))
					pro<-pos
				}
				#Else dans ce cas garde pro comme il etait						

				if (is(dim1,"GLnothing")) class(dim1)<-"missing"

				#Execution du subscript
				depth <- x@depth
				x <- as(x,"matrix")
				if (!l2$named) 
					dimnames(x)<-NULL
				x[dim1,pro]<-value
				GLmulti(x,depth)
		}else
			stop("You can only use one or two subscript (exactly like a matrix)")
			#stop("Vous ne pouvez utiliser qu'un ou deux subscript (exactement comme une matrice)")					
})
##########################
#GLnoone
GLnoone <- -999
class(GLnoone) <- "GLnothing"
##########################################################################################
#Definition des classes servant a contenir et gerer l'information du phi moyen par group
#
#1-GLmultiPhiGroupSingle
#2-GLmultiPhiGroup
#GLPhiGroup
#GLPrivExtPHISINGLE
#GLPrivExtPHI 
#GLapplyPhi
#GLOverGroup2


#class d'inspiration : setClass("GLmultiMatrix",representation("array",depth="integer"))
setClass("GLmultiPhiGroup",representation(group="GLgroup",grindex="list"), contains="GLmultiMatrix") #Plusieurs depths
setClass("GLmultiPhiGroupSingle",representation(group="GLgroup",grindex="list"), contains="matrix")  #Une depth

#FONCTION VIRTUELLE ASSOCIeE (jusqu'ici pas necessaire)

setMethod("group","GLmultiPhiGroup",function(x) x@group)

setMethod("group","GLmultiPhiGroupSingle",function(x) x@group)

setMethod("Dim","GLmultiPhiGroup",function(object, ...){	
	param <- list(...)
	if (!is.null(param$drop) && param$drop==F){	
		#Dans ce cas, group seulement
		z=dim(object)	
		c(z[length(z)],length(object@group))
	}else{
		#Dans ce cas, ces les groups par group
		z=dim(object)	
		c(z[length(z)],length(object@group),length(object@group))
	}
})

setMethod("Dim","GLmultiPhiGroupSingle",function(object, ...){
	param <- list(...)
	if (!is.null(param$drop) && param$drop==F){	
		#Dans ce cas, group seulement
		c(length(object@group))
	}else{	
		#Dans ce cas, ce sont les groups par group
		c(length(object@group),length(object@group))
	}
})

setMethod("show","GLmultiPhiGroup",function(object){
		depth 		<- object@depth
		xgrindex <- object@grindex
		xgroupe   <- object@group		
		object 	<- as(object,"array")		
		lapply(1:length(depth),function(x,obj,depth,xgroupe,xgrindex){
		 	cat("depth : ",depth[x],"\n")
		
			m <- GLapplyGroup(obj[,,x],xgrindex,gen.phiMean,check=0,named=F)
			dimnames(m) <- list(names(xgroupe),names(xgroupe))
			prmatrix(m)

			cat("\n")
		},obj=object,depth=depth,xgroupe=xgroupe,xgrindex=xgrindex)
})

setMethod("show","GLmultiPhiGroupSingle",function(object){
	#Affiche les donnees concernant une genealogie	
	m <- GLapplyGroup(as(object,"matrix",strict=T),object@grindex,gen.phiMean,check=0,named=F)
	dimnames(m) <- list(names(object@group),names(object@group))
	prmatrix(m)
})

#EXTRACTION
setMethod("[","GLmultiPhiGroupSingle",function(x,...,drop) GLPrivExtPHISINGLE(x=x,...,drop=drop) )

setMethod("Ops","GLmultiPhiGroupSingle",function(e1,e2=NULL){
	stop("\nYou can not use operations on a GLmultiPhiGroupSingle object\nUse the empty [] operator to convert into matrix. (ex: z[]+1)")
	#stop("\nVous ne pouvez pas faire d'operation sur un GLmultiPhiGroupSingle\nUtiliser l'operateur [] vide pour convertir en matrice EX: z[]+1")
})

setReplaceMethod("[","GLmultiPhiGroupSingle",function(x,...,value){
	stop("\nYou can not directly modify a GLmultiPhiGroupSingle object")
})

setMethod("Math","GLmultiPhiGroupSingle",function(x){
	callGeneric(x[])
})

setMethod("Math2","GLmultiPhiGroupSingle",function(x, digits){
	callGeneric(x[],digits)
})

setMethod("Summary","GLmultiPhiGroupSingle",function(x, ..., na.rm = F){
	callGeneric(x[], ..., na.rm)
})

#OBJET 2: Plusieurs depths
setMethod("[","GLmultiPhiGroup",function(x,...,drop) GLPrivExtPHI(x=x,...,drop=drop) )

#Operateur mathematique (Pas les operateurs car sa marche pas sur l'autre sens 1 + obj)
setMethod("Ops","GLmultiPhiGroup",function(e1,e2=NULL){
	stop("\nYou can not use operations on a GLmultiPhiGroup object\nUse the empty [] operator to convert into matrix. (ex: z[]+1)")
	#stop("\nVous ne pouvez pas faire d'operation sur un GLmultiPhiGroup\nUtiliser l'operateur [] vide pour convertir en matrice EX: z[]+1")
})

setReplaceMethod("[","GLmultiPhiGroup",function(x,...,value){
	stop("\nYou can not directly modify a GLmultiPhiGroup object")
})

setMethod("Math","GLmultiPhiGroup",function(x){
	callGeneric(x[])
})

setMethod("Math2","GLmultiPhiGroup",function(x, digits){
	callGeneric(x[],digits)
})
setMethod("Summary","GLmultiPhiGroup",function(x, ..., na.rm = F){
	callGeneric(x[], ..., na.rm)
})


# Classes ajoutees par JFL car inexistantes dans ce qu'on a recu de la version SPlus
# mais qui sont necessaire pour certaines fonctions STAT.

setClass("GLmultiList"		,representation(liste = "list"), contains="list")

setMethod("show", "GLmultiList", function(object){
	cat("\n")
	cat("GLmultiList\n\n")
	nbDeb <- as.numeric(strsplit(names(object)[1]," ")[[1]][2])
	nbgen <- paste(unlist(lapply(c(nbDeb:(nbDeb+length(object)-1)), function(i){i})), collapse=" ")
	cat("Generations:", nbgen,"\n\n")
	x<-lapply(c(1:length(object)), function(i){
		cat(names(object)[i],"\n")
		print(object[[i]])
		cat("\n")
	})
})

