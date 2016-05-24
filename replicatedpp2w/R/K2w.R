K2w <-
function(pplist=NULL, dataKijk=NULL,nijk=NULL,  r, r0=NULL, rmax=NULL,  tratA, tratB=NULL,wt=NULL, nsim=999,correction="trans",...)
{
       
      if(!is.null(pplist)){
         dataKijk<- sapply(pplist, function(x) Kest(x, r=r,correction=correction,...)[[correction]])
          nijk <- sapply(pplist, function(x) x$n)
      }
      if(is.null(r0)) r0<-min(r)
      if(is.null(rmax)) rmax<-max(r)
      
        nameA<-deparse(substitute(tratA))
	tratA <- factor(tratA)
	l.tratA <- levels(tratA)
	nsumA <- tapply(nijk,tratA,sum)
       
  
  #null objects to store results in the one-way case
   nsumB <- NULL; nsumAB <- NULL; tratAB <- NULL
   KrepB <-NULL; KrepAB <-NULL
   K0j <-  NULL; K0ij <-  NULL
  nameB<-NULL
  
  
if(!is.null(tratB)){         ##### in case of TWO WAY
        nameB<-deparse(substitute(tratB) )
	tratB <- factor(tratB)
	l.tratB <- levels(tratB)
        

	#combinacion de factores / combination of factors
	tratAB <- factor(apply(cbind(as.character(tratA),as.character(tratB)),1,paste, collapse=""))
        l.tratAB <- levels(tratAB)
	
	nsumB <- tapply(nijk,tratB,sum)
	nsumAB <- tapply(nijk,tratAB,sum)
	}                          #####  end in case of TWO WAY
	
	nsum<- sum(nsumA) # numero total de puntos / total number of points

	###
	###  Medias por level para cada tratamiento y para la combinacion (interaccion)
	###  Weighted averages per level for each  tratment and for each combination (interaction)
	###
	
	# primero hacemos un apply para multiplicar todas las filas de Kijk por nijk (cada Kijk por su nijk)
	# a la data.frame resultante le hacemos otro apply (por columnas porque esta transpuesta)
	# que a cada columna le aplique con tapply la suma segun los levels del factor "trata"
	# por ultimo, dividimos cada una de las sumas por el nsum correspondiente
	# Es decir, se trata del weighted average (formula 8.11) de Diggle 2003: 123, i.e.
	# K_hat(t) = Sigma ni * Ki_hat(t) / Sigma ni
	
	KrepA <- t(apply(apply((apply(dataKijk, 1, function(x) x*nijk)),2,tapply, tratA, sum, na.rm=TRUE),2, function(x) x/nsumA))
  
  
	if(!is.null(tratB)){
		KrepB <- t(apply(apply((apply(dataKijk, 1, function(x) x*nijk)),2,tapply, tratB, sum, na.rm=TRUE),2, function(x) x/nsumB))
		KrepAB <- t(apply(apply((apply(dataKijk, 1, function(x) x*nijk)),2,tapply, tratAB, sum, na.rm=TRUE),2, function(x) x/nsumAB))
	}
	
	###
	### Gran media, i.e., weigthed average of Ki_hat, with weights proportional to the total number
	###                        of events in each group, i.e., ni = Sigma nij (Diggle 2003: 125). [ni= nsumA]
	
	# Todas las K0 son la misma !!
	K0i <-  apply(apply(KrepA, 1, function(x) x*nsumA),2,sum, na.rm=TRUE)/nsum
	
	if(!is.null(tratB)){         ##### TWO WAY
		K0j <-  apply(apply(KrepB, 1, function(x) x*nsumB),2,sum, na.rm=TRUE)/nsum
		K0ij <-  apply(apply(KrepAB, 1, function(x) x*nsumAB),2,sum, na.rm=TRUE)/nsum
	}                            ##### TWO WAY


	
	###
	### Calculo del sumatorio de cuadrados // Computing sum of squares
	###
		
	ok <- (r >= r0 & r <= rmax)        # distancias sobre las que se va a integrar //  r distances for which summation will be made
	if(r[ok][1]==0){
                ok[1]<- FALSE
		# la primera distancia la eliminamos por que 0^-2 =Inf // first distance is not analyzed
                warning("first distance not anlyzed: 0^-2 =Inf", call. =FALSE)  
	}
   wt<-eval(wt)
  if(is.null(wt)) wt <-  r^-2        #weight de Diggle 2003:127 para funciones K  // weight for K functions (Diggle 2003:127 )
	
  btss.i <-NULL
  btss.j <-NULL
  btss.ij <-NULL
  
	# Factor A
	integrand.i <- apply(KrepA, 2, function(x) wt*(x-K0i)^2)
	integrand.i <- apply(integrand.i[ok,], 2, sum, na.rm=TRUE)               # sumamos solo en el rango de r // sum only within the selected r range 
	btss.i <- sum(integrand.i*nsumA)
	
	if(!is.null(tratB)){         ##### TWO WAY
	#Factor B
	integrand.j <- apply(KrepB, 2, function(x) wt*(x-K0j)^2)
	integrand.j <- apply(integrand.j[ok,], 2, sum, na.rm=TRUE)               # sumamos solo en el rango de r
	btss.j <- sum(integrand.j*nsumB)
	
		
	# Interacction AB
	integrand.ij <- NULL
	for (i in 1: length(l.tratA)){
		for(j in 1: length(l.tratB)){
			chose<-paste(l.tratA[i],l.tratB[j],sep="",colapse="")==names(data.frame(KrepAB))
			integrand.ij <- cbind(integrand.ij, wt*(KrepAB[,chose]-KrepA[,i]-KrepB[,j]+K0ij)^2)
		}
	}
	integrand.ij <- apply(integrand.ij[ok,], 2, sum, na.rm=TRUE) # sumamos solo en el rango de r
	btss.ij <- sum(integrand.ij*nsumAB)
  }                           ##### TWO WAY  

	#################	
	### Calculo de los residuales // Computing residuals
	#################
	
		
	Rik <- Rjk <- Rijk <- dataKijk *0 # matrices de ceros
	
	# Residuos para A // residuals for A
	for (i in 1: length(l.tratA)){
		Rik[,tratA==l.tratA[i]] <- (apply(dataKijk[,tratA==l.tratA[i]],2, function(x) (x-KrepA[,i])))
	}
	    # idem para B y AB si el analisis es 2way // idem for B and AB in the 2way case
	    if(!is.null(tratB)){         ##### TWO WAY
	          # Residuales para B
	          for (i in 1: length(l.tratB)){
	 	       Rjk[,tratB==l.tratB[i]] <- (apply(dataKijk[,tratB==l.tratB[i]],2, function(x) (x-KrepB[,i])))
		  }
	          # Residuales para AB	
	          # definimos el residual como Kij-Ki.- K.j + K0 //  Residuals are defined as Kij-Ki.- K.j + K0
	         for (i in 1: length(l.tratA)){
		     for(j in 1: length(l.tratB)){
		         Rijk[,tratA==l.tratA[i]& tratB==l.tratB[j]] <- (apply(dataKijk[,tratA==l.tratA[i]& tratB==l.tratB[j]],2, function(x) (x-KrepA[,i]-KrepB[,j]+K0ij)))
		     }
	         }
	    }                           ##### TWO WAY  
	
	# multiplicamos todos los residuales pr nijk^0.5 //  multiply all resduals by nijk^0.5
	Rik <- t(apply(Rik, 1, function(x) x*nijk^0.5))
	    # idem para B y AB si el analisis es 2way
   	    if(!is.null(tratB)){         ##### TWO WAY
		Rjk <- t(apply(Rjk, 1, function(x) x*nijk^0.5))
		Rijk <- t(apply(Rijk, 1, function(x) x*nijk^0.5))
	     }                            ##### TWO WAY  
        
	
	######################
	### Resampling
	#######################
	
	btss.i.res <- NULL
	btss.j.res <- NULL
	btss.ij.res <- NULL
	
	KA.res<-NULL   # tablas para almacenar las K medias resampleadas y obtener envueltas bootstrapeadas para dibujar
	KB.res<-NULL
	KAB.res<-NULL
  
       #declaracion de productos que seran nulos si el analisis es solo one way
       Rjk.res <-NULL
       Rijk.res <-NULL
       K0j.res <-  NULL
       K0ij.res <-  NULL
          
	#  Comienzan las simulaciones // start simulations
		
	nKi <- length(tratA) # cuantas replicas hay que sacar en cad simulacion (tantas como replicas originales) // number of replicates per simulation
	for (ns in 1:nsim){
		progressreport(ns,nsim)
		
		# definimos un vector de muestras aleatorias // define a vector of random samples
		resamp <- sample(nKi, nKi, replace=TRUE) # elegimos la muestra de residuales
		
		#seleccionamos los residuales aleatorios // random sample of residuals
		Rik.res <- Rik[,resamp]
		    # idem para B y AB si el analisis es 2way
		    if(!is.null(tratB)){         ##### TWO WAY
			Rjk.res <- Rjk[,resamp]
			Rijk.res <- Rijk[,resamp]
		    }
		    
		# los multiplicamos por nijk^-0.5
		Rik.res <- t(apply(Rik.res, 1, function(x) x*nijk^-0.5))
		     # idem para B y AB si el analisis es 2way
		     if(!is.null(tratB)){         ##### TWO WAY
			Rjk.res <- t(apply(Rjk.res, 1, function(x) x*nijk^-0.5))
			Rijk.res <- t(apply(Rijk.res, 1, function(x) x*nijk^-0.5))
		     }		
		
		# reconstruimos las nuevas Kijk, (un set distinto para cada test) // reconstruct simulated Kijk (a different set for each test)
		dataKik.res <- apply(Rik.res, 2, function(x) x+K0i)
		      # idem para B y AB si el analisis es 2way
		      if(!is.null(tratB)){         ##### TWO WAY
			 dataKjk.res <- apply(Rjk.res, 2, function(x) x+K0j)
			 
			 dataKijk.res <- Rijk.res + KrepA[,tratA] + KrepB[,tratB] -K0ij    
		      }
		      
        	# Obtenemos las nuevas medias para cada nivel de cada factor y de la interaccion // new weighted averages for each factor level and combination of factor levels
		KrepA.res <- t(apply(apply((apply(dataKik.res, 1, function(x) x*nijk)),2,tapply, tratA, sum, na.rm=TRUE),2, function(x) x/nsumA))
		      # idem para B y AB si el analisis es 2way
		      if(!is.null(tratB)){         ##### TWO WAY
			 KrepB.res <- t(apply(apply((apply(dataKjk.res, 1, function(x) x*nijk)),2,tapply, tratB, sum, na.rm=TRUE),2, function(x) x/nsumB))
			 KrepAB.res <- t(apply(apply((apply(dataKijk.res, 1, function(x) x*nijk)),2,tapply, tratAB, sum, na.rm=TRUE),2, function(x) x/nsumAB))
			   # medias resampleadas para restar de KrepAB.res
			   KrepAB.A.res <- t(apply(apply((apply(dataKijk.res, 1, function(x) x*nijk)),2,tapply, tratA, sum, na.rm=TRUE),2, function(x) x/nsumA))
			   KrepAB.B.res <- t(apply(apply((apply(dataKijk.res, 1, function(x) x*nijk)),2,tapply, tratB, sum, na.rm=TRUE),2, function(x) x/nsumB))
		      }
		      
		# Calculamos las nuevas medias globales //new global averages after resampling
		K0i.res <-  apply(apply(KrepA.res, 1, function(x) x*nsumA),2,sum, na.rm=TRUE)/nsum
		      # idem para B y AB si el analisis es 2way
		      if(!is.null(tratB)){         ##### TWO WAY
			  K0j.res <-  apply(apply(KrepB.res, 1, function(x) x*nsumB),2,sum, na.rm=TRUE)/nsum
			  K0ij.res <-  apply(apply(KrepAB.res, 1, function(x) x*nsumAB),2,sum, na.rm=TRUE)/nsum
		      }
		
		# calculamos la suma de cuadrados para cada factor e interaccion //  sum of squares
		integrand.i.res <- apply(KrepA.res, 2, function(x) wt*(x-K0i.res)^2)
		integrand.i.res<- apply(integrand.i.res[ok,], 2, sum, na.rm=TRUE)
		btss.i.res<- c(btss.i.res, sum(integrand.i.res*nsumA) )# este es el btss simuladoempirico
		      # idem para B y AB si el analisis es 2way
		      if(!is.null(tratB)){         ##### TWO WAY
			   integrand.j.res <- apply(KrepB.res, 2, function(x) wt*(x-K0j.res)^2)
			   integrand.j.res<- apply(integrand.j.res[ok,], 2, sum, na.rm=TRUE)
			   btss.j.res<- c(btss.j.res, sum(integrand.j.res*nsumB) ) # este es el btss simuladoempirico

			   integrand.ij.res <- NULL
			   for (i in 1: length(l.tratA)){
				for(j in 1: length(l.tratB)){
					chose<-paste(l.tratA[i],l.tratB[j],sep="",colapse="")==names(data.frame(KrepAB.res))
					integrand.ij.res <- cbind(integrand.ij.res, wt*(KrepAB.res[,chose]-KrepAB.A.res[,i]-KrepAB.B.res[,j]+K0ij.res)^2)
				}
			   }
			  integrand.ij.res<- apply(integrand.ij.res[ok,], 2, sum, na.rm=TRUE)
			  btss.ij.res<- c(btss.ij.res, sum(integrand.ij.res*nsumAB) ) # este es el btss simuladoempirico
                      }            
	   
                # Almacenamos los valores de las K medias por grupoo resampleadas // store average K per resampled group
		KA.res <- cbind(KA.res, KrepA.res)
		      # idem para B y AB si el analisis es 2way
		      if(!is.null(tratB)){         ##### TWO WAY
			   KB.res <- cbind(KB.res, KrepB.res)
			   KAB.res <- cbind(KAB.res,KrepAB.res)
		      }		
	}
	
	##############
	### Escribir resultados // write results
	##############
	
	result<- list(btss.i=btss.i, btss.i.res=btss.i.res,btss.ij=btss.ij, btss.ij.res=btss.ij.res,
                btss.j=btss.j, btss.j.res=btss.j.res, KrepA=KrepA, KrepB=KrepB, KrepAB=KrepAB,
                K0i=K0i, K0j=K0j, K0ij=K0ij, Rik=Rik, Rjk=Rjk, Rijk=Rijk, nsumA=nsumA, nsumB=nsumB, nsumAB=nsumAB,
                wt=wt, tratA=tratA, tratB=tratB, tratAB=tratAB, dataKijk=dataKijk, nijk=nijk, r=r, r0=r0,
             #   K0i.res=K0i.res, K0j.res=K0j.res, K0ij.res=K0ij.res,
		KA.res=KA.res, KB.res=KB.res, KAB.res=KAB.res,
               nameA=nameA, nameB=nameB)
	class(result)=c("k2w", class(result))
	return(result)
}
