ToProcess <-
function(N,n,parMode,pI,eq){

#require (MASS)
cat(paste("\nInputs " , N," ",n ," ",  parMode,"  ",pI, "  ",eq,"\n")) 
if ( is.na(N) || is.na(n))
cat ("N and n values must be an integer greater than zero (try again)")
else 
   if (is.numeric(N)==FALSE || is.numeric(n)==FALSE) 
   cat ("N and n values must be an integer greater than zero (try again)")
   else
   if (length(pI)==0)
   cat ("Must select an indicador for each block variables")  
   else 
   {
   #Input retrieve
   #------------------------------------------------#
   if (is.na(pI)) pI<- 2
   if (eq==1) eq<- TRUE
   else eq<-FALSE

   #Setting parameters to functions
   #------------------------------------------------#
   r.s <- matrix(c(0,0,0,1,
		   0,0,0,1,
		   0,0,0,1,
		   1,1,1,0),4,4,byrow=TRUE)

   r.ie <- matrix(c(0,0,1,			
		    0,1,0,
		    0,0,0),3,3,byrow=TRUE)
		    
   modo <- c(parMode,parMode,parMode,parMode)
   bv <- c(pI,pI,pI,pI)
   path.coef <- c(0.5,0.4,0.3,0,0.3,0,0.3,0,0)
   #------------------------------------------------#  
   if (eq==TRUE) #equal
   {
    if (pI==2){ #Case 1: two indicators per construct, all outer weights equal
    wei.ex <- matrix (c(0.63,0.63,
                        0.63,0.63,
                        0.63,0.63), 3,2,byrow=TRUE)#cambio 07.11
    wei.en <- matrix(c(0.63,0.63),1,2,byrow=TRUE)  
    }
    if (pI==4){#Case 2: four indicators per construct, all outer weights equal
    wei.ex <- matrix (c(0.42,0.42,0.42,0.42,
                        0.42,0.42,0.42,0.42,
                        0.42,0.42,0.42,0.42), 3,4,byrow=TRUE)
    wei.en <- matrix(c(0.42,0.42,0.42,0.42),1,4,byrow=TRUE)
    }
    if (pI==6){#Case 3: six indicators per construct, all outer weights equal
    wei.ex <- matrix (c(0.35,0.35,0.35,0.35,0.35,0.35,
                                  0.35,0.35,0.35,0.35,0.35,0.35,
                                  0.35,0.35,0.35,0.35,0.35,0.35),             
                                  3,6,byrow=TRUE)
    wei.en <- matrix(c(0.35,0.35,0.35,0.35,0.35,0.35),1,6,byrow=TRUE)
    }
    if (pI==8){#Case 4: eight indicators per construct, all outer weights equal
    wei.ex <- matrix (c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,
                        0.3,0.3,0.3,0.3,0.3,0.3, 0.3,0.3,
                        0.3,0.3,0.3,0.3,0.3,0.3, 0.3,0.3),             
                        3,8,byrow=TRUE)
    wei.en <- matrix(c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3),1,8,byrow=TRUE)
    }
    
   }else{#Different
   
    if (pI==2) {#Case 1: two indicators per construct, all differents
    wei.ex <- matrix(c(0.8,0.4,
                      0.4,0.8,
                      0.1,0.9),3,2,byrow=TRUE)
    wei.en <- matrix(c(0.4,0.8),1,2,byrow=TRUE)
    }
    if (pI==4){#Case 2: four indicators per construct, all differents
    wei.ex <- matrix(c(0.2,0.3,0.5,0.7,
		                   0.2,0.4,0.6,0.5,
                       0.3,0.5,0.7,0.2),3,4,byrow=TRUE)
    wei.en <- matrix(c(0.2,0.3,0.5,0.5),1,4,byrow=TRUE)
   
    }
    if (pI==6){#Case 3: six indicators per construct, all differents
     wei.ex <- matrix(c(0.5,0.3,0.4,0.3,0.5,0.1,
	                   0.2,0.4,0.6,0.4,0.2,0.3,
	                   0.3,0.6,0.2,0.3,0.4,0.2),
                                 3,6,byrow=TRUE)
     wei.en <- matrix(c(0.5,0.3,0.4,0.3,0.5,0.1),1,6,byrow=TRUE)
     }
    if (pI==8){#Case 4: eight indicators per construct, all differents
     wei.ex <- matrix(c(0.3,0.3,0.4,0.3,0.4,0.3,0.2,0.3,
		                   0.3,0.3,0.4,0.3,0.2,0.3,0.4,0.2,
                       0.4,0.5,0.4,0.3,0.2,0.1,0.3,0.2),
                       3,8,byrow=TRUE)
     wei.en <- matrix(c(0.3,0.3,0.4,0.3,0.4,0.3,0.2,0.3),1,8,byrow=TRUE)
    }
}

#Call to functions simulators
#-------------------------------------------------------------------#
intpar <- IntPar(r.s,r.ie,modo,bv)
xex <- ExoMVs(N,n,intpar$ind.ex)
yex <- ExoLVs(N,n,bv,intpar$bex,xex$x.ex,wei.ex) 		
nlie <- NIEffects(N,n,yex$y.ex,intpar$rie)
yexcor <- ExoLVsCor(N,n,intpar$bex,intpar$rie,yex$y.ex,nlie$a.nle,nlie$a.ie)
err <- ErrEnLV(N,n,intpar$ben,path.coef,yexcor$y.ex.cor)

if (as.numeric(err$vis) == 0)
 {		
 yen <- EnLVs(N,n,intpar$ben,path.coef,err$elv,yexcor$y.ex.tot)
 xen <- EnMVs(N,n,intpar$ind.en,wei.en,yen$y.en)
 x <- XexXen(N,n,intpar$ind.ex,intpar$ind.en,xex$x.ex,xen$x.en)
 x <- x$x

#Output to text file and csv file
#-------------------------------------------------------------------#
outputfile <- tclvalue(tkgetSaveFile(initialdir="~"))
   if (outputfile != "") 
   { 
    cat ("Generating data...")
    write.table(x, file=paste(outputfile,".txt"), sep=",", col.names=TRUE, row.names=FALSE)
    write.csv(x, file = paste(outputfile,".csv"))
    msg<-paste("csv and txt data were successfully saved ")
    }else msg<-paste("The data were not successfully saved")
  tkmessageBox(title="Data Output",message=msg)
  }
 } 
}#Fin ToProcess

