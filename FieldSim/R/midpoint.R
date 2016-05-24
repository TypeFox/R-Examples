### midpoint.R  (2006-16)
###
###   Fractional Brownian field simulation by the midpoint displacement method
###
### Copyright 2006-16 Alexandre Brouste and Sophie Lambert-Lacroix
###
###

midpoint <- function(process){

##    INPUT VARIABLES
#########################
##  process   : an object of class process (namely, a fBm with finer or visualization grid)

##    OUTPUT VARIABLES
##########################
## process with the simulation in the corresponding slot values

if(missing(process)){
	cat("Error from midpoint.R: parameter process is missing\n")
	return(NULL)
}

if(!isS4(process)){
	cat("Error from midpoint.R: parameter process is not of type process\n")
	return(NULL)
}else if(!(class(process)[1]=="process")){
	cat("Error from midpoint.R: parameter process is not of type process\n")
	return(NULL)
}

name<-process@name
nameofgrid<-whichgrid(process@manifold)

if (name=="fBm"){
    
    nameofmanifold<-process@manifold@name
    
    if (nameofmanifold=="line"){
    cat("Error from midpoint.R: mdipoint for this process is not yet available\n")
	return(NULL)
    }else if (nameofmanifold=="plane"){
    
    if (nameofgrid=="visualization"){
        
        N<-sqrt(length(process@manifold@atlas))
        nblevel<-log(N-1)/log(2)
    
        H<-process@parameter
        
        ##  TEST ON INPUT VARIABLES
        ##############################
        
        #if ((is.numeric(H)==FALSE)||(H>=1)||(H<=0)){
        #    stop("Message from midpoint.R: H is not of valid type")}
        
        #if ((is.numeric(nblevel)==FALSE)||(round(nblevel)-nblevel!=0)||(nblevel<0)){
        #    stop("Message from midpoint.R: nblevel is not of valid type")}
        
        Z<-matrix(0,2,2)
        Z[1,2] <- rnorm(1)
        Z[2,1] <- rnorm(1)
        Z[2,2] <- rnorm(1)*sqrt(2^H)
        
        level<-1
        while (level<=nblevel){
            Y<- matrix(0,2^(level)+1,2^(level)+1)
            
            
            #--------Build Y----------------------
            
            for (l in 1:(2^(level)+1)){ #les colonnes
                for (m in 1:(2^(level)+1)){ #les lignes
                    
                    if (((m/2-floor(m/2))!=0) & ((l/2-floor(l/2))!=0)){
                        Y[m,l]<-Z[((m-1)/2+1),((l-1)/2+1)]
                    }
                    
                }}
            
            #--------------Simulation of the other points --------------
            
            #////////variance of the center point  //////
            varC<-(1-1/4*2^(H)-1/8*2^(2*H))*2^(-2*level*H+H)
            
            #////////variance of the remaining points //////
            varA<-(1-2^(2*H-2))/2^(2*level*H)
            
            
            for (m in 1:2^(level-1)){
                for (l in 1:2^(level-1)){
                    
                    pc_x<-2*l
                    pc_y<-2*m
                    
                    
                    ###################################
                    #Simulation of the center point
                    ###################################
                    
                    Y[pc_x,pc_y]<-sqrt(varC)*rnorm(1) + 1/4*(Y[pc_x-1,pc_y-1]+Y[pc_x-1,pc_y+1]+Y[pc_x+1,pc_y+1]+Y[pc_x+1,pc_y-1])
                    
                    ##############################################
                    #Simulation for the point at right
                    ##############################################
                    
                    Y[pc_x+1,pc_y]<-sqrt(varA)*rnorm(1) + 1/2*(Y[pc_x+1,pc_y-1]+Y[pc_x+1,pc_y+1])
                    
                    ##############################################
                    #Simulation for the point above
                    ##############################################
                    
                    Y[pc_x,pc_y+1]<-sqrt(varA)*rnorm(1) + 1/2*(Y[pc_x-1,pc_y+1]+Y[pc_x+1,pc_y+1])
                    
                    ##############################################
                    #Simulation for the point below
                    #only for m=1
                    ##############################################
                    if (m==1){
                        
                        Y[pc_x,pc_y-1]<-sqrt(varA)*rnorm(1) + 1/2*(Y[pc_x-1,pc_y-1]+Y[pc_x+1,pc_y-1])
                        
                    }
                    
                    ##############################################
                    #Simulation for the point at left 
                    #only for l=1
                    ##############################################
                    if (l==1){
                        
                        Y[pc_x-1,pc_y]<-sqrt(varA)*rnorm(1) + 1/2*(Y[pc_x-1,pc_y-1]+Y[pc_x-1,pc_y+1])
                        
                    }
                    
                }
            }
            
            level<-level+1
            Z<-Y
        }
        
        nameProcess<-deparse(substitute(process))
        process@values<-as.vector(Z)
        assign (nameProcess,process,envir = parent.frame())
        return(invisible(1))
        
    }else{
       
       cat("Error from midpoint.R: mdipoint for this process is not yet available\n")
       return(NULL)
    }
    
        
    }else{
        
        cat("Error from midpoint.R: mdipoint for this process is not yet available\n")
        return(NULL)
        
    }
    
    
    
}else{
    cat("Error from midpoint.R: parameter process is not of type process\n")
	return(NULL)
}





}
