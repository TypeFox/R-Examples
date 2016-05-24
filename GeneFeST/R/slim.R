BBB2 <- new.env()

GeneFeST <- function(input,GROUP=FALSE,nb.pilot=20,pilot.runtime=500,main.runtime=5000,type=1,only.pilots=FALSE,h.average.P=0.2,h.step.width=1,mcmc.diag=FALSE,h=TRUE){


discard <- 10 
BBB2$prior.odds   <- h.average.P  #prior.odds
BBB2$h.step.width <- h.step.width
       
name         <- input 
pauses       <- 0
modification <- FALSE


	if(is.character(name)){
	
		 ## CHECKING THE GROUPS
        	 if(GROUP[1]==FALSE){
            		stop("Please verify the SNP groups")
         	 }else{
		  #parse the GROUP file
			
			 GROUP.file    <- read.table(GROUP,sep="\t",colClasses=c("character","numeric","numeric"))
			 GROUP.names   <- GROUP.file[,1]
                         GROUP.regions <- as.matrix(GROUP.file[,2:3])
                         GROUP         <- apply(GROUP.regions,1,function(x){ret <- x[1]:x[2];return(x)})
		         GROUP         <- unlist(GROUP)



           		 BBB2$GROUP      <- GROUP
          		 BBB2$N.REGIONS  <- length(unique(GROUP))
         	 }



		myheader <- scan(name,what="",nlines=5)

		#read out loci number
		locnum <- as.integer(sub("\\[loci\\]=", "", myheader[1]))
		#read out population number
		popnum <- as.integer(sub("\\[populations\\]=", "", myheader[2]))



		#import file
		if(myheader[3] == "[pop]=1" )  {
			print("Format's fine")
		} else { stop("Something's wrong with the format")  }

		modus <<- 0 #0 dominant 1 codominant
		if(length(scan(name,what="",skip=5,nlines=1)) > 4 ){
			modus <-1
			print("You've put in a codominant set of data")
		} else {
			print("You've put in a dominant set of data, call the other function")
			break
		}

		 population <- import_file(modus, name, locnum, popnum) # extract dataset from file	
		 population <- lapply(population,function(x)            # delete redundant data
		{ 
		return(x[,-(1)])
		})
		

		hapcount <- population[[1]][1,2]
         		
        

        ### INPUT is an R-object
	} else {

		modus <- 1

		
		 popnum     <- length(name$LISTE)
		 locnum     <- nrow(name$LISTE[[1]])
		 population <- name$LISTE
                 BBB2$GROUP      <- name$FUNC         
                 BBB2$N.REGIONS  <- length(unique(BBB2$GROUP))

	
	}
	
	if(pauses==1){
		print("Import finished, press any key to start initialisation")
		readline()
	}
		
        nullzeilen   <- which(population[[1]][,1]==0)
        einhaplotype <- which(population[[1]][,2]==1)
	delete       <- unique(c(nullzeilen,einhaplotype))
	
       

	if(length(delete)>0){    
         BBB2$GROUP      <- GROUP[-delete] 
         BBB2$N.REGIONS  <- length(unique(GROUP))
         for(xx in 1:popnum){
           population[[xx]] <- population[[xx]][-delete,]
         }
	}

 	 BBB2$sample_size <- sapply(population,function(x){return(x[,1])})	       
         locnum           <- locnum - length(delete)
         BBB2$locnum      <- locnum    
	 BBB2$hapcount    <- population[[1]][,2] # ist ab sofort nen array # brauch nur eine Pop, weil gleich		
         BBB2$popnum      <- popnum  
	 BBB2$d_alpha     <- numeric(locnum)      
         BBB2$d_beta      <- rnorm(popnum,-2,1.8)       
	 BBB2$acc_alpha   <- numeric(locnum)
         BBB2$acc_beta    <- numeric(popnum)

	 BBB2$acc_freq_ancestral <- numeric(locnum)

	
	BBB2$freq_pop <- lapply(population,function(x) #delete redundant data
	{ 
		return(x[,-(1:2)])
	})
	    
	BBB2$summe <- matrix(0,locnum,dim(population[[1]])[2]-2)
	lapply(BBB2$freq_pop,function(x){BBB2$summe <- BBB2$summe + x})
	sums  <- BBB2$summe + 1
	
	#Dirichlet
        rueck       <- rep(NA,dim(population[[1]])[2]-2)
        diri        <- apply(sums,1,function(x){
		      x  <- x[!is.na(x)] 	
                      dd <- my_rdirichlet(1,x)
                      rueck[1:length(dd)] <- dd
                      return(rueck)

                      })
	
	BBB2$freq_locus <- t(diri)
	
	BBB2$m1_prior_alpha     <-  0    
	BBB2$m2_prior_alpha     <-  0    
	BBB2$sd_prior_alpha     <-  1    
        BBB2$sd_prior_beta      <-  1     
        BBB2$mean_prior_beta    <- -1     
        BBB2$e_ancestral 	<- rep(0.2,locnum)   
	BBB2$var_prop_alpha     <- rep(1,locnum)
	BBB2$var_prop_beta      <- rep(1,popnum) # ""
	nb_pilot_alpha          <- 0 # " "	
	BBB2$mean_alpha 	<- array(0,c(locnum)) 
	BBB2$var_alpha  	<- array(5,c(locnum)) 
	BBB2$m2         	<- array(0,c(locnum))
	
	BBB2$alpha_included     <- array(FALSE,locnum) 
	
	
	if(pauses==1){
		print("About to start the pilot runs, press any key to do so")
		readline()
	}

	
	p <- 0
	cat("Pilotruns\n")


        #test
        nb_pilot      <-  nb.pilot          
        pilot_runtime <-  pilot.runtime   
        e_f 	      <-  0.05

        BBB2$GLOBAL_INIT1  <- matrix(0,locnum,dim(BBB2$freq_locus)[2])
	BBB2$GLOBAL_INIT2  <- rep(NA,2)
        BBB2$GLOBAL_INIT3  <- numeric(locnum)

        BBB2$PILOT     <- TRUE

        PILOT_TOTAL   <-  nb_pilot * pilot_runtime 
	
       
## PROGRESS #########################
 progr <- progressBar()
#####################################

MCMC.matrix <- NULL
MCMC.list   <- list()


	for(k in 1:nb_pilot){

		if(mcmc.diag){
		MCMC.matrix <- NULL
		}

		for(i in 1:pilot_runtime){ 
	                 
	                
    			update_d_betaco()
    			
		
			if(type==1){				
			X_update_d_alphaco_i_new()
			}
			if(type==2){
			X_update_d_alphaco_i_new2()
			}
                        if(type==3){
			X_update_d_alphaco_i_new3()
			}
                        
                        update_freq_codominant()  
		                              		         
			if(2*(k+1)>=nb_pilot){ 
					BBB2$mean_alpha <- BBB2$mean_alpha + BBB2$d_alpha
					BBB2$m2         <- BBB2$m2 + BBB2$d_alpha^2   					
			}		
		
		if(mcmc.diag){

		 if(type==1){	
		 MCMC.matrix <- rbind(MCMC.matrix,tapply(BBB2$d_alpha,BBB2$GROUP,unique))
                 }

		 if(type==2){	
		 MCMC.matrix <- rbind(MCMC.matrix,BBB2$d_alpha)
                 }		

		}

		}#end of one pilot run 
		
		if(mcmc.diag){
		MCMC.list[[k]] <- as.mcmc(MCMC.matrix)	
		}

			if(2*(k+1)>=nb_pilot){ 
				nb_pilot_alpha <- nb_pilot_alpha + 1 
                           
     			}
						
      				 # alpha var updates                             
     				 check <- (BBB2$acc_alpha/pilot_runtime)>0.4
      				 BBB2$var_prop_alpha[check] <- BBB2$var_prop_alpha[check]*1.2
     				 check <- (BBB2$acc_alpha/pilot_runtime)<0.2
      				 BBB2$var_prop_alpha[check] <- BBB2$var_prop_alpha[check]/1.2
		  		 
				  
     				 # var beta updates
      				 check <- BBB2$acc_beta/(pilot_runtime)>0.4
     				 BBB2$var_prop_beta[check] <- BBB2$var_prop_beta[check]*1.1
     				 check <- BBB2$acc_beta/(pilot_runtime)<0.2
     				 BBB2$var_prop_beta[check] <- BBB2$var_prop_beta[check]/1.1
		  		 
      				 
				 # frquencies var update
      				 check <- (BBB2$acc_freq_ancestral/pilot_runtime)>0.4
     				 BBB2$e_ancestral[check] <- BBB2$e_ancestral[check]*1.2
      				 check <- (BBB2$acc_freq_ancestral/pilot_runtime)<0.2
      				 BBB2$e_ancestral[check] <- BBB2$e_ancestral[check]/1.2
				 	

      				 BBB2$acc_alpha          <- numeric(locnum)
      				 BBB2$acc_beta  	 <- numeric(popnum)
      				 BBB2$acc_freq_ancestral <- numeric(locnum)
			
		
  # PROGRESS #######################################################
    progr <- progressBar(k,nb_pilot, progr)
  ###################################################################
	
}# end of pilot runs

if(mcmc.diag){
  return(as.mcmc(MCMC.list))	
}

        #################################################################
        # end of pilot runs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#################################################################


          BBB2$mean_alpha <- BBB2$mean_alpha/(nb_pilot_alpha*pilot_runtime) 
          BBB2$var_alpha  <- BBB2$m2/(nb_pilot_alpha*pilot_runtime)  - BBB2$mean_alpha^2 

output                 <- new("BAYESRETURN")

# TYPE 1
if(type==1 || type==3){
output@pilot.alpha     <- as.numeric(tapply(as.numeric(BBB2$mean_alpha),BBB2$GROUP,unique))
output@pilot.beta      <- as.numeric(BBB2$d_beta)
output@pilot.var_alpha <- as.numeric(tapply(as.numeric(BBB2$var_alpha),BBB2$GROUP,unique))
# calculate FST
fst                   <- outer(output@pilot.alpha,output@pilot.beta,"+")
fst                   <- 1/(1+exp(-fst))
fst                   <- rowSums(fst) 
fst                   <- fst/length(output@pilot.beta)
output@pilot.fst      <- fst 
output@pilot.P.norm   <- 1 - dnorm(output@pilot.alpha,mean(output@pilot.alpha),sd(output@pilot.alpha))
output@pilot.Q.norm   <- 1 - calcQ(output@pilot.P.norm)

#if(mcmc.diag){
#output@mcmc.diag    <- as.mcmc(MCMC.list)
#}

}

#TYPE 2
if(type==2){
num                    <- 1:length(BBB2$mean_alpha)
extrem.ids             <- as.numeric(tapply(num,BBB2$GROUP,function(x){
			  extrem    <- which.max(abs(BBB2$mean_alpha[x]));
                          extrem    <- x[extrem]
                          return(extrem)
                          }))
output@pilot.alpha     <- as.numeric(BBB2$mean_alpha[extrem.ids])
output@pilot.beta      <- as.numeric(BBB2$d_beta)
output@pilot.var_alpha <- as.numeric(BBB2$var_alpha[extrem.ids])
# calculate FST
fst                   <- outer(output@pilot.alpha,output@pilot.beta,"+")
fst                   <- 1/(1+exp(-fst))
fst                   <- rowSums(fst) 
fst                   <- fst/length(output@pilot.beta)
output@pilot.fst      <- fst 
output@pilot.P.norm   <- 1 - dnorm(output@pilot.alpha,mean(output@pilot.alpha),sd(output@pilot.alpha))
output@pilot.Q.norm   <- 1 - calcQ(output@pilot.P.norm)
}

if(is.character(name)){
output@region.names <- GROUP.names
}

#if(mcmc.diag){
#output@mcmc.diag    <- as.mcmc(MCMC.list)
#}

if(only.pilots){
return(output)
}


       ################################
       ## MAIN                        #
       ################################
       
        BBB2$alpha_updates <- 0 
        BBB2$d_alpha       <- numeric(locnum) 
        BBB2$PILOT         <- FALSE
	BBB2$acc.ratio     <- 0.25
	

	#  MAIN LOOP MCMC loop
	main_runtime <- main.runtime   
	post_alpha   <- numeric(locnum)
	post_fst     <- numeric(locnum)
	nb_alpha     <- numeric(locnum)
	cur_fst      <- numeric(locnum)
	cur_out      <- 0
        
        cat("\n")
	cat("Main loop\n")
  
## PROGRESS #########################
 progr <- progressBar()
#####################################

### Approx
BBB2$acc.ratio <- 0
BBB2$tempting  <- 1
###


for(i in 1:main_runtime){
    

cat(i,"of",main_runtime)
cat("\n")
 
### Approx #################################
BBB2$acc.ratio   <- BBB2$acc.ratio + sum(tapply(BBB2$alpha_included,BBB2$GROUP,unique))/BBB2$N.REGIONS
if(h){
if(i%%100 == 0){
BBB2$acc.ratio <- BBB2$acc.ratio/100
print("acc.ratio")
print(BBB2$acc.ratio)
driveP()
print("tempering")
print(BBB2$tempting)
}
}else{
BBB2$tempting  <- h.average.P
}
#### ######################################

    update_d_betaco()
   
    if(type==1){
     X_update_d_alphaco_i_new()
    }
    if(type==2){
     X_update_d_alphaco_i_new2()
    }
    if(type==3){
     X_update_d_alphaco_i_new3()
    }	
     
    update_freq_codominant()

    if(type==1){		
    X_jump_model_codominant()
    #X_jump_model_codominant3()
    }
    if(type==2){
    X_jump_model_codominant2()
    }	
    if(type==3){
    X_jump_model_codominant3()
    }	
    
   
			        
               if(i>discard){ 
               			 
                cur_out <- cur_out + 1
               ## Output

                   drinne                        <- which(BBB2$alpha_included)
 	           if(length(drinne)>0){
 	                   nb_alpha[drinne]      <- nb_alpha[drinne] + 1
                   }
                   
                     post_alpha            <- post_alpha + BBB2$d_alpha 
                     val                   <- outer(BBB2$d_alpha,BBB2$d_beta,"+")
                     val                   <- 1/(1+exp(-val))
                     val                   <- rowSums(val) 
                     cur_fst               <- val/popnum
                     post_fst              <- post_fst  + cur_fst
             }

   # PROGRESS #######################################################
    progr <- progressBar(i,main_runtime, progr)
   ###################################################################
         
  
	} # End of Main Loop
	
        outalpha      <- post_alpha/cur_out
        outfst        <- post_fst/cur_out
        outnb_alpha   <- nb_alpha/cur_out
	names(outfst) <- rownames(population[[1]])
        
	# TYPE 1
	if(type==1||type==3){
	output@post.alpha     <- as.numeric(tapply(as.numeric(outalpha),BBB2$GROUP,unique))
	output@post.beta      <- as.numeric(BBB2$d_beta)
	output@post.fst       <- as.numeric(tapply(outfst,BBB2$GROUP,unique))
        output@post.P         <- as.numeric(tapply(as.numeric(outnb_alpha),BBB2$GROUP,unique))
        output@post.Q         <- 1-calcQ(output@post.P)
        }
	# TYPE 2
	if(type==2){
	num            <- 1:length(outalpha)
        extrem.ids     <- as.numeric(tapply(num,BBB2$GROUP,function(x){
			  extrem    <- which.max(abs(outalpha[x]));
                          extrem    <- x[extrem]
                          return(extrem)
                          }))
	output@post.alpha     <- as.numeric(outalpha[extrem.ids])
        output@post.beta      <- as.numeric(BBB2$d_beta)
	output@post.fst       <- as.numeric(outfst[extrem.ids])
        output@post.P         <- as.numeric(tapply(as.numeric(outnb_alpha),BBB2$GROUP,unique))
        output@post.Q         <- 1-calcQ(output@post.P)
	}
		
	return(output) 
}


