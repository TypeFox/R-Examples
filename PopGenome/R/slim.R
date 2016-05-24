#BayeScan R slim version
#BAYESRETURN Krebs
#18.05.2011

#HEADER
#setwd("/home/berend/Desktop/ma_code")

#dependencies
#source("import_export_slim.R") #data import
#source("update_slim.R")
#source("math.R")
#source("BAYESRETURN.R")

#library("gregmisc") #provides rdirichlet

#pilot run statics
#m1_prior_alpha <- 0
#m2_prior_alpha <- 0
#sd_prior_alpha <- 1

#e_ancestral <- 0.2
#e_freq      <- 0.2

#var_prop_alpha <- 1.0
#var_prop_beta <- 0.7

# var_prop_a_p <- 0.2

# pilot runs
# nb_pilot      <<- 2 # for over pilot runs
# pilot_runtime <<- 100


BBB <- new.env()


BayeScanR <- function(input,nb.pilot=10,pilot.runtime=2500,main.runtime=100000,discard=50000){

	# modification <<- modification
        

name         <- input 
pauses       <- 0
modification <- FALSE


	# Muss nen neuer Parser her !!!
	if(is.character(name)){
	
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
		
        ### INPUT is a R-object
	} else {
		modus <- 1

		#hapcount   <<- name[[1]][,2] # ist ab sofort nen array # brauch nur eine Pop, weil gleich
		 popnum     <- length(name$LISTE)
		 locnum     <- nrow(name$LISTE[[1]])
		 population <- name$LISTE
             # mod  # GROUP      <<- name$FUNC
             # mod  # GROUP2     <<- name$FUNC2
             # mod  # N.REGIONS  <<- length(unique(GROUP))

	
	}
	
	if(pauses==1){
		print("Import finished, press any key to start initialisation")
		readline()
	}
		
        nullzeilen   <- which(population[[1]][,1]==0)
        einhaplotype <- which(population[[1]][,2]==1)
	delete       <- unique(c(nullzeilen,einhaplotype))
	
       

	if(length(delete)>0){
         # GROUP2 eventuell noch aendern ....
         # mod # GROUP      <<- GROUP[-delete] 
         # mod # N.REGIONS  <<- length(unique(GROUP))
         for(xx in 1:popnum){
           population[[xx]] <- population[[xx]][-delete,]
         }
	}
	
 	 BBB$sample_size <- sapply(population,function(x){return(x[,1])})
	
         # print(sample_size)

         locnum      <- locnum - length(delete)
         BBB$locnum  <- locnum
       

	 BBB$hapcount    <- population[[1]][,2] # ist ab sofort nen array # brauch nur eine Pop, weil gleich
	
    
	#INITIALISATION

         BBB$popnum     <- popnum
     
	 BBB$d_alpha    <- numeric(locnum)
        
         BBB$d_beta     <- rnorm(popnum,-2,1.8) # passt start.cpp Ln 225
        
	 BBB$acc_alpha  <- numeric(locnum)

         BBB$acc_beta   <- numeric(popnum)

	 BBB$acc_freq_ancestral <- numeric(locnum)

	#add up frequencies
	BBB$freq_pop <- lapply(population,function(x) #delete redundant data
	{ 
		return(x[,-(1:2)])
	})
	
        
	BBB$summe <- matrix(0,locnum,dim(population[[1]])[2]-2)
	lapply(BBB$freq_pop,function(x){BBB$summe <- BBB$summe + x})
	sums  <- BBB$summe + 1
	
 
	#Dirichlet
        rueck       <- rep(NA,dim(population[[1]])[2]-2)
        diri        <- apply(sums,1,function(x){
		      x  <- x[!is.na(x)] 	
                      dd <- my_rdirichlet(1,x)
                      rueck[1:length(dd)] <- dd
                      return(rueck)

                      })
	
	BBB$freq_locus <- t(diri)
	

	# log_likelihood <<- allelecount_loglikelihood() #
	
        
	
	#d_old_alpha    <<- d_alpha
	
	

	
	#alpha_back <<- list()
	
	#BwPrp <- 0
	#f     <- runif(1)  # 0
	
	BBB$m1_prior_alpha <- 0 # passt start.cpp Ln219
	BBB$m2_prior_alpha <- 0 # passt ""
	BBB$sd_prior_alpha <- 1 # passt ""
        BBB$sd_prior_beta          <-  1 #3
        BBB$mean_prior_beta        <- -1  #-5.15

#	e_ancestral <<- 0.2 # ""
        BBB$e_ancestral <- rep(0.2,locnum)

	# e_freq <<- 0.2 # ""

        #var_prop_alpha <<- 1.0 # ""
        
	BBB$var_prop_alpha <- rep(1,locnum)
	BBB$var_prop_beta  <- rep(0.7,popnum) # ""

	#var_prop_a_p   <<- 0.2 # ""
	nb_pilot_alpha    <- 0 # " "
	
	BBB$mean_alpha 	  <- array(0,c(locnum)) # passt start.cpp Ln 306
	BBB$var_alpha  	  <- array(5,c(locnum)) # hm bei berend vorher 0, passt auch in start.cpp 5
	BBB$m2         	  <- array(0,c(locnum)) # passt
	
	BBB$alpha_included    <- array(FALSE,locnum) # array(1,c(locnum)) # in C ist es false !!!
	



        #nb_alpha_included <<- 0
	
	#test 
	#mean_test <<- 0
	#loglog    <<- 0
	#BwBw      <<- 0
	#staba <<- 1
	###
	
	if(pauses==1){
		print("About to start the pilot runs, press any key to do so")
		readline()
	}

	# outputalpha <<- NULL
        # outputbeta  <<- NULL
	# logout <<- NULL	
  
	p <- 0
	cat("Pilotruns\n")


        #test
        nb_pilot      <-  nb.pilot        #20    
        pilot_runtime <-  pilot.runtime   #5000 
        e_f 	      <-  0.05

        BBB$GLOBAL_INIT1  <- matrix(0,locnum,dim(BBB$freq_locus)[2])
	BBB$GLOBAL_INIT2  <- rep(NA,2)
        BBB$GLOBAL_INIT3  <- numeric(locnum)

        BBB$PILOT     <- TRUE

        PILOT_TOTAL   <-  nb_pilot * pilot_runtime 
	
        # rm(population)
        # print(is(population))

## PROGRESS #########################
 progr <- progressBar()
#####################################



	for(k in 1:nb_pilot){

  #print("var_prop_alpha")
  #print(var_prop_alpha)
  #print("var_prop_beta")
  #print(var_prop_beta)
  #print("e_ancestral")
  #print(e_ancestral)
  
  
		for(i in 1:pilot_runtime){ # in C pilot_length
	                 
	               # cat(i,"\n")
                       
			update_d_betaco()         # alpha update		
			#cat("bet finished \n")
                        update_d_alphaco_i(FALSE) # beta update	-> wenn nur fst berechnet werden soll anscheinend nicht notwendig
			#cat("alpha finished \n")
                        update_freq_codominant()  # freq update
		        #cat("freq finished \n")
		
			if(2*(k+1)>=nb_pilot){ 
					BBB$mean_alpha <- BBB$mean_alpha + BBB$d_alpha
					BBB$m2         <- BBB$m2 + BBB$d_alpha^2			
			}		

			
		}#end of one pilot run (innere Schleife)
			

			if(2*(k+1)>=nb_pilot){ 
				nb_pilot_alpha <- nb_pilot_alpha + 1 
     			 }
			
      				 # alpha var updates
     				 check <- (BBB$acc_alpha/pilot_runtime)>0.4
      				 BBB$var_prop_alpha[check] <- BBB$var_prop_alpha[check]*1.2
     				 check <- (BBB$acc_alpha/pilot_runtime)<0.2
      				 BBB$var_prop_alpha[check] <- BBB$var_prop_alpha[check]/1.2
		  
		
     				 # var beta updates
      				 check <- BBB$acc_beta/(pilot_runtime)>0.4
     				 BBB$var_prop_beta[check] <- BBB$var_prop_beta[check]*1.2
     				 check <- BBB$acc_beta/(pilot_runtime)<0.2
     				 BBB$var_prop_beta[check] <- BBB$var_prop_beta[check]/1.2
		  
		
      				# frquencie var update
      				check <- (BBB$acc_freq_ancestral/pilot_runtime)>0.4
     				BBB$e_ancestral[check] <- BBB$e_ancestral[check]*1.1
      				check <- (BBB$acc_freq_ancestral/pilot_runtime)<0.2
      				BBB$e_ancestral[check] <- BBB$e_ancestral[check]/1.1
		

      				BBB$acc_alpha          <- numeric(locnum)
      				BBB$acc_beta  	       <- numeric(popnum)
      				BBB$acc_freq_ancestral <- numeric(locnum)
			
		
  # PROGRESS #######################################################
    progr <- progressBar(k,nb_pilot, progr)
  ###################################################################
	
  

}# end of pilot runs 

  	
        #################################################################
        # end of pilot runs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#################################################################


          BBB$mean_alpha <- BBB$mean_alpha/(nb_pilot_alpha*pilot_runtime) 
          BBB$var_alpha  <- BBB$m2/(nb_pilot_alpha*pilot_runtime)  - BBB$mean_alpha^2 
#print(BBB$mean_alpha)
#print(BBB$var_alpha)

       ################################
       ## MAIN                        #
       ################################
       
        BBB$alpha_updates <- 0 
        BBB$d_alpha       <- numeric(locnum) 
        BBB$PILOT         <- FALSE

	#  MAIN LOOP MCMC loop
	main_runtime <- main.runtime    # 150000 # tot_nr_of_iter
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

	

for(i in 1:main_runtime){
    
   #cat(i,"\n")
   #print("alpha")
   #print(d_alpha)
   #print("beta")
   #print(d_beta)
   #print(BBB$d_alpha)

    # update beta  
    update_d_betaco()
    #------------	

    # update alpha
    # wenn nur fst berechnet werden soll anscheinend nicht notwendig	
    
     id   <- which(BBB$alpha_included)
    if(length(id)!=0){
     update_d_alphaco_i(id)
     BBB$alpha_updates <- BBB$alpha_updates + length(id)
    }
    # ------------
    		
    update_freq_codominant()		
    jump_model_codominant()
	

    #if(!(i%%1000)){
    #   log_likelihood <<- allelecount_loglikelihood()
    #}	

               
               if(i>discard){ # 50000 discard
                cur_out <- cur_out + 1
                ## Output
                     drinne                      <- which(BBB$alpha_included)
 	           if(length(drinne)>0){
 	                   nb_alpha[drinne]      <- nb_alpha[drinne] + 1
                   }
                   
                     post_alpha            <- post_alpha + BBB$d_alpha 
                     val                   <- outer(BBB$d_alpha,BBB$d_beta,"+")
                     val                   <- 1/(1+exp(-val))
                     val                   <- rowSums(val) 
                     cur_fst               <- val/popnum
                     post_fst              <- post_fst  + cur_fst
              }

   # PROGRESS #######################################################
    progr <- progressBar(i,main_runtime, progr)
   ###################################################################
   #print(nb_alpha)
	} # End of Main Loop
	
        outalpha      <- post_alpha/cur_out
        outfst        <- post_fst/cur_out
        outnb_alpha   <- nb_alpha/cur_out
	names(outfst) <- rownames(population[[1]])
        


	output           <- new("BAYESRETURN")
	output@alpha     <- as.numeric(outalpha)
	output@beta      <- as.numeric(BBB$d_beta)
	output@a_inc     <- as.numeric(BBB$alpha_included)
	output@var_alpha <- as.numeric(BBB$var_alpha)
	output@fst       <- outfst
	output@P         <- as.numeric(outnb_alpha)
	
	
	return(output) 
}


