sim.bdsky.stt <- function (n, lambdasky, deathsky, timesky, sampprobsky, omegasky=rep(0,length(timesky)), rho = 0, timestop = 0, model="BD", N=0, trackinfecteds=FALSE )
{	# detail when one wants to print all extant lineages for serially sampled tree, with rho=0 (no sampling at present)
    detail <- 0
    # define the model being used
    if (model=="BD"){
        SIS=SIR=SIRS=FALSE
    } else if (model=="SIS"){
        SIR=SIRS=FALSE
        SIS=TRUE
    } else if (model=="SIR"){
        SIS=SIRS=FALSE
        SIR=TRUE
    } else if (model=="SIRS"){
        SIS=SIR=FALSE
        SIRS=TRUE
    }
    psampsky <- sampprobsky
    # recalculate death rate so it is mu instead of delta
    musky <- deathsky * (1 - psampsky)
    # recalculate sampling rate
    psisky <- deathsky * (psampsky)
    extant <- 0
    if (rho > 0 && timestop == 0) {
        print("rho > 0 and timestop = 0 is not a valid stopping condition")
        return()
    }
    if (n == 0 && timestop == 0) {
        print("n = 0 and timestop = 0 is not a valid stopping condition")
        return()
    }
    if (n > 0 && timestop > 0) {
        print("n > 0 and timestop > 0 are two stopping condition. Condition n is ignored.")
    }
    if (SIR & n > N){
    	print("n > N cannot be ran")
    	return()
    }
    if (SIS || SIR) {
        for(i in (1:length(omegasky))){
            if (omegasky[i]!=0){
                print("omega is not a parameter for SIS/SIR models")
                return()
            }
        }
    }
    timesky <- c(timesky, 10^100)
    timesky <- timesky[-1]
    extincttree = 1
    while (extincttree == 1) {
        edge <- c(-1, -2)
        leaves <- c(-2)
        sampled <- vector()
        timecreation <- c(0, 0)
        extinct <- vector()
        time <- 0
        maxspecies <- -2
        edge.length <- c(0)
        extincttree = 0
        stop = 0
        susc=N-1
        recovered=0
        totalsusceptibles=N-1
        totalrecovereds=0
        totalinfecteds=1
        totalsampleds=0
        timeevent=0
        cumulativeinfected=1
        lambda <- lambdasky[1]
        mu <- musky[1]
        psi <- psisky[1]
        omega <- omegasky[1]
        timecut <- timesky[1]
        shift <- 2
        while (stop == 0) {
            if (length(leaves) == 0) {
                phy2 = 0
                extincttree = 1
                print(paste ("extinct tree at time ", time ,sep=""))
                #print("extinct tree")
                if (timestop > 0) {
                  phy2 <- 0
                  return(phy2)
                }
                stop = 1
            }
            else {
            	# SAMPLE EXPONENTIAL WAITING TIME TILL EVENT AT LEAVES HAPPENS
                timestep <- rexp(1, (length(leaves) * (lambda + mu + psi)))
                if (SIR || SIS || SIRS) {
                	timestep <- rexp(1, (length(leaves) * (lambda/N*susc + mu + psi) + omega*recovered))
                }
                time = time + timestep
                #print(paste ("time ", time, " lambda=",lambda," #infecteds=",length(leaves), " #sampleds=", totalsampleds[length(totalsampleds)], sep=""))
                if (timestop > 0 && time > timestop) {
                  stop = 1
                }
                else {
                  if (time < timecut) {
                  	if ( SIS || SIR || SIRS ){
                    		species <- sample(leaves, 1)
                    		del <- which(leaves == species)
                    		# SAMPLE RANDOM NUMBER TO DETERMINE EVENT THAT HAPPENS
                    		specevent <- runif(1, 0, 1)
                    		# identifies the leave in edgelist (always second column, and counting goes columnwise) - number
                            #   of nodes in first column in edge list == number of lengths in edge.length
                    		edgespecevent <- which(edge == species) - length(edge.length)
                    		##BIRTH
                    		if (((length(leaves)*lambda/N*susc)/(length(leaves)*(lambda/N*susc + mu + psi) + omega*recovered)) > specevent) {
                      			edge.length[edgespecevent] <- time - timecreation[-species]
                      			edge <- rbind(edge, c(species, maxspecies - 1))
                      			edge <- rbind(edge, c(species, maxspecies - 2))
                      			edge.length <- c(edge.length, 0, 0)
                      			leaves <- c(leaves, maxspecies - 1, maxspecies - 2)
                      			maxspecies <- maxspecies - 2
                      			leaves <- leaves[-del]
                      			timecreation <- c(timecreation, time, time)
                      			
                                susc=susc-1
                                timeevent=c(timeevent,time)
                                totalinfecteds=c(totalinfecteds,totalinfecteds[length(totalinfecteds)]+1)
                                totalsampleds=c(totalsampleds,totalsampleds[length(totalsampleds)])
                                
                                if (SIR || SIRS) {
                                totalrecovereds=c(totalrecovereds,totalrecovereds[length(totalrecovereds)])
                                totalsusceptibles=c(totalsusceptibles,totalsusceptibles[length(totalsusceptibles)]-1)
                                }
                                
                                cumulativeinfected=c(cumulativeinfected,cumulativeinfected[length(cumulativeinfected)]+1)
                                
                   			}
                    		## SAMPLING
                   	 		else if (((length(leaves)*(lambda/N*susc + psi))/(length(leaves)*(lambda/N*susc + mu + psi) + omega*recovered)) > specevent) {
                    	  		sampled <- c(sampled, leaves[del])
                      			leaves <- leaves[-del]
                     			edge.length[edgespecevent] <- time - timecreation[-species]
                      			if (length(sampled) == n && timestop == 0) {
                        			stop = 1
                      			}
                                
                                # in SIS model, the S+I=constant, so move the sampled individual back to the susceptibles
                                if (SIS) {
                      				susc=susc+1
                      			}
                                # in SIRS model, the S+I+R=constant, so move the sampled individual to the recovereds
                                if (SIRS) {
                                    recovered=recovered+1
                                }
                                
                                timeevent=c(timeevent,time)
                                totalinfecteds=c(totalinfecteds,totalinfecteds[length(totalinfecteds)]-1)
                                totalsampleds=c(totalsampleds,totalsampleds[length(totalsampleds)]+1)
                                
                                if (SIR || SIRS) {
                                    totalrecovereds=c(totalrecovereds,totalrecovereds[length(totalrecovereds)]+1)
                                    totalsusceptibles=c(totalsusceptibles,totalsusceptibles[length(totalsusceptibles)])
                                }
                                
                                cumulativeinfected=c(cumulativeinfected,cumulativeinfected[length(cumulativeinfected)])
                                
                    		}
                    		## DEATH
                            else if (((length(leaves)*(lambda/N*susc + psi + mu))/(length(leaves)*(lambda/N*susc + mu + psi) + omega*recovered)) > specevent) {
                      			extinct <- c(extinct, leaves[del])
                      			leaves <- leaves[-del]
                      			edge.length[edgespecevent] <- time - timecreation[-species]
                                
                                if (SIS) {
                      				susc=susc+1
                      			}
                                if (SIRS) {
                                    recovered=recovered+1
                                }
                                
                                timeevent=c(timeevent,time)
                                totalinfecteds=c(totalinfecteds,totalinfecteds[length(totalinfecteds)]-1)
                                totalsampleds=c(totalsampleds,totalsampleds[length(totalsampleds)])

                                if (SIR || SIRS) {
                                    totalrecovereds=c(totalrecovereds,totalrecovereds[length(totalrecovereds)]+1)
                                    totalsusceptibles=c(totalsusceptibles,totalsusceptibles[length(totalsusceptibles)])
                                }
                                
                                cumulativeinfected=c(cumulativeinfected,cumulativeinfected[length(cumulativeinfected)])
                                
                    		}
                            ## LOOSE IMMUNITY
                            else {
                                if (SIRS) {
                                    recovered=recovered-1
                                    susc=susc+1
                                    timeevent=c(timeevent,time)
                                    totalrecovereds=c(totalrecovereds,totalrecovereds[length(totalrecovereds)]-1)
                                    totalinfecteds=c(totalinfecteds,totalinfecteds[length(totalinfecteds)])
                                    totalsusceptibles=c(totalsusceptibles,totalsusceptibles[length(totalsusceptibles)]+1)
                                    totalsampleds=c(totalsampleds,totalsampleds[length(totalsampleds)])
                                    
                                    cumulativeinfected=c(cumulativeinfected,cumulativeinfected[length(cumulativeinfected)])
                                }
                            }
                    }
                    
                    else {
                    	species <- sample(leaves, 1)
                    	del <- which(leaves == species)
                    	# SAMPLE RANDOM NUMBER TO DETERMINE EVENT THAT HAPPENS
                    	specevent <- runif(1, 0, 1)
                    	# identifies the leave in edgelist (always second column, and counting goes columnwise) - number
                        #	of nodes in first column in edge list == number of lengths in edge.length
                    	edgespecevent <- which(edge == species) - length(edge.length)
                    	##BIRTH
                    	if ((lambda/(lambda + mu + psi)) > specevent) {
                      		edge.length[edgespecevent] <- time - timecreation[-species]
                      		edge <- rbind(edge, c(species, maxspecies - 1))
                      		edge <- rbind(edge, c(species, maxspecies - 2))
                      		edge.length <- c(edge.length, 0, 0)
                      		leaves <- c(leaves, maxspecies - 1, maxspecies - 2)
                      		maxspecies <- maxspecies - 2
                      		leaves <- leaves[-del]
                      		timecreation <- c(timecreation, time, time)
                            
                            timeevent=c(timeevent,time)
                            totalinfecteds=c(totalinfecteds,totalinfecteds[length(totalinfecteds)]+1)
                            totalsampleds=c(totalsampleds,totalsampleds[length(totalsampleds)])

                            cumulativeinfected=c(cumulativeinfected,cumulativeinfected[length(cumulativeinfected)]+1)
                            
                   		}
                    	## SAMPLING
                   	 	else if (((lambda + psi)/(lambda + mu + psi)) > specevent) {
                            sampled <- c(sampled, leaves[del])
                      		leaves <- leaves[-del]
                     		edge.length[edgespecevent] <- time - timecreation[-species]
                      		if (length(sampled) == n && timestop == 0) {
                        		stop = 1
                      		}
                            
                            timeevent=c(timeevent,time)
                            totalinfecteds=c(totalinfecteds,totalinfecteds[length(totalinfecteds)]-1)
                            totalsampleds=c(totalsampleds,totalsampleds[length(totalsampleds)]+1)
                            
                            cumulativeinfected=c(cumulativeinfected,cumulativeinfected[length(cumulativeinfected)])
                    	}
                    	## DEATH
                    	else {
                      		extinct <- c(extinct, leaves[del])
                      		leaves <- leaves[-del]
                      		edge.length[edgespecevent] <- time - timecreation[-species]
                            
                            timeevent=c(timeevent,time)
                            totalinfecteds=c(totalinfecteds,totalinfecteds[length(totalinfecteds)]-1)
                            totalsampleds=c(totalsampleds,totalsampleds[length(totalsampleds)])
                            
                            cumulativeinfected=c(cumulativeinfected,cumulativeinfected[length(cumulativeinfected)])
                            
                    	}
                    }
                  }
                  ## CHANGE IN SKYLINE INTERVAL PARAMETERS
                  else {
                    time <- timecut
                    timecut <- timesky[shift]
                    lambda <- lambdasky[shift]
                    mu <- musky[shift]
                    psi <- psisky[shift]
                    omega <- omegasky[shift]
                    shift <- shift + 1
                  }
                }
            }
        }
    }
    # AFTER SIMULATION, SET ALL CURRENT NON-SAMPLED LEAVES TO EXTINCT AFTER TIME ELAPSED
	if (rho == 0) {
    	extant <- (length(leaves))
    	while (length(leaves) > 0) {
			del <- 1
			extinct <- c(extinct, leaves[del])
			k = which(edge == leaves[del]) - length(edge.length)
			edge.length[k] <- time - timecreation[-leaves[del]]
			leaves <- leaves[-del]
       	}
    } 
    else {
        time <- timestop
        todiesample = runif(length(leaves), 0, 1)
        todie <- which(todiesample > rho)
        extinct <- c(extinct, leaves[todie])
        if (length(todie) == length(leaves)) {
            print("no sampled extant individuals")
        }
        for (j in (1:length(leaves))) {
            k = which(edge == leaves[j]) - length(edge.length)
            edge.length[k] <- time - timecreation[-leaves[j]]
        }
        if (length(todie) > 0) {
            leaves <- leaves[-todie]
        }
        
        timeevent=c(timeevent,time)
        totalinfecteds=c(totalinfecteds,totalinfecteds[length(totalinfecteds)]-length(leaves))
        totalsampleds=c(totalsampleds,totalsampleds[length(totalsampleds)]+length(leaves))
        
        if (SIR || SIRS) {
            totalrecovereds=c(totalrecovereds,totalrecovereds[length(totalrecovereds)]+length(leaves))
            totalsusceptibles=c(totalsusceptibles,totalsusceptibles[length(totalsusceptibles)]-0)
        }
        
        cumulativeinfected=c(cumulativeinfected,cumulativeinfected[length(cumulativeinfected)])
        
        sampled <- c(sampled, leaves)
        leaves <- vector()
    }
    if (length(sampled) > 0) {
        if (length(extinct > 0)) {
            for (j in 1:length(extinct)) {
                del <- which(edge == extinct[j]) - length(edge.length)
                surpress <- edge[del, 1]
                edge.length <- edge.length[-del]
                edge <- edge[-del, ]
                del2 <- which(edge[, 1] == surpress)
                modify <- which(edge[, 2] == surpress)
                edge[modify, 2] <- edge[del2, 2]
                edge.length[modify] <- edge.length[modify] + 
                  edge.length[del2]
                edge.length <- edge.length[-del2]
                edge <- edge[-del2, ]
            }
        }
        
        nodes <- (length(sampled)) * 2
        leaf = 1
        interior = length(sampled) + 1
        edgetemp <- edge
        if (length(sampled) == 1) {
            print("only one sampled individual")
            return(edge.length)
        }
        temp <- unique(c(edge))
        temp <- sort(temp, decreasing = TRUE)
        for (j in temp) {
            if (sum(match(sampled, j, 0)) == 0 || j == -1) {
                posvalue <- interior
                interior <- interior + 1
            }
            else {
                posvalue <- leaf
                leaf <- leaf + 1
            }
            replacel <- which(edge == j)
            if (length(replacel) > 0) {
                for (k in 1:length(replacel)) {
                  if ((replacel[k] - 1) < length(edge.length)) {
                    edge[replacel[k], 1] <- posvalue
                  }
                  else {
                    edge[(replacel[k] - length(edge.length)), 
                      2] <- posvalue
                  }
                }
            }
        }
        
        timeevent=timeevent[length(timeevent)]-timeevent
        timecuts=c(timesky-time)
        timecuts=-(timecuts[timecuts<1])
        susceptiblesinfectedsrecoveredsvstimelist=list("timesky"=timecuts, "eventtimes"=timeevent, "infecteds"=totalinfecteds, "cumulativeinfecteds"=cumulativeinfected, "cumulativesampleds"=totalsampleds)
        
        if (SIR || SIRS) {
            susceptiblesinfectedsrecoveredsvstimelist=list("timesky"=timecuts, "eventtimes"=timeevent, "infecteds"=totalinfecteds, "cumulativeinfecteds"=cumulativeinfected, "cumulativesampleds"=totalsampleds, "susceptibles"=totalsusceptibles, "recovereds"=totalrecovereds)
        }
        
        phy <- list(edge = edge)
        phy$tip.label <- paste("t", sample(length(sampled)), sep = "")
        phy$edge.length <- edge.length
        phy$Nnode <- length(sampled)
        class(phy) <- "phylo"
        br<-sort(getx(phy,sersampling=1),decreasing=TRUE)
        phy$root.edge<-br[1]-br[2]
        phy<-collapse.singles(phy)
        phy2 <- phy
    } 
    else {
        print("no sampled extant or extinct individuals after time age")
        phy2 <- 0
    }
    if (detail == 1) {
        phy2 <- c(phy2, extant)
    }
    
    ## RETURN TREE AND POSSIBLY THE TRACK LIST OF INFECTEDS OVER TIME
    if (trackinfecteds){
         list(phy2,susceptiblesinfectedsrecoveredsvstimelist)
        
    }
    else {
        list(phy2)
    }

}

