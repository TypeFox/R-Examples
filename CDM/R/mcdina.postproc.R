
#################################################
# mcdina information criteria
mcdina.calc.ic <- function( dev, weights , itemstat , pi.k , G , I ,
			zeroprob.skillclasses ,  reduced.skillspace , Z ){
	ic <- list( "deviance" = dev , "n" = sum(weights) , "loglik" = -dev/2 )
	ic$G <- G
	ic$itempars <- sum( itemstat$N.pars)
	ic$traitpars <- G*(nrow(pi.k)-1 - length( zeroprob.skillclasses ) )
	if ( reduced.skillspace ){
		ic$traitpars <- G * ncol(Z) 
				}	
	ic$np <- ic$itempars + ic$traitpars	
	ic$Nskillclasses <- nrow(pi.k) - length( zeroprob.skillclasses )
	# AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (consistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
    # corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	return(ic)
		}
#######################################################
# standard errors
mcdina.calc.se.delta <- function( delta , n.ik , probs , lr_list , lc_list , 
			itemstat , I , G , itempars , lr_counts , CC ){	
	se.delta <- delta
	for (ii in 1:I){
	for (gg in 1:G){
		lc.ii <- lc_list[[ii]]
		lr.ii <- lr_list[[ii]]
		se.delta[ii,,,gg] <- sqrt( delta[ii,,,gg] * ( 1 - delta[ii,,,gg] ) / 
			matrix( lr_counts[ii,,gg] , nrow= CC  , ncol= CC , byrow=TRUE) )
					} # end gg 
				} # end ii
		return(se.delta)
			}
###############################################################
# collect item parameters
mcdina.collect.itempars <- function( I , lc , itempars , itemstat , dat ,
	G , CC , delta , se.delta , group0_unique  ){
	item <- NULL
	for (ii in 1:I){
	# ii <- 1
	lc.ii <- lc[ lc$item == ii , ]
	ip.ii <- itempars[ii]
	itemstat.ii <- itemstat[ii,]
	if ( ip.ii == "gr" ){ G1 <- G } else {G1 <- 1 }
	for (gg in 1:G1){  #  gg <- 1
	delta.ii <- delta[ ii ,,,gg ]
	se.delta.ii <- se.delta[ ii ,,,gg ]
		for (cc in 1:itemstat.ii$N.lr ){ # cc <- 1
			lc.ii.cc <- lc.ii[ lc.ii$lr_index == cc  , ]
			lc.ii.cc <- lc.ii.cc[1,]
			item.cc <- data.frame( "item" = colnames( dat )[ii] , "itemnr" = lc.ii.cc$item )
			item.cc$lr <- lc.ii.cc$lr
			item.cc$lr_level <- lc.ii.cc$lr_level
			item.cc$Q <- lc.ii.cc$Q
			item.cc$lr_index <- lc.ii.cc$lr_index
			item.cc$max.cat <- lc.ii.cc$max.cat
			item.cc$partype <- itemstat.ii$partype
			item.cc$group <- group0_unique[gg]
			if ( item.cc$partype != "gr" ){ item.cc$group <- NA }
			d1 <- t( delta.ii[,cc] )
			colnames(d1) <- paste0( "Cat" , 1:CC )
			if ( itemstat.ii$N.cat < CC ){
				d1[ 1, seq(itemstat.ii$N.cat + 1 , CC ) ] <- NA 
					}
			item.cc <- cbind( item.cc , d1 )
			d1 <- t( se.delta.ii[,cc] )
			if ( itemstat.ii$N.cat < CC ){
				d1[ 1, seq(itemstat.ii$N.cat + 1 , CC ) ] <- NA 
					}
			colnames(d1) <- paste0( "se.Cat" , 1:CC )
			item.cc <- cbind( item.cc , d1 )
			item <- rbind( item , item.cc )
					}
		}	
		}
	return(item)
	}
#########################################################
# skill probabilities
mcdina.skill.patt <- function( q.matrix , skillclasses , G , pi.k ,
	group0_unique){	
	maxK <- max( q.matrix[ , -c(1:2) ] )
	K <- ncol(skillclasses)
	skill.patt <- matrix( NA , nrow=K , ncol=(maxK+1)*G )
	skill.patt <- as.data.frame( skill.patt)
	zz <- 1
	for (kk in 0:maxK){ # kk <- 0
	for (gg in 1:G){  # gg <- 1
	for (ss in 1:K){   #ss <- 1
		skill.patt[ ss , zz ] <- sum( pi.k[  skillclasses[ ,ss] == kk  ,gg]  )
		ind <- which( skillclasses[ ,ss] == kk )
		if ( length(ind) == 0 ){ 
				skill.patt[ss,zz] <- NA 
						}
				}
		colnames(skill.patt)[zz] <- paste0("skill.prob" , kk , ".group", group0_unique[gg] )					
		zz <- zz+1
				}
				}
	
				
	rownames(skill.patt) <- colnames(q.matrix)[ - c(1:2) ]
	return(skill.patt)				
			}