


#******************************************************************************************************************
# Function R2Conquest for Rasch Model
##NS export(rasch.conquest)
rasch.conquest <- function( dat , path.conquest , conquest.name = "console" ,
                            converge = 0.001 , deviancechange = 0.0001 , 
                            iter = 800 , nodes = 20 , 
							minnode=-6 , maxnode = 6 , show.conquestoutput = FALSE , 
                            name = "rasch" ,   pid = 1:(nrow(dat) )  , 
                            wgt = NULL , X = NULL , set.constraints = NULL , 
                            model = "item" , regression = NULL , itemcodes = seq( 0 , max(dat,na.rm=TRUE) ) , 
                            constraints = NULL  , digits = 5 , onlysyntax = FALSE , qmatrix = NULL ,
                            import.regression = NULL , anchor.regression = NULL , anchor.covariance = NULL , 
							pv = TRUE , designmatrix = NULL ,
                            only.calibration = FALSE , init_parameters = NULL , n_plausible = 10 ,
							persons.elim = TRUE , est.wle = TRUE  , save.bat = TRUE , use.bat = FALSE ,
							read.output=TRUE , ignore.pid = FALSE
                                ){
        # INPUT:
        # pv ... estimation of plausible values?
#        pv <- TRUE
		s1 <- Sys.time()
        if ( is.null(wgt) ){ wgt <- 1 + 0 * pid }
        if (onlysyntax == FALSE){ 
            cat("---------------------------------------------------------------------------------------------------------- \n")
            cat("Marginal Maximum Likelihood Estimation \n")
            cat("Estimation of the Rasch Model in ConQuest\n")
            cat("---------------------------------------------------------------------------------------------------------- \n")
            utils::flush.console()
                }
        # shQuote() for path.conquest
        path.conquest <- shQuote( path.conquest )
        # pv estimation?
        if (only.calibration){ pv <- FALSE }
        
        # calculate score
        score1 <- rowMeans( dat , na.rm=T )
        # character vector for items in constraints matrix
        if (! is.null( constraints) ){ 
                constraints[,1] <- paste( constraints[,1] ) 
                constraints <- constraints[ constraints[,1] %in% colnames(dat) , ]
                                }
        # data preparation (eliminate items)
        items.elim <- which( colSums( 1 - is.na(dat) ) == 0 )
        # initial warnings (no responses on some items)
        if ( length(items.elim) > 0 & is.null(constraints) ){ 
                    cat( "\nSome items have been removed because they have no responses!\n" ) 
                    utils::flush.console()
                    dat <- dat[ , - items.elim ]
                    }
        # data preparation (eliminate persons)
		if ( ! is.null(X) ){ persons.elim <- FALSE }
		
		if ( persons.elim ){ 
			persons.elim <- which( rowSums( 1 - is.na(dat) ) == 0 )        
			if ( length(persons.elim) > 0 ){ 
						cat( "\nSome persons have been removed because they have no responses!\n" ) 
						utils::flush.console()
						dat <- dat[ - persons.elim ,  ]
						pid <- pid[ - persons.elim ]
						wgt <- wgt[ - persons.elim ]
						}
				}
        # number of items
        I <- ncol(dat)
        p1 <- path.conquest
        # remove files
        h1 <- c( "analysis.bat" , "rasch.cqc" ,  "rasch.nam" , "rasch.shw" , "rasch.itn" , "rasch.wle" , "rasch.des" , 
                        "rasch.par" , "rasch.constrpar" , "rasch.regr" , "rasch.cov" , "rasch.pv" , "rasch.dat" , "rasch.log" ,
						"rasch.designmatrix"	)
        h1 <- gsub( "rasch" , name , h1 )
        f1 <- getwd()
        h1 <- intersect( list.files( f1 ) , h1 )
        base::file.remove(h1 )  
        # generate label file for ConQuest
        base::writeLines( c( "===> item" , paste( 1:I , colnames(dat) ) ) , 
                        paste( name ,".nam" ,sep="") )

										
        # which items are being constrained
        if ( ! is.null( constraints) ){
            constraints$itemnr <- match( constraints[,1] , colnames(dat) )
            constraints <- constraints[ order( constraints$itemnr) , ]                       
            c1 <- paste( constraints[,3] ,  constraints[,2] ,       paste(  " /* " , constraints[,1] , "*/" ) )
            utils::write.table( c1  ,  paste(  name ,".constrpar" ,sep="") , quote=F , row.names=F , col.names=F )
                    }											
        # extract decimals for person ID
        PP <- round( max( floor( log(  pid , 10 ) + 1 ) ) )
        # extract decimals for weights
        WW <- round( max( floor( log( wgt , 10 )) + 1 ) + digits )
        # do covariates exist? If, then add constant column
        if (! is.null(X)){ 
            X <- data.frame( "one" = 1 , X )
            XNcols <- apply( abs( X ) , 2 , FUN = function(ll){ floor(max( log( ll , 10 ) )+1)+1 } )
            rdcols <- digits * ( 1 - colMeans( X == round(X) ) )
            XNcols <- XNcols + rdcols + 1
            # possible correction
            XNcols <- ifelse( floor(XNcols ) != XNcols , floor(XNcols)+1 , XNcols )		
			
            # write data
            dat2 <- data.frame( pid , wgt ,  X , dat )
            i1 <- write.fwf2( dat = dat2  , format.full = c( PP , WW , XNcols ,  rep(1,I)) , 
                   format.round  = c( 0 , digits , rdcols , rep(1,I) ) , savename = name )
                        }				
        if (is.null(X)){ 
            dat2 <- data.frame( pid , wgt ,   dat )
            i1 <- write.fwf2( dat = dat2  , 
					format.full = c( PP , WW ,   rep(1,I)) , 
                    format.round  = c( 0 , digits ,  rep(1,I) ) , 
					savename = name )
                    }
        # statement data input
        state1 <-   paste0( "format " )
		if ( ! ignore.pid ){
			state1 <- paste0( state1 , 	"pid " , 1 , "-" , PP )
							}
            state1 <- paste0( state1 , " wgt " , PP+1 , "-" , PP+WW ,  
                            " responses " , PP+1+WW , "-" , PP + I + WW ,";" )

        if ( ! is.null(X) ){
            i1a <- i1[ seq( 1 , 2 + ncol(X) ) , ]
            state2 <- paste( paste( i1a$variable  , " " , i1a$begin , "-" , i1a$end , sep="" ) , collapse = " " )
            state1 <- paste( "format " , state2 , paste( "responses " , i1[ 3 + ncol(X) , "begin" ] , "-" , i1[ nrow(i1) , "end" ] , ";" , sep="") , sep=" " )
                }
										
        # input designmatrix
        if ( ! is.null(designmatrix) ){        
            utils::write.table( designmatrix , paste( name , ".des" , sep="") , col.names=F , row.names=F , quote=F )
            desformat <- readLines(  paste( name , ".des" , sep="")  )    
            desformat <- gsub( "0" , "-0" , desformat )
            desformat <- c( ncol(designmatrix) , desformat )
            base::writeLines( desformat , paste( name , ".des" , sep="") )         
                                    }
        # regression statement
        if (is.null(X)){ regrstate <- "" } else { 
                    regrstate <- paste( "regression " , paste( colnames(X)[-1] , collapse=" " ) , " ;" , sep="") }
        if ( ! is.null(regression) ){ regrstate <- paste( "regression " , regression , " ;" , sep="") }
        # constraints statement
        constrstate <- ""
        if( ! is.null(constraints)  ){ constrstate <-  "set constraints = none;" } else
            {   if( is.null(constraints) & is.null(X) ) { 
                    constrstate <- "set constraints = cases;" } else
                    { constrstate <- "set constraints = items;" }
            }  
        if ( ! is.null( set.constraints )){ 
                constrstate <- paste( "set constraints =" , set.constraints , ";" , sep="") 
                        }
        # generate score statement
        if (! is.null(qmatrix)){ 
#                q1 <- qmatrix[  match( colnames(dat) , qmatrix[,1]) , ]
				q1 <- qmatrix		
                DD <- ncol(q1)
                h1 <- paste( "score (" , paste( itemcodes , collapse= " " ) , ")" , sep="" )
                for (dd in 1:DD){
                    # dd <- 1
                    h1.dd <- outer( q1[,dd] , itemcodes ,  FUN = function(x,y){ x * y } ) 
                    h1.dd <- apply( h1.dd , 1 , FUN = function(ll){ paste( ll , collapse= " " ) } )
                    h1 <- paste( h1 , "(" , h1.dd , ")" , sep="")
                        }
                scorestate <- paste( h1 , " !items(" , 1:I , ");   /* " , q1[,1] , " */" , sep="" )
                    } else { 
                                DD <- 1
                                scorestate <- NULL 
                                }
						
						
        # change number of nodes in monte carlo integration
        if (DD > 2 & nodes < 100 ){ nodes <- 2500 }
        # generate ConQuest Syntax for Rasch model
        cqc <- c( 
                paste( "title " , name , " " , Sys.time() , " ;" ,  sep="") , 
                paste( "let name = " , name , ";" , sep="") , 
                "datafile %name%.dat ;" , 
                state1 ,
                "caseweight wgt;" ,
                paste( "codes " , paste( itemcodes,collapse=",") ,";",sep="") , 
                ifelse( is.null(designmatrix) , "labels << %name%.nam;"  , "" )  ,          "" ,
                scorestate , ""     ,
                # initial item parameters
                ifelse( ! is.null(init_parameters) , paste( "import init_parameters << " , init_parameters , " ; " , sep="")  , "") , 
                # import designmatrix
                ifelse( is.null(designmatrix)  , "" , "import designmatrix << %name%.des ;" ) , 
                # anchor item parameters
                ifelse( is.null(constraints)  , "" , "import anchor_parameters << %name%.constrpar ;" ) , 
                ifelse( is.null(import.regression) , "" , paste( "import init_reg_coefficients << " , import.regression , " ;" , sep="") ) ,
                ifelse( is.null(anchor.regression) , "" , paste( "import anchor_reg_coefficients << " , anchor.regression , " ;" , sep="") ) ,
                ifelse( is.null(anchor.covariance) , "" , paste( "import anchor_covariance << " , anchor.covariance , " ;" , sep="") ) ,
                constrstate , 
                paste( "model " , model , ";" ,sep="") ,
                regrstate ,
                paste( "estimate !iter=" , iter , ", nodes = " , nodes , 
				    ", minnode = " , minnode ,  ", maxnode = " , maxnode , 
                    ifelse( DD > 2 , ", method = montecarlo " , "" ) , 
                    ", converge =" , gsub( " " , "" , formatC(converge,8)) ,
                    ", deviancechange = " , gsub( " " , "" , formatC(deviancechange,8)) , "; " ,  sep="" ) ,
                "export parameters              >> %name%.par;" , 
                "export reg_coefficients        >> %name%.regr;" ,
                "export covariance              >> %name%.cov;" , 
                "export logfile                 >> %name%.log;" , 
				"export designmatrix >> %name%.designmatrix ;" ,			
                ifelse( ( ! only.calibration ) &  est.wle , "show cases ! estimates = wle >> %name%.wle ;", "" ) , 
                ifelse( pv  , paste( "set n_plausible=" , n_plausible , ";",sep="") , "" ) ,
                ifelse( pv  , "show cases !estimates=latent >> %name%.pv;" , "" ) , 
                "show >> %name%.shw ;" ,
                ifelse(!only.calibration , "itanal >> %name%.itn;","") ,
                "quit;"               )	
        # write ConQuest syntax
        base::writeLines( cqc , paste( name , ".cqc" , sep="") )
        # write bat file
		if ( save.bat ){ 
			base::writeLines( paste( path.conquest , "\\" , conquest.name , " " ,  
					paste( name , ".cqc" , sep="") , sep = "") , "analysis.bat" )
						}
        if ( onlysyntax ){ 
				cat( paste( "Conquest Input Syntaxes are in " , getwd() ,  "\n", sep="") )
				res <- NULL 
						} 
                else {
        # link to conquest console
		if ( use.bat ){
#			system( "analysis.bat" , show.output.on.console = show.conquestoutput , invisible = FALSE )
			base::system( "analysis.bat"  )			
					} else {
        base::system(paste(path.conquest,"\\" , conquest.name , ".exe ", name,".cqc",sep=""),
					show.output.on.console = show.conquestoutput , invisible = FALSE)
							}
							
						


		#######################################################
		# READING CONQUEST OUTPUT
		if (read.output){ 					
			if ( ( ! only.calibration ) & est.wle  ){
				# read WLEs and add pid
				wle <- utils::read.table( paste( name , ".wle" , sep="")  )		
				v1 <- c("score" , "max" , "wle" , "se.wle" )
				if (DD > 1 ){ 
					v1 <- rep( v1 , each = DD )
					v1 <- paste( rep( paste( "dim" , 1:DD , sep="")  , 4 ) , v1 , sep="")
							}
				if ( ignore.pid ){			
					colnames(wle) <- c("case" ,  v1 )
					wle <- data.frame( "case" = wle[,"case"] , "pid" = wle[,"case"] , 
								wle[,-1] )						
								} else {
					colnames(wle) <- c("case" , "pid" , v1 )								
								}
				wle <- data.frame(  wle )
			
						} else { wle <- NA }
												
			# read shw file
			shw <- base::readLines( paste( name , ".shw" , sep="") )
			deviance <- as.numeric( substring( shw[ grep( "Final Deviance:" , shw ) ]  , 17 ) )
			numbiter <- as.numeric( substring( shw[ grep( "The number of iterations:" , shw ) ]  , 26 ) )
			ind1 <- grep("TERM 1: item" , shw ) + 6
			ind2 <- grep("An asterisk next to a" , shw )[1] - 2
						
			# extract trait mean and trait variance
#			m3 <- read.table( paste( name , ".regr" , sep="") , header=F)
			m3 <- utils::read.table( paste( name , ".regr" , sep="") )			
			mean.trait <- m3[ 1, 3]
			trait.variance <- as.numeric( substring( shw[ grep( "Variance " , shw ) ] , 25 ) )
			wle.rel <-  as.numeric( substring( shw[ grep( " WLE Person separation RELIABILITY:" , shw ) ] , 37 ) )
			itemdiff <- as.numeric( substring( shw[ ind1:ind2 ] , 20 , 27 ) )

			# read pv file
			if ( pv  ) {                
					if (DD == 1){ 
						pv1 <- .read.pv( pvfile = paste( name , ".pv" , sep="") , npv=n_plausible ) 						
								} else {
									pv1 <- .read.multidimpv( pvfile = paste( name , ".pv" , sep="") , ndim = DD , 
												npv = n_plausible  )
										}
				if ( mean( is.na( pv1$pid ) ) == 1 ){
							pv1$pid <- pv1$case 
									}										
					pv1 <- pv1[ , - grep( "case" , colnames(pv1) ) ]			
					wle <- merge( x = wle , y = pv1 , by = "pid" , all=T  )
								 }
								 
								 
								 
				 
			#*****
			# reliability estimation
			reliability <- NULL
			if ( ! only.calibration & DD==1){ 
				# WLE reliability
				reliability$wle.reliability <- 1 - mean( wle$se.wle^2 ) / stats::var( wle$wle )
				# EAP reliability
				reliability$eap.reliability <- 1 - mean( wle$se.eap^2 ) / ( stats::var( wle$eap ) +  mean( wle$se.eap^2 ) )		
									}
			# calculate expected probability
			if (DD == 1 & pv ){ 
				thetagrid <- seq( -10 , 10 , .01 )
				mean.trait <- mean( pv1$PV1)
				sd.trait <- stats::sd( pv1$PV1 )
				dens.thetagrid <- stats::dnorm( thetagrid , mean = mean.trait  , sd = sd.trait )
				dens.thetagrid <- dens.thetagrid / sum( dens.thetagrid )
				p.exp <- sapply( itemdiff , FUN = function(bii) { 
					   stats::weighted.mean( stats::plogis( thetagrid - bii ) , dens.thetagrid )
					} )    }
		if ( DD > 1 |  ( !pv )) { p.exp <- rep(NA,I)  }
		if ( DD == 1 & pv  ) { 
	#		emp.discrim <- round( item.discrim( dat , wle$wle ),3) } 
			emp.discrim <- stats::cor( dat , wle$wle , use="pairwise.complete.obs") }
					else { emp.discrim <- NA  }			

					
			diffs <- as.numeric( substring( shw[ ind1:ind2 ] , 20 , 27 ) ) 
			diffs <- diffs - mean(diffs)	
					
		item <- data.frame( "item" = colnames(dat) ,  
						"N" = colSums( 1 -is.na(dat) , na.rm=T) ,
						"p" = round(colMeans( dat , na.rm=T ),3) , 
						"p.exp" = p.exp , 
						"itemdiff" = itemdiff ,
						"itemdiff.cent" = diffs , 
						"se.itemdiff" = as.numeric( substring( shw[ ind1:ind2 ] , 29 , 35 ) ) ,
						"emp.discrim" = round(emp.discrim,3) ,
						"outfit" = as.numeric( substring( shw[ ind1:ind2 ] , 38 , 43 ) ) ,
						"t.outfit" = as.numeric( substring( shw[ ind1:ind2 ] , 58 , 62 ) ) ,
						"infit" = as.numeric( substring( shw[ ind1:ind2 ] , 64 , 69 ) ) ,
						"t.infit" = as.numeric( substring( shw[ ind1:ind2 ] , 84 , 88 ) ) 
								)
			res <- list( "person" = wle , "item" = item , "deviance" = deviance , "numbiter" = numbiter , "mean.trait" = mean.trait , 
					"sd.trait" = sqrt(trait.variance) , "wle.rel" = wle.rel , "itemmean" = mean( item$itemdiff) ,
						"reliability" = reliability )  
			s2 <- Sys.time()
			res$sys.time <- list( "start" = s1  , "end" = s2 , "timediff" = s2-s1 )
	
			
			if ( DD == 1){
					res$shw.itemparameter <- read.show( showfile = paste( name , ".shw" , sep="") )
					res$shw.regrparameter <- read.show.regression( showfile = paste( name , ".shw" , sep="") )        
					res$shw.pimap <- read.pimap( showfile = paste( name , ".shw" , sep="") )        
							}
			class(res) <- "rasch.conquest" 
					}			
			return(res)
			}
        }
#**************************************************************************************************************************


#...............................................................#
# function for reading conquest files of plausible values
.read.pv <- function( pvfile , npv = 5 ){
    pv <- readLines( pvfile )   # read file with plausible values
    nl <- length(pv)            # length of pv file
    nseq <- npv + 3             # number of lines corresponding to one person
    N <- nl/nseq                # number of persons 
    # all plausible values
    pv.all <- as.numeric( substring( pv[  - c( nseq * ( 1:N  ) - 1 , nseq * ( 1:N  ) - 0 , nseq * ( 1:N -1 ) +1  ) ] , 11 ) )
    pv1 <- matrix( pv.all , ncol=npv , byrow=T )
    colnames(pv1) <- paste("PV" , 1:npv, sep="")

    pid <- as.numeric( substring( pv[ nseq * ( 1:N  ) - npv -2 ] , 6 , 40 ) )
    pv1 <- data.frame(  "case" = 1:N , 
                        "pid" = pid , 
                        "eap" = as.numeric( pv[ nseq * ( 1:N  ) - 1  ] ) , 
                        "se.eap" = as.numeric( pv[ nseq * ( 1:N  )  ] ) ,
                        pv1
                            )
    return(pv1)
    }
#................................................................#



#...............................................................#
# function for reading conquest files of plausible values
.read.multidimpv <- function( pvfile , ndim , npv = 5  ){
    # INPUT:
    # pvfile    ... name of file with plausible values
    # ndim      ... number of dimensions for which plausible values were drawn
    pv <- scan( pvfile )
    pv1 <- readLines( pvfile )
    a1 <- pv1[1]
    a1 <- gsub( " " , "" , a1 )
	pid.pres <- 1
#    if ( a1 == "1" ){ pid.pres <- 0 } else { pid.pres <- 1 }
    nl <- length(pv1)            # length of pv file
    nseq <- npv + 3             # number of lines corresponding to one person
    N <- nl/nseq                # number of persons 

    # numbers for each person
    L <- npv * ( ndim + 1 ) + 2 * ndim +  1   + pid.pres
    # arrange matrix
    dfr <- matrix( pv , ncol = L , byrow=T)
    dfr <- cbind( dfr[,1] , seq(1,N) , dfr[,-1] )
    # which columns should be deleted from dfr
	
#    dfr <- data.frame(dfr[ , - c( 1 + seq( 2 , 1 + (ndim+1)*npv , ndim+1 )) ])
	del <- c( 3 ,  ( 0:(npv-1))*(ndim+1) +1+3 )
	dfr <- data.frame(dfr[ , - del ] )
    colnames(dfr)[1:2] <- c( "case" , "pid" )
    cdim <- matrix( sapply( 1:npv , FUN = function(pp){ paste( "dim" , 1:ndim , "pv" , pp , sep="") } ) , ncol= 1)[,1]
    colnames(dfr)[ seq( 1 , length(cdim) ) + 2 ] <- cdim
    cdim <- c( paste( "dim" , 1:ndim , "eap" , sep="") , paste( "dim" , 1:ndim , "seeap" , sep="") )
    colnames(dfr)[ seq( 1 + ncol(dfr) - 2* ndim  , ncol(dfr) )  ] <- cdim
    return(dfr)
    }
#................................................................#
R2conquest <- rasch.conquest


read.pv <- .read.pv
read.multidimpv <- .read.multidimpv

#*******************************************************
# Summary for rasch.conquest object                         *
##NS S3method(summary,rasch.conquest)
summary.rasch.conquest <- function( object , ... ){
    # object      ... object from rasch.conquest                #
    cat("---------------------------------------------------------------------------------------------------------- \n")
    cat("Marginal Maximum Likelihood Estimation (Normal Distribution) \n")
    cat("Estimation of the Rasch Model in ConQuest\n" )
    cat("---------------------------------------------------------------------------------------------------------- \n")
    cat(paste( "Rasch Model with ", nrow(object$item) , " Items (Average Number of Items per Person: " , 
                 round( mean( object$person$max ) , 1 ) ,  " )\n"  , sep="") )
    cat( "Deviance = " , round( object$deviance , 2 ) , "\n" )
    cat( "Trait Distribution: " , 
              "Mean=" , round(object$mean.trait,3) , " SD=" , round( object$sd.trait , 3) , "  WLE Reliability =" , round( object$wle.rel , 3) , "\n") 
    cat("Mean of Item Difficulty Parameters: ", round( object$itemmean , 3 ) , "\n" )
    cat("---------------------------------------------------------------------------------------------------------- \n")
	cat( paste( "WLE Reliability:" , round( object$reliability$wle.reliability,3 )) , "\n")
	cat( paste( "EAP Reliability:" , round( object$reliability$eap.reliability,3 )) , "\n")
    cat("---------------------------------------------------------------------------------------------------------- \n")
    cat("Item Parameter \n")
    print( object$item[ , c("N" , "p" , "itemdiff" , "emp.discrim" , "outfit" , "infit" )] )                
                }
#*******************************************************
summary.R2conquest <- summary.rasch.conquest



##   #*******************************************************
##   # plot rasch.conquest object                   
##   ##NS S3method(plot,rasch.conquest)
##   plot.rasch.conquest <- function( x , plottitle = "Person-Item-Map" , 
##       n.gridpoints = 41 , itemcex = 1 , ...){
##       # calculate empirical WLE distribution
##       object <- x
##       trait.grid <- seq( -4 , 4 ,len= n.gridpoints )
##       h <- diff(trait.grid)[1]
##       breaks.grid <- c( - 1000 ,  ( trait.grid + h / 2 )[ - n.gridpoints ] , 1000 )
##       freq <- hist( object$person$wle , breaks = breaks.grid , plot=F)
##       object$trait.distr <- cbind( trait.grid , freq$counts / sum( freq$counts)  )
##       object$fixed.a <- 1 + 0 * 1:(nrow(object$item))
##   #    class(object) <- "rasch.mml"
##       plot.rasch.mml( object , plottitle = plottitle , itemcex = itemcex)
##       }
##   #******************************************************




#-----------------------------------------------------------------------------------------------------------#
# Read parameter for one term from a ConQuest output shw file                                                     #
##NS export(read.show.term)
read.show.term <- function( showfile , term ){
        shw <- readLines( showfile )
        ind <- intersect( grep( "TERM" , shw ) , grep( term , shw ) )[1]
        i1 <- ind+6
        i2 <- grep( "--------------" , shw )
        i2 <- i2[ i2 > i1 ][1] - 1
        table1 <- shw[i1:i2]
        h1 <- regexpr( "ESTIMAT" , shw[i1-2] )[1]        
        h2 <- h1 - 21
        dfr1 <- data.frame(     "term" = shw[ind] ,
                        "parameter" = substring( table1 , 6, h1- 1 )  , 
                        "parameter2" = gsub( " " , "" , substring( table1 , 6, h1- 1 ) ) , 
                        "itemdiff" = as.numeric( substring( table1 , 21+h2,27+h2 ) ) , 
                        "constrained" = 1* ( substring( table1 , 28+h2,28+h2 ) == "*" ) , 
                        "se.itemdiff" = as.numeric( substring( table1 , first=31+h2 , last=35+h2 )) , 
                                "outfit" = as.numeric( substring( table1 , first=39+h2 , last=43+h2 )) , 
                                "t.outfit" = as.numeric( substring( table1 , first=58+h2 , last=62+h2 )) , 
                                "infit" = as.numeric( substring( table1 , first=65+h2 , last=69+h2 ) )  ,
                                "t.infit" = as.numeric( substring( table1 , first=84+h2 , last=88+h2 ) )               
                                        )
        p2 <- paste(dfr1$parameter)
        p2 <- unlist( lapply( strsplit( p2 , split=" " ) , FUN = function(ll){ 
                    paste( ll[ ll != "" ] , collapse="-") } ) )
        dfr1$parameter2 <- p2
        return(dfr1)
            }
#-----------------------------------------------------------------------------------------------------------#



#.............................................................................
# read all terms from a ConQuest shw file
##NS export(read.show)
read.show <- function( showfile){
    shw <- readLines( showfile )
    # liste alle verschiedenen auftretenden Terme auf
    terms.shw <- shw[ grep( "TERM" , shw ) ]
    terms.shw <- strsplit( shw[ grep( "TERM" , shw ) ] , split =": ")
    terms.shw <- unlist( lapply( terms.shw , FUN = function(ll){ paste(ll[1],sep=":") } ) )
    dfr <- NULL
    for (tt in terms.shw ){
        dfr1 <- read.show.term( showfile = showfile , term = tt )
        dfr <- rbind( dfr , dfr1 )
            }
    return(dfr)
    }
#.............................................................................






#......................................................................................
# read regression coefficients and residual variance from a ConQuest output showfile
##NS export(read.show.regression)
read.show.regression <- function( showfile ){
    shw <- readLines( showfile )
    # Regression Coefficients
    i1 <- grep("Regression Variable" , shw )+2
    i2 <- grep("--------------------------" , shw )
    i2 <- i2[ i2>i1][1] - 1
    dfr1 <- data.frame( "parameter" = gsub( " " , "" , substring( shw[i1:i2] , 1 , 25 ) ) , 
                "estimate" = as.numeric( substring( shw[i1:i2] , 26 , 32 ) ) ,
                "se" = as.numeric( substring( shw[i1:i2] , 35 , 40 ) ) 
                    )
    # Variance
    i1 <- grep("Variance" , shw )
    dfr2 <- data.frame( "parameter" = "Variance" , 
                "estimate" = as.numeric( substring( shw[i1] , 26 , 33 ) ) ,
                "se" = NA
                    )
    dfr <- rbind( dfr1 , dfr2 )
    dfr$par.index <- seq( 1 , nrow(dfr) )
    return(dfr)
        }
#......................................................................................




#---------------------------------------------------------------------------
# read person-item map from ConQuest output shw file
##NS export(read.pimap)
read.pimap <- function( showfile ){
    shw <- readLines( showfile )
    ind1 <- grep( "MAP OF LATENT " , shw )[1] + 5
    ind2 <- grep( "===========================" , shw )
    ind2 <- ( ind2[ ind2 > ind1 ] )[1] - 1
    ZZ <- ind2 - ind1 
    dfr <- data.frame( matrix(NA , nrow=ZZ , ncol=6 ) )
    colnames(dfr) <- c( "row" , "first.bar" , "second.bar" , 
                "ylab" ,  "X" , "items" )
    dfr$row <- seq( 1 , ZZ )
    for (zz in 1:ZZ){
        # zz <- 25
        shw.zz <- unlist( strsplit( shw[ind1+zz] , split="") )
        a1 <- grep( "\\|" , shw.zz )
        dfr$first.bar[zz] <- a1[1]
        dfr$second.bar[zz] <- a1[2]
        dfr$X[zz] <- sum( shw.zz[ seq(1,a1[1]) ] == "X" )
        items.zz <- paste( shw.zz[ seq(a1[1]+1,a1[2]-1) ] , collapse="" )
        items.zz <- unlist( strsplit( items.zz , split=" " ) )
        items.zz <- items.zz[ items.zz != "" ]
        dfr$items[zz] <- paste( items.zz , collapse="-" )
        dfr$ylab[zz] <- as.numeric(shw.zz[4])
            }
    return(dfr)
        }
#---------------------------------------------------------------------------


