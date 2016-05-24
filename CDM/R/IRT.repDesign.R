
########################################################
IRT.repDesign <- function( data , wgt = NULL, jktype = "JK_TIMSS", 
     jkzone = NULL, jkrep = NULL, jkfac = NULL, fayfac = 1 , 
     wgtrep = "W_FSTR" , ngr=100 , Nboot=200 , seed = .Random.seed ){
		#---------------------------------
		# BIFIEsurvey must be loaded
		base::requireNamespace("BIFIEsurvey")
		
		data <- as.data.frame(data)	
		
		#*************************
		# apply Jackknife (Bootstrap) routine		
		if ( jktype != "BOOT"){
			bdat2 <- BIFIEsurvey::BIFIE.data.jack(data, wgt = wgt, jktype = jktype, 
						pv_vars = NULL, jkzone = jkzone , jkrep = jkrep, 
						jkfac = jkfac, fayfac = 1 , 
						wgtrep = wgtrep , pvpre = NULL , ngr=ngr ,
						seed = seed , cdata=FALSE)
							} else {
			bdat2 <- BIFIEsurvey::BIFIE.data.boot( data , wgt=wgt ,  pv_vars = NULL ,
						Nboot = Nboot , seed = seed , cdata=FALSE)					
								}
		cat("+++ Generated IRT.repDesign object\n")			
		# output
		res <- list( "wgt" = bdat2$wgt , "wgtrep" = bdat2$wgtrep , 
                "fayfac" = bdat2$fayfac , "RR" = ncol( bdat2$wgtrep) )			
		class(res) <- "IRT.repDesign"			
		
		return(res)			
					}
########################################################################					