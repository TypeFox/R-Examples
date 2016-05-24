

epp <- function(breedingDat, polygonsDat, eppDat, maxlag = 3) { 

	# bricks
	epd = data.frame(z = paste(eppDat@male, eppDat@female), epp = 1, stringsAsFactors = FALSE)
  	if( missing(polygonsDat) )   polygonsDat = DirichletPolygons(breedingDat)
	
  	nb  = poly2nb(polygonsDat, row.names = polygonsDat$ID, queen = TRUE) 
  	hnb = higherNeighborsDataFrame(nb, maxlag = maxlag)
  	b   = data.frame(breedingDat@data, id = breedingDat@id, male = breedingDat@male, female = breedingDat@female, stringsAsFactors = FALSE)
	b$k = NULL
		
	# pre new() validity
	if( length(setdiff(polygonsDat@data[, 1], breedingDat@id ) ) > 0  )
      stop( "the 1st column of ", dQuote("polygonsDat"), " should be identical with ",  dQuote("breedingDat"), " id." )
            
    if( length(intersect(breedingDat@male, eppDat@male ) ) < 1 )
      stop("no extra-pair males found in breedingDat.")
    
	noepp = intersect(epd$z, paste(breedingDat@male, breedingDat@female) )
	if( length( noepp) > 0 ){ 
       warning("extra-pair partners cannot be social partners. The following pairs in eppDat are disregarded:\n", paste( sQuote(noepp), collapse = ",") ) }
    
    # build up epp set
    d = merge(hnb, b, by = "id") 
    d = merge(d, b, by.x = 'id_neigh', by.y = 'id',  all.x = TRUE, suffixes= c("_MALE","_FEMALE") )
    d$z = paste(d$male_MALE, d$female_FEMALE)    
    d = merge(d, epd, by = "z", all.x = TRUE)
	d[is.na(d$epp), "epp"] = 0
    d$z = NULL
    
    # fix names
    names(d) [which(names(d) == "male_MALE")] = "male"
    names(d) [which(names(d) == "female_FEMALE")] = "female"
    d$male_FEMALE = NULL; d$female_MALE = NULL    
	
	names(d) [which(names(d) == "id")] = "id_MALE"
	names(d) [which(names(d) == "id_neigh")] = "id_FEMALE"
	
    d = d[, union(c("id_FEMALE", "id_MALE", "rank", "male", "female", "epp"), names(d)) ]
    
	
	# post-merge validity
	eppInSet = apply(unique((d[d$epp == 1, c('male', 'female')] )), 1, paste, collapse = " ")
	lostEpPairs = setdiff(eppInSet,  epd$z)	
	
	if( length(lostEpPairs) > 0  ) {
	warning("something wicked happened merging datasets; some extra-pair partners are not in the final dataset:\n", paste( sQuote(lostEpPairs), collapse = ",") )
    }
	
	
	
    # new
	new("epp", breedingDat = breedingDat, polygonsDat = polygonsDat, eppDat = eppDat, maxlag = maxlag, EPP = d)
	
	
	
	
	}

if (!isGeneric("plot"))
  setGeneric("plot", function(x, y, ...)
    standardGeneric("plot"))

setMethod("plot", signature(x = "epp", y = "missing"),
          function(x, zoom, maxlag = 3, zoom.col = 'grey', ...) {
			
			p = x@polygonsDat
			b = x@breedingDat
			emat = x@eppDat
			e = x@EPP	
				
			if( !missing(zoom)) { 
				set = unique( c(zoom, 
					e[e$id_FEMALE%in%zoom & e$rank <= maxlag, 'id_MALE'], 
					e[e$id_MALE%in%zoom & e$rank <= maxlag, 'id_FEMALE']) 
					)
				
				p = p[p$ID%in%set, ]	
				
				bset = which(b@id%in%set)
				b = b[bset, ]
				b@male = b@male[bset]
				b@female = b@female[bset]
				b@id = b@id[bset]
        
        emat = e[ (e$id_FEMALE%in%set | e$id_MALE%in%set) & e$epp == 1, c("male", "female")]
				emat = eppMatrix(emat)
        
			}
				
		    plot(p, ...)
			if(!missing(zoom) )
				plot(p[p$ID == zoom, ], col = zoom.col, add = TRUE)
			plot(b, emat, add = TRUE, ...)
       
          })


if (!isGeneric("barplot")) {
    setGeneric("barplot", function(height,...)
      standardGeneric("barplot"))
   }  
    

setMethod("barplot", signature(height = "epp"),
          function(height, relativeValues = FALSE, ...) {

		  p = table(height@EPP[,c('rank', 'epp')])
            
            if(relativeValues == FALSE) {
                p = p[,2]
                plot(p, type = 'h', axes = FALSE, ylab ='No. of EPP events', xlab = 'Distance', ...)
                axis(1, at = 1:max(height@EPP$rank), labels = 1:max(height@EPP$rank))
                axis(2, at = 0:(max(p)), labels = 0:(max(p)))
              }
            
            if(relativeValues == TRUE) {
                p[,1] = p[,1]+p[,2]
                p = apply(p, MARGIN = 2, FUN = function(x) x/sum(x))
                plot(p[,2], type = 'h', axes = FALSE, ylab ='', xlab = '', ...)
                par(new = TRUE)
                plot(p[,1], type = 'l', axes = FALSE, ylab ='Proportion of EPP events', xlab = 'Distance', lty = 2, ...)
                axis(1, at = 1:max(height@EPP$rank), labels = 1:max(height@EPP$rank))
                axis(2, labels = (0:10)/10, at = (0:10)/10)  
            }
            
          })

if (!isGeneric("as.data.frame")) {
  setGeneric("as.data.frame", function(x)
    standardGeneric("as.data.frame"))
}	
	
setMethod('as.data.frame', signature(x='epp'), 
          function(x) {
            return(x@EPP)
          } )



	
	
	
	
	
	
	




















