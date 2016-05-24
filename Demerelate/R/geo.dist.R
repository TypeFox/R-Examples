geo.dist <- function(pop1, pop2, onlypairs=FALSE, value)
  { 

	# Calculates distances between samples omitting any population information.
	# Input:
  
  # coord.x/y may be "relative" or "latitude" and "longitude" in decimal degrees
  # 	ind1 	pop coord.x coord.y
	# ind1 	EGB 	2 	3	 
	# ind2 	EGB 	2 	3
	# ind3 	EGB 	2 	3
	# .    	. 	. 	. 	.  
  
  # Returns:
  # matrix.share 
	# 	ind1 	ind2 	ind3 	ind4
	# ind1 	4 	0 	2 	3	 
	# ind2 	4 	0 	2 	3
	# ind3 	4 	0 	2 	3
	# .    	. 	. 	. 	.


	# Erstelle matrix	
	matrix.share <- matrix(numeric(0),nrow=length(pop1[,1]),ncol=length(pop2[,1]))
	row.names(matrix.share) <- pop1[,1]
	colnames(matrix.share) <- pop2[,1]
	pop.size <- length(pop1[,1])+length(pop2[,1])-2

	if (onlypairs==FALSE)
		
		{

			
	for (i in 1:(length(matrix.share)))

		{


		# Berechnet aus vortlaufender Nummerierung der Matrix anhand der Reihenanzahl die Spalten und Reihen Position um Namen abzugreifen
			col.position <- ceiling(i/length(row.names(matrix.share)))  #rundet auf
			row.position <- (i-(length(row.names(matrix.share))*(ceiling(i/length(row.names(matrix.share)))-1)))
			
			# Berechne Allelvergleiche nach TRUE/FALSE fuer Vergleich a,b,c,d und gibt als as.numeric 1/0
			if (col.position<row.position)
			{

	
  	# Nach Pythagoras
  	if (value=="relative") {	matrix.share[row.position,col.position] <- 
          sqrt(
            (pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),3]-pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),3])^2
            +
            (pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),4]-pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),4])^2
                )
                        }
    
    if (value=="decimal") {   matrix.share[row.position,col.position] <-
          6378137*((
                      2*asin(sqrt((sin((
                                          180/pi*(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),3])-
                                          180/pi*(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),4]))/2)^2)+cos(
                                          180/pi*(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),3]))*cos(
                                          180/pi*(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),4]))*(sin((
                                          180/pi*(pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),3])-
                                          180/pi*(pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),4]))/2)^2))
                             )
                      )
                   )
}
	
		}
		}

		
}
	if (onlypairs==TRUE)
	{
	     
	     for (i in 1:length(colnames(matrix.share)))
        	{

			# Nach Pythagoras
			if (value=="relative") 
				{
         matrix.share[i,i] <- sqrt(
            (pop1[which(pop1[,1]==row.names(matrix.share)[i]),3]-pop2[which(pop2[,1]==colnames(matrix.share)[i]),3])^2
            +
            (pop1[which(pop1[,1]==row.names(matrix.share)[i]),4]-pop2[which(pop2[,1]==colnames(matrix.share)[i]),4])^2
                )
				}
      
      if (value=="decimal") 
        {
        matrix.share[row.position,col.position] <-
            6378137*(
                      (2*asin(sqrt(
                                    (sin((
                                  180/pi*(pop1[which(pop1[,1]==row.names(matrix.share)[i]),3])-
                                  180/pi*(pop1[which(pop1[,1]==row.names(matrix.share)[i]),4]))/2)^2)+
                                    cos(
                                  180/pi*(pop1[which(pop1[,1]==row.names(matrix.share)[i]),3]))*
                                    cos(
                                  180/pi*(pop1[which(pop1[,1]==row.names(matrix.share)[i]),4]))*
                                    (sin((
                                  180/pi*(pop2[which(pop2[,1]==colnames(matrix.share)[i]),3])-
                                  180/pi*(pop2[which(pop2[,1]==colnames(matrix.share)[i]),4]))/2)^2)
                                   )
                              )
                       )
                      )

			
        	}

		
		}
	}

		return(matrix.share)

	}
