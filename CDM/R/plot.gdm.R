

##########################################################
# S3 plot method for gdm function
plot.gdm <- function( x , perstype="EAP" , group=1 , barwidth=.1 , histcol=1 ,
        cexcor=3 , pchpers=16 , cexpers=.7 , ... ){
	object <- x
    theta.k <- object$theta.k
    pi.k <- object$pi.k[,group]
    Ndim <- ncol(theta.k)
    mean.trait <- object$mean.trait[,group]	
    sd.trait <- object$sd.trait[,group]
    cor.trait <- object$correlation.trait[[group]]
    # extract person parameters
    person <- object$person[ object$group == group , ]
    person <- person[ , grep( paste0(perstype,".") , colnames(person) ) ]                            

    # define plot grid
    plotgrid <- as.data.frame( expand.grid( 1:Ndim , 1:Ndim )[,c(2,1) ] )
    plotgrid$type <- ""
    plotgrid[ plotgrid[,1] == plotgrid[,2] , "type" ] <- "hist"
    plotgrid[ plotgrid[,1] < plotgrid[,2] , "type" ] <- "cornumber"
    plotgrid[ plotgrid[,1] > plotgrid[,2] , "type" ] <- "scatterEAP"
    PG <- nrow(plotgrid)

    graphics::par( mfrow=c(Ndim,Ndim) )
    
    for (pp in 1:PG){
        if ( paste(plotgrid$type)[pp] == "cornumber"){
        plot.cornumber.gdm( cor.trait , dim1=plotgrid[pp,1] , dim2=plotgrid[pp,2] , cexcor)
                                            }
        if ( paste(plotgrid$type)[pp] == "scatterEAP"){
            plot.pers.gdm( person , dim1=plotgrid[pp,1] , dim2=plotgrid[pp,2] , pchpers , cexpers ,
					perstype )
                                }                                                        
        if ( paste(plotgrid$type)[pp] == "hist"){
            plot.hist.gdm( theta.k , pi.k , object , dim=plotgrid[pp,1] , group , barwidth , histcol ,
						mean.trait , sd.trait )
                    }
            }
    graphics::par( mfrow=c(1,1))
            }
########################################################################

#################################################################
# plot person parameters
plot.pers.gdm <- function( person , dim1 , dim2 , pchpers , cexpers , perstype ){
    dd1 <- dim1
    dd2 <- dim2
    graphics::plot( person[,dd1] , person[,dd2] , xlab=paste0(perstype , " Dim",dd1) ,
                ylab=paste0(perstype , " Dim",dd2) ,pch = pchpers , cex= cexpers)
                }
#################################################################


#########################################################################
# correlations between dimensions
plot.cornumber.gdm <- function( cor.trait , dim1 , dim2 , cexcor){
    dd1 <- dim1 
    dd2 <- dim2 
    graphics::plot( c(0,1) , c(0,1) , type="n" , axes=FALSE , xlab="" , ylab="")
    #text( .5 , .60 , paste0( "Cor(Dim" , dd1 , ",Dim" , dd2 ,")="   ) , cex= cexcor)
    graphics::text( .5 , .50 , paste0( round( cor.trait[dd1,dd2] ,3)  ) , cex= cexcor)
                            }
#########################################################################


#######################################################################
# histogram plot
plot.hist.gdm <- function( theta.k , pi.k , object , dim , group , barwidth , histcol ,
				mean.trait , sd.trait  ){
    dd <- dim
    a1 <- stats::aggregate( pi.k , list( theta.k[ , dd] ) , sum )
    mainpl <- paste0("Dim" , dd , " | M=" , round( mean.trait[dd] , 3 ) ,
                " | SD=" ,round( sd.trait[dd] , 3 ) )
    graphics::plot( a1[,1] , a1[,2] , type="n" , xlab= paste0("theta (Dim" , dd , ")" ),
                ylab="Probability" , main=mainpl)
    AA <- nrow(a1)            
    for ( aa in 1:AA){
        #aa <- 1
        graphics::rect(xleft = a1[aa,1] - barwidth/2, ybottom = 0 , 
                    xright=a1[aa,1] + barwidth/2, ytop = a1[aa,2] , col= histcol)
                }
        }
#######################################################################		