################################################################################
# summary method for objects of class "din"                                    #
################################################################################
plot.din <- 
function(x, items=c(1:ncol(x$data)), pattern="", 
  uncertainty=0.1, top.n.skill.classes=6, pdf.file="", 
  hide.obs = FALSE, display.nr=1:4, ask = TRUE, ...){

# Call: generic
# Input: 
#	x: object of class din
#	items:	an index vector giving the items to be visualized in the first plot
#	pattern: an optional character specifying a response pattern of an examinee, 
#     whose attributes are then analyzed in a seperate grafic.
#	uncertainty: a numeric between 0 and 0.5 giving the uncertanity bounds for 
#     deriving the observed skill occurrence probabilities in plot 2 
#	and the simplified deterministic attribute profiles in plot 5.
# top.n.skill.classes: a numeric, specifying the number of skill classes, 
#    starting with the most frequent, to be plotted in plot 3. Default value is 6.
#	pdf.file: an optional character string. If specified the graphics obtained 
#           from the function plot.din are provided in a pdf file.
# hide.idi: an optional logical value. If set to \code{TRUE}, the IDI curve in 
#           first graphic is not displayed.
# hide.obs: an optional logical value. If set to \code{TRUE}, the polygonal 
#           chain for observed frequencies of skill class probabilities in the 
#           second graphic is not displayed.
# display.nr: an optional numeric or numeric vector. If specified, only the 
#             plots in display.nr are shown. 
# ask: request a user input before the next figure is drawn; cf. ask option in ?par manual page
# Output: none

################################################################################
# check consistency of input                                                   #
################################################################################

hide.idi <- FALSE

	# define pattern in case of vector inputz
	if ( ( pattern[1] != "" ) & ( length( pattern) > 1 ) ){
		pattern <- paste( pattern , collapse="")
		}


	suc <- which(unique(x$pattern[,"pattern"]) == pattern) #subject under control
	if(pdf.file!="") try(pdf(file=pdf.file, ...))

	try(if(uncertainty<0||uncertainty>.5|| #uncertainty >= 0, <= .5
	  top.n.skill.classes<0||top.n.skill.classes>2^length(x$skill.patt) )  #top.n.skill.classes >= 0, <=2^K
		warning("check your plot parameter specifications. See Help-files."))   
	
	old.par <- graphics::par(no.readonly = TRUE)
	on.exit( graphics::par(old.par))
			
################################################################################
# Plot 1                                                                       #
################################################################################

if(1 %in% display.nr){
  
	# extract information
	errors <- rbind(x$guess[,1],x$slip[,1])[,items]
	colnames(errors) <- items
	errors[2,] <- 1 - errors[2,]
	
# 	if(!is.null(colnames(x$data)[items])){ 
# 		colnames(errors) <- colnames(x$data)[items]
# 	}else{
# 		colnames(errors) <- items
# 	}

	# generate plot
	graphics::barplot(errors, ylim=c(0,1.19), beside=TRUE, col=c("gray","darkred"), 
	        xlab="Item index", ylab="Probability", cex.lab=1.3)
	if(!hide.idi){
#    if(any(apply(errors, 2, function(x) 1-x[1] - x[2] < 0 ))){
    if( FALSE ){
#		  warning(paste("Item discrimination index < 0 for item",
#				which(apply(errors, 2,  function(x) 1-x[1]-x[2] < 0 )),"\n"))
    }else{
      graphics::legend("topright",c("Guessing probability","Non-Slipping probability"), 
             lty=c(1,1), pch=c(NA,NA), lwd=c(2,2), col=c("gray","darkred"), bg = "gray97")
#          lines(seq(2,2+3*(ncol(errors)-1),3),apply(errors, 2, function(x) 1-x[1]-x[2] ), lty=1)
#		  points(seq(2,2+3*(ncol(errors)-1),3),apply(errors, 2, function(x) 1-x[1]-x[2] ), pch=19, cex=1.5)
# 		  legend("topright",c("guessing parameter","slipping parameter", "item discrimination index"), 
# 	      lty=c(1,1,1), pch=c(NA,NA,19), lwd=c(2,2,2), col=c("gray","darkred", "black"), bg = "gray97")
		}
  }else{
	  graphics::legend("topright",c("guessing parameter","slipping parameter"), 
	    lty=c(1,1), lwd=c(2,2), col=c("gray","darkred"), bg = "gray97")
	}
  
	if(pdf.file=="" & ask) graphics::par(ask=TRUE)
	if(1 == display.nr[length(display.nr)]) graphics::par(ask=FALSE)
}
	   
################################################################################
# Plot 2                                                                       #
################################################################################

if(2 %in% display.nr){
  
  # extract information
	skill.patterns <- x$skill.patt[length(x$skill.patt):1,]
  
  ind <- match( apply(x$item.patt.split, 1, paste, collapse ="") , 
                unique(x$pattern)[,"pattern"] )
  EAP <- ifelse( unique(x$pattern)[ ind, grep("post.attr", colnames(x$pattern) ) ] > 
    0.5 + uncertainty, 1, NA )
  
	master <- colSums( apply( EAP , 2, function(y) y*x$item.patt.freq), na.rm=TRUE )
  master <- ( master/ nrow(x$data) )[length(x$skill.patt):1]

	# generate plot
	graphics::par(yaxt="n")
	graphics::barplot(skill.patterns, horiz=TRUE, 
#			ylim=c(0,length(skill.patterns)*1.2+0.9),
			ylim=c(0,length(skill.patterns)*1.2),			
			xlim=c(0,1), 
			xlab="Skill mastery probability", axes=F, cex.lab=1.3, col="gray")
	graphics::axis(1,at=seq(0,1,0.2))
	graphics::axis(3 ,at=seq(0,1,0.2))
	
	if ( is.null( attributes(x$q.matrix)$skill.labels ) ){
		attr( x$q.matrix , "skill.labels") <- colnames(x$q.matrix )
						}
	
	graphics::text(attributes(x$q.matrix)$skill.labels[length(x$skill.patt):1], 
       x=c(rep(0.01,length(skill.patterns))), y=seq(0.7,
        0.7+1.2*(length(skill.patterns)-1),1.2), col="black", pos=4, cex=1.3)
	
  if(!hide.obs){
#			legend("topright",c("Marginal skill probability", "Percentage of masters (EAP)"),
#           lty=c(1,1), pch=c(NA,19), lwd=c(2,2), col=c("gray", "black"), bg = "gray97")
#    points(x=master, y=seq(0.7,0.7+1.2*(length(skill.patterns)-1),1.2),pch=19, cex=1.3)
#    lines(x=master, y=seq(0.7,0.7+1.2*(length(skill.patterns)-1),1.2),lty=1)
  }
	
	if(pdf.file=="" & ask) graphics::par(ask=TRUE)
	if(2 == display.nr[length(display.nr)]) graphics::par(ask=FALSE)
	
	graphics::par(yaxt="s")
}
	  
################################################################################
# Plot 3                                                                       #
################################################################################

if(3 %in% display.nr){
  
  # extract information
	patt.fq <- x$attribute.patt[,1]
	main.effects <- which(rownames(x$attribute.patt)%in%
	  rownames(x$attribute.patt[order(x$attribute.patt[,1], decreasing = TRUE),][
	  1:min(top.n.skill.classes, 2^length(x$skill.patt)), ])
	    )
	
	# generate plot
	graphics::par(xaxt="n"); graphics::par(mar=c(6,4,4,2)+0.1)
	graphics::plot(c(0:(length(patt.fq)+1)),c(0,t(patt.fq),0),type="h", 
	   ylab="Skill class probability",
		xlab="", ylim=c(0,max(patt.fq)+.02),cex.lab=1.3, col=c(NA ,rep("black",length(patt.fq)), NA))
	graphics::par(xaxt="s")
	graphics::axis(1, at=main.effects, las=2, labels=rownames(x$attribute.patt)[main.effects], cex.axis=.8)
	
# print(patt.fq)	
	eps <- .2
	PP <- length(patt.fq)
	if (PP<65){
	for (pp in 1:PP){
#		pp <- 1
		graphics::rect( xleft= pp-eps , ybottom=0 , xright=pp+eps , ytop=patt.fq[pp] , col="black" )
					} }
  graphics::par(mar=c(5,4,4,2)+0.1)
  
	if(pdf.file=="" & ask) graphics::par(ask=TRUE)
	if(3 == display.nr[length(display.nr)]) graphics::par(ask=FALSE)
	
}

################################################################################
# Plot 4                                                                       #
################################################################################

if(4 %in% display.nr){   
	if(pattern!=""){
    if(length(suc) == 0)
      warning("The specified pattern was not achieved.")
         
		# if a pattern is specified extract information
		post.skill <- as.matrix(unique(x$pattern)[suc ,grep("post.attr", colnames(x$pattern))])[nrow(x$skill.patt):1]
		names(post.skill) <- colnames(x$q.matrix)[nrow(x$skill.patt):1]
		
		# generate plot
    graphics::par(mar=c(5,4,4,2)+0.1)
    graphics::par(mgp=c(3.5,1,0))
    graphics::par(yaxt="n")
    graphics::barplot(post.skill, horiz=TRUE, xlab=paste("Skill probabilities conditional on response pattern\n",
			pattern ,sep=""), 
			xlim=c(0,1), axes=FALSE, cex.lab=1.3, col="gray")

		graphics::axis(1,at=seq(0,1,0.2))
		graphics::abline(v=c(.5-uncertainty,.5+uncertainty),lty=1,col="darkred" , lwd=2)
		graphics::axis(3, at=c( (.5-uncertainty)/2, .5, .5+uncertainty+(1-(.5+uncertainty))/2 ), tick=F, 
			labels=c("not mastered", "unclassified", "mastered"), cex.axis=1.3, mgp=c(3,0,0))
			
	if ( is.null( attributes(x$q.matrix)$skill.labels ) ){
		attr( x$q.matrix , "skill.labels") <- colnames(x$q.matrix )
						}			
			
		graphics::text(attributes(x$q.matrix)$skill.labels[length(x$skill.patt):1], 
         x=c(rep(0.01,length(row.names(x$skill.patt)))),
         y=seq(0.7,0.7+1.2*(length(row.names(x$skill.patt))-1),1.2),col="black", pos=4, cex=1.3)
		graphics::par(yaxt="s")
    graphics::par(mar=c(5,4,4,2)+0.1)
    graphics::par(mgp=c(3,1,0))
    
    if(4 == display.nr[length(display.nr)]) graphics::par(ask=FALSE)  
    
	}
	
}
	# reset open plot parameter
	if(pdf.file!="") try(dev.off())     
	graphics::par(old.par)
	invisible()
}


