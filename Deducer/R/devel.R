# TODO: Add comment
# 
# Author: Ian
###############################################################################





addComponent <- function(container, component, top, right, bottom, left, topType = "REL", 
		rightType = "REL", bottomType = "REL", leftType = "REL"){
	AnchorConstraint <- J("org.rosuda.JGR.layout.AnchorConstraint")
	top <- as.integer(top)
	right <- as.integer(right)
	bottom <- as.integer(bottom)
	left <- as.integer(left)
	d <- getSize(container)
	ht <- d[2]
	wd <- d[1]
	if(topType=="ABS")
		top <- as.integer(round(ht*top/1000))
	if(rightType=="ABS")
		right <- as.integer(round(wd-wd*right/1000))
	if(bottomType=="ABS")
		bottom <- as.integer(round(ht-ht*bottom/1000))
	if(leftType=="ABS")
		left <- as.integer(round(wd*left/1000))
	
	topType <-	if(topType=="REL") AnchorConstraint$ANCHOR_REL else if(topType=="ABS") AnchorConstraint$ANCHOR_ABS else AnchorConstraint$ANCHOR_NONE
	rightType <-	if(rightType=="REL") AnchorConstraint$ANCHOR_REL else if(rightType=="ABS") AnchorConstraint$ANCHOR_ABS else AnchorConstraint$ANCHOR_NONE
	bottomType <-	if(bottomType=="REL") AnchorConstraint$ANCHOR_REL else if(bottomType=="ABS") AnchorConstraint$ANCHOR_ABS else AnchorConstraint$ANCHOR_NONE
	leftType <-	if(leftType=="REL") AnchorConstraint$ANCHOR_REL else if(leftType=="ABS") AnchorConstraint$ANCHOR_ABS else AnchorConstraint$ANCHOR_NONE
	container$add(component,new(AnchorConstraint,top,right,bottom,left,
					topType,rightType,bottomType,leftType))
	
	#container$add(component,top, right, bottom, left, topType, rightType, bottomType, leftType)
	
}

getSize <- function(component){
	dim <- component$getSize()
	c(dim$getWidth(),dim$getHeight())
}

setSize <- function(component,width,height){
	component$setSize(as.integer(width),as.integer(height))
	component$setPreferredSize(new(J("java.awt.Dimension"),as.integer(width),as.integer(height)))
}

execute <-function(cmd){
	if(!.jgr)
		DeducerMain$execute(cmd)
	else
		DeducerMain$executeAndContinue(cmd)
}