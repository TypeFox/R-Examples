
spinControl <- function(base, dev = rgl.cur(), continue=FALSE,
			speed=30, scale=100){

	slider <- tclVar(speed)
	getZooms <- function() {
	 	old <- rgl.cur()
	    	on.exit(rgl.set(old))
	    	result <- numeric(max(dev))
	    	for (device in dev) {
	    	    rgl.set(device)
	 	    result[device] <- par3d("zoom")
	 	}
	 	result
	}
	zooms <- getZooms()
		
	scale <- tclVar(scale)
	continuous <- tclVar(as.numeric(continue))
	
	buttonPress <- NULL
	direction <- NULL
	
	lastTime <- proc.time()[3]
	timeDiff <- 0

	rotateUp <- function(){
		angle <- timeDiff*as.numeric(tclObj(slider))*pi/180
		par3d(userMatrix = rotationMatrix(-angle, 1,0,0) %*% par3d("userMatrix"))
	}

	rotateLeft <- function(){
		angle <- timeDiff*as.numeric(tclObj(slider))*pi/180
		par3d(userMatrix = rotationMatrix(-angle, 0,1,0) %*% par3d("userMatrix"))
	}

	rotateRight <- function(){
		angle <- timeDiff*as.numeric(tclObj(slider))*pi/180
		par3d(userMatrix = rotationMatrix(angle, 0,1,0) %*% par3d("userMatrix"))
	}

	rotateSpin <- function(){
		angle <- timeDiff*as.numeric(tclObj(slider))*pi/180
		par3d(userMatrix = rotationMatrix(-angle, 0,0,1) %*% par3d("userMatrix"))
	}

	rotateDown <- function(){
		angle <- timeDiff*as.numeric(tclObj(slider))*pi/180
		par3d(userMatrix = rotationMatrix(angle, 1,0,0) %*% par3d("userMatrix"))
	}	

	rotate <- function(){
	    	old <- rgl.cur()
	    	on.exit(rgl.set(old))	
	    	if (buttonPress) {
		    if ((currentTime <- proc.time()[3]) > lastTime) {
			timeDiff <<- currentTime - lastTime
			lastTime <<- currentTime
			for (device in dev) {
			    rgl.set(device)
			    if (direction == "up")
				    rotateUp()
			    if (direction == "left")
				    rotateLeft()
			    if (direction == "spin")
				    rotateSpin()
			    if (direction == "right")
				    rotateRight()
			    if (direction == "down")
				    rotateDown()
			}
		    }
		    tcl("after",5,rotate)
	    }

	}

	# rotation button callback functions
	# note that "..." argument is necessary
	upButtonPress <- function(...){
		buttonPress <<- TRUE
		lastTime <<- proc.time()[3]
		direction <<- "up"
		rotate()
	}

	leftButtonPress <- function(...){
		buttonPress <<- TRUE
		lastTime <<- proc.time()[3]
		direction <<- "left"
		rotate()
	}

	spinButtonPress <- function(...){
		buttonPress <<- TRUE
		lastTime <<- proc.time()[3]
		direction <<- "spin"
		rotate()
	}

	rightButtonPress <- function(...){
		buttonPress <<- TRUE
		lastTime <<- proc.time()[3]
		direction <<- "right"
		rotate()
	}

	downButtonPress <- function(...){
		buttonPress <<- TRUE
		lastTime <<- proc.time()[3]
		direction <<- "down"
		rotate()
	}

        onIdle <- function(...){
        	buttonPress <<- TRUE
        	rotate()
        	buttonPress <<- FALSE
		if (as.numeric(tclObj(continuous)))
		    tcl("after", "idle", onIdle)
	}
        	
	buttonRelease <- function(...){
		buttonPress <<- FALSE
		if (as.numeric(tclObj(continuous)))
		    tcl("after", "idle", onIdle)
	}


	resetAxis <- function(...){
	    	old <- rgl.cur()
	    	on.exit(rgl.set(old))
	    	for (device in dev) {
		    rgl.set(device)
		    par3d(userMatrix = diag(4))
		}
	}

	setScale <- function(...){
	    	old <- rgl.cur()
	    	on.exit(rgl.set(old))		    
	    	scale <- as.numeric(tclObj(scale))
	    	for (device in dev) {
		    rgl.set(device)
		    par3d(zoom = 10^((scale - 100)/50)*zooms[device])
		}
	}

	spec.frm <- tkframe(base)
	first.frm <- tkframe(spec.frm)
	second.frm <- tkframe(spec.frm)
	third.frm <- tkframe(spec.frm)
	fourth.frm <- tkframe(spec.frm)

	# rotations buttons
	upButton <- tkbutton(first.frm, text="  ^  ")
	leftButton <- tkbutton(first.frm, text="  <  ")
	spinButton <- tkbutton(first.frm, text="  O  ")
	rightButton <- tkbutton(first.frm, text="  >  ")
	downButton <- tkbutton(first.frm, text="  v  ")


	tkgrid(tklabel(first.frm, text=" "), upButton,	tklabel(first.frm, text=" "))
	tkgrid(leftButton,spinButton,rightButton)
	tkgrid(tklabel(first.frm, text=" "), downButton, tklabel(first.frm, text=" "))

	tkbind(upButton, "<Button-1>", upButtonPress)
	tkbind(leftButton, "<Button-1>", leftButtonPress)
	tkbind(spinButton, "<Button-1>", spinButtonPress)
	tkbind(rightButton, "<Button-1>", rightButtonPress)
	tkbind(downButton, "<Button-1>", downButtonPress)

	tkbind(upButton,"<ButtonRelease-1>", buttonRelease)
	tkbind(leftButton,"<ButtonRelease-1>", buttonRelease)
	tkbind(spinButton,"<ButtonRelease-1>", buttonRelease)
	tkbind(rightButton,"<ButtonRelease-1>", buttonRelease)
	tkbind(downButton,"<ButtonRelease-1>", buttonRelease)

	# control buttons
	frameAxis <- tkframe(second.frm,borderwidth=2)
	tkpack(frameAxis)

	buttonAxis <- tkbutton(frameAxis,text="     Reset Axis     ",command=resetAxis)
	contBox <- tkcheckbutton(frameAxis,text="Continue rotating", variable=continuous)
	
	tkpack(contBox, buttonAxis,fill="x")

	#control scale frame
	frameControl <- tkframe(third.frm,borderwidth=4)
	tkpack(frameControl)
	frameControlSpeed <- tkframe(frameControl)
	sliderSpeed <- tkscale(frameControlSpeed, showvalue=FALSE,
			orient="horiz", from=0, to=400,resolution=1,variable=slider)
	tkpack(tklabel(frameControlSpeed,text="Speed:"),sliderSpeed,side="left")

	frameControlScale <- tkframe(frameControl)
	sliderScale <- tkscale(frameControlScale,showvalue=FALSE,orient="horiz",
			from=0, to=200,resolution=1, variable=scale, command=setScale)
	tkpack(tklabel(frameControlScale,text="Scale:"),sliderScale,side="left")
	tkpack(frameControlSpeed,frameControlScale)


	tkpack(first.frm, second.frm, third.frm, fourth.frm)
	tkpack(spec.frm,expand=TRUE)

	spec.frm
}


spin3d <- function(dev = rgl.cur(), ...){

	base <- tktoplevel()
	tkwm.title(base, "Spin")
	
	spin <- spinControl(base, dev, ...)		

	quit <- tkbutton(base,text="Quit", command=function()tkdestroy(base))

	tkpack(spin, quit)

}

