reggeom <-
    function(matrx=matrix(c(
  0,5780,-1156,3468,3468,3468,-867,4335,0,0,-612,4080,5440,2652,3468,3420,3468,
  0,  0,4624,3468,3468, 0,3468,0,3468,4624,2448,1020,1360,3264,3264,3456,3456,
  0,  0,   0,4624,   0, 0,0,0,0,0,0,0,0,0,0,0,0),
		   nrow = 17,
		   ncol = 3),
      collab = c("U","V","W"),
      rowlab = c("o","x1","x2", "y", letters[2:8], "k","m","p","q","r","s"),
      colors = NULL,
      glyphs = NULL,
      erase = NULL,
      lines = matrix(c(1,6,8,1,11,7,1,1,5,6,6,15,17,8,5,9,1,9,10,
		     6,8,2,11,7,3,4,5,4,4,15,17,5,5,9,7,9,10,3),
		   nrow = 19,
		   ncol = 2),

      linecolors=c("red", "yellow", "yellow", "yellow", "yellow", "yellow",
		   "orchid", "green", "green", "red", "skyblue", "skyblue",
		   "skyblue", "white", "white", "white",
		   "slateblue", "slateblue", "slateblue"),

      resources=c("*showLines: True",
		  "*showAxes: False",
		  "*showPoints: False",
		  "*XGobi*PlotWindow.height: 500",
		  "*XGobi*PlotWindow.width: 500",
		  "*XGobi*VarPanel.width: 50"),

      title="Regression Geometry",

      vgroups	 = c(1,1,1),
      std	 = "msd",
##    dev	 = 1.5, #default is 2
      nlinkable	 = NULL,
      subset	 = NULL,
      display	 = NULL)
{
  xgobi(matrx=matrx,
	collab=collab,
	rowlab=rowlab,
	colors=colors,
	glyphs=glyphs,
	erase=erase,
	lines=lines,
	linecolors=linecolors,
	resources=resources,
	title=title,
	vgroups=vgroups,
	std=std,
	## dev=dev,
	nlinkable=nlinkable,
	subset=subset,
	display=display)
}

