###### Menu functions ##### 

plot.MenuItem <- function(MenuItem, x=0, y=0, pos = NULL, Plot=TRUE){

  	# if(!all(MenuItem$active))
		# return(NULL)

	if(is.null(MenuItem$TextParams))
		MenuItem$TextParams = list()
 
 
	XBufCof = MenuItem$RegionParams$XBuf
	YBufCof = MenuItem$RegionParams$YBuf
	Width = MenuItem$RegionParams$Width
	Height = MenuItem$RegionParams$Height
 
	if(is.null(XBufCof)){
		XBufCof = 1	
	}
	if(is.null(YBufCof)){
		YBufCof = 1
	}
	
	XBuf = XBufCof*do.call(strwidth, args=c(" ", MenuItem$TextParams))
	YBuf = YBufCof*do.call(strheight, args=c(" ", MenuItem$TextParams))/2
	
	if(is.null(Width)){
		Width = do.call(strwidth, args=c(MenuItem$label, MenuItem$TextParams)) + 2*XBuf
	}	
	if(is.null(Height)){
		Height = do.call(strheight, args=c(MenuItem$label, MenuItem$TextParams)) + 2*YBuf
	}	


	
	MenuItem$RegionParams$XBufCof = XBufCof
	MenuItem$RegionParams$YBufCof = YBufCof
	MenuItem$RegionParams$Width = Width
	MenuItem$RegionParams$Height = Height
	MenuItem$Region = c(x1=NA, x2=NA, y1=NA, y2=NA)

 
	if(Plot){
		if(is.null(pos)){
			pos = 1
		}
		
		
		switch(pos
			,{PlotX2 = x+Width;PlotY2 = y-Height; TextX=x+XBuf; TextY = y-Height/2}
			,{PlotX2 = x+Width;PlotY2 = y-Height; TextX=x+XBuf; TextY = y-Height/2}
			,{PlotX2 = x+Width;PlotY2 = y-Height; TextX=x+XBuf; TextY = y-Height/2}
			,{PlotX2 = x+Width;PlotY2 = y-Height; TextX=x+XBuf; TextY = y-Height/2}
		)
		
		xx = range(c(x, PlotX2))
		yy = range(c(y, PlotY2))
			
		
		do.call(rect, args=as.list(c(xleft=xx[1], ybottom=yy[1], xright=xx[2], ytop=yy[2], MenuItem$RecParams)))
		do.call(text, args=as.list(c(x=unname(TextX), y=unname(TextY), labels = MenuItem$label, pos = 4, offset = 0, MenuItem$TextParams)))
		
		MenuItem$Region = c(x1=xx[1], x2=xx[2], y1=yy[1], y2=yy[2])
	}

	return(MenuItem)
}
 
plot.MenuLine <- function(MenuLine, x=0, y=0, pos = NULL, Plot=TRUE){
  
  	# if(!all(MenuLine$active))
		# return(NA)
  
	CumWidth = 0
	MaxHeight = 0
	ind = rep(TRUE, length(MenuLine$MenuItems))
	for(i in seq_along(MenuLine$MenuItems)){
		MenuLine$MenuItems[[i]] = plot.MenuItem(MenuLine$MenuItems[[i]], x=0, y=0, pos = 1, Plot=FALSE)
		if(all((MenuLine$MenuItems[[i]]$active))){
			CumWidth = CumWidth + MenuLine$MenuItems[[i]]$RegionParams$Width	
			MaxHeight = max(c(MaxHeight,MenuLine$MenuItems[[i]]$RegionParams$Height))	
		} else {
			ind[i] = FALSE
		}
	}  
    MenuLine$MenuItems = MenuLine$MenuItems[ind]
  
  
	MenuLine$RegionParams$Width = max(c(CumWidth, MenuLine$RegionParams$Width) , na.rm = TRUE)
	MenuLine$RegionParams$Height = max(c(MaxHeight, MenuLine$RegionParams$Height) , na.rm = TRUE)

  
	if(Plot){
		if(is.null(pos)){
			pos = 1
		}
		
		
		LR = switch(pos
			,c(x=x, y=y)
			,c(x=x-MenuLine$RegionParams$Width, y=y)
			,c(x=x, y=y+MenuLine$RegionParams$Height)
			,c(x=x-MenuLine$RegionParams$Width, y=y+MenuLine$RegionParams$Height)
		)
		
		MenuLine$Region = c(x1=unname(LR["x"]), x2=unname(LR["x"]+MenuLine$RegionParams$Width), y1=unname(LR["y"]-MenuLine$RegionParams$Height), y2=unname(LR["y"]))
		

		do.call(rect, args=as.list(c(
			xleft=unname(MenuLine$Region["x1"]),
			ybottom=unname(MenuLine$Region["y1"]), 
			xright=unname(MenuLine$Region["x2"]),
			ytop=unname(MenuLine$Region["y2"]), 
			MenuLine$RecParams))
		)
			
		# res = list(Region=c(x1=NA, x2=x, y1=NA, y2=y))
		for(i in seq_along(MenuLine$MenuItems)){
			MenuLine$MenuItems[[i]]$RegionParams$Height = MenuLine$RegionParams$Height
			MenuLine$MenuItems[[i]] = plot.MenuItem(MenuLine$MenuItems[[i]], x=LR["x"], y=LR["y"], pos = 1, Plot=TRUE)
			LR["x"] = LR["x"] + MenuLine$MenuItems[[i]]$RegionParams$Width
		}  
	}
  
	return(MenuLine)
  
}
 
plot.Menu <- function(Menu, x=Menu$Params$x, y=Menu$Params$y, pos=Menu$Params$pos, Plot=Menu$Params$Plot){
 
	if(is.null(Menu))
		return(NULL) 
		
	if(!all(Menu$active))
		return(NULL)
	
	### pasirupinam nustatymais
	if(is.null(x)){
		x=0	
	}
	if(is.null(y)){
		y=0	
	}
	if(is.null(pos)){
		pos=NA	
	}
	if(is.null(Plot)){
		Plot=TRUE
	}	
	
	### jei isina uz rezius
	if(x<par("usr")[1])
		x=par("usr")[1]
	if(x>par("usr")[2])
		x=par("usr")[2]
	if(y<par("usr")[3])
		y=par("usr")[3]
	if(y>par("usr")[4])
		y=par("usr")[4]
		
		
	Menu$Params$x = x
	Menu$Params$y = y
	Menu$Params$pos = pos
	Menu$Params$Plot = Plot
	
	
	### susirandam ismatavimus
 	CumHeight = 0
	MaxWidth = 0
	ind = rep(TRUE, length(Menu$MenuLines))
	for(i in seq_along(Menu$MenuLines)){
		Menu$MenuLines[[i]] = plot.MenuLine(Menu$MenuLines[[i]], x=0, y=0, pos = 1, Plot=FALSE)
		if(all((Menu$MenuLines[[i]]$active))){
			CumHeight = CumHeight +  Menu$MenuLines[[i]]$RegionParams$Height
			MaxWidth = max(c(MaxWidth, Menu$MenuLines[[i]]$RegionParams$Width))	
		} else {
			ind[i] = FALSE
		}
	}  
	Menu$MenuLines = Menu$MenuLines[ind]
  
	Menu$RegionParams$Width = max(c(MaxWidth, Menu$RegionParams$Width) , na.rm = TRUE)
	Menu$RegionParams$Height = max(c(CumHeight, Menu$RegionParams$Height) , na.rm = TRUE)

	### pasirupiname pos
	if(is.na(pos)){
		pos = 1
		PlotX2 = x+Menu$RegionParams$Width
		PlotY2 = y-Menu$RegionParams$Height
		
		if(PlotX2>par("usr")[2]){
			pos = pos+1
			PlotX2 = x-Menu$RegionParams$Width
		}
		if(PlotY2<par("usr")[3]){
			pos = pos + 2
			PlotY2 = y+Menu$RegionParams$Height
		}
	}
	Menu$Params$pos = pos
	
 	if(Plot){
		if(is.null(pos)){
			pos = 1
		}
		
		LR = switch(pos
			,c(x=x, y=y)
			,c(x=x-Menu$RegionParams$Width, y=y)
			,c(x=x,  y=y+Menu$RegionParams$Height)
			,c(x=x-Menu$RegionParams$Width, y=y+Menu$RegionParams$Height)
		)	
		
		Menu$Region = c(x1=unname(LR["x"]), x2=unname(LR["x"]+Menu$RegionParams$Width), y1=unname(LR["y"]-Menu$RegionParams$Height), y2=unname(LR["y"]))
		
		do.call(rect, args=as.list(c(
			xleft=unname(Menu$Region["x1"]),
			ybottom=unname(Menu$Region["y1"]), 
			xright=unname(Menu$Region["x2"]),
			ytop=unname(Menu$Region["y2"]), 
			Menu$RecParams))
		)
		
			
		# res = list(Region=c(x1=NA, x2=x, y1=NA, y2=y))
		for(i in seq_along(Menu$MenuLines)){
			Menu$MenuLines[[i]]$RegionParams$Width = Menu$RegionParams$Width
			Menu$MenuLines[[i]] = plot.MenuLine(Menu$MenuLines[[i]], x=unname(LR["x"]), y=unname(LR["y"]), pos = 1, Plot=TRUE)	
			LR["y"] = LR["y"] - Menu$MenuLines[[i]]$RegionParams$Height
		}  
	}
  
	return(Menu) 
  
 }
 
is.InRegion <- function(x, y, Region){
	if(is.null(Region))
		return(FALSE)
	
	unname((Region[1] <= x & x <= Region[2]) & (Region[3] <= y & y <= Region[4]))
}
 
GetMenuItem <- function(x, y, MenuList){

	ans = NULL

	if(is.null(MenuList))
		return(ans)
		
	for(i in seq_along(MenuList)){
		if(is.InRegion(x, y, MenuList[[i]]$Region)){
			# cat("Menu!("%.%i%.%",")
			for(j in seq_along(MenuList[[i]]$MenuLines)){
				# MenuLine = MenuList[[i]]$MenuLines[[j]]
				if(is.InRegion(x, y, MenuList[[i]]$MenuLines[[j]]$Region)){
					# cat(j%.%",")	
					for(k in seq_along(MenuList[[i]]$MenuLines[[j]]$MenuItems)){
						MenuItem = MenuList[[i]]$MenuLines[[j]]$MenuItems[[k]]
						if(is.InRegion(x, y, MenuItem$Region)){
							ans = MenuItem
							break;
						}							
					}
					break;
				}
			}				
			break;
		}
	} 
	
	return(ans)
}
  
ClearBottomMenuPlot <- function(g){
	reg = g$Menu$ActiveMenuList$BottomMenu$Region
	rect(reg[1], reg[3], reg[2], reg[4], col=rgb(240,240,240,255,maxColorValue=255), border = NA)
}

### BottomMenu info eilutes tvarkymas
MenuMainMessage <- function(Msg, MenuLine='Info', g=g){
	
	g$Menu$MenuList$BottomMenu$MenuLines[[MenuLine]]$MenuItems[[1]]$label <- Msg
	g$Menu$MenuList$BottomMenu$MenuLines[['Info']]$active <- nchar(Msg)>0
	g
}

MenuClearExtraItems <- function(MenuLine='Info', g=g){
	g$Menu$MenuList$BottomMenu$MenuLines[[MenuLine]]$MenuItems <- g$Menu$MenuList$BottomMenu$MenuLines[[MenuLine]]$MenuItems[1]
	g$Menu$MenuList$BottomMenu$MenuLines[[MenuLine]]$active <- FALSE
	g
}

MenuAddItems <- function(Items, MenuLine='Info', g=g){
	g$Menu$MenuList$BottomMenu$MenuLines[[MenuLine]]$MenuItems <- c(g$Menu$MenuList$BottomMenu$MenuLines[[MenuLine]]$MenuItems, Items)
	if(length(Items)>0){
		g$Menu$MenuList$BottomMenu$MenuLines[[MenuLine]]$active = TRUE
	}
	g
}

MenuSelectEge <- function(eProgIds=NULL, g=g){	# taisytina
	cat("\nSelected edge: "%.% eProgIds %.%"\n")
	SelectEdge(eProgIds=eProgIds, g=g)
}

MenuSelectVertex <- function(vProgIds=NULL, g=g){ # taisytina	
	cat("\nSelected vertex: "%.% vProgIds %.%"\n")
	SelectVertex( vProgIds=vProgIds, g=g)
}

MenuSelectGroup <- function(id=NULL, oProgId=NULL, g=g){ # taisytina
	g = PutActiveObject(id=id, oProgId=oProgId, type="G", g=g)
	g
}

