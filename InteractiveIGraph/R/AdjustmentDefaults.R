##########################################################################################
################################## Default'tai ###########################################
PlotAdjustment.default <- function(g){
	# browser(skipCalls=2)
	### graziname sarasa su sluoksniuotais grafais
	ans = list()
	
#1	### pirmas sluoksninis - visi labai pilki
	gi = g
	V(gi)$color = "#EEEEEE"
	V(gi)$label.color="#CCCCCC"
	V(gi)$frame.color = NA
	E(gi)$color = "#EEEEEE"
	E(gi)$arrow.mode = 0
	E(gi)$label=NA
	E(gi)[!active]$lty = 2
	
	VO = GetViewObject(gi)
	if(VO$type!="none"){
		gi$PlotParam$main = paste("ViewObject: ", GetObjectLabel(oProgIds=VO$ProgId, type=VO$type, g=gi), sep="")
	}
	
	### issirasom 
	ans[[length(ans) + 1]] = gi

#2	### antras sluoksninis - tik aktyvus
	gi = induced.subgraph(g, vids=V(g)[active])
	gi = delete.edges(gi, edges=E(gi)[!active])
	# braisom taip kaip yra, tik sugraibom id
	
	### sugraibom grupes
	# browser()
	groups = list()
	for(i in seq_along(gi$groups)){
		if(all(gi$groups[[i]]$active)){		
			vids = as.numeric(V(gi)[ProgId %in% gi$groups[[i]]$MemberProgIds])
			if(length(vids)>0){
				if(gi$groups[[i]]$closed){	# jei uzdarytoje grupeje yra ir neuzdaryto kintamuju, tai grupes bloka irgi itraukim i pazymetus
					vids = c(vids, as.numeric(V(gi)[ProgId==gi$groups[[i]]$ProgId]))
				}
				groups[[length(groups)+1]] = vids
			}
		}	
	}
	if(length(groups)>0) gi$PlotParam$mark.groups = groups
	
	
	### issirasom 
	ans[[length(ans) + 1]] = gi
	
#3	### Trecias sluoksninis - tik pazymeti
	gg = g
	
	gg = induced.subgraph(gg, vids=V(gg)[active | selected])
	gg = delete.edges(gg, edges=E(gg)[!(active | selected)])	
	
	# pazykim tas krastines, kurios jungiasi su pazymetom virsunem
	gi = SelectEdge(vids=V(gg)[selected], g=gg)
	
		

	### normalus taskai
	V(gi)[selected]$label.color = "#FF0000"
	V(gi)[selected]$frame.color = "#FF0000"
	E(gi)[selected]$width <- 3	
	E(gi)[!selected]$color <- NA
	E(gi)$label = NA

	### Neaktyvus ir nepasymeti pasliapiami
	V(gi)[!active]$color = NA
	V(gi)[!active & !selected]$frame.color = NA
	V(gi)[!active & !selected]$label = NA
	E(gi)[selected & adj(V(gi)[!active])]$color = "darkgrey"
	E(gi)[E(gg)[selected]]$color = "red"
	E(gi)[!active]$lty = 2
	
	### jeigu pazymeta tik viena krastine, tai rodomi ir rodykliu lablai
	if(length(V(gi)[selected])==1){	
		E(gi)[selected]$label <- as.character(E(gi)[selected]$ProgId)
	}
	### Jei mode nurodo rodytas visas Ege, taip ir darom aktyviam objektui	
	if(gi$mode$AllEdges & GetActiveObject(g)$type=="V"){
		vid = V(g)[ProgId==GetActiveObject(g)$ProgId]
		eids = E(gi)[adj(V(g)[vid])]
		E(gi)[eids]$label <- as.character(E(gi)[eids]$ProgId)
	}
	
	# ### issirasom 
	ans[[length(ans) + 1]] = gi
	
	ans
}

ExtraInfo.default <- function(type="V", ProgId=0, g=g){	# taisytina
	
	msg = ""
	
	if(type=="V"){
		msg = ", V..."		
	}
	if(type=="E"){
		msg = ", E..."		
	}
	if(type=="G"){
		msg = ", G..."		
	}
		
	Item = list(list(label=msg, RegionParams=list(XBufCof = 0, YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8)))
	g <- MenuAddItems(Items=Item,  MenuLine='Info', g=g)	
	
	return(g)
}

BottomMenuAdjustment.default <- function(g){ # taisytina

	g = SelectEdge(vids=V(g)[selected], g=g)

	### issivalom, bedroji tvarka
	g <- MenuMainMessage(Msg = "", MenuLine='Info', g=g)
	g <- MenuClearExtraItems(MenuLine='A', g=g)
	g <- MenuClearExtraItems(MenuLine='G', g=g)

	
	MainInfo = ""
	ActiveObject = GetActiveObject(g)
	oProgId = ActiveObject$ProgId
	type = ActiveObject$type
	if(type != "none"){

		### pagrindinio info rodymas
		# virsunes
		if(type=="V"){
			MainInfo = "Active object: "%.% GetObjectLabel(oProgIds=oProgId, type=type, g=g)
			g <- do.call(g$functions$ExtraInfo, args=list(type="V", ProgId=oProgId, g=g))			
		}
		# krastai
		if(type=="E"){
			MainInfo = "Active object: "%.% GetObjectLabel(oProgIds=oProgId, type=type, g=g)				
			g <- do.call(g$functions$ExtraInfo, args=list(type="E", ProgId=oProgId, g=g))			
		}	
		# grupe
		if(type=="G"){
			MainInfo = "Active object: "%.% GetObjectLabel(oProgIds=oProgId, type=type, g=g)		
			g <- do.call(g$functions$ExtraInfo, args=list(type="G", ProgId=oProgId, g=g))
		}	
		g <- MenuMainMessage(Msg = MainInfo, MenuLine='Info', g=g)


		### Papildomu pasirinkimu(Menu isvedimas). Pagal aktyvu elementa.	
		# virsunes
		if(type=="V" | type=="G"){
			Items = list()
			vid = V(g)[ProgId==oProgId]
			for(ProgId in E(g)[adj(V(g)[vid])]$ProgId){
				Items[[as.character(ProgId)]] = list(label=as.character(ProgId), command="g <<- MenuSelectEge(eProgIds="%.%ProgId%.%", g=g)", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8, col=2))
			}					
			g <- MenuMainMessage(Msg = "E: ", MenuLine='A', g=g)
			g <- MenuAddItems(Items,  MenuLine='A', g=g)				
		}
		if(type=="E"){		
			### Todo: paradyti, kad rodytu rodykles krypti
			Items = list()
			eid = E(g)[ProgId==oProgId]
			from = V(g)[from(E(g)[eid])]$ProgId
			to = V(g)[to(E(g)[eid])]$ProgId
			
			separator = "--"
			
	
			if(!is.null(E(g)[eid]$arrow.mode)){
				if(E(g)[eid]$arrow.mode==1 | E(g)[eid]$arrow.mode=="<" | E(g)[eid]$arrow.mode=="<-") separator = "<-"
				if(E(g)[eid]$arrow.mode==2 | E(g)[eid]$arrow.mode==">" | E(g)[eid]$arrow.mode=="->") separator = "->"
				if(E(g)[eid]$arrow.mode==3 | E(g)[eid]$arrow.mode=="<>" | E(g)[eid]$arrow.mode=="<->") separator = "<->"
			} else {
				if(is.directed(g)){
					separator = "->"
				}
			}
			
			Items[[as.character(from)]] = list(label=as.character(from), command="g <<- MenuSelectVertex(vProgIds="%.%from%.%", g=g)", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8, col=2))
			Items[["separator"]] = list(label=separator, RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8, col=1))
			Items[[as.character(to)]] = list(label=as.character(to), command="g <<- MenuSelectVertex(vProgIds="%.%to%.%", g=g)", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8, col=2))
			g <- MenuMainMessage(Msg = "V: ", MenuLine='A', g=g)
			g <- MenuAddItems(Items,  MenuLine='A', g=g)				
		}		
		if(type=="G"){	
		}

	}
	
	### Sugraibom grupiu pasirinkimus
	# 1. rodom bendras grupes visiems elemnetams
	# 2. jei nera bendru grupiu, rodom tas grupes kurios yra issemtos pazymejimu
	# 3. jei ir to nera rodom grupe(s), kuriom priklauso daugiausiai pazymetu elementu 
	vProgIds = V(g)[selected]$ProgId
	if(length(vProgIds)>0 & length(g$groups)>0){
		gbs1 = gbs2 = rep(FALSE, length(g$groups))
		gnvs = rep(0, length(g$groups))		
		for(i in seq_along(g$groups)){
			VinG = vProgIds %in% g$groups[[i]]$MemberProgIds
			
			# 1. 
			if(all(VinG) & length(VinG)>0) gbs1[i] = TRUE
			
			# 2.
			if(all(g$groups[[i]]$MemberProgIds %in% vProgIds ) & length(g$groups[[i]]$MemberProgIds)>0) gbs2[i] = TRUE
		
			# 3.
			gnvs[i] = sum(VinG)			
		}		
		gids = vector("numeric", 0)
		if(any(gbs1)){
			gids = which(gbs1)
		} else {
			if(any(gbs2)){
				gids = which(gbs2)			
			} else {
				gbs3 = gnvs==max(gnvs)
				gbs3[gnvs==0] = FALSE
				if(any(gbs3)){
					gids = which(gbs3)	
				} 		
			}
		}		 
		Items = list()
		for(i in gids){		
			Items[[as.character(i)]] = list(label="<"%.%GetObjectAttribute(ids=i, Attr="name", type="G", g=g)%.%">", command="g <<- MenuSelectGroup(id="%.%i%.%", g=g)", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8, col=2))
		}					
		g <- MenuAddItems(Items,  MenuLine='G', g=g)			
	}
	
	g
}

