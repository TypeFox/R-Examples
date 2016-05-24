
GroupFix <- function(gid=NULL, gProgId=NULL, g=g){
	if(length(gProgId)>0) gid = which(names(g$groups)==as.character(gProgId))

	group = g$groups[[gid]]
	
	if(length(group$active)==0) group$active = TRUE
	if(length(group$closed)==0) group$closed = FALSE
	if(length(group$name)==0) group$name = "group_"%.%gid
	if(length(group$label)==0) group$label = group$name
	if(length(group$ProgId)==0) group$ProgId = max(V(g)$ProgId)+1

	### pakartojam:
	if(any(is.na(group$active))) group$active = TRUE
	if(any(is.na(group$closed))) group$closed = FALSE
	if(any(is.na(group$name))) group$name = "group_"%.%gid
	if(any(is.na(group$label))) group$label = group$name
	if(any(is.na(group$ProgId))) group$ProgId = max(V(g)$ProgId)+1
	
	if(group$closed) group$active = TRUE	### uztikrinam defolta, jei grupe uzdaryta, tai sis atributas labiau taikytinas verrtexui
	
	if(group$ProgId %in% V(g)$ProgId){ ### patikrinam ar sis prog id tikrai geras
		warning("Wrong group's 'ProgId' provided. It allready taken. New id was generated.")
		group$ProgId = max(V(g)$ProgId)+1
	}
	
	if(is.null(group$MemberProgIds))  group$MemberProgIds = vector("numeric", 0)
	group$MemberProgIds = group$MemberProgIds[!is.na(group$MemberProgIds)]
	
	if(is.null(group$StartX) | is.null(group$StartY)){
		if(length(group$MemberProgIds)>0){
			group$StartX = mean(V(g)[ProgId %in% group$MemberProgIds]$x)
			group$StartY = mean(V(g)[ProgId %in% group$MemberProgIds]$y)
		} else {
			group$StartX = 0
			group$StartY = 0
		} 	
	}
	
	group
}

GroupCreate <- function(g, GroupsInfoList=NULL, vProgIds=vector(), ...){
	### susirandam dabartini (busima) indeksa
	Ng = length(g$groups)+1

	VertexAttributes = list(...)
	VActive = TRUE
	
	### priklausomai nuo to ar tai tik vProgIds ar visas info sarasas elgiamas skirtingai. GroupsInfoList dominuoja
	if(!is.null(GroupsInfoList)){
		GroupsParam = c("ProgId", "MemberProgIds", "active", "name", "label", "closed", "StartX", "StartY")
		ind = names(GroupsInfoList) %in% GroupsParam
		g$groups[[Ng]] = GroupsInfoList[ind]
		
		VertexParam = c("x", "y", "size", "size2", "color", "frame.color", "shape", "label", "label.family", "label.font", "label.cex", "label.dist", "label.degree", "label.color")
		ind = names(GroupsInfoList) %in% VertexParam
		
		VertexAttributes = GroupsInfoList[ind]
		### issivalom nuo NULL ir NA
		for(i in rev(seq_along(VertexAttributes))){
			if(is.null(VertexAttributes[[i]])){
				VertexAttributes[[i]] = NULL
			} else {
				if(is.na(VertexAttributes[[i]])){
					VertexAttributes[[i]] = NULL
				}
			}
		}		
		### specialiai del aktyvumo jis yra bieju parametras. Jei grupe uzdaryta, tai tampa verteko parametru
		if(!is.null(g$groups[[Ng]]$active) & !is.null(g$groups[[Ng]]$closed)){
			if(!is.na(g$groups[[Ng]]$active) & !is.na(g$groups[[Ng]]$closed)){
				if(g$groups[[Ng]]$closed){
					VActive = g$groups[[Ng]]$active
				}	
			}
		}		
	} else {
		g$groups[[Ng]] = list(MemberProgIds = vProgIds)
	}
	 # browser()
	g$groups[[Ng]] = GroupFix(gid = Ng, g = g)
	names(g$groups)[Ng] = g$groups[[Ng]]$ProgId
	

	### pridedam vertex'a
	G = g$groups[[Ng]]
	attr = list(name=G$name, label=G$label, ProgId=G$ProgId, x=G$StartX, y=G$StartY
		, active=(G$active & VActive), hidden=!G$closed, selected=FALSE
		, ProgType = "groupe", shape = "crectangle", size=30, size2=15,  color='grey'
		, frame.color = "black")
	attr[names(VertexAttributes)] = VertexAttributes # overadinam rekmes	
	g <- add.vertices(g, 1, attr=attr)		
	
	### pasirupiname, kad visi Vertexai butu sutvarkyti,
	if(G$closed){
		g = GroupClose(gProgId=G$ProgId, g=g, passive=TRUE)
	} else {
		g = GroupOpen(gProgId=G$ProgId, g=g, passive=TRUE)
	}

	if(class(g)[1]!="InteractiveIGraph") class(g) <- c("InteractiveIGraph", class(g))
	g
}

GroupOpen <- function(gid=NULL, gProgId=NULL, g, passive=FALSE){
	if(length(gProgId)>0) gid = which(names(g$groups)==as.character(gProgId))

	G = g$groups[[gid]]	
	gvid = as.numeric(V(g)[ProgId==G$ProgId])[1]	
	vids = as.numeric(V(g)[(ProgId %in% G$MemberProgIds) & (block == G$name)])

	
	
	V(g)[vids]$hidden = FALSE
	V(g)[vids]$block = ""
	
	V(g)[vids]$x = V(g)[vids]$x - (G$StartX - V(g)[gvid]$x)
	V(g)[vids]$y = V(g)[vids]$y - (G$StartY - V(g)[gvid]$y)
	
	g$groups[[gid]]$closed = FALSE
	
	### susitvarkom su grupem
	if(!passive){
		for(i in seq_along(g$groups)){
			# gi = g$groups[[i]]
			# visos virusnes perkialaimos i tas grupes, kurioms priklauso blokas
			if(G$ProgId %in% g$groups[[i]]$MemberProgIds){	
				g <- GroupMembersAdd(gid=i, VProgIds=V(g)[vids]$ProgId, g=g)
			}
		
			### Jei visos priklause vir-grupei, bet blokas ne, tai virsunes irgi pasalinamos is vir-grupei(isskyrus aisku savaja)
				# tipo, jei blokas buvo ismestas is grupes
			if(all(V(g)[vids]$ProgId %in% g$groups[[i]]$MemberProgIds) & !(G$ProgId %in% g$groups[[i]]$MemberProgIds) & gid!=i){
				g <- GroupMembersDelete(gid=i, VProgIds=V(g)[vids]$ProgId, g=g)
			}

			### o pats blokas turi dingti is grupes sarasu. 
			g <- GroupMembersDelete(gid=i, VProgIds=G$ProgId, g=g)
		}
		g <- PutActiveObject(id=gid, type="G", g=g)
	}
	
	g <- delete.edges(g, E(g)[adj(V(g)[gvid])])
	V(g)[gvid]$hidden = TRUE
	
	
	### perreinam per grupes ir jei kuri nors tucia tai ja trinam
	# DeleteGroups = rep(FALSE, length(g$groups))
	# for(i in seq_along(g$groups)){
		# DeleteGroups[i] = length(g$groups[[i]]$MemberProgIds)==0	
	# }
	# g$groups[DeleteGroups] <- NULL
	

	g <- XY.norm(g)	

	if(class(g)[1]!="InteractiveIGraph") class(g) <- c("InteractiveIGraph", class(g))	
	g	
}

GroupClose <- function(gid=NULL, gProgId=NULL, g, passive=FALSE){
	

	if(length(gProgId)>0) gid = which(names(g$groups)==as.character(gProgId))

	vids = as.numeric(V(g)[ProgId %in% g$groups[[gid]]$MemberProgIds & block==""])
	if(length(vids)>0){	
		V(g)[vids]$hidden = TRUE
		V(g)[vids]$block = g$groups[[gid]]$name
		
		if(!passive){	
			g$groups[[gid]]$StartX = mean(V(g)[vids]$x)
			g$groups[[gid]]$StartY = mean(V(g)[vids]$y)
		}
	} else {
		if(!passive){	
			g$groups[[gid]]$StartX = 0
			g$groups[[gid]]$StartY = 0
		}
	}
	
	g$groups[[gid]]$closed = TRUE
	
	### pasirupinam vertexo atributais
	gvid = as.numeric(V(g)[ProgId==g$groups[[gid]]$ProgId])[1]
	V(g)[gvid]$hidden = FALSE
	if(!passive){
		V(g)[gvid]$active = TRUE
		V(g)[gvid]$selected = TRUE	
		V(g)[gvid]$x = g$groups[[gid]]$StartX
		V(g)[gvid]$y = g$groups[[gid]]$StartY			
	}			
		
	### pridedame krastines 
	G = g$groups[[gid]]
	MemberVids = as.numeric(V(g)[ProgId %in% G$MemberProgIds])
	eids = E(g)[adj(MemberVids)]
	if(length(eids)>0){
		EfromDF = as.data.frame(get.edgelist(g, names=FALSE)[eids,], stringsAsFactors = FALSE)
		FromB = EfromDF[,1] %in% MemberVids

		EAList = list.edge.attributes(g)
		for(A in EAList){
			EfromDF[A] = get.edge.attribute(g, name=A, index=eids)
		}
		EfromDF$name = "B."%.%EfromDF$name 
		EfromDF$ProgId = max(E(g)$ProgId) + seq_along(EfromDF$ProgId)
		EfromDF$ProgType = "block" 

		EfromDF[FromB,1] = gvid
		EfromDF[!FromB,2] = gvid
		rem = duplicated(EfromDF[,1:2])
		EfromDF = EfromDF[!rem,]

		for(i in 1:nrow(EfromDF)){
			g = add.edges(g, as.vector(as.matrix(EfromDF[i,1:2])), attr=as.list(EfromDF[i,-(1:2)]))
		}
	}
		
	### pereinam per grupes, jei yra grupe, kurioje visi kintamieji priklaus blokui, tai ir blokas jai priklauso
	for(i in seq_along(g$groups)){
		if(all(vids %in% g$groups[[i]]$MemberProgIds) & i!=gid & length(vids)>0){	
			g <- GroupMembersAdd(gid=i, VProgIds=V(g)[gvid]$ProgId, g=g)
		}
	}
	if(class(g)[1]!="InteractiveIGraph") class(g) <- c("InteractiveIGraph", class(g))	
	g
}

GroupDelete <- function(gid=NULL, gProgId=NULL, g){
	
	if(length(gProgId)>0) gid = which(names(g$groups)==as.character(gProgId))
	
	### visu pirma atidarom
	if(all(g$groups[[gid]]$closed)){
		g <- GroupOpen(gid=gid, g=g)
	}

	### istrinam vertexa
	vid = as.numeric(V(g)[ProgId==g$groups[[gid]]$ProgId])[1]
	g <- delete.vertices(g, vid)
	g <- XY.norm(g)
	
	### istrinam grupe
	g$groups[[gid]] = NULL
	
	if(class(g)[1]!="InteractiveIGraph") class(g) <- c("InteractiveIGraph", class(g))	
	g

}


GroupMembersAdd <- function(gid=NULL, gProgId=NULL, VProgIds, g=g){
	if(length(gProgId)>0) gid = which(names(g$groups)==as.character(gProgId))
	g$groups[[gid]]$MemberProgIds = union(g$groups[[gid]]$MemberProgIds, VProgIds)
	g
}

GroupMembersDelete <- function(gid=NULL, gProgId=NULL, VProgIds, g=g){
	if(length(gProgId)>0) gid = which(names(g$groups)==as.character(gProgId))
	g$groups[[gid]]$MemberProgIds = setdiff(g$groups[[gid]]$MemberProgIds, VProgIds)
	g
}

