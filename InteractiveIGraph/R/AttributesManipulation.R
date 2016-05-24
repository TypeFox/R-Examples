################ Atributu keitimas ##############

AttributeCreate <- function(name, otype=c("V", "E", "G"), g){
	otype = match.arg(otype)
	
	if(otype=="V"){
		g <- set.vertex.attribute(graph=g, name=name, value=NA)	
	}
	if(otype=="E"){
		g <- set.edge.attribute(graph=g, name=name, value=NA)	
	}
	
	g 
}

AttributeIsMeaningful <- function(name, otype=c("V", "E", "G"), g){
	otype = match.arg(otype)
	res = switch(EXPR=otype
		,V = name %in% list.vertex.attributes(g)
		,E = name %in% list.edge.attributes(g)
		,G = TRUE # is principo tai sarasa, taigi leidziam betkokias nesamones
	)
	res
}

AttributeSet <- function(name, value, oids=NULL, oProgIds=NULL, otype=c("V", "E", "G"), g, add=TRUE, quiet=FALSE){
	if(length(oProgIds)>0){
		oids = switch(EXPR=otype
			, none = vector()
			, V = as.numeric(V(g)[ProgId %in% oProgIds])
			, E = as.numeric(E(g)[ProgId %in% oProgIds])
			, G = which(names(g$groups) %in% oProgIds)	
		)
	}
	
	
	otype = match.arg(otype)
	if(!AttributeIsMeaningful(name=name, otype=otype, g=g)){
		if(add){
			g <- AttributeCreate(name=name, otype=otype, g=g)
		} else {
			if(!quiet){
				type = switch(EXPR=otype, V="vertex", E="edge", G="group")
				g <- MsgToLogObj(Msg="Attribute `"%.%name%.%"` has no mening for `"%.%type%.%"`. Use add=TRUE or create one by your self.", g=g, add=FALSE)
			}
			return(g)
		}	
	}
	
	
	if(otype=="V"){
		g <- set.vertex.attribute(graph=g, name=name, index=oids, value=value)	
	}
	if(otype=="E"){
		g <- set.edge.attribute(graph=g, name=name, index=oids, value=value)	
	}
	if(otype=="G"){
		if (length(g$groups) >= oids){
			if(any(g$groups[[oids]]$closed)){
				oids = as.numeric(V(g)[ProgId==g$groups[[oids]]$ProgId])[1]
				g <- Recall(name, value, oids,  otype=c("V"), g=g, add=FALSE, quiet=FALSE)
			} else {
				if(!all(name %in% c("ProgId", "ProgType"))){ #sito tai tikrai neleidziame keisti
					g$groups[[oids]][[name]] = value
				}
			}
		}
	}	

	g
}

AttributeAllSet <- function(name, value, g){
	### visu prima keiciam active (jei tik galim
	### jei ne tai keiciam visus paselectintus, t.y. verteces arba edge (tame tarpe ir bloko virsune)
	### jei ir to nera, bet yra view, tai keiciam jo parametrus, kitu atveju nieko nekeiciam

	### keiciam tik tada, jei atributas turi tam prasme


	# selected V or E (iskaitant bloko virsune)
	
	AO = GetActiveObject(g=g)
	if(AO$type!="none"){		
		g <- AttributeSet(name=name, value=value, oProgIds=AO$ProgId,  otype=AO$type, g=g, add=FALSE, quiet=FALSE)	
	} else { 
		if(length(V(g)[selected])>0 || length(E(g)[selected])>0){
			VMeaningful = AttributeIsMeaningful(name, otype="V", g=g)
			EMeaningful = AttributeIsMeaningful(name, otype="E", g=g)
			
			if(VMeaningful|EMeaningful){
				g <- AttributeSet(name=name, value=value, oids=as.numeric(V(g)[selected]),  otype="V", g=g, add=FALSE, quiet=TRUE)		
				g <- AttributeSet(name=name, value=value, oids=as.numeric(E(g)[selected]),  otype="E", g=g, add=FALSE, quiet=TRUE)		
			} else {
				g <- MsgToLogObj(Msg="Attribute `"%.%name%.%"` has no mening nor for vertex nor for edge. Use AttributeCreate to create one.", g=g, add=FALSE)
			}
		} else {
			VO = GetViewObject(g=g)
			if(VO$type!="none"){
				g <- AttributeSet(name=name, value=value, oProgIds=VO$ProgId,  otype=VO$type, g=g, add=FALSE, quiet=FALSE)	
			}
		}		
	}
	g
}

AttributeChange <- function(AttrCommand, g){

	#AttrCommand = "color='red'"

	eq = tryCatch({parse(text=AttrCommand)[[1]]}, 
		error = function(TechnicalMessage) {
			NULL
		}, warning = function(TechnicalMessage) {
			NULL
		})
	
	name = NULL
	value = NULL
	if(is.call(eq)){
		if(as.character(eq[[1]])=="=" && length(eq)==3){
			name = as.character(eq[[2]])
			value = eq[[3]]
		}
	}

	if(!is.null(name) && !is.null(value)){
		g <- AttributeAllSet(name=name, value=value, g=g)	
	} else {
		g <- MsgToLogObj(Msg="Spelling error... try next time.", g=g, add=FALSE)
	}
	
	g
}

