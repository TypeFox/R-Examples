


menuFunctions<-(function(){
			menus<-list()
			
			getMenus<-function(){
				return(menus)
			}
			
			setMenus<-function(newMenus){
				menus<<-newMenus
			}
			
			addMenu<-function(name, pos=length(menus)+1){
				newMenu<-list(name=name,items = list() )
				position<-pos
				if(length(menus)==0)
					menus <<-list(newMenu)
				else if(pos==1)
					menus <<- c(list(newMenu),menus)
				else if(pos<=length(menus))
					menus <<- c(menus[1:(pos-1)],list(newMenu),menus[pos:length(menus)])
				else if(pos>length(menus))
					menus <<- c(menus,list(newMenu))
				names(menus)[pos] <<- name 
			}
			
			addMenuItem<-function(name, pos=NULL, command, menuName, silent=TRUE){
				m.ind<-which(sapply(menus,function(x) x$name==menuName))
				m<-menus[[m.ind]]$items
				if(is.null(pos)) pos <-length(m)+1
				newItem<-list(name=name,command=command,silent=silent)
				if(length(m)==0){
					menus[[m.ind]]$items<<-list(newItem)
					names(menus[[m.ind]]$items)<<-name
					return(NULL)
				}
				if(pos==1)
					m <- c(list(newItem),m)
				else if(pos<=length(m))
					m <- c(m[1:(pos-1)],list(newItem),m[pos:length(m)])
				else if(pos>length(m))
					m <- c(m,list(newItem))
				menus[[m.ind]]$items <<- m
				names(menus[[m.ind]]$items)[pos] <<- name 
			}
			return(list(deducer.getmenus=getMenus,
							deducer.setmenus=setMenus,
							deducer.addmenu=addMenu,
							deducer.addMenuItem=addMenuItem))
})()

deducer.getMenus <- menuFunctions[["deducer.getmenus"]]
deducer.setMenus <- menuFunctions[["deducer.setmenus"]]
deducer.addMenu <- menuFunctions[["deducer.addmenu"]]
deducer.addMenuItem <- menuFunctions[["deducer.addMenuItem"]]



































