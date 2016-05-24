#' Function to convert phylogenies from the class 'phylo' to the class 'ggphy'
#'
#' @param phylo an object of the class "phylo"
#' @param tip.dates a vector containing the sample dates of the tip in "Date" format, the dates must be ordered like the tips
#' @param branch.unit the unit of the branch. Either "year", "month", "day" or "subst". If a time unit is provided, together with tip.dates, then the x-axis of the phylogeny will be in the Date format
#' @param verbose if \code{TRUE} additional information is provided at execution
#' @author Anton Camacho
phylo2ggphy<-function(phylo,tip.dates=NULL,branch.unit=NULL,verbose=FALSE){

	phy<-phylo
	has.node.label<-(!is.null(phy$node.label))

	N.tips<-length(phy$tip.label)
    if(!is.null(tip.dates) & length(tip.dates)!=N.tips){stop("tip.dates must be the same length as the number of tips")}

	edge<-as.data.frame(phy$edge)
	names(edge)<-c("beg","end")

	#if edge.length is not provided, set all to 1 by default
	if(is.null(edge$length<-phy$edge.length)){edge$length<-rep(1,nrow(edge))}


#find root
	ind<-which(!edge$beg%in%edge$end)
	phy.root<-unique(edge$beg[ind])

	if(length(phy.root)!=1){
        cat(length(phy.root),"root(s) found!!\n")
        stop("Algorithm cannot handle more than one root in the phylo at the moment!")
	}

	if(!is.rooted(phy)){

		sel <- which(edge$beg == phy.root)[1]
		outgroup <- edge$end[sel]

		if (outgroup > N.tips){
                    outgroup <- prop.part(phy)[outgroup - N.tips]
                }

		phy <- root(phy, outgroup = unlist(outgroup), resolve.root = TRUE)

		edge<-as.data.frame(phy$edge)
		names(edge)<-c("beg","end")

		#if edge.length is not provided, set all to 1 by default
		if(is.null(edge$length<-phy$edge.length)){edge$length<-rep(1,nrow(edge))}


		if(has.node.label){
#add NA to the root label
            ind<-phy.root-N.tips
            tmp<-rep(NA,phy$Nnode)
            tmp[-ind]<-phy$node.label
            phy$node.label<-tmp
		}
	}

	Nt<-N.tips
	Nn<-phy$Nnode
	Ne<-nrow(phy$edge)
	Nnt<-Nn+Nt

	tip<-1:Nt
	node<-(Nt+1):Nnt

#get the y coord of nodes: start from tips and browse tree to the root
	y<-rep(NA,Nnt)
	order.tip<-edge$end[edge$end<=Nt]
	y[order.tip]<-tip
	cur.tip<-tip
	removed<-c()

	while(length(removed)<(Nnt-1)){

#how many nodes have two children among cur.tip
		ind<-which(edge$end%in%cur.tip)
		n.children<-table(edge$beg[ind])
		ind<-which(n.children>1)

		new.tip<-as.numeric(names(ind))
		old.tip<-edge$end[edge$beg%in%new.tip]

		y.new.tip<-sapply(new.tip,function(i){
			children<-edge$end[edge$beg==i]
			return(mean(y[children]))
		})

		y[new.tip]<-y.new.tip
		removed<-c(removed,old.tip)
		cur.tip<-cur.tip[!cur.tip%in%old.tip]
		cur.tip<-c(cur.tip,new.tip)

	}

#get x coord: start from root and browse tree until tips
	x<-rep(NA,Nnt)
	x[phy.root]<-0
	cur.node<-phy.root
	visited<-c(cur.node)

	while(length(visited)<Nnt){
		df1<-data.frame(beg=cur.node)
		df2<-merge(df1,edge)
		x[df2$end]<-x[df2$beg]+df2$length
		cur.node<-df2$end
		visited<-c(visited,cur.node)
	}


    if(!is.null(tip.dates) & !is.null(branch.unit))
    	if(branch.unit%in%c("year","month","day")){
#tip.date<-as.Date(extract.string(phy$tip.label,"_",2))
        df<-data.frame(tip,date=tip.dates,age=x[tip])
        time.unit=switch(branch.unit,year=365.25,month=30.5,day=1)
        root.date<-mean(df$date-df$age*time.unit)
#check that the new tip dates do not vary much from sample dates
		if(verbose){
	        cat("X axis is converted into date\n")
			new.tip.date<-root.date+df$age*time.unit
			print(paste("reconstruction of nodes ages in Date format led to change of the tip age ranging:",range(new.tip.date-df$date)))
		}

#x coord as date
        x<-root.date+x*time.unit
    }


#build data frame for nodes and tips
    if(has.node.label){
        df.node<-data.frame(node,x=x[node],y=y[node],label=phy$node.label)
    }else{
        df.node<-data.frame(node,x=x[node],y=y[node])
    }

    df.tip<-data.frame(tip,x=x[tip],y=y[tip],label=phy$tip.label)

#build segment data frame
#horizontal segments
    edge.h<-edge
    edge.h$direction<-"H"
    edge.h$x.beg<-x[edge.h$beg]
    edge.h$x.end<-x[edge.h$end]
    edge.h$y.beg<-y[edge.h$end]
    edge.h$y.end<-edge.h$y.beg

#vertical segments
    edge.v<-edge
    edge.v$direction<-"V"
    edge.v$x.beg<-x[edge.v$beg]
    edge.v$x.end<-edge.v$x.beg
    edge.v$y.beg<-y[edge.v$beg]
    edge.v$y.end<-y[edge.v$end]

#combine
    df.edge<-rbind(edge.h,edge.v)
    df.edge<-df.edge[,-which(names(df.edge)=="length")]

#if has.node.label==T put label on edges
    if(has.node.label){
        tmp<-data.frame(end=node,label=phy$node.label)
        df.edge<-merge(df.edge,tmp,all.x=TRUE)
    }

    return(list(df.tip,df.node,df.edge))
}
