library(clv)

# ----------- PLOT TO PNG function ---------------------------------------------

plotToPng <- function(data, file.name="plot.png", ... )
{
	png( file=file.name, width=1000, height=1000 )
	plot( data , ... )
	dev.off()
}

pairsToPng <- function(x, file.name="pairsplot.png", ... )
{
	png( file=file.name, width=1000, height=1000 )
	pairs( x , ... )
	dev.off()
}

matplotToPng <- function(x, y=NULL, file.name="matplot.png", ... )
{
	png( file=file.name, width=1000, height=1000 )
	matplot( x , y , ... )
	dev.off()
}

boxplotToPng <- function(x, file.name="matplot.png", ... )
{
	png( file=file.name, width=1000, height=1000 )
	boxplot( x , ... )
	dev.off()
}

# ---------- WRAPERS: AGNES, PAM --------------------------------------

pam.wrapp <-function(data)
{
	return( data$clustering )
}

identity <- function(data) { return(data) }

agnes.single <- function(data, clust.num)
{
	return( cutree( agnes(data,method="single"), clust.num ) )
}

agnes.complete <- function(data, clust.num)
{
	return( cutree( agnes(data,method="complete"), clust.num ) )
}

agnes.average <- function(data, clust.num)
{
	return( cutree( agnes(data,method="average"), clust.num ) )
}

agnes.ward <-function(data, clust.num)
{
	return( cutree( agnes(data,method="ward"), clust.num ) )
}

agnes.weighted <-function(data, clust.num)
{
	return( cutree( agnes(data,method="weighted"), clust.num ) )
}

diana.norm <-function(data, clust.num)
{
	return( cutree( diana(data), clust.num ) )
}
# ------------------- STD INTERNAL - AGNES ------------------------

agnes.std.internal <- function( data, iter, intra, inter, method.name )
{
	dim.mult = length(intra)*length(inter)	
	
	dunn = as.list(iter)
	db = as.list(iter)

	mx.dunn.plot = matrix(0, length(iter), dim.mult)
	mx.db.plot = matrix(0, length(iter), dim.mult)
	vect.conn = rep(0, times=length(iter))

	k = 1
	agn.result = agnes(data,method=method.name)
	
	# iter - vector with numbers of cluster to which data will be partitioned
	for( i in iter )
	{
		clust <- cutree( agn.result, i )
		tmp = cls.scatt.data(data,clust)
		
		# gather particular indecies for each iteration (cluster number)
		# and put it in the matrix 
		dunn[[k]] = clv.Dunn(tmp, intra, inter)
		dunn.v = as.vector(dunn[[k]])
		for( t in 1:dim.mult) mx.dunn.plot[k,t] = dunn.v[t]
		
		db[[k]] = clv.Davies.Bouldin(tmp, intra, inter)
		db.v = as.vector(db[[k]])
		for( t in 1:dim.mult) mx.db.plot[k,t] = db.v[t]

		vect.conn[k] = connectivity(data, clust, 10)
		
		k = k+1
	}
	return( list(Dunn.matrix=mx.dunn.plot, DB.matrix=mx.db.plot, Conn.vect=vect.conn) )
}

# ------------------- STD INTERNAL - PAM ------------------------

pam.std.internal <- function( data, iter, intra, inter )
{
	dim.mult = length(intra)*length(inter)	
	
	dunn = as.list(iter)
	db = as.list(iter)

	mx.dunn.plot = matrix(0, length(iter), dim.mult)
	mx.db.plot = matrix(0, length(iter), dim.mult)

	k = 1
	# iter - vector with numbers of cluster to which data will be partitioned
	for( i in iter )
	{
		clust <- pam(data, i)$clustering
		tmp = cls.scatt.data(data,clust)
		
		# gather particular indecies for each iteration (cluster number)
		# and put it in the matrix 		
		dunn[[k]] = clv.Dunn(tmp, intra, inter)
		dunn.v = as.vector(dunn[[k]])
		for( t in 1:dim.mult) mx.dunn.plot[k,t] = dunn.v[t]
		
		db[[k]] = clv.Davies.Bouldin(tmp, intra, inter)
		db.v = as.vector(db[[k]])
		for( t in 1:dim.mult) mx.db.plot[k,t] = db.v[t]
		
		k = k+1
	}
	return( list(Dunn.matrix=mx.dunn.plot, DB.matrix=mx.db.plot) )
}

# ----------------- Cluster Stablility -------------------------------------

cls.stabb <- function( data, clust.num, sample.num , ratio, method, wrapp  )
{
	rand = 0
	jacc = 0
	flkm = 0
	dot.pr  = 0
	sim.ind = 0

	dot.pr.vec = rep(0, times=sample.num)
	sim.ind.vec = rep(0, times=sample.num)

	obj.num = dim(data)[1]

	for( j in 1:sample.num )
	{
		smp1 = sort( sample( 1:obj.num, ratio*obj.num ) )
		smp2 = sort( sample( 1:obj.num, ratio*obj.num ) )

		d1 = data[smp1,]
		cls1 = wrapp( method(d1,clust.num) )

		d2 = data[smp2,]
		cls2 = wrapp( method(d2,clust.num) )

		clsm1 = t(rbind(smp1,cls1))
		clsm2 = t(rbind(smp2,cls2))

		m = cls.set.section(clsm1, clsm2)
		cls1 = as.integer(m[,2])
		cls2 = as.integer(m[,3])
		cnf.mx = confusion.matrix(cls1,cls2)
		std.ms = std.ext(cls1,cls2)
		
		# external measures - compare partitioning
		rd = clv.Rand(std.ms)
		dt = dot.product(cls1,cls2)
		si = similarity.index(cnf.mx)

		rand = rand + rd/sample.num

		if( !is.nan(dt) ) dot.pr = dot.pr + dt/sample.num 
		dot.pr.vec[j] = dt

		sim.ind = sim.ind + si/sample.num
		sim.ind.vec[j] = si
	}

	print(paste("Sample number = ", j, " for clust num = ", clust.num))	
	sim.ind.vec = ((clust.num+1)/clust.num)*sim.ind.vec - 1/clust.num
	sim.ind = ((clust.num+1)/clust.num)*sim.ind - 1/clust.num
	print("--------------")

	return( list( ind.means = c(rand, dot.pr, sim.ind), dt=dot.pr.vec, si=sim.ind.vec ) )
}


# ------------------ stability for Agnes --------------------------------------- 

agnes.stability <- function( data, iter, sample.num, ratio, method.name )
{
	ind.num = 3 # 1 x standard, cor & sim ind
	res.mean.mx = matrix(0,0,ind.num)
	res.dt.mx = matrix(0,sample.num,0)
	res.si.mx = matrix(0,sample.num,0)
	for( i in iter )
	{
		res.ls = cls.stabb(data, clust.num=i, sample.num=sample.num, 
						   ratio=ratio, method=method.name, wrapp=identity ) 
		res.mean.mx = rbind( res.mean.mx, res.ls$ind.means )

		res.dt.mx = cbind( res.dt.mx, res.ls$dt )
		res.si.mx = cbind( res.si.mx, res.ls$si )
	}

	names = paste("C", iter, sep="")
	r1 = as.data.frame(res.mean.mx)
	rownames(r1) = names
	r2 = as.data.frame(res.dt.mx)
	colnames(r2) = names
	r3 = as.data.frame(res.si.mx)
	colnames(r3) = names
	
	return(
		list(
			mean.mx=r1,
			dt.mx = r2,
			si.mx = r3
			)
		)
}

# ------------------ stability for PAM -------------------------------------------

pam.stability <- function( data, iter, sample.num, ratio )
{
	ind.num = 3 # 3 x standard, cor & sim ind
	res.mean.mx = matrix(0,0,ind.num)
	res.dt.mx = matrix(0,sample.num,0)
	res.si.mx = matrix(0,sample.num,0)
	for( i in iter )
	{
		res.ls = cls.stabb(data, clust.num=i, sample.num=sample.num, ratio=ratio, method=pam, wrapp=pam.wrapp)
		res.mean.mx = rbind( res.mean.mx, res.ls$ind.means )

		res.dt.mx = cbind( res.dt.mx, res.ls$dt )
		res.si.mx = cbind( res.si.mx, res.ls$si )
	}

	names = paste("C", iter, sep="")
	r1 = as.data.frame(res.mean.mx)
	rownames(r1) = names
	r2 = as.data.frame(res.dt.mx)
	colnames(r2) = names
	r3 = as.data.frame(res.si.mx)
	colnames(r3) = names
	
	return(
		list(
			mean.mx=r1,
			dt.mx = r2,
			si.mx = r3
			)
		)
}

# ----------------- Additional data functions: eg. circles ---------------------------------

circles <- function(obj.num, circ.radius)
{
	norm.vect = rnorm(obj.num,sd=1)
	norm.rad = rnorm(obj.num,sd=circ.radius/50)

	vect1 = c()
	vect2 = c()
	i = 1

	if( length(obj.num) != length(circ.radius) ) error("Bad vector length!")

	for( rad in circ.radius )
	{
		norm.vect = rnorm(obj.num[i],sd=1)
		norm.rad = rnorm(obj.num[i],sd=circ.radius/20)

		vect1 = c( vect1, ((rad + norm.rad) * cos( pi * norm.vect )) )
		vect2 = c( vect2, ((rad + norm.rad) * sin( pi * norm.vect )) )
		i = i+1
	}

	r = rbind(vect1, vect2)
	return(t(r))
}


# ------------------------- MAIN ----------------------------------------------------

main <- function()
{
	data(iris)
	data.new = iris[,1:4]
	#data.new = rFace(1000)
	
	pl_num = 1
	plotToPng(data.new, file.name=paste(pl_num,"_dataplot.png")) #, xlab="V_1", ylab="V_2")
	pl_num = 1 + pl_num

	#data.new = data.new[data.new[,5] < 11 & data.new[,5] > -11 & data.new[,2] < 23 & data.new[,6] < 3,]
	data = as.data.frame(data.new)
	n <- dim(data)[1]
	m <- dim(data)[2]

	mean.vec <- mean(data)
	sd.vec <- sd(data)
	
	mean.mx <- matrix(rep(mean.vec,times=n),n,m,byrow=T)
	sd.mx <- matrix(rep(sd.vec,times=n),n,m,byrow=T)

	data <- (data - mean.mx)/sd.mx

	# subset size 
	rt = 0.75
	# number of samples in each iteration 
	smp = 5
	# vector with interesting numbers of clusters
	iter = c(2,3,4,5,6,8,10,12)
	
	
	intraclust = c("complete") #,"average","centroid")
	interclust = c("single", "average") #, "complete", "centroid", "aveToCent", "hausdorff")

	plotToPng(data, file.name=paste(pl_num,"_dataplot.png"))
	pl_num = 1 + pl_num

	ls = list(agn.sin=T, agn.ave=T, agn.com=T, agn.wrd=T, agn.wgt=T, diana=T, pam=T)

	if( ls$agn.sin == TRUE )
	{
		# results and plots for Agnes "single linkage"
		method = "AGNES"
		type = "single"
		print(paste(method," and cluster stablitity:"))

		agn.stab = agnes.stability( data, iter, sample.num=smp, ratio=rt, method.name=agnes.single )
		matplotToPng(x=iter, y=agn.stab$mean.mx, file.name=paste(pl_num,"_matplot.png",sep=""), 
			     type="b", main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="")
		pl_num = 1 + pl_num
		print(agn.stab)
		boxplotToPng(agn.stab$dt, file.name=paste(pl_num, "_boxplot.png",sep=""), 
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Correlation similarity index")
		pl_num = 1 + pl_num
		boxplotToPng(agn.stab$si, file.name=paste(pl_num, "_boxplot.png",sep=""),
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Similarity index")
		pl_num = 1 + pl_num
	}
	
	if( ls$agn.ave == TRUE )
	{
		# results and plots for Agnes "complete linkage"
		method = "AGNES"
		type = "average"
		print(paste(method," and cluster stablitity:"))

		agn.stab = agnes.stability( data, iter, sample.num=smp, ratio=rt, method.name=agnes.average )
		matplotToPng(x=iter, y=agn.stab$mean.mx, file.name=paste(pl_num,"_matplot.png",sep=""), 
			     type="b", main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="")
		pl_num = 1 + pl_num
		print(agn.stab)
		boxplotToPng(agn.stab$dt, file.name=paste(pl_num, "_boxplot.png",sep=""), 
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Correlation similarity index")
		pl_num = 1 + pl_num
		boxplotToPng(agn.stab$si, file.name=paste(pl_num, "_boxplot.png",sep=""),
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Similarity index")
		pl_num = 1 + pl_num
	}
	if( ls$agn.com == TRUE )
	{
		# results and plots for Agnes "average linkage" 
		method = "AGNES"
		type = "complete"

		print(paste(method," and cluster stablitity:"))
		print(type)

		agn.stab = agnes.stability( data, iter, sample.num=smp, ratio=rt, method.name=agnes.complete )
		matplotToPng(x=iter, y=agn.stab$mean.mx, file.name=paste(pl_num,"_matplot.png",sep=""), 
			     type="b", main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="")
		pl_num = 1 + pl_num
		print(agn.stab)
		boxplotToPng(agn.stab$dt, file.name=paste(pl_num, "_boxplot.png",sep=""), 
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Correlation similarity index")
		pl_num = 1 + pl_num
		boxplotToPng(agn.stab$si, file.name=paste(pl_num, "_boxplot.png",sep=""),
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Similarity index")
		pl_num = 1 + pl_num
	}

	if( ls$agn.wrd == TRUE )
	{
		method = "AGNES"
		type = "ward"

		print(paste(method," and cluster stablitity:"))
		print(type)

		agn.stab = agnes.stability( data, iter, sample.num=smp, ratio=rt, method.name=agnes.ward )
		matplotToPng(x=iter, y=agn.stab$mean.mx, file.name=paste(pl_num,"_matplot.png",sep=""), 
			     type="b", main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="")
		pl_num = 1 + pl_num
		print(agn.stab)
		boxplotToPng(agn.stab$dt, file.name=paste(pl_num, "_boxplot.png",sep=""), 
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Correlation similarity index")
		pl_num = 1 + pl_num
		boxplotToPng(agn.stab$si, file.name=paste(pl_num, "_boxplot.png",sep=""),
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Similarity index")
		pl_num = 1 + pl_num
	}

	if( ls$agn.wgt == TRUE )
	{
		method = "AGNES"
		type = "weighted"

		print(paste(method," and cluster stablitity:"))
		print(type)

		agn.stab = agnes.stability( data, iter, sample.num=smp, ratio=rt, method.name=agnes.weighted )
		matplotToPng(x=iter, y=agn.stab$mean.mx, file.name=paste(pl_num,"_matplot.png",sep=""), 
			     type="b", main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="")
		pl_num = 1 + pl_num
		print(agn.stab)
		boxplotToPng(agn.stab$dt, file.name=paste(pl_num, "_boxplot.png",sep=""), 
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Correlation similarity index")
		pl_num = 1 + pl_num
		boxplotToPng(agn.stab$si, file.name=paste(pl_num, "_boxplot.png",sep=""),
				main=paste(method, " (", type, ")", sep=""), xlab="Number of clusters", ylab="Similarity index")
		pl_num = 1 + pl_num
	}

	if( ls$diana == TRUE )
	{
		method = "DIANA"
		print(paste(method," and cluster stablitity:"))

		agn.stab = agnes.stability( data, iter, sample.num=smp, ratio=rt, method.name=diana.norm )
		matplotToPng(x=iter, y=agn.stab$mean.mx, file.name=paste(pl_num,"_matplot.png",sep=""), 
			     type="b", main=paste(method, sep=""), xlab="Number of clusters", ylab="")
		pl_num = 1 + pl_num
		print(agn.stab)
		boxplotToPng(agn.stab$dt, file.name=paste(pl_num, "_boxplot.png",sep=""), 
				main=paste(method, sep=""), xlab="Number of clusters", ylab="Correlation similarity index")
		pl_num = 1 + pl_num
		boxplotToPng(agn.stab$si, file.name=paste(pl_num, "_boxplot.png",sep=""),
				main=paste(method, sep=""), xlab="Number of clusters", ylab="Similarity index")
		pl_num = 1 + pl_num
	}

	if( ls$pam == TRUE )
	{
		# plots and prints for PAM alghorithm
		method = "PAM" 
	
		print("PAM and cluster stablitity:")
		
		pam.stab = pam.stability( data, iter, sample.num=smp, ratio=rt )
		matplotToPng(x=iter, y=pam.stab$mean.mx, file.name=paste(pl_num,"_matplot.png",sep=""), 
				     type="b", main=paste(method, sep=""), xlab="Number of clusters", ylab="")
		pl_num = 1 + pl_num
		print(pam.stab)
		boxplotToPng(pam.stab$dt, file.name=paste(pl_num, "_boxplot.png",sep=""), 
					main=paste(method, sep=""), xlab="Number of clusters", ylab="Correlation similarity index")
		pl_num = 1 + pl_num
		boxplotToPng(pam.stab$si, file.name=paste(pl_num, "_boxplot.png",sep=""), 
					main=paste(method, sep=""), xlab="Number of clusters", ylab="Similarity index")
		pl_num = 1 + pl_num
	}

	print(paste(pl_num, "of \"*.png\" plots generated in folder:",getwd()))
}

# exec main function
main()