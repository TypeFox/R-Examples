
#-------------------------------------------------------------------------------------------------#

# external measures constant values 

dt.sim.ind.const = 1
optas.sim.ind.const = 2
rand.const = 3
jaccard.const = 4

supp.cls.stab.sim.ind.vec.const = c("dot.pr", "sim.ind", "rand", "jaccard")

cls.stab.names.const = list(
	dot.product = "Cluster stability - dot product",
	opt.assgn = "Cluster stability - optimal assignment",
	rand.ind = "Cluster stability - Rand index",
	rand.ind = "Cluster stability - Jaccard index"
)

# cluster methods constant values

agnes.num.const	 = 1
diana.num.const  = 2
#mona.num.const   = 3
hclust.num.const = 3
kmeans.num.const = 4
pam.num.const    = 5
clara.num.const  = 6
#model.num.const  = 8
#sota.num.const   = 9

supp.cls.methods.list.const = list(
	list(
		num = agnes.num.const,
		alg = agnes.clust,
		wrp = agnes.wrap,
		hrr = TRUE,
		sup = TRUE
	),
	list(
		num = diana.num.const,
		alg = diana.clust,
		wrp = diana.wrap,
		hrr = TRUE,
		sup = TRUE
	),
#	list(
#		num = mona.num.const,
#		alg = mona.clust,
#		wrp = mona.wrap,
#		hrr = TRUE,
#		sup = FALSE
#	), 
	list(
		num = hclust.num.const,
		alg = hclust.clust,
		wrp = hclust.wrap,
		hrr = TRUE,
		sup = TRUE
	),
	list(
		num = kmeans.num.const,
		alg = kmeans.clust,
		wrp = kmeans.wrap,
		hrr = FALSE,
		sup = TRUE
	),
	list(
		num = pam.num.const,
		alg = pam.clust,
		wrp = pam.wrap,
		hrr = FALSE,
		sup = TRUE
	),
	list(
		num = clara.num.const,
		alg = clara.clust,
		wrp = clara.wrap,
		hrr = FALSE,
		sup = TRUE
#	), 
#	list(
#		num = model.num.const,
#		alg = FALSE, #model.clust,
#		wrp = FALSE,  # model.wrap,
#		hrr = FALSE,
#		sup = FALSE
#	), 
#	list(
#		num = sota.num.const,
#		alg = FALSE, #sota.clust,
#		wrp = FALSE,  #sota.wrap
#		hrr = FALSE,
#		sup = FALSE
	)
)

#supp.cls.methods.alg.const = list( agnes.clust, diana.clust, mona.clust, hclust.clust, 
#			        kmeans.num.const, pam.num.const, clara.num.const, model.num.const, sota.num.const )
supp.cls.methods.vec.const = c("agnes", "diana", "hclust", "kmeans", "pam", "clara")

# agnes/hclust cluster method constant values

single.num.const   = 1
average.num.const  = 2
complete.num.const = 3
ward.num.const     = 4
weighted.num.const = 5

hierarhical.method.types.vec.const = c("single", "average", "complete", "ward", "weighted")

#-------------------------------------------------------------------------------------------------#

# classification/prediction algorithms constants

knn.num.const = 1
pred.method.types.vec.const = c("knn")

supp.pred.methods.list.const = list(
	list(
		num = knn.num.const,
		alg = knn.pred,
		sup = TRUE
	)
)

#-------------------------------------------------------------------------------------------------#
