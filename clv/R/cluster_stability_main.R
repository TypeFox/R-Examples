
#-------------------------------------------------------------------------------#
#----------------- USER INTERFACES FOR CLUSTER STABILITY APPROACH --------------#

# cluster stability -> similarity index approach

cls.stab.sim.ind <- function( data, cl.num, 
				rep.num=10, 
				subset.ratio=0.75, 
				clust.method=c("agnes","pam"),
				method.type=c("single","average"),
				sim.ind.type=c("dot.pr","sim.ind"), 
				fast=TRUE, ... )
{
	# check input arguments
	data = data.validity(data, "data")
	cl.num = cls.num.vect.validity(cl.num, dim(data)[1], "cl.num")

  rep.num = check.rep.num(rep.num)
  subset.ratio = check.subset.ratio(data, subset.ratio)
  cl.num = cut.cl.num(data, cl.num, subset.ratio)

	cls.method.type.bool = check.avail.methods(clust.method, "clust.method", supp.cls.methods.vec.const)
	sim.ind.type.bool = check.avail.methods(sim.ind.type, "sim.ind.type", supp.cls.stab.sim.ind.vec.const )
	method.type.bool = check.avail.methods(method.type, "method.type", hierarhical.method.types.vec.const )
	
	if( !is.logical(fast) ) fast = TRUE

	iter = 1
	result.list = vector( "list", length=length( cls.method.type.bool[cls.method.type.bool] ) )

	for( method.num in 1:length(cls.method.type.bool) )
	{
		if( cls.method.type.bool[method.num] && supp.cls.methods.list.const[[method.num]]$sup )
		{
			if( method.num != agnes.num.const && method.num != hclust.num.const )
			{
				clust.alg.pver <- function(data, clust.num)
				{
					return( supp.cls.methods.list.const[[method.num]]$wrp
                  (
                    supp.cls.methods.list.const[[method.num]]$alg(data, clust.num, NULL, ...),
                    clust.num
                  )
                )
				}

				result.list[[iter]] = clust.stab.sim.ind.pver( 
					data=data, cl.num=cl.num, 
					sample.num=rep.num, ratio=subset.ratio,
					clust.method=clust.alg.pver,
					sim.ind=sim.ind.type.bool
				)
				names(result.list)[iter] = supp.cls.methods.vec.const[method.num]
				iter = iter + 1
			}
		
			else
			{
				# be careful - this line depends on the fact that 
				# first is computed "agnes" (if choosen) and always after "agnes", "hclust" (see constants variables)
				if( method.num == hclust.num.const ) method.type.bool = method.type.bool[1:4]
				
				if( any(method.type.bool) )
				{
					for( i in 1:length(method.type.bool))
					{
						if( method.type.bool[i] == TRUE )
						{
              clust.alg.hver <- function(data)
              {
                return( supp.cls.methods.list.const[[method.num]]$alg(data, 0, hierarhical.method.types.vec.const[i], ...) )
              }

              if( fast )
              {
							  result.list[[iter]] = clust.stab.sim.ind.hver( 
									  data=data, cl.num=cl.num, 
									  sample.num=rep.num, ratio=subset.ratio,
									  clust.method=clust.alg.hver,
									  clust.wrap=supp.cls.methods.list.const[[method.num]]$wrp,
									  sim.ind=sim.ind.type.bool
								  )
              }
              else
              {
                clust.alg.pver <- function(data, clust.num)
                {
                  return( supp.cls.methods.list.const[[method.num]]$wrp
                          (
                            clust.alg.hver(data),
                            clust.num
                          )
                        )
                }

                result.list[[iter]] = clust.stab.sim.ind.pver( 
                    data=data, cl.num=cl.num, 
                    sample.num=rep.num, ratio=subset.ratio,
                    clust.method=clust.alg.pver,
                    sim.ind=sim.ind.type.bool
                  )
              }
							names(result.list)[iter] = paste(
									supp.cls.methods.vec.const[method.num],
									hierarhical.method.types.vec.const[i],
									sep=".")
		
							iter = iter + 1
						}
					}
				}
			}
		}
	}

	return(result.list)
}

# cluster stability -> optimal assignment approach
 
cls.stab.opt.assign <- function( data, cl.num, 
				 rep.num=10, 
				 subset.ratio=0.75, 
				 clust.method=c("agnes","pam"),
				 method.type=c("single","average"),
				 fast=TRUE, ... )
{
	# check input arguments
	data = data.validity(data, "data")
	cl.num = cls.num.vect.validity(cl.num, dim(data)[1], "cl.num")
	
	rep.num = check.rep.num(rep.num)
	subset.ratio = check.subset.ratio(data, subset.ratio)
  cl.num = cut.cl.num(data, cl.num, subset.ratio) 

	cls.method.type.bool = check.avail.methods(clust.method, "clust.method", supp.cls.methods.vec.const )
	method.type.bool = check.avail.methods(method.type, "method.type", hierarhical.method.types.vec.const )

	if( !is.logical(fast) ) fast = TRUE

	iter = 1
	result.list = vector( "list", length=length( cls.method.type.bool[cls.method.type.bool] ) )

	for( method.num in 1:length(cls.method.type.bool) )
	{
		if( cls.method.type.bool[method.num] && supp.cls.methods.list.const[[method.num]]$sup )
		{
			if( method.num != agnes.num.const && method.num != hclust.num.const )
			{
				clust.alg.pver <- function(data, clust.num)
				{
					return( supp.cls.methods.list.const[[method.num]]$wrp
                  ( 
                     supp.cls.methods.list.const[[method.num]]$alg(data, clust.num, NULL, ...),
                     clust.num
                  )
                )
				}
				
				result.list[[iter]] = clust.stab.opt.assign.pver(
					data=data, cl.num=cl.num, 
					sample.num=rep.num, ratio=subset.ratio,
					clust.method=clust.alg.pver
					)
				names(result.list)[iter] = supp.cls.methods.vec.const[method.num]
				iter = iter + 1	
			}
			else
			{
				# be careful - this line depends on the fact that 
				# first is computed "agnes" (if choosen) and always after "agnes", "hclust" (see constants variables)
				if( method.num == hclust.num.const ) method.type.bool = method.type.bool[1:4]
	
				if( any(method.type.bool) )
				{
					for( i in 1:length(method.type.bool))
					{
						if( method.type.bool[i] == TRUE )
						{
              clust.alg.hver <- function(data)
              {
                return( supp.cls.methods.list.const[[method.num]]$alg(data, 0, hierarhical.method.types.vec.const[i], ...) )
              }

              if( fast )
              {
							  result.list[[iter]] = clust.stab.opt.assign.hver( 
									  data=data, cl.num=cl.num, 
									  sample.num=rep.num, ratio=subset.ratio,
									  clust.method=clust.alg.hver,
									  clust.wrap=supp.cls.methods.list.const[[method.num]]$wrp
								  )
              }
              else 
              {
                clust.alg.pver <- function(data, clust.num)
                {
                  return( supp.cls.methods.list.const[[method.num]]$wrp
                           ( 
                              clust.alg.hver(data),
                              clust.num
                           ) 
                        )
                }
 
                result.list[[iter]] = clust.stab.opt.assign.pver( 
                    data=data, cl.num=cl.num, 
                    sample.num=rep.num, ratio=subset.ratio,
                    clust.method=clust.alg.pver
                  )
              }
							names(result.list)[iter] = paste(
									supp.cls.methods.vec.const[method.num],
									hierarhical.method.types.vec.const[i],
									sep=".")
		
							iter = iter + 1
						}
					}
				}
			}
		}
	}

	return(result.list)
}

