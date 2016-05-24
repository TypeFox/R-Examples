#--------------------------------------------------------------------------------------------------#
#------------------- CLUSTER STABILITY -> DISTANCE BETWEEN CLUSTERING RESULTS ---------------------#

# some check functions

check.pred.wrapper <- function(pred.method)
{
  return( list(pred.method=pred.method) )
}

#---------------------------------#
# ALGORITHM
#---------------------------------#

get.perm.sample <- function(prev.perm, clust.num)
{
  head = 1:2
  
  if(is.vector(prev.perm)) 
  {
    if(prev.perm[1] == clust.num)
    {
      if(prev.perm[2] >= (clust.num-1))
      { 
        head = NULL
        tail = NULL
      }
      else
      {
        head[1] = clust.num
        head[2] = prev.perm[2] + 1
      }
    }
    else
    {
      if(prev.perm[2] >= clust.num)
      {
        head[1] = prev.perm[1] + 1
        head[2] = 1
      }
      else
      {
        head[1] = prev.perm[1]
        if( (prev.perm[2] + 1) == prev.perm[1] ) head[2] = prev.perm[2] + 2
        else head[2] = prev.perm[2] + 1
      }
    }
  }

  if(clust.num > 2) tail = sample(1:clust.num)[-head]
  else tail = NULL

  return(c(head,tail))
}

stab.assoc.factor <- function(cnf.mx, clust.num)
{
  factor = 0
  perm.num = clust.num *(clust.num-1)

  iter = 0
  perm.smp.vec = NULL
  smp.tmp = 1:clust.num

  while( iter < perm.num )
  {
    # subsampling algorithm,- to be investigated which approach to use
    # perm.smp.vec = get.perm.sample(perm.smp.vec, clust.num)
    perm.smp.vec = sample(smp.tmp)
    factor = factor + ( similarity.index.int(cnf.mx, perm.smp.vec) / perm.num )
    iter = iter + 1
  }

  return(factor)
}

# -------------------------------------------------------------------------------------------
# CLUSTER STABILITY ALGORITHM - transfer by prediction
# version for all - hierarchical and aglomerative algorithms (for hierarcical is not optimal)

clust.stab.predict.pver.internal <- function( data, clust.num, sample.num, ratio, clust.method, pred.method )
{
# prepare results
	result = 0
	si.res = 0
	factor.res = 0

	obj.num = dim(data)[1]
  
	for( j in 1:sample.num )
	{
		smp = sort( sample( 1:obj.num, ratio * obj.num ) )

		base.data = as.matrix(data[smp,])
		cls.base = clust.method( base.data, clust.num )

		rest.data = as.matrix(data[-smp,])
		cls.rest = clust.method(rest.data, clust.num)

		cls.pred = pred.method( base.data, cls.base, rest.data )

    print("CLUST STAB")
    print(cls.rest)
    print(cls.pred)
    # JEST PROBLEM jesli cnf.mx nie jest kwardratowa!!!!!!!!
		cnf.mx = confusion.matrix(cls.rest, cls.pred)
    # claculate fi_s  
		si = similarity.index( cnf.mx )
    # claculate R_s
		factor = stab.assoc.factor(cnf.mx, clust.num)

    # calculate S(fi_s)
		si.res = ( si / sample.num ) + si.res
    # calculate S(R_s)  
		factor.res = ( factor / sample.num ) + factor.res
	}

  # normalize the result
  tmp = 1/(1-factor.res)
	result = ( tmp * si.res ) + 1 - tmp
	return( result )
}

# ------- cluster stability - similarity index -> ver. for partitionings algorithms ------ #

clust.stab.predict.pver <- function( data, cl.num, sample.num, ratio, clust.method, pred.method )
{
	# result = vector("list", length=ind.num )
	result = as.data.frame(matrix(0, 1, length(cl.num)))

	iter = 1
	for( cls.num in cl.num )
	{
		result[[iter]] = clust.stab.predict.pver.internal(data, clust.num=cls.num, sample.num=sample.num, ratio=(1-ratio), clust.method=clust.method, pred.method=pred.method )
		iter = iter + 1
	}
  
	colnames(result) = paste("C", cl.num, sep="")

	return(result)
}

# ------- cluster stability - similarity index -> ver. for hierarhical algorithms ------ #

clust.stab.predict.hver <- function(data, cl.num, sample.num, ratio, clust.method, clust.wrap, pred.method)
{
  cl.num.len = length(cl.num)
	obj.num = dim(data)[1]
	
  si.vec  = rep(0, times=cl.num.len) 
  fac.vec = rep(0, times=cl.num.len) 
  
	for( j in 1:sample.num )
	{
		smp = sort( sample( 1:obj.num, (1-ratio) * obj.num ) )
		base.data = data[smp,]

		clust.tree = clust.method(base.data)

		rest.data = data[-smp,]
		rest.tree = clust.method(rest.data)

		iter = 1
		for( clust.num in cl.num )
		{
			base.cls = clust.wrap( clust.tree, clust.num )
			rest.cls = clust.wrap( rest.tree, clust.num )
			rest.pred = pred.method( base.data, base.cls, rest.data )

			cnf.mx = confusion.matrix(rest.cls, rest.pred)
      # claculate fi_s
			si = similarity.index( cnf.mx )
      # calculate R_s
			factor = stab.assoc.factor(cnf.mx, clust.num)

      # calculate S(fi_s)
			si.vec[iter]  = si.vec[iter]  + si/sample.num
      # calculate S(R_s)
      fac.vec[iter] = fac.vec[iter] + factor/sample.num
			iter = iter + 1 
		}
	}

  result  = as.data.frame(matrix(0,1,length(cl.num)))

  # normalize all results
  for(i in 1:cl.num.len)
  {
    tmp = 1/(1-fac.vec[i])
    result[1, i] = tmp * si.vec[i] + 1 - tmp
  }

	colnames(result) = paste("C", cl.num, sep="")
	return(result)
}

# cluster stability -> similarity index approach

cls.stab.pred <- function( data, cl.num, 
								rep.num=10, 
								subset.ratio=0.75, 
								clust.method=c("agnes","pam"),
								method.type=c("single","average"),
								pred.type=c("knn"),
								fast=TRUE, ... )
{
# check input arguments
	data = data.validity(data, "data")
	cl.num = cls.num.vect.validity(cl.num, dim(data)[1], "cl.num")

	rep.num = check.rep.num(rep.num)
  subset.ratio = check.subset.ratio(data, subset.ratio)
  cl.num = cut.cl.num(data, cl.num, subset.ratio)

	cls.method.type.bool = check.avail.methods(clust.method, "clust.method", supp.cls.methods.vec.const)
	method.type.bool = check.avail.methods(method.type, "method.type", hierarhical.method.types.vec.const )
	pred.method.type.bool = check.avail.methods(pred.type, "pred.type", pred.method.types.vec.const )

	pred.alg.num = length( pred.method.type.bool[pred.method.type.bool] )
	result.list = vector( "list", length=pred.alg.num )

	for( pred.method.num in 1:length(pred.method.type.bool) )
	{
		if( pred.method.type.bool[pred.method.num] )
		{
			iter = 1
			cls.result.list = vector( "list", length=length( cls.method.type.bool[cls.method.type.bool] ) )

			for( method.num in 1:length(cls.method.type.bool) )
			{
				if( cls.method.type.bool[method.num] && supp.cls.methods.list.const[[method.num]]$sup )
				{
# if algorithm is hierarhical and user wants fast computation set proper function
					if( supp.cls.methods.list.const[[method.num]]$hrr && fast ) algorithm = clust.stab.predict.hver
					else algorithm = clust.stab.predict.pver

					if( method.num != agnes.num.const && method.num != hclust.num.const )
					{
						clust.alg.pver <- function(data, clust.num) 
						{ 
              return( supp.cls.methods.list.const[[method.num]]$wrp(
                        supp.cls.methods.list.const[[method.num]]$alg(
                          data, clust.num=clust.num, 
                          method.type=NULL, ...
                        )
                      )
                    )
						}
						
            cls.result.list[[iter]] = clust.stab.predict.pver( 
														data=data, cl.num=cl.num, 
														sample.num=rep.num, ratio=subset.ratio,
														clust.method=clust.alg.pver,
														pred.method=supp.pred.methods.list.const[[pred.method.num]]$wrp
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
										return( supp.cls.methods.list.const[[method.num]]$alg(data, 0, method.type=hierarhical.method.types.vec.const[i], ...) ) 
									}
                  
                  if( fast )
                  {
									  cls.result.list[[iter]] = clust.stab.predict.hver( 
																	  data=data, cl.num=cl.num, 
																	  sample.num=rep.num, ratio=subset.ratio,
																	  clust.method=clust.alg.hver,
																	  clust.wrap=supp.cls.methods.list.const[[method.num]]$wrp, 
				                            pred.method=supp.pred.methods.list.const[[pred.method.num]]$wrp
																  )
                  }
                  else
                  {
                    clust.alg.pver <- function(data, clust.num) 
                    { 
                      return( supp.cls.methods.list.const[[method.num]]$wrp(clust.alg.hver(data), clust.num ) )
                    }
                    
                    cls.result.list[[iter]] = clust.stab.predict.pver( 
                                    data=data, cl.num=cl.num, 
                                    sample.num=rep.num, ratio=subset.ratio,
                                    clust.method=clust.alg.hver,
                                    pred.method=supp.pred.methods.list.const[[pred.method.num]]$wrp
                                  )
                  }
									names(cls.result.list)[iter] = paste(
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
			
			result.list[[i]] = cls.result.list
		}
	}
	
	return(result.list)
}

cls.stab.pred.usr <- function( 
                data, cl.num,
                clust.alg,
                pred.alg, 
                rep.num=10, 
                subset.ratio=0.75
                )
{
# check input arguments
  data = data.validity(data, "data")
  cl.num = cls.num.vect.validity(cl.num, dim(data)[1], "cl.num")

  rep.num = check.rep.num(rep.num)
  subset.ratio = check.subset.ratio(data, subset.ratio)
  cl.num = cut.cl.num.pred(data, cl.num, subset.ratio)

  alg.cls = check.clust.wrappers(clust.alg)
  alg.prd = check.pred.wrapper(pred.alg)

  # run algorithms
  result = NULL

  if( alg.cls$fast )
  {
    result  = clust.stab.predict.hver( 
                data=data, cl.num=cl.num, 
                sample.num=rep.num, ratio=subset.ratio,
                clust.method=alg.cls$clust.method,
                clust.wrap=alg.cls$clust.wrap, 
                pred.method=alg.prd$pred.method
              )
   }
   else
   {
     result = clust.stab.predict.pver( 
                data=data, cl.num=cl.num, 
                sample.num=rep.num, ratio=subset.ratio,
                clust.method=alg.cls$clust.method,
                pred.method=alg.prd$pred.method
              )
   }

   return(result)
}
