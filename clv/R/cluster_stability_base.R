
#-------------------------------------------------------------------------------------------#
#------------------- CLUSTER STABILITY -> STABLE POINTS APPROACH ---------------------------#

# cluster stability (optimal assignment/stable objects) approach for hierarchical algorithms 
# fast version

clust.stab.opt.assign.hver <- function(data, cl.num, sample.num, ratio, clust.method, clust.wrap)
{
  obj.num = dim(data)[1]

  cls.stabb.vec = rep( 0, length(cl.num) )

  obj.assgn.freq.mx = matrix( 0, length(cl.num), obj.num )
  obj.cls.freq.vec = rep(0, obj.num) 

  clust.tree = clust.method(data)
  for( j in 1:sample.num )
  { 
    smp = sort( sample( 1:obj.num, ratio*obj.num ) )
    subset = as.matrix(data[smp,])
    subset.tree = clust.method(subset)

    for( i in 1:length(smp) ) obj.cls.freq.vec[ smp[i] ] = obj.cls.freq.vec[ smp[i] ] + 1

    iter = 1
    for( clust.num in cl.num )
    {
      data.cls = clust.wrap( clust.tree, clust.num )[smp]
      subset.cls = clust.wrap( subset.tree, clust.num )

      section.matrix = cbind(smp,data.cls, subset.cls)
      cnf.mx = confusion.matrix(data.cls, subset.cls)

      tmp1 = 2
      tmp2 = 3
# be careful - trasposition possible
      if( dim(cnf.mx)[1] > dim(cnf.mx)[2] )
      {
# situation when rows means second clustering collumns means first clustering 
        cnf.mx = t(cnf.mx)
        tmp1 = 3
        tmp2 = 2
      }
      
      opt.assignment = ( .Call( "clv_optimalAssignment" , cnf.mx, PACKAGE="clv" ) + 1 )

# this part should be speed up in C language
      for( i in 1:dim(section.matrix)[1] )
      {
        if( opt.assignment[ section.matrix[i,tmp1] ] == section.matrix[i,tmp2] ) 
        obj.assgn.freq.mx[ iter, section.matrix[i,1] ] = obj.assgn.freq.mx[ iter, section.matrix[i,1] ] + 1
      }
      iter = iter + 1
    }
  }

  tmp.mx = matrix( 0, length(cl.num), max(cl.num) )
  
  for( iter in 1:length(cl.num) )
  {
    data.cls = clust.wrap(clust.tree, cl.num[iter])
    cls.size = cluster.size(data.cls, as.integer(cl.num[iter]))

    for( obj.id in 1:obj.num )
    if( obj.cls.freq.vec[obj.id] != 0 && cls.size[ data.cls[obj.id] ] != 0 )
    {
      tmp.mx[ iter, data.cls[obj.id] ] = 
      tmp.mx[ iter, data.cls[obj.id] ] + ( obj.assgn.freq.mx[iter,obj.id] / obj.cls.freq.vec[obj.id] ) / cls.size[ data.cls[obj.id] ]
    }
  }
  
  result = rep(0,length(cl.num))
  for( i in 1:length(cl.num))
  result[i] = min(tmp.mx[i, 1:cl.num[i]])
  
  return(result)
}

# cluster stability (optimal assignment/stable objects approach) for partitioning algorithms
# (can be applied also for hierarchical algorithms but it is not optimal - see: clust.stab.opt.assign.hver)

clust.stab.opt.assign.pver <- function(data, cl.num, sample.num, ratio, clust.method)
{
  obj.num = dim(data)[1]

  cls.stabb.vec = rep( 0, length(cl.num) )
  iter = 1

  for( cls.num in cl.num )
  {
    data.cls = clust.method(data, cls.num)
    obj.assgn.freq.vec = rep( 0, obj.num )
    obj.cls.freq.vec = rep( 0, obj.num )

    for( j in 1:sample.num )
    {
      smp = sort( sample( 1:obj.num, ratio*obj.num ) )

      subset = as.matrix(data[smp,])
      subset.cls = clust.method(subset, cls.num)

      clsm1 = cbind(1:obj.num, data.cls)
      clsm2 = cbind(smp,       subset.cls)

      section.matrix = cls.set.section(clsm1, clsm2)
      cls1 = as.integer(section.matrix[,2])
      cls2 = as.integer(section.matrix[,3])
      cnf.mx = confusion.matrix(cls1,cls2)

      tmp1 = 2
      tmp2 = 3

# be careful - trasposition possible
      if( dim(cnf.mx)[1] > dim(cnf.mx)[2] )
      {
# situation when rows means second clustering collumns means first clustering 
        cnf.mx = t(cnf.mx)
        tmp1 = 3
        tmp2 = 2
      }

      opt.assignment = ( .Call( "clv_optimalAssignment" , cnf.mx, PACKAGE="clv" ) + 1 )

# this part should be speed up in C language
      for( i in 1:dim(section.matrix)[1] )
      {
        obj.cls.freq.vec[ section.matrix[i,1] ] = obj.cls.freq.vec[ section.matrix[i,1] ] + 1
        if( opt.assignment[ section.matrix[i,tmp1] ] == section.matrix[i,tmp2] ) 
          obj.assgn.freq.vec[ section.matrix[i,1] ] = obj.assgn.freq.vec[ section.matrix[i,1] ] + 1
      }
    }

# we have everything what is needed to find the less stable cluster
    cls.att = cls.attrib(data, data.cls)
    tmp.vec = rep(0,length(cls.att$cluster.size))
    for( k in 1:obj.num )
    {
      if( obj.cls.freq.vec[k] != 0 && cls.att$cluster.size[ data.cls[k] ] != 0 )
      {
        tmp.vec[ data.cls[k] ] = tmp.vec[ data.cls[k] ] + (obj.assgn.freq.vec[k]/obj.cls.freq.vec[k])/cls.att$cluster.size[ data.cls[k] ]
      }
    }
    
    cls.stabb.vec[iter] = min(tmp.vec)
    iter = iter + 1
  }
  
  return(cls.stabb.vec)
}

#--------------------------------------------------------------------------------------------------#
#------------------- CLUSTER STABILITY -> DISTANCE BETWEEN CLUSTERING RESULTS ---------------------#

# version for all - hierarchical and aglomerative algorithms (for hierarcical is not optimal)
# arguments explanation:
# measure - vector with information about choosen measures, the sequence is: dt, sim.ind, rand, jacc

clust.stab.sim.ind.pver.internal <- function( data, clust.num, sample.num, ratio, clust.method, measure )
{
# prepare result matrix
  result = matrix(0, sample.num, length(measure[measure == TRUE]))

  obj.num = dim(data)[1]

  for( j in 1:sample.num )
  {
    smp1 = sort( sample( 1:obj.num, ratio*obj.num ) )
    smp2 = sort( sample( 1:obj.num, ratio*obj.num ) )

    d1 = as.matrix(data[smp1,])
    cls1 = clust.method(d1, clust.num)

    d2 = as.matrix(data[smp2,])
    cls2 = clust.method(d2, clust.num)

    clsm1 = cbind(smp1,cls1)
    clsm2 = cbind(smp2,cls2)

    m = cls.set.section(clsm1, clsm2)

    # protect situation if the section is empty
    if(dim(m)[1] < 1)
    {
      for(iter in 1:4) 
        if( measure[iter] == TRUE ) 
          result[j, iter] = NaN
    }
    else
    {
      cls1 = as.integer(m[,2])
      cls2 = as.integer(m[,3])

      iter = 1
      if( measure[dt.sim.ind.const] == TRUE ) 
      {
# remember - result might be NaN
        result[j, iter] = dot.product(cls1,cls2)
        iter = iter + 1
      }

      if( measure[optas.sim.ind.const] == TRUE )
      {
        si = similarity.index( confusion.matrix(cls1,cls2) )
        result[j, iter] = ((clust.num+1)/clust.num)*si - 1/clust.num
        iter = iter + 1
      }
  
# external measures - compare partitionings
      if( measure[rand.const] == TRUE || measure[jaccard.const] == TRUE ) 
      {
        std.ms = std.ext(cls1,cls2)
        if( measure[rand.const] == TRUE ) 
        {
          result[j, iter] = clv.Rand(std.ms)
          iter = iter + 1
        }

        if( measure[jaccard.const] == TRUE ) 
        {
          result[j, iter] = clv.Jaccard(std.ms)
          iter = iter + 1
        }
      }
    } 
  }

  return( result )
}

# ------- cluster stability - similarity index -> ver. for partitionings algorithms ------ #

clust.stab.sim.ind.pver <- function( data, cl.num, sample.num, ratio, clust.method, sim.ind )
{
  ind.num = length( sim.ind[sim.ind == TRUE] )

  result = vector("list", length=ind.num )
  for( i in 1:ind.num ) result[[i]] = as.data.frame(matrix(0,sample.num,length(cl.num)))

  tmp = 1
  for( cls.num in cl.num )
  {
#print(paste( i, "-ta liczba klastrow =", cls.num))
    cls.stab.mx = clust.stab.sim.ind.pver.internal(data, clust.num=cls.num, sample.num=sample.num, 
    ratio=ratio, clust.method=clust.method, measure=sim.ind)

    for( index.type in 1:ind.num ) result[[index.type]][,tmp] = cls.stab.mx[,index.type]
    tmp = tmp + 1
  }

  names = paste("C", cl.num, sep="")

  iter = 1
  for( i in 1:length(sim.ind) )
  {
    if(sim.ind[i] == TRUE)
    {
      colnames(result[[iter]]) = names
      names(result)[iter] = supp.cls.stab.sim.ind.vec.const[i]
      iter = iter + 1
    }
  }

  return(result)
}

# ------- cluster stability - similarity index -> ver. for hierarhical algorithms ------ #

clust.stab.sim.ind.hver <- function( data, cl.num, sample.num, ratio, clust.method, clust.wrap, sim.ind )
{
  obj.num = dim(data)[1]
 
  ind.num = length( sim.ind[sim.ind == TRUE] )

  result = vector("list", length=ind.num )
  for( i in 1:ind.num ) result[[i]] = as.data.frame(matrix(0,sample.num,length(cl.num)))

  clust.tree = clust.method(data)
  for( j in 1:sample.num )
  { 
    smp = sort( sample( 1:obj.num, ratio*obj.num ) )
    subset = as.matrix(data[smp,])
    subset.tree = clust.method(subset)

    iter = 1
    for( clust.num in cl.num )
    {
      data.cls = clust.wrap( clust.tree, clust.num )[smp]
      subset.cls = clust.wrap( subset.tree, clust.num )

      tmp = 1
      if( sim.ind[dt.sim.ind.const] == TRUE )
      {
# remember - result might be NaN 
        result[[tmp]][j, iter] = dot.product(data.cls, subset.cls)
        tmp = tmp + 1 
      }

      if( sim.ind[optas.sim.ind.const] == TRUE )
      {
        si = similarity.index( confusion.matrix(data.cls, subset.cls) )
        result[[tmp]][j, iter] = ((clust.num+1)/clust.num)*si - 1/clust.num
        tmp = tmp + 1
      }
# external measures - compare partitionings
      if( sim.ind[rand.const] == TRUE || sim.ind[jaccard.const] == TRUE ) 
      {
        std.ms = std.ext(data.cls, subset.cls)
        if( sim.ind[rand.const] == TRUE ) 
        {
          result[[tmp]][j, iter] = clv.Rand(std.ms)
          tmp = tmp + 1
        }
        if( sim.ind[jaccard.const] == TRUE ) result[[tmp]][j, iter] = clv.Jaccard(std.ms)
      }
      iter = iter + 1 
    }
  }

  names = paste("C", cl.num, sep="")

  iter = 1
  for( i in 1:length(sim.ind) )
  {
    if(sim.ind[i] == TRUE)
    {
      colnames(result[[iter]]) = names
      names(result)[iter] = supp.cls.stab.sim.ind.vec.const[i]
      iter = iter + 1
    }
  } 

  return(result)
}
