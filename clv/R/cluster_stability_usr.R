
# ##################################################### #
# Cluster stability with user's wrappers
# ##################################################### #

# some globals

please.see.doc.txt = "Please see DOCUMENTATION for more details.\n"
exec.example.txt   = "clust.method(data = matrix(as.numeric(1:100),20,5), clust.num=5)"

print.see.doc.txt <- function() { cat(please.see.doc.txt) }
cls.alg.is.not.good <-function() stop("Cluster algorithm is not well deined!")

# check input data with this function
check.res <- function(res, data.size, msg1, msg2)
{
  if( (data.class(res) == "try-error") )
  {
    msg1()
    cat("\nretunrs the following error:")
    print(res)
    msg2()
    stop("Error!")
  }

  if( !is.integer(res) )
  {
    msg1()
    cat("\nis not integer vector but:")
    print(res)
    msg2()
    cls.alg.is.not.good()
  }

  if( length(res) != data.size )
  {
    msg1()
    cat("\nlength of the integer vector is not equal to the data size.")
    cat("\ndata size is:")
    cat(data.size)
    cat("\nand vector length is:")
    print(length(res))
    msg2()
    cls.alg.is.not.good()
  }
  # these numbers depends on what it in the msg1 (normally there is an information 
  # about number of clusters you should have at most)
  if( min(res) < 1 || max(res) > 5 )
  {
    msg1()
    cat("\nis integer vector where some elements are not filling following requirement:\n")
    cat("{ foreach {\"elem\" in \"res\" }  1 <= elem <= clust.num }\n")
    cat("instead we have this vector:\n")
    print(res)
    msg2()
    cls.alg.is.not.good()
  }
}

check.clust.wrappers <- function(clust.alg)
{
  print.std.error.msg <- NULL

  if( !is.function(clust.alg) && class(clust.alg) != "cls.alg" )
  {
      stop( clv.paste( "\"clust.alg\" parameter is nor a function neither \"clustalg\" type.\n", please.see.doc.txt ) )
  }

  result <- NULL
  
  if( is.function(clust.alg) )
  {
    arg.names<- names(formals(clust.alg))
    if( length(arg.names) != 2 || arg.names[1] != "data" || arg.names[2] != "clust.num" )
      stop( clv.paste( "\"clust.alg\" is a function with incorrect argument count or with wrong parameter names.\n", please.see.doc.txt ) ) 

    # check if algorithm works more less
    data.size = 20
    dt = matrix(as.numeric(1:100),data.size,5)
    cln = 5
    res = try( clust.alg(data=dt, clust.num=cln) )

    print.std.error.msg <-function ()
    {
      cat("\nFunction \"clust.method\" is not correct.")
      cat("\nAn output comming from execution:\n")
      cat(exec.example.txt)
    }
    check.res(res, data.size, print.std.error.msg, print.see.doc.txt)

    result = list(clust.method=clust.alg, fast=FALSE)
  }
  else
  {
    if( clust.alg$fast == TRUE ) 
    {
      result = clust.alg
    }
    else 
    {
      tmp <- function(data, clust.num) clust.alg$clust.wrap( clust.alg$clust.method(data), clust.num )
      result = list(clust.method=tmp, fast=FALSE)
    }
  }

  return(result)
}



# helper VISIBLE TO THE END USER function to create well defined object for hierarhical cluster algorithms
cls.alg <-function(clust.method, clust.wrap, fast=TRUE)
{
  if( !is.function(clust.method) )
    stop( clv.paste( "\n\"clust.method\" parameter is not a function type.\n", please.see.doc.txt ) )
  if( !is.function(clust.wrap) )
    stop( clv.paste( "\n\"clust.wrap\" parameter is not a function type.\n", please.see.doc.txt ) )
  
  if( !is.logical(fast) ) fast = TRUE
  
  arg.names<- names(formals(clust.method))
  if( length(arg.names) != 1 || arg.names[1] != "data" )
    stop( clv.paste( "\n\"clust.method\" is a function with incorrect count of arguments or with wrong parameter names.\n", please.see.doc.txt ) )

  arg.names<- names(formals(clust.wrap))
  if( length(arg.names) != 2 || arg.names[1] != "clust.res" || arg.names[2] != "clust.num" )
    stop( clv.paste( "\n\"clust.wrap\" is a function with incorrect count of arguments or with wrong parameter names.\n", please.see.doc.txt ) )

  # check if algorithm works more less
  data.size = 20
  dt = matrix(as.numeric(1:100),data.size,5)
  cln = 5
  res = try( clust.wrap( clust.res=clust.method( data=dt ), clust.num=cln ) )

  print.std.error.msg <-function()
  {
    cat("\nCombination of functions \"clust.method\" and \"clust.wrap\" is not correct.")
    cat("\nAn output comming from execution:")
    cat(exec.example.txt)
  }

  check.res(res, data.size, print.std.error.msg, print.see.doc.txt)

  result = list(clust.method=clust.method, clust.wrap=clust.wrap, fast=fast)
  class(result) <- "cls.alg"

  return(result)
}


# cluster stability ->similarity index approach
 
cls.stab.sim.ind.usr <- function( 
         data, cl.num,
         clust.alg,
         sim.ind.type=c("dot.pr","sim.ind"),
         rep.num=10, 
         subset.ratio=0.75
         )
{
  # check input arguments
  data = data.validity(data, "data")
  cl.num = cls.num.vect.validity(cl.num, dim(data)[1], "cl.num")

  rep.num = check.rep.num(rep.num)
  subset.ratio = check.subset.ratio(data, subset.ratio)
  cl.num = cut.cl.num(data, cl.num, subset.ratio)

  sim.ind.type.bool = check.avail.methods(sim.ind.type, "sim.ind.type", supp.cls.stab.sim.ind.vec.const )

  alg = check.clust.wrappers(clust.alg)

  # run algorithms
  result = NULL

  if( alg$fast )
  {
    result = clust.stab.sim.ind.hver( 
      data=data, cl.num=cl.num, 
      sample.num=rep.num, ratio=subset.ratio,
      clust.method=alg$clust.method,
      clust.wrap=alg$clust.wrap,
      sim.ind=sim.ind.type.bool
    )
  }
  else
  {
    result = clust.stab.sim.ind.pver(
      data=data, cl.num=cl.num, 
      sample.num=rep.num, ratio=subset.ratio,
      clust.method=alg$clust.method,
      sim.ind=sim.ind.type.bool
    )
  }

  return(result)
}

# cluster stability -> optimal assignment approach
cls.stab.opt.assign.usr <- function( 
         data, cl.num,
         clust.alg,
         rep.num=10, 
         subset.ratio=0.75
         )
{
  # check input arguments
  data = data.validity(data, "data")
  cl.num = cls.num.vect.validity(cl.num, dim(data)[1], "cl.num")

  rep.num = check.rep.num(rep.num)
  subset.ratio = check.subset.ratio(data, subset.ratio)
  
  cl.num = cut.cl.num(data, cl.num, subset.ratio)

  alg = check.clust.wrappers(clust.alg)

  result = NULL

  if( alg$fast )
  {
    result = clust.stab.opt.assign.hver( 
      data=data, cl.num=cl.num, 
      sample.num=rep.num, ratio=subset.ratio,
      clust.method=alg$clust.method,
      clust.wrap=alg$clust.wrap
    )
  }
  else
  {
    result = clust.stab.opt.assign.pver(
      data=data, cl.num=cl.num, 
      sample.num=rep.num, ratio=subset.ratio,
      clust.method=alg$clust.method
    )
  }

  return(result)
}
