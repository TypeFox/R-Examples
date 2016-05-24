#<<BEGIN>>
mcdata <- function(data, type=c("V","U","VU","0"), nsv=ndvar(), nsu=ndunc(),nvariates=1,outm="each")
#TITLE Build mcnode Objects from Data or other mcnode Objects
#NAME mcnode
#DESCRIPTION
# Creates a \samp{mcnode} object from a vector, an array or a \samp{mcnode}.
#KEYWORDS methods
#INPUTS
#{data}<<The numeric/logical vector/matrix/array of data or the \samp{mcnode} object.>>
#[INPUTS]
#{type}<<The type of node to be built. By default, a \samp{"V"} node.>>
#{nsv}<<The variability dimension (\samp{type="V"} or \samp{type="VU"}) of the node.
#By default: the current value in \code{\link{mc.control}}>>
#{nsu}<<The uncertainty dimension (\samp{type="U"} or \samp{type="VU"}) of the node.
#By default: the current value in \code{\link{mc.control}}>>
#{nvariates}<<The number of variates.
#By default: 1>>
#{outm}<<The output of the \samp{mcnode} for multivariates nodes. May be "each" (default)
#if output should be provided for each variates considered independently, "none" for no output
#or a vector of name of function(s) (as a character string) that will be applied on the variates dimension
#before any output (ex: \samp{"mean"}, \samp{"median"}, \samp{c("min", "max")}). The function should have no other arguments
#and send one value per vector of values (ex. do not use \samp{"range"}).
#Note that the \samp{outm} attribute may be changed at any time using the \code{\link{outm}} function.>>
#VALUE
#An \samp{mcnode} object.
#DETAILS
#A \samp{mcnode} object is the basic element of a \code{\link{mc}} object.
#It is an array of dimension \samp{(nsv x nsu x nvariates)}. Four types of \samp{mcnode} exists:</>
#{*}<<\samp{"V" mcnode}, for "Variability", are arrays of dimension \samp{(nsv x 1 x nvariates)}. The alea in the data should
#reflect variability of the parameter.>>
#{*}<<\samp{"U" mcnode}, for "Uncertainty", are arrays of dimension \samp{c(1 x nsu x nvariates)}. The alea in the data should
#reflect uncertainty of the parameter.>>
#{*}<<\samp{"VU" mcnode}, for "Variability and Uncertainty", are arrays of dimension \samp{(nsv x nsu x nvariates)}. The alea in the data
#reflects separated variability (in rows) and uncertainty (in columns) of the parameter.>>
#{*}<<\samp{"0" mcnode}, for "Neither Variability or Uncertainty", are arrays of dimension \samp{(1 x 1 x nvariates)}. No alea is
#considered for these nodes. \samp{"0" mcnode} are not necessary in the univariate context (use scalar instead) but
#may be useful for operations on multivariate nodes.>>
#
#Multivariate nodes (i.e. \samp{nvariates != 1}) should be used for multivariate distributions implemented in \samp{mc2d}
#(\code{\link{rmultinomial}}, \code{\link{rmultinormal}}, \code{\link{rempiricalD}} and \code{\link{rdirichlet}}).
#
#For security, recycling rules are limited to fill the array using \samp{data}. The general rules is that recycling is
#only permitted to fill a dimension from 1 to the final size of the dimension.</>
#If the final dimension of the node is \samp{(nsv x nsu x nvariates)} (with \samp{nsv = 1} and \samp{nsu = 1} for \samp{"0"} nodes,
#\samp{nsu = 1} for \samp{"V"} nodes
#and \samp{nsv = 1} for \samp{"U"} nodes), \samp{mcdata} accepts :
#{*}<<Vectors of length \samp{1} (recycled on all dimensions), vectors of length \samp{(nsv * nsu)} (filling first the dimension of variability, then
#the dimension of uncertainty then recycling on nvariates),
#or vectors of length \samp{(nsv * nsu * nvariates)} (filling first the dimension of variability, then the uncertainty,
#then the variates).>>
#{*}<<Matrixes of dimensions \samp{(nsv x nsu)}, recycling on variates.>>
#{*}<<Arrays of dimensions \samp{(nsv x nsu x nvariates)} or \samp{(nsv x nsu x 1)}, recycling on variates.>>
#{*}<<For \samp{data} as \samp{mcnode}, recycling is dealt to proper fill the array:>>
#{#}<<a \samp{"V"} node accepts a \samp{"0"} node of dimension \samp{(1 x 1 x nvariates)} (recycling on variability)
#or of dimension \samp{(1 x 1 x 1)} (recycling on variability and variates),
#or a \samp{"V"} node of dimension \samp{(nsv x 1 x nvariates)}
#or \samp{(nsv x 1 x 1)} (recycling on variates),>>
#{#}<<a \samp{"U"} node accepts a \samp{"0"} node of dimension \samp{(1 x 1 x nvariates)} (recycling on uncertainty)
#or of dimension \samp{(1 x 1 x 1)} (recycling on uncertainty and variates),
#or a \samp{"U"} node of dimension \samp{(1 x nsu x nvariates)},
#or \samp{(1 x nsu x 1)} (recycling on variates),>>
#{#}<<a \samp{"VU"} node accepts a \samp{"0"} node of dimension \samp{(1 x 1 x nvariates)} (recycling on varaiability and uncertainty)
#or of dimension \samp{(1 x 1 x 1)} (recycling on variability, uncertainty and variates),
#a \samp{"U"} node of dimension \samp{(1 x nsu x nvariates)}(recycling "by row" on the variability dimension),
#or of dimension \samp{(1 x nsu x 1)}(recycled "by row" on the variability dimension then on variates),
#a \samp{"V"} node of dimension \samp{(nsv x 1 x nvariates)}(recycling on the uncertainty dimension)
#or of dimension \samp{(nsv x 1 x 1)}(recycled on the uncertainty dimension then on variates),
#and a \samp{"VU"} node of dimension \samp{(nsv x nsu x nvariates)} or of dimension \samp{(nsv x nsu x 1)} (recycling on variates).>>
#{#}<<a \samp{"0"} node accepts a \samp{"0"} node of dimension \samp{(1 x 1 x nvariates)} or
#\samp{(1 x 1 x 1)} (recycling on variates).>>
#
#
#SEE ALSO
#\code{\link{mcstoc}} to build a stochastic \samp{mcnode} object, \code{\link{mcprobtree}} to build a
#stochastic node fro a probability tree.</>
#\code{\link{Ops.mcnode}} for operations on \samp{mcnode} objects.</>
#\code{\link{mc}} to build a Monte-Carlo object.</> </>
#Informations about an mcnode: \code{\link{is.mcnode}}, \code{\link{dimmcnode}}, \code{\link{typemcnode}}.</>
#To build a correlation structure between \samp{mcnode}: \code{\link{cornode}}.</>
##To apply a function on a \samp{mcnode} object: \code{\link{mcapply}}</>
#To study \samp{mcnode} objects: \code{\link{print.mcnode}}, \code{\link{summary.mcnode}}, \code{\link{plot.mcnode}},
#\code{\link{converg}}, \code{\link{hist.mcnode}}</>
#To modify \samp{mcnode} objects:
##\code{\link{subset.mc}}
#\code{\link{NA.mcnode}}</>
##To transform \samp{mc} objects in vector or matrix: \code{\link{unmcnode}}.
#EXAMPLE
#oldvar <- ndvar()
#oldunc <- ndunc()
#ndvar(3)
#ndunc(5)
#
#(x0 <- mcdata(100,type="0"))
#mcdata(matrix(100),type="0")
#
#(xV <- mcdata(1:ndvar(),type="V"))
#mcdata(matrix(1:ndvar(),ncol=1),type="V")
#
#(xU <- mcdata(10*1:ndunc(),type="U"))
#mcdata(matrix(10*1:ndunc(),nrow=1),type="U")
#
#(xVU <- mcdata(1:(ndvar()*ndunc()),type="VU"))
#mcdata(matrix(1:(ndvar()*ndunc()),ncol=5,nrow=3),type="VU")
#
###Do not use
#\dontrun{
#mcdata(matrix(1:5,nrow=1),type="VU")
#}
###use instead
#mcdata(mcdata(matrix(1:ndunc(),nrow=1),type="U"),"VU")
###or
#mcdata(matrix(1:ndunc(),nrow=1),type="U") + mcdata(0,"VU")
#
#mcdata(x0,type="0")
#
#mcdata(x0,type="V")
#mcdata(xV,type="V")
#
#mcdata(x0,type="U")
#mcdata(xU,type="U")
#
#mcdata(x0,type="VU")
#mcdata(xU,type="VU")
#mcdata(xV,type="VU")
#
###Multivariates
#(x0M <- mcdata(1:2,type="0",nvariates=2))
#mcdata(1,type="0",nvariates=2)
#
#(xVM <- mcdata(1:(2*ndvar()),type="V",nvariates=2))
#mcdata(1:ndvar(),type="V",nvariates=2)
#mcdata(array(1:(2*ndvar()),dim=c(3,1,2)),type="V",nvariates=2)
#
#mcdata(1,type="V",nvariates=2)
#mcdata(x0,type="V",nvariates=2)
#mcdata(x0M,type="V",nvariates=2)
#mcdata(xV,type="V",nvariates=2)
#mcdata(xVM,type="V",nvariates=2)
#
#(xUM <- mcdata(10*(1:(2*ndunc())),type="U",nvariates=2))
#mcdata(array(10*(1:(2*ndunc())),dim=c(1,5,2)),type="U",nvariates=2)
#
#mcdata(1,type="U",nvariates=2)
#mcdata(x0,type="U",nvariates=2)
#mcdata(x0M,type="U",nvariates=2)
#mcdata(xU,type="U",nvariates=2)
#mcdata(xUM,type="U",nvariates=2)
#
#(xVUM <- mcdata(1:(ndvar()*ndunc()),type="VU",nvariates=2))
#mcdata(array(1:(ndvar()*ndunc()),dim=c(3,5,2)),type="VU",nvariates=2)
#
#mcdata(1,type="VU",nvariates=2)
#mcdata(x0,type="VU",nvariates=2)
#mcdata(x0M,type="VU",nvariates=2)
#mcdata(xV,type="VU",nvariates=2)
#mcdata(xVM,type="VU",nvariates=2)
#mcdata(xU,type="VU",nvariates=2)
#mcdata(xUM,type="VU",nvariates=2)
#mcdata(xVU,type="VU",nvariates=2)
#mcdata(xVUM,type="VU",nvariates=2)
#
#ndvar(oldvar)
#ndunc(oldunc)

#CREATED 08-01-25
#--------------------------------------------
{
  type <- match.arg(type)

    if(!is.character(outm)  || (outm != "none" && outm != "each" && !all(sapply(outm,exists,mode="function"))))
      stop("outm should be 'none','each' or the name a valid function")

  dimf <- switch(type,"V"=c(nsv,1,nvariates),"U"=c(1,nsu,nvariates),"VU"=c(nsv,nsu,nvariates),"0"=c(1,1,nvariates))

  if(inherits(data,"mcnode") ){             # data = mcnode
    typem <- attr(data,"type")
    dimd <- dim(data)
    err <- paste("The output node dimension is not compatible with the input node dimension. Should be of dim: ",paste(dimf,collapse=" "))
    if(dimd[3] != 1 && dimd[3]!=dimf[3]) stop(err)
    if(typem=="VU") {
      if(dimd[1]!=dimf[1] || dimd[2]!=dimf[2]) stop(err)
      if(type!="VU") stop("Incompatible node type. A 'VU' mcnode makes only a 'VU'mcnode.")
      data <- array(data,dim=dimf)}                                             # recycle on nvariates
    else if(typem=="U") {
      if(dimd[2]!=dimf[2]) stop(err)
      if(type!="VU" && type!="U") stop("Incompatible node type. A 'U' mcnode makes a 'U' or a 'VU' mcnode.")
      data <- array(apply(data,3,rep,each=dimf[1]),dim=dimf)}                   # recycle on nsu and nvariates
    else if(typem=="V") {
      if(dimd[1]!=dimf[1]) stop(err)
      if(type!="VU" && type!="V") stop("Incompatible node type. A 'V' mcnode makes a 'V' or a 'VU' mcnode.")
      data <- array(apply(data,3,rep,times=dimf[2]),dim=dimf)}                  # recycle on nsv and nvariates
    else if(typem=="0") {
      data <- array(apply(data,3,rep,each=dimf[1]*dimf[2]),dim=dimf)}           # recycle on nvariates
  }
  else {                                    # data non mcnode
      if(!is.numeric(data) && !is.logical(data)) stop("data should be numeric or logical")
      if(is.vector(data)){
        l <- length(data)
        if(l != 1 && l != dimf[1]*dimf[2] && l!= prod(dimf)) stop("The vector size is not compatible
          with the node dimension. Length should be 1 or n=",dimf[1]*dimf[2]," or n=",prod(dimf))
        data <- array(data,dim=dimf)}

      else if(is.array(data)){
        dimd <- dim(data)
        l <- length(dimd)
        if(l > 3) stop("Maximum accepted dim of arrays is 3.")
        if(l == 2) dimd <- c(dimd,0)  # just to simplify the following tests
        if( (type=="VU" && dimd[1]== nsv && dimd[2]==nsu && dimd[3] %in% c(0,1,nvariates)) ||
            (type=="V" &&  dimd[1]== nsv && dimd[2]==1 && dimd[3] %in% c(0,1,nvariates)) ||
            (type=="U" &&  dimd[1]== 1 && dimd[2]==nsu && dimd[3] %in% c(0,1,nvariates)) ||
            (type=="0" &&  dimd[1]== 1 && dimd[1]==1 && dimd[3] %in% c(0,1,nvariates)) )
        data <- array(data,dim=dimf)
        else stop("The array size is not compatible with the node dimension. Should be of dim: ",paste(dimf,collapse=" "))
        }

      else stop("data should be a vector, a matrix, an array or a mcnode")}

  class(data) <- "mcnode"
  attr(data,"type") <- type
  attr(data,"outm") <- outm
  return(data)
}
