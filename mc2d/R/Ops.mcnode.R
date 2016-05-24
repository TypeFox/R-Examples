#<<BEGIN>>
Ops.mcnode <- function(e1,e2)
#TITLE Operations on mcnode Objects
#DESCRIPTION
# This function alters the way operations are performed on \samp{mcnode} objects for a better consistancy of the theory.
#KEYWORDS utilities
#INPUTS
#{e1}<<An \samp{mcnode} object, a vector or an array.>>
#{e2}<<An optionnal \samp{mcnode} object, a vector or a matrix with at least one of both objects as an \samp{mcnode}.>>
#VALUE
#The results as a \samp{mcnode} object.
#DETAILS
#This method will be used for any of the Group \code{\link{Ops}} functions.
#
#The rules are as following (illustrated with a \samp{+} function and ignoring the \samp{nvariates} dimension):
#{*}<<\samp{0 + 0 = 0};>>
#{*}<<\samp{0 + V = V}: classical recycling of the scalar;>>
#{*}<<\samp{0 + U = U}: classical recycling of the scalar;>>
#{*}<<\samp{0 + VU = VU}: classical recycling of the scalar;>>
#{*}<<\samp{V + V = V}: if both of the same \samp{(nsv)} dimension;>>
#{*}<<\samp{V + U = VU}: the \samp{U} object will be recycled "by row". The \samp{V} object will be recycled classically "by column";>>
#{*}<<\samp{V + VU = VU}: if the dimension of the \samp{V} is \samp{(nsv)} and the dimension of the \samp{VU} is \samp{(nsv x nsu)}. The \samp{V} object will be recycled classically "by column";>>
#{*}<<\samp{U + U = U}: if both of the same \samp{(nsu)} dimension;>>
#{*}<<\samp{U + VU = VU}: if the dimension of the \samp{U} is \samp{(nsu)} and the dimension of the \samp{VU} is \samp{(nsv x nsu)}. The \samp{U} object will be recycled "by row";>>
#{*}<<\samp{VU + VU = VU}: if the dimension of the \samp{VU} nodes is \samp{(nsu x nsv)};>>
#
#A vector or an array may be combined with an \samp{mcnode} of size \samp{(nsv x nsu)} if an \samp{mcnode} of this dimension
#may be built from this vector/array using the \samp{mcdata} function. See \code{\link{mcdata}} for the rules.
#
#The \samp{outm} attribute is transferred as following: \samp{each + each = each}; \samp{none + other = other};
#\samp{other1 + other2 = other1}. The \samp{outm} attribute of the resulting node may be changed using the \code{\link{outm}} function.
#
#For multivariate nodes, a recycling on the \samp{nvariates} dimension is done if a \samp{(nsu x nsv x nvariates)} node
#is combined with a \samp{(nsu x nsv x 1)} node.

#SEE ALSO
#\code{\link{mcdata}}, \code{\link{mcstoc}}
#EXAMPLE
#oldvar <- ndvar()
#oldunc <- ndunc()
#ndvar(30)
#ndunc(20)
#
### Given
#x0 <- mcdata(3,type="0")
#xV <- mcdata(1:ndvar(),type="V")
#xU <- mcdata(1:ndunc(),type="U")
#xVU <- mcdata(1:(ndunc()*ndvar()),type="VU")
#x0M <- mcdata(c(5,10),type="0",nvariates=2)
#xVM <- mcdata(1:(2*ndvar()),type="V",nvariates=2)
#xUM <- mcdata(1:(2*ndunc()),type="U",nvariates=2)
#xVUM <- mcdata(1:(2*(ndunc()*ndvar())),type="VU",nvariates=2)
#
### All possible combinations
### "0"
#-x0
#x0 + 3
#
### "V"
#-xV
#3 + xV
#xV * (1:ndvar())
#xV * x0
#xV - xV
#
### "U"
#-xU
#xU + 3
#(1:ndunc()) * xU
#xU * x0
#xU - xU
#
### Watch out the resulting type
#xV + xU
#xU + xV
#
### "VU"
#-xVU
#3 + xVU
#(1:(ndunc()*ndvar())) * xVU
#xVU + xV
#x0 + xVU
#xU + xVU
#xVU - xVU
#
### Some Multivariates
#x0M+3
#xVM * (1:ndvar())
#xVM - xV
#xUM - xU
#xVUM - xU

#CREATED 08-01-25
#--------------------------------------------
{
  err <- "Incompatible mcnode dimensions"

  if(missing(e2)) {    # only e1 like -X
    type <- attr(e1,"type")
    outm <- attr(e1,"outm")
    dimf <- dim(e1)
  }   

  else {
      if(!inherits(e1,"mcnode")) {   # only e2 mcnode
          dimf <- dim(e2)
          type <- attr(e2,"type")
          outm <- attr(e2,"outm")
        
          if(!is.numeric(e1) && !is.logical(e1)) stop("e1 should be numeric or logical")

          if(is.vector(e1)){
            l <- length(e1)
            if(l != 1 && l != dimf[1]*dimf[2] && l!= prod(dimf)) stop("The vector size is not compatible
              with the node dimension. Length should be 1 or n=",dimf[1]*dimf[2]," or n=",prod(dimf))
            }
    
          else if(is.array(e1)){
            dim1 <- dim(e1)
            l <- length(dim1)
            if(l > 3) stop("Maximum accepted dim of arrays is 3.")
            if(l == 2) dim1 <- c(dim1,0)  # just to simplify the following tests
            if(! (type=="VU" && dim1[1]== dimf[1] && dim1[2]==dimf[2] && dim1[3] %in% c(0,1,dimf[3])) ||
                 (type=="V" &&  dim1[1]== dimf[1] && dim1[2]==1 && dim1[3] %in% c(0,1,dimf[3])) ||
                 (type=="U" &&  dim1[1]== 1 && dim1[2]==dimf[2] && dim1[3] %in% c(0,1,dimf[3])) ||
                 (type=="0" &&  dim1[1]== 1 && dim1[1]==1 && dim1[3] %in% c(0,1,dimf[3])) )
                  stop("The array size is not compatible with the node dimension. Should be of dim: ",paste(dimf,collapse=" "))
            }
          else stop("data should be a vector, a matrix, an array or a mcnode")
      }
  
      else {
          if(!inherits(e2,"mcnode")) {   # only e1 mcnode
              dimf <- dim(e1)
              type <- attr(e1,"type")
              outm <- attr(e1,"outm")
            
              if(!is.numeric(e2) && !is.logical(e2)) stop("e2 should be numeric or logical")
    
              if(is.vector(e2)){
                l <- length(e2)
                if(l != 1 && l != dimf[1]*dimf[2] && l!= prod(dimf)) stop("The vector size is not compatible
                  with the node dimension. Length should be 1 or n=",dimf[1]*dimf[2]," or n=",prod(dimf))
                }
        
              else if(is.array(e2)){
                dim2 <- dim(e2)
                l <- length(dim2)
                if(l > 3) stop("Maximum accepted dim of arrays is 3.")
                if(l == 2) dim2 <- c(dim2,0)  # just to simplify the following tests
                if(! (type=="VU" && dim2[1]== dimf[1] && dim2[2]==dimf[2] && dim2[3] %in% c(0,1,dimf[3])) ||
                     (type=="V" &&  dim2[1]== dimf[1] && dim2[2]==1 && dim2[3] %in% c(0,1,dimf[3])) ||
                     (type=="U" &&  dim2[1]== 1 && dim2[2]==dimf[2] && dim2[3] %in% c(0,1,dimf[3])) ||
                     (type=="0" &&  dim2[1]== 1 && dim2[1]==1 && dim2[3] %in% c(0,1,dimf[3])) )
                      stop("The array size is not compatible with the node dimension. Should be of dim: ",paste(dimf,collapse=" "))
                }
              else stop("data should be a vector, a matrix, an array or a mcnode")
          }
      
          else {                            # e1 and e2 mcnode
              dim1 <- dim(e1)
              dim2 <- dim(e2)
              if(!(dim1[3] == 1 || dim2[3] == 1 || dim1[3] == dim2[3])) stop(err)    # Controle dimension 3
              dimf <- pmax(dim1,dim2)
        
              outm1 <- attr(e1,"outm")
              outm2 <- attr(e2,"outm")
              if(is.null(outm1) && is.null(outm2)) outm <- "each"
              else if(is.null(outm1)) outm <- outm2
              else if(is.null(outm2)) outm <- outm1
              else if(outm1 == "each" || outm2 == "each") outm <- "each"
                else if(outm1=="none") outm <- outm2
                  else outm <- outm1
        
              type1 <- attr(e1,"type")
              type2 <- attr(e2,"type")
              if(type1==type2){ 
                type <- type1 
                if(any(dim1[1:2] != dim2[1:2])) stop(err)        # 0+0, U+U, V+V, VU + VU, gere les deux premieres dimensions 
              }                                                                                # pas de pb pour la troisieme (recycle)
        
              else if(type1=="0"){ 
                type <- type2                                                                    # 0 + others. gere la troisieme dim si necessaire
                if(dim1[3] != 1) e1 <- rep(e1, each=dimf[1]*dimf[2])
                }
              
              else if(type2=="0"){ 
                type <- type1 
                if(dim2[3] != 1) e2 <- rep(e2, each=dimf[1]*dimf[2])
                }
        
              else {                                                                                 #U+V, V+U, V+VU and U+VU
                type <- "VU"
        
                if(type1=="U") {                                                                    #U+V and U+VU
                  if(type2=="VU" && dim1[2]!= dim2[2]) stop(err)
                  e1 <- apply(e1,3,matrix,ncol=dimf[2],nrow=dimf[1],byrow=TRUE)                    # recycling dim 1,2
                  }
                
                else if(type2=="U") {                                                               #V+U and VU+U dim 1,2
                  if(type1=="VU" && dim1[2]!= dim2[2]) stop(err)
                  e2 <- apply(e2,3,matrix,ncol=dimf[2],nrow=dimf[1],byrow=TRUE)                    # recycling dim 1,2
                  }
        
                if(type1=="V") {
                  if(type2=="VU" && dim1[1]!= dim2[1]) stop(err)
                  if(dim1[3]!=1 && dimf[2]!=1) e1 <- apply(e1,3,matrix,ncol=dimf[2],nrow=dimf[1])  # necessary recycling of V for arrays
                  }
                
                else if(type2=="V"){
                  if(type1=="VU" && dim1[1]!= dim2[1]) stop(err)
                  if(dim2[3]!=1 && dimf[2]!=1) e2 <- apply(e2,3,matrix,ncol=dimf[2],nrow=dimf[1])
                  }
              }
          } # fin e1 e2 mcnode
    } # fin else

    e1 <- as.vector(e1)
    e2 <- as.vector(e2)

    } # fin else
    

    res <- NextMethod(.Generic)
    res <- array(res,dim=dimf)
    class(res) <- "mcnode"
    attr(res,"type") <- type
    attr(res,"outm") <- outm
    return(res)
}

