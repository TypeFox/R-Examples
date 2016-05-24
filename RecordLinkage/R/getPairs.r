# getPairs.r: functions to view and edit record pairs

# moved to getPairs-methods.r:
#getPairs <- function(rpairs,max.weight=Inf,min.weight=-Inf,
#					single.rows=FALSE, show="all",
#					sort=!is.null(rpairs$Wdata))


setGeneric(
  name = "editMatch",
  def = function(rpairs) standardGeneric("editMatch")
)

setMethod(
  f = "editMatch",
  signature = "RecLinkData",
  def = function (rpairs)
  {
    if (rpairs$type=="deduplication")
    {
        data1=rpairs$data
        data2=data1
    } else
    {
        data1=rpairs$data1
        data2=rpairs$data2
    }
    # get number of columns depending on type (linkage or dedup)
    numCol <- switch(rpairs$type, deduplication = ncol(rpairs$data),
           linkage = ncol(rpairs$data1))
    p=data.frame(data1[rpairs$pairs[,1],],
                     data2[rpairs$pairs[,2],],
                     matrix("",nrow=nrow(rpairs$pairs),
                        ncol=numCol))

    # Transformation of "is_match" allows displaying empty cells in the table   
    p=matrix(as.matrix(t(p))[TRUE],nrow=nrow(p)*3,byrow=TRUE)
    # unlist(lapply) instead of sapply to avoid matrix result
    p=data.frame(p,is_match=unlist(lapply(rpairs$pairs$is_match,function(x) c(x,"",""))))
    colnames(p)=c(colnames(data1),"is_match")
    p=edit(p)
    is_match=p[seq(1,nrow(p)-2,3),"is_match"]
    is_match=as.integer(levels(is_match)[as.integer(is_match)])
    rpairs$pairs$is_match <- is_match
    return(rpairs)
  }
)


setMethod(
  f = "editMatch",
  signature = "RLBigData",
  def = function (rpairs)
  {
    if (is(rpairs, "RLBigDataDedup"))
    {
        data1=rpairs@data
        data2=data1
    } else if (is(rpairs, "RLBigDataLinkage"))
    {
        data1=rpairs@data1
        data2=rpairs@data2
    } else stop("No support for class %s", class(rpairs))
    # get number of columns depending on type (linkage or dedup)
    numCol <- ncol(data1)
    p=data.frame(data1[rpairs@pairs[,1],],
                     data2[rpairs@pairs[,2],],
                     matrix("",nrow=nrow(rpairs@pairs),
                        ncol=numCol))

    # Transformation of "is_match" allows displaying empty cells in the table   
    p=matrix(as.matrix(t(p))[TRUE],nrow=nrow(p)*3,byrow=TRUE)
    # unlist(lapply) instead of sapply to avoid matrix result
    p=data.frame(p,is_match=unlist(lapply(as.ram(rpairs@pairs$is_match),function(x) c(x,"",""))))
    colnames(p)=c(colnames(data1),"is_match")
    p=edit(p)
    is_match=p[seq(1,nrow(p)-2,3),"is_match"]
    is_match=as.integer(levels(is_match)[as.integer(is_match)])
    rpairs@pairs$is_match <- ff(is_match)
    return(rpairs)
  }
)

