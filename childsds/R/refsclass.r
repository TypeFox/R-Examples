## Class refs
#' Class of references
#'
#' Container for reference tables 
#'
#' @slot name of the reference group
#' @slot ref List of references, each reference refers to one item and contains age, l, m, and s values for both genders
#' @slot citation information about the sources of the references
#' @exportClass refs
#' @examples
#' data(kiggs.ref)
#' print(kiggs.ref)
#' data(ukwho.ref)
#' print(ukwho.ref)
#' data(who.ref)
#' print(who.ref)
setClass(
  Class = "refs",
  representation=representation(
    name = "character",
    ref = "list",
    citation = "list"
    )
  )


#' Shows \code{items} and additional information of an \code{refs} object
#' @param object \code{refs} object
#' @rdname show-method
#' @importFrom methods show
setMethod("show","refs",
          function(object){
            cat("*** Class of Reference Values ***\n")
            cat("*** References and Information:\n"); print(object@citation)
            cat("*** Items:\n");
            print(lapply(object@ref,function(liste){
  data.frame(sex=c("male","female"),minage=c(min(liste[["male"]]$age,na.rm=T),min(liste[["female"]]$age,na.rm=T)),maxage=c(max(liste[["male"]]$age,na.rm=T),max(liste[["female"]]$age,na.rm=T)))
})
)})

#' calculates standard deviation scores according to a given reference
#' @param object an \code{\linkS4class{refs}} object 
#' @param ... further arguments to \code{sds}
#' @exportMethod sds

setGeneric(
  name="sds",
  def=function(object,...){standardGeneric("sds")}
  )

#' calculates standard deviation scores according to a given reference 
#' @param object an \code{\linkS4class{refs}} object 
#' @param ... further arguments to \code{sds}
#' @exportMethod sds.df

setGeneric(
  name="sds.df",
  def=function(object,...){standardGeneric("sds.df")}
  )


#' calculates standard deviation scores according to a given reference
#' @param object \code{\linkS4class{refs}} object
#' @param item reference contained in \code{object}
#' @param sex gender
#' @param age age
#' @param value numeric value which should transform into SDscore
#' @rdname sds-methods
#' @examples
#' data(who.ref)
#' sds(who.ref,item="height",sex="female",age=0,value=c(50,51,52,53))
#' sds(who.ref,item="weight",sex="female",age=0,value=c(2.5,3,3.5))
#' sds(who.ref,item="bmi",sex="male",age=c(10,11,12),value=c(14.5,15,19.5))
setMethod("sds","refs",
          definition=function(object,item,sex,age,value){
            m <- approx(object@ref[[item]][[sex]]$age,object@ref[[item]][[sex]]$m,xout=age,rule=1)$y
            l <- approx(object@ref[[item]][[sex]]$age,object@ref[[item]][[sex]]$l,xout=age,rule=1)$y
            s <- approx(object@ref[[item]][[sex]]$age,object@ref[[item]][[sex]]$s,xout=age,rule=1)$y
            ((value/m)**l-1)/(l*s)
          }
          )

#' calculates standard deviation scores according to a given reference and adds the respective column
#' @param object \code{\linkS4class{refs}} object
#' @param df dataframe 
#' @param item reference contained in \code{object}
#' @param sex name of the column containing sex
#' @param age name of the column containing age
#' @param value name of the column containing numeric value which should transform into SDscore
#' @param male coding for male
#' @param female coding for female
#' @rdname sdsdf-methods
#' @examples
#' data(who.ref)
#' x <- data.frame(height=c(50,100,60,54),
#'                 sex=c("m","f","f","m"),
#'                 age=c(0,2.9,0.6,0.2))
#' sds.df(who.ref,df=x,item="height",value="height",sex="sex",age="age",male="m",female="f")
setMethod("sds.df","refs",
          definition=function(object,df,item,value="value",age="age",sex="sex",male=0,female=1){
#           df <- df[!is.na(df[,value]),]
            tmp.list <- split(df,df[,sex])
            tmp.m <- tmp.list[[as.character(male)]]
            if(!is.null(tmp.m)){
              tmp.m$sds <- sds(object,item=item,sex="male",age=tmp.m[,age],value=tmp.m[,value])
            }
            tmp.f <- tmp.list[[as.character(female)]]
            if(!is.null(tmp.f)){
              tmp.f$sds <- sds(object,item=item,sex="female",age=tmp.f[,age],value=tmp.f[,value])
            }
            res <- rbind(tmp.f,tmp.m)
            nc <- ncol(res)
            names(res)[nc] <- paste(item,"sds",object@name,sep="")
            return(res)
          }
          )

#' give references like articles or uris
#' @param object an object
#' @param ... further arguments to \code{articles}
setGeneric(
  name="articles",
  def=function(object,...){standardGeneric("articles")}
  )


#' Gives back references like articles or uris from an refs object, if \code{item} is given only the respective information will be extracted, if omitted information for all items will be returned
#' @param object object of class refs
#' @param item reference contained in \code{object}
#' @param ... further arguments to \code{articles}
#' @rdname articles-methods
setMethod("articles","refs",
          definition=function(object,item=NA){
            if(!is.na(item)){
              text <- object@citation[[item]]
              print(text)
            }else{
              return(object@citation)
            }
          }
          )
