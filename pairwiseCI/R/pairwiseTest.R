"pairwiseTest" <-
function(formula, data, by=NULL, method="t.test",
 control=NULL, ...)

{

#alternative<-match.arg(alternative, choices=c("two.sided", "less", "greater"))

# # # check the arguments

# first argument must be a formula
# check the structure of the formula: length=3

if(class(formula) != "formula")
 {stop(" first argument must be a two.sided formula of the structure 'response ~ treatment'  ")}

if(length(formula) != 3)
 {stop(" argument formula must have the structure 'response ~ treatment' ")}

if(length(formula[[3]]) != 1)
 {stop(" argument formula must have the structure 'response ~ treatment' ")}

# # #

if (all(class(data) != "data.frame"))
 {stop("argument data must be of class 'data.frame' ")}

if( "Prop.test"==method )
 {PWTINT<-"pairwiseTestProp"}
  else
  {PWTINT<-"pairwiseTestCont"}


# # #

args<-list(...)
args$formula <- formula
#args$alternative <- alternative
args$method <- method
args$control <- control

if(any(c("t.test", "ratio.t.test")==method) & is.null(args$var.equal))
 {args$var.equal<- FALSE}


if(is.null(by)==TRUE)
{
args$data <- data
out<-do.call(PWTINT,args)
output <- list(
 byout=list(out),
 bynames="",
# alternative=alternative,
 method=method,
 list(...) )
}


if(!is.null(by))
 {

 # check the by argument:

  if(!is.character(by))
   {stop("argument 'by' should be a character vector specifying names of columns of factors of the data.frame specified in 'data'")}
  else
   {
    for(i in 1:length(by))
     { if(all(names(data)!=by[i]))
        {stop(paste("column", by[i], "can not be found in data"))}
     }
   }

# create a list of by-factors to be used by "split":

  by.list <- list()
  for(i in seq(along.with=by))
   {by.list[[i]]<-as.factor(data[[by[i]]])}

  data.by <- split(data, f=by.list, drop=TRUE)

  bynames <- names(data.by)
 
  byout<-list()

  for(e in 1:length(data.by))
   {args$data <- data.by[[e]]
    byout[[bynames[e]]] <- do.call(PWTINT, args)
   }

 output <- list(
 byout=byout,
 bynames=bynames,
# alternative=alternative,
 method=method,
 control=control,
 by=by,
 list(...) )
 }


class(output) <- "pairwiseTest"
return(output)

}


