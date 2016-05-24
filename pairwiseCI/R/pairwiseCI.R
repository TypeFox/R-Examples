"pairwiseCI" <-
function(formula, data, by=NULL, 
alternative="two.sided", conf.level=0.95,
 method="Param.diff", control=NULL,...)

{

# # # # #

method<-match.arg(method, choices=c("Param.diff", "Param.ratio", "Lognorm.diff", "Lognorm.ratio","HL.diff", "HL.ratio",
 "Median.diff" ,"Median.ratio", "Prop.diff", "Prop.ratio", "Prop.or", "Negbin.ratio", "Poisson.ratio", "Quasipoisson.ratio", "Quasibin.ratio", "ODbin.ratio", "Betabin.ratio", "np.re"))


alternative<-match.arg(alternative, choices=c("two.sided", "less", "greater"))

# # # check the arguments

# first argument must be a formula
# check the structure of the formula: length=3

if(class(formula) != "formula")
 {stop(" first argument must be a two.sided formula of the structure 'response ~ treatment'  ")}

if(length(formula[[3]]) != 1)
 {stop(" argument formula must have the structure 'response ~ treatment' ")}

if (all(class(data) != "data.frame"))
 {stop("argument data must be of class 'data.frame' ")}


lhs<-formula[[2]]

if(method %in% c("Prop.diff", "Prop.ratio", "Prop.or", "Quasibin.ratio", "ODbin.ratio", "Betabin.ratio"))
 {
 if(as.character(formula[[2]])[1]!="cbind")
  {stop("For methods for proportions please specify formula with structure: cbind(success, failure) ~ treatment")}

 if(length(as.character(formula[[2]]))!=3)
  {stop("For methods for proportions please specify formula with structure: cbind(success, failure) ~ treatment")}

 clhs<-as.character(lhs)[c(2,3)]

 if(!any(c(is.integer(data[[clhs[1]]]), is.numeric(data[[clhs[1]]]) )) )
  {stop(paste("Response variable", clhs[1], "is not an integer variable"))}

 if(!any(c(is.integer(data[[clhs[2]]]), is.numeric(data[[clhs[2]]]) )) )
  {stop(paste("Response variable", clhs[2], "is not an integer variable"))}
 
 }

if(method %in% c("Param.diff","Param.ratio","Lognorm.diff","Lognorm.ratio", "HL.diff","HL.ratio", "Median.diff","Median.ratio","Negbin.ratio", "Poisson.ratio", "Quasipoisson.ratio", "np.re"))
 {
 if(class(formula[[2]])[1]!="name")
  {stop("For the specified method, please specify formula with structure: response ~ treatment")}
 }


# # #

if(conf.level <= 0 || conf.level >= 1 || length(conf.level) > 1 || is.numeric(conf.level) == FALSE )
 {stop("conf.level must be a single number between 0 and 1")}


args<-list(...)
args$formula <- formula
args$alternative <- alternative
args$conf.level <- conf.level
args$method <- method
args$control <- control

if(any(c("Prop.diff", "Prop.or", "Prop.ratio", "Quasibin.ratio", "ODbin.ratio", "Betabin.ratio")==method))
 {PWCINT<-"pairwiseCIProp"}
 else
 {PWCINT<-"pairwiseCICont"}


if(any(c("Param.diff", "Param.ratio")==method) & is.null(args$var.equal))
 {args$var.equal<- FALSE}


if(is.null(by)==TRUE)
{
args$data <- data
out<-do.call(PWCINT,args)

output <- list(byout=list(out), bynames="", alternative=alternative, conf.level=conf.level, method=method )

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

  by.list <- list()
  for(i in seq(along.with=by))
   {by.list[[i]]<-data[[by[i]]]}

  data.by <- split(data, f=by.list, drop=TRUE)

  bynames <- names(data.by)
 
  byout<-list()

  for(e in 1:length(data.by))
   {args$data <- data.by[[e]]
    byout[[e]] <- do.call(PWCINT, args)
   }
 names(byout)<-bynames
 output <- list(byout=byout, bynames=bynames, alternative=alternative, conf.level=conf.level, method=method, var.equal=args$var.equal )
 }


class(output) <- "pairwiseCI"
return(output)

}


