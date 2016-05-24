`summary.pairwiseCI` <-
function( object, digits=4, ...)
{
byout <- object$byout
bynames <- object$bynames
method <- object$method

if(length(bynames)==1){METHOD<-object[[1]]$method}else{METHOD<-object$byout[[1]]$method}

if(length(bynames)==1)
 {
   estimate <- as.data.frame(matrix(round( byout[[1]]$estimate, digits=digits), ncol=1))
   conf.int <- as.data.frame(round( cbind(byout[[1]]$lower, byout[[1]]$upper), digits=digits))

   colnames(estimate) <- c("estimate")
   rownames(estimate) <- byout[[1]]$compnames

   colnames(conf.int) <- c("lower", "upper")
   rownames(conf.int) <- byout[[1]]$compnames

   out<-list(estimate=estimate, conf.int=conf.int)
 }
else
 {
  
  if(length(byout) != length(bynames))
   {stop("INTERNAL: bynames and byout of different length!! ")}

   out<-list()
   for (i in 1:length(byout))
    {
     estimate <- as.data.frame(matrix(round( byout[[i]]$estimate, digits=digits), ncol=1))
     conf.int <- as.data.frame(round( cbind(byout[[i]]$lower, byout[[i]]$upper), digits=digits))

     colnames(estimate) <- c("estimate")
     rownames(estimate) <- byout[[i]]$compnames

     colnames(conf.int) <- c("lower", "upper")
     rownames(conf.int) <- byout[[i]]$compnames

     out[[i]]<-list(estimate=estimate, conf.int=conf.int)   
    }
   names(out)<-bynames
 }

attr(out, "methodname") <- METHOD
attr(out, "conf.level") <- object$conf.level

class(out)<-"summary.pairwiseCI"
return(out)
}





#####################################



`as.data.frame.pairwiseCI` <-
function( x, ...)
{
byout <- x$byout
bynames <- x$bynames
method <- x$method

METHOD<-byout[[1]]$method

if(length(bynames)==1)
 {
   dat<-data.frame(byout[[1]]$estimate, byout[[1]]$lower, byout[[1]]$upper, byout[[1]]$compnames)

   colnames(dat)<-c("estimate", "lower", "upper", "comparison")
   rownames(dat)<-NULL
 }
else
 {
  
  if(length(byout) != length(bynames))
   {stop("INTERNAL: bynames and byout of different length!! ")}

   bynames<-x$bynames

   dat<-data.frame()
   for (i in 1:length(byout))
    {
     estimate <- byout[[i]]$estimate
     lower <- byout[[i]]$lower
     upper <- byout[[i]]$upper
     compnames <- byout[[i]]$compnames

     by<-rep( bynames[i], times=length(estimate) )

     dati<-data.frame( estimate, lower, upper, compnames, by)

     dat<-rbind(dat, dati)  
    }
   colnames(dat)<-c("estimate","lower","upper", "comparison", "by")
   rownames(dat)<-NULL
 }

attr(dat, "methodname") <- METHOD
attr(dat, "conf.level") <- x$conf.level

class(dat)<-"data.frame"
return(dat)
}


###############################





