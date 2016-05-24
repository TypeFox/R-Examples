

#--------------
#seriation framework for dendser
#--------------


seriation_method_dist_dser <- function(x, control) {
   
    if(!is.null(control$hclust)) 
        h<- control$hclust
    
    else if(is.null(control$method)) h<- hclust(x) 
    else h<- hclust(x, control$method)
    
     foo<- as.list(formals(DendSer)  )
     foo <- modifyList(foo, control)
     foo$h<-h
     if (is.null(control$ser_weight)) foo$ser_weight<- x   
    do.call("DendSer", foo)
    }
    
seriation::set_seriation_method("dist", "Dendser", seriation_method_dist_dser,
"Dendrogram seriation")

seriation_method_leafsort <- function(x, control) {
   
    h<- control$hclust
    if (is.null(h)) stop("'control' must have a component hclust")
    
     
     foo<- as.list(formals(DendSer)  )
     foo <- modifyList(foo, control)
     foo$h<-h
     foo$cost <- costLS
     if (is.null(control$ser_weight)) foo$ser_weight<- x[,1]   
    do.call("DendSer", foo)
    }

seriation::set_seriation_method("matrix", "Dendser", seriation_method_leafsort,
"Dendrogram seriation with leaf sort")


#--------------
# examples
#--------------



   
# all.equal(get_order(seriate(d,method="Dendser", control=list(hclust=h,cost=costARc))),
# DendSer(h,d,cost=costARc))



# all.equal(get_order(seriate(d,method="Dendser", control=list(hclust=h,cost=costLS,ser_weight=1:15))),
# DendSer(h,1:15,cost=costLS))

#--------------
# Seriation criteria as cost functions for dendser
#--------------


crit2cost <- function(crit){
	merit <- seriation::get_criterion_method("dist", crit)$merit
	function(dm,o,...){
	  if (merit)
		- seriation::criterion(as.dist(dm),o, crit)
	else seriation::criterion(as.dist(dm),o,crit)		
	}
}


#--------------
# examples
#--------------

# DendSer(h,d,cost=function(x,o,...) criterion(as.dist(x),o,method="AR_deviations"))
# DendSer(h,d,cost=crit2cost("AR_deviations"))
# DendSer(h,d,cost=crit2cost("ME"))



#--------------
# DendSer cost as seriation criteria
#--------------

criterion_method_dist_ARc <- function(x, order, ...) {
	x <- as.matrix(x)
	if (is.null(order)) order <- 1:nrow(x)
    costARc(x,order,...)
}


seriation::set_criterion_method("dist", "ARc", criterion_method_dist_ARc, 
    "AR cost", FALSE)

criterion_method_dist_BAR <- function(x, order, ...) {
	x <- as.matrix(x)
	if (is.null(order)) order <- 1:nrow(x)
    costBAR(x,order,...)
}


seriation::set_criterion_method("dist", "BAR", criterion_method_dist_BAR, 
    "Banded AR cost", FALSE)

criterion_method_dist_LPL <- function(x, order, ...) {
	x <- as.matrix(x)
	if (is.null(order)) order <- 1:nrow(x)
    costLPL(x,order,...)
}


seriation::set_criterion_method("dist", "LPL", criterion_method_dist_LPL, 
    "Lazy path cost", FALSE)

#--------------
# Examples
#--------------

# show_criterion_methods("dist")

# criterion(d)
