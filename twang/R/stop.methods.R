
stop.methods <- matrix(c("es.mean","ks.mean","es.max","ks.max","ks.max.direct","es.max.direct"), nr = 1)
names(stop.methods) <- c("es.stat.mean","ks.stat.mean", "es.stat.max", "ks.stat.max",
"ks.max.direct", "es.max.direct")


ks.mean <- list(metric=ksStat,
                rule.summary=mean,
                direct=FALSE,
                estimand = NULL,
                na.action="level",
                name=NULL)
                
class(ks.mean) <- "stop.method"
                
es.mean <- list(metric=esStat,
                rule.summary=mean,
                direct=FALSE,
                estimand = NULL,
                na.action="level",
                name=NULL)
                
class(es.mean) <- "stop.method"
                
ks.max <- list(metric=ksStat,
				rule.summary=max,
				direct=FALSE,
				estimand = NULL,
				na.action="level",
				name=NULL)
				
class(ks.max) <- "stop.method"
				
es.max <- list(metric=esStat,
				rule.summary=max,
				direct=FALSE,
				estimand = NULL,
				na.action="level",
				name=NULL)
				
class(es.max) <- "stop.method"
				
ks.max.direct <- list(metric=ksStat,
				rule.summary=max,
				direct=TRUE,
				estimand= NULL,
				na.action="level",
				name=NULL)
				
class(ks.max.direct) <- "stop.method"
				
es.max.direct <- list(metric=esStat,
				rule.summary=max,
				direct=TRUE,
				estimand= NULL,
				na.action="level",
				name=NULL)
				
class(es.max.direct) <- "stop.method"