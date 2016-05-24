`print.SCSnp` <-
function(x, ...)
{
args<-list(...)

alt<-x$alternative
cl<-x$conf.level
N<-x$N

switch(alt,
two.sided={ALT="Two-sided simultaneous "; BD<-"intervals "},
less={ALT="Upper simultaneous "; BD<-"bounds "},
greater={ALT="Lower simultaneous "; BD<-"bounds "}
)

Text<-paste(ALT, round(cl*100,2), "% credible ", BD, "\n", "based on ",N, " simulation runs","\n", sep="")
cat(Text)

dat<-data.frame(x$estimate, x$conf.int)
names(dat)<-c("Estimate","Lower","Upper")

pargs<-args
pargs$x<-dat

if(is.null(args$digits))
pargs$digits<-4

do.call("print.data.frame", pargs)

invisible(x)
}

