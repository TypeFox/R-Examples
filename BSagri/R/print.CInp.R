`print.CInp` <-
function(x, ...)
{
args<-list(...)

alt<-x$alternative
cl<-x$conf.level
N<-x$N

switch(alt,
two.sided={ALT="Two-sided local "; BD<-"intervals "},
less={ALT="Upper local "; BD<-"bounds "},
greater={ALT="Lower local "; BD<-"bounds "}
)

Text<-paste(ALT, round(cl*100,2), "% confidence ", BD, "\n", "based on ",N, " simulation runs","\n", sep="")
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

