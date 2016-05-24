plotexmap<-function(sp,mt,
xaxt="n",
lift=0.1,leaflift=0.1,ylim=NULL,
leafcolors=NULL
)
{
if (is.null(leafcolors)) lc<-mt$colot
c2s<-colo2scem(sp,mt,lc)

plotvecs(sp$bigvecs,sp$bigdepths,
lift=lift,xaxt="n",
ylim=ylim,
#ylim=c(horilines[length(horilines)],horilines[1]),   #hseq[1]),
leafcolors=c2s,leaflift=leaflift)                        #log="y")

}

