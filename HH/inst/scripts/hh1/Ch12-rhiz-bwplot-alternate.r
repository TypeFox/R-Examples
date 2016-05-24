data(rhiz.alfalfa)
data(rhiz.clover)

## horizontal boxplots conditioned on strain
if.R(r={
## rhiz-bwplot.r       # rhiz, (Fig 11.2 11.3 11.4) alternate
## horizontal boxplots conditioned on strain

print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(comb ~ Npg | strain, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = TRUE,  # right
bwplot(comb ~ Npg | strain, data=rhiz.clover,
       main="clover", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Nitrogen per gram", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhnpg.eps"))

print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(comb ~ nitro | strain, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = TRUE,  # right
bwplot(comb ~ nitro | strain, data=rhiz.clover,
       main="clover", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Nitrogen", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhn.eps"))

print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(comb ~ weight | strain, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = TRUE,  # right
bwplot(comb ~ weight | strain, data=rhiz.clover,
       main="clover", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Weight", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhw.eps"))
}, s={
## rhiz-bwplot.s       # rhiz, (Fig 11.2 11.3 11.4) alternate
## horizontal boxplots conditioned on strain

print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(comb ~ Npg | strain, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = FALSE,  # right
bwplot(comb ~ Npg | strain, data=rhiz.clover,
       main="clover", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
mtext("Nitrogen per gram", outer=TRUE, side=3)
## export.eps(hh("tway/figure/rhnpg.eps"))

print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(comb ~ nitro | strain, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = FALSE,  # right
bwplot(comb ~ nitro | strain, data=rhiz.clover,
       main="clover", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
mtext("Nitrogen", outer=TRUE, side=3)
## export.eps(hh("tway/figure/rhn.eps"))

print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(comb ~ weight | strain, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = FALSE,  # right
bwplot(comb ~ weight | strain, data=rhiz.clover,
       main="clover", layout=c(1,6),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
mtext("Weight", outer=TRUE, side=3)
## export.eps(hh("tway/figure/rhw.eps"))
})




if.R(r={
## rhiz-bwplot.t.r     # rhiz, (Fig 11.2 11.3 11.4) alternate
## vertical boxplots conditioned on strain

print(position = c(0, .47, 1, .94), more = TRUE,  # top
bwplot(Npg ~ comb | strain, data=rhiz.alfalfa,
         main="alfalfa", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, x=list(at=c(1,2),labels=c("alf","alf+clov"))))
)
print(position = c(0, 0, 1, .47), more = TRUE,  # bottom
bwplot(Npg ~ comb | strain, data=rhiz.clover,
         main="clover", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, x=list(at=c(1,2),labels=c("alf","alf+clov"))))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Nitrogen per gram", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhnpg.t.eps"))


print(position = c(0, .47, 1, .94), more = TRUE,  # top
bwplot(nitro ~ comb | strain, data=rhiz.alfalfa,
         main="alfalfa", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, x=list(at=c(1,2),labels=c("alf","alf+clov"))))
)
print(position = c(0, 0, 1, .47), more = TRUE,  # bottom
bwplot(nitro ~ comb | strain, data=rhiz.clover,
         main="clover", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, x=list(at=c(1,2),labels=c("alf","alf+clov"))))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Nitrogen", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhn.t.eps"))


print(position = c(0, .47, 1, .94), more = TRUE,  # top
bwplot(weight ~ comb | strain, data=rhiz.alfalfa,
         main="alfalfa", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, x=list(at=c(1,2),labels=c("alf","alf+clov"))))
)
print(position = c(0, 0, 1, .47), more = TRUE,  # bottom
bwplot(weight ~ comb | strain, data=rhiz.clover,
         main="clover", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, x=list(at=c(1,2),labels=c("alf","alf+clov"))))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Weight", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhw.t.eps"))
}, s={
## rhiz-bwplot.t.s     # rhiz, (Fig 11.2 11.3 11.4) alternate
## vertical boxplots conditioned on strain


print(split = c(1,2,1,2), more = TRUE,  # top
t(bwplot(comb ~ Npg | strain, data=rhiz.alfalfa,
         main="alfalfa", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, y=list(at=c(1,2),labels=c("alf","alf+clov"))))
))
print(split = c(1,1,1,2), more = FALSE,  # bottom
t(bwplot(comb ~ Npg | strain, data=rhiz.clover,
         main="clover", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, y=list(at=c(1,2),labels=c("alf","alf+clov"))))
))
mtext("Nitrogen per gram", outer=TRUE, side=3)
## export.eps(hh("tway/figure/rhnpg.t.eps"))


print(split = c(1,2,1,2), more = TRUE,  # top
t(bwplot(comb ~ nitro | strain, data=rhiz.alfalfa,
         main="alfalfa", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, y=list(at=c(1,2),labels=c("alf","alf+clov"))))
))
print(split = c(1,1,1,2), more = FALSE,  # bottom
t(bwplot(comb ~ nitro | strain, data=rhiz.clover,
         main="clover", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, y=list(at=c(1,2),labels=c("alf","alf+clov"))))
))
mtext("Nitrogen", outer=TRUE, side=3)
## export.eps(hh("tway/figure/rhn.t.eps"))


print(split = c(1,2,1,2), more = TRUE,  # top
t(bwplot(comb ~ weight | strain, data=rhiz.alfalfa,
         main="alfalfa", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, y=list(at=c(1,2),labels=c("alf","alf+clov"))))
))
print(split = c(1,1,1,2), more = FALSE,  # bottom
t(bwplot(comb ~ weight | strain, data=rhiz.clover,
         main="clover", layout=c(6,1),
         par.strip.text=list(cex=.9),
         scales=list(cex=.9, y=list(at=c(1,2),labels=c("alf","alf+clov"))))
))
mtext("Weight", outer=TRUE, side=3)
## export.eps(hh("tway/figure/rhw.t.eps"))
})




if.R(r={
## rhiz-bwplot.i.r     # rhiz, (Fig 11.2 11.3 11.4) alternate
## horizontal boxplots conditioned on combination

print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(strain ~ Npg | comb, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = TRUE,  # right
bwplot(strain ~ Npg | comb, data=rhiz.clover,
       main="clover", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Nitrogen per gram", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhnpg.i.eps"))


print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(strain ~ nitro | comb, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = TRUE,  # right
bwplot(strain ~ nitro | comb, data=rhiz.clover,
       main="clover", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Nitrogen", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhn.i.eps"))


print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(strain ~ weight | comb, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = TRUE,  # right
bwplot(strain ~ weight | comb, data=rhiz.clover,
       main="clover", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(position=c(0, .95, 1, 1), more=FALSE,
      xyplot(0 ~ 0, panel=function(...){},
             main=list("Weight", cex=1.2),
             ylab=expression(""),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
## export.eps(hh("tway/figure/rhw.i.eps"))
}, s={
## rhiz-bwplot.i.s     # rhiz, (Fig 11.2 11.3 11.4) alternate
## horizontal boxplots conditioned on combination

print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(strain ~ Npg | comb, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = FALSE,  # right
bwplot(strain ~ Npg | comb, data=rhiz.clover,
       main="clover", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
mtext("Nitrogen per gram", outer=TRUE, side=3)
## export.eps(hh("tway/figure/rhnpg.i.eps"))


print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(strain ~ nitro | comb, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = FALSE,  # right
bwplot(strain ~ nitro | comb, data=rhiz.clover,
       main="clover", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
mtext("Nitrogen", outer=TRUE, side=3)
## export.eps(hh("tway/figure/rhn.i.eps"))


print(split = c(1,1,2,1), more = TRUE,  # left
bwplot(strain ~ weight | comb, data=rhiz.alfalfa,
       main="alfalfa", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
print(split = c(2,1,2,1), more = FALSE,  # right
bwplot(strain ~ weight | comb, data=rhiz.clover,
       main="clover", layout=c(1,2),
       par.strip.text=list(cex=.9),
       scales=list(cex=.9))
)
mtext("Weight", outer=TRUE, side=3)
## export.eps(hh("tway/figure/rhw.i.eps"))
})
