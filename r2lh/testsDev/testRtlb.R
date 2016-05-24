
cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++                       Begin testRtlb                    +++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
source("../R/rtlbDisplay.R")
source("../R/rtlb.R")



cat("####################
#   Latex output   #
####################\n")

### Styles de display
ds <- list(lim=4,wide=c("f1","o2"),long=c("f2","y2"))
ds2 <- list(wide=c("f1","o2"),long=c("f2","y2"))

### Logical

rtlb(f1~f1b)
rtlb(f1~f2)
rtlb(f1~f3)

rtlb(f1~o1)
rtlb(f1~o2)
rtlb(f1~o3)

rtlb(f1~i1)
rtlb(f1~i2)
rtlb(f1~i3)
rtlb(f1~i2,limDiscreteX=3)

rtlb(f1~n2)
rtlb(f1~n3)

rtlb(f1b~df,displayStyle=ds)


### Factor
rtlb(f2~f1)
rtlb(f2~f2)
rtlb(f2~f3)

rtlb(f2~o1)
rtlb(f2~o2)
rtlb(f2~o3)

rtlb(f2~i1)
rtlb(f2~i2)
rtlb(f2~i3)
rtlb(f2~i2,limDiscreteX=3)
rtlb(f2~i2,displayStyle=3)

rtlb(f2~n2)
rtlb(f2~n3)

rtlb(f2~df)


### Ordered
rtlb(o2~f1)
rtlb(o2~f2)
rtlb(o2~f3)

rtlb(o2~o1)
rtlb(o2~o2)
rtlb(o2~o3)

rtlb(o2~i1)
rtlb(o2~i2)
rtlb(o2~i3)
rtlb(o2~i2,limDiscreteX=3)
rtlb(o2~i2,displayStyle=3)

rtlb(o2~n2)
rtlb(o2~n3)

rtlb(o2~df)


### Discrete
rtlb(i2~f1)
rtlb(i2~f2)
rtlb(i2~f3)

rtlb(i2~o1)
rtlb(i2~o2)
rtlb(i2~o3)

rtlb(i2~i1)
rtlb(i2~i3)
rtlb(i2~i3,limDiscreteX=3)
rtlb(i2~i3,displayStyle=3)

rtlb(i2~n2)
rtlb(i2~n3)

rtlb(i2~df)


### Continuous
rtlb(n1~f1)
rtlb(n1~f2)
rtlb(n1~f3)

rtlb(n1~o1)
rtlb(n1~o2)
rtlb(n1~o3)

rtlb(n1~i1)
rtlb(n1~i2)
rtlb(n1~i3)
rtlb(n1~i2,limDiscreteX=3)
rtlb(n1~i2,displayStyle=3)

rtlb(n1~n2)
rtlb(n1~n3)

rtlb(n1~df)

cat("---------------------------------------------------------------
---                       Begin testRtlb                    ---
---------------------------------------------------------------\n")
