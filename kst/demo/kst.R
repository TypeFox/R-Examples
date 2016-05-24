### create input
input <- set(set("a"), set("a","b"), set("a","c"), set("d","e"), 
   set("a","b","d","e"), set("a","c","d","e"), set("a","b","c","d","e"))
input

### convert to knowledge structure
kst <- kstructure(input)
kst

### compute domain of knowledge structure
kdomain(kst)

### compute notions of knowledge structure
knotions(kst)

### compute atoms of knowledge structure
katoms(kst, items=set("a","b","c"))

### compute trace of knowledge structure
ktrace(kst, items=set("c","d","e"))

### is knowledge structure well-graded?
kstructure_is_wellgraded(kst)

### compute inner and outer fringe of whole knowledge structure
kfringe(kst, state=NULL, fringe="inner")
kfringe(kst, state=NULL, fringe="outer")

### compute inner and outer fringe of a particular knowledge state
kfringe(kst, state=set("a", "b"), fringe="inner")
kfringe(kst, state=set("a", "b"), fringe="outer")

### convert to relation
as.relation(kst)

### assess individuals
rp <- data.frame(a=c(1,1,0,1,1,1,1,0,0,0),b=c(0,1,0,1,0,1,0,1,0,0),
   c=c(0,0,0,0,1,1,1,0,1,0),d=c(0,0,1,1,1,1,0,0,0,1), e=c(0,0,1,1,1,1,0,0,0,0))
kassess(kst, rpatterns=rp)

### validate knowledge structure
kvalidate(kst, rpatterns=rp, method="gamma")
kvalidate(kst, rpatterns=rp, method="percent")
kvalidate(kst, rpatterns=rp, method="VC")
kvalidate(kst, rpatterns=rp, method="DA")

### compute closure under union
closure(kst, operation="union")

### compute discriminative reduction
reduction(kst, operation="discrimination")

### is structure a space?
kstructure_is_kspace(kst)

### convert to knowledge space
ksp <- kspace(kst)
ksp

### is structure a space?
kstructure_is_kspace(ksp)

### base of knowledge space
kbase(ksp)

### compute learning paths in knowledge structure
lp <- lpath(kst)
lp

### are learning paths gradations?
lpath_is_gradation(lp)
