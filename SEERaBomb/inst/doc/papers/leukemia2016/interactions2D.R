# interactions2D.R (Fig S2: show NHL and liver need an interaction term in the model)
library(SEERaBomb)
load("~/data/SEER/mrgd/cancDef.RData")
levels(canc$cancer)
load("~/data/SEER/mrgd/popsa.RData")
zap=c("cervixCIS","benign")
# mk2D cannot handle cervixCIS because it ends abruptly in 1995. Benigns 
# starting in 2004 are also a problem. The problem is that a plane of zeros ends
# up with too much of a smooth landing from the hills, i.e. tall cliffs aren't
# well approximated by splines. So for CMML and MDS special handing of smaller
# PY matrices takes place, but such code is ugly and hard to think about, so we
# don't want to clutter SEERaBomb source with additional special cases unless
# the second cancers are very important to us. So cervicalCIS and benigns are
# being left out. Note that benigns cause an additional problem, which is that
# their seqnums go 60 to 88. Handling this would be additional work that
# probably isn't worth it from a heme cancer centric perspective.
canc=canc%>%filter(!cancer%in%zap)
canc$cancer=factor(canc$cancer)
pm=seerSet(canc,popsa,Sex="male")
pm=mk2D(pm,txt="interaction") # pooled race males, all male cancers not zapped get fitted
plot2D(pm)
