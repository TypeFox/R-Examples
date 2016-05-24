require(multilevelPSA)
require(party)
require(ggplot2)
require(mice)
data(pisana)
data(pisa.psa.cols)

table(pisana$CNT, pisana$PUBPRIV, useNA='ifany')

##### Canada
can <- pisana[pisana$CNT == 'CAN',]
can$treat <- can$PUBPRIV == 'Private'
can <- can[,c('treat',pisa.psa.cols)]

# Classification Trees
psranges.can.tree <- psrange(can, can$treat, treat ~ ., type='ctree', nboot=20)
plot(psranges.can.tree, labels=c('Public', 'Private')) + ggtitle('Classification Trees')
#ggsave('~/Dropbox/School/Dissertation/Figures/PSRanges-tree-PISA.pdf', width=9, height=12)

# Logistic Regression
psranges.can.lr <- psrange(can, can$treat, treat ~ ., type='logistic', nboot=20)
plot(psranges.can.lr, labels=c('Public', 'Private')) + ggtitle('Logistic Regression')
#ggsave('~/Dropbox/School/Dissertation/Figures/PSRanges-lr-PISA.pdf', width=9, height=12)

