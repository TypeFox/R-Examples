bb <- mlb2004
bb$runmargin = (bb$R - bb$OR) / (bb$G)
glm(cbind(W,L)~runmargin,data=bb,family='binomial') -> glm.bb
bb$winP <- bb$W/bb$G
bb$glmPredWinP <- predict(glm.bb,
    newdata=data.frame(runmargin=bb$runmargin),type='response')
lm(winP~runmargin,data=bb) -> lm.bb
summary(lm.bb)
bb$lmPredWinP <- predict(lm.bb,
    newdata=data.frame(runmargin=bb$runmargin),type='response')
# observations 8 and 27 have largest residuals
bb[c(8,27,1:2,29:30),c("Team","winP","glmPredWinP","lmPredWinP")]
