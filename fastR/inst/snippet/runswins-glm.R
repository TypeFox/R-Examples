bb <- mlb2004
bb$runmargin = (bb$R - bb$OR) / (bb$G)

# data frame has summarized data for each team, so different syntax here:
glm(cbind(W,L)~runmargin,data=bb,family='binomial') -> glm.bb
###hop:3-9
summary(glm.bb)
