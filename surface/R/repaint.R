repaint <-
function(otree,regshifts,stem=TRUE){
	subtrees<-branches<-regshifts[order(otree@times[match(names(regshifts),otree@nodes)])]
	if(stem==FALSE)branches<-branches[which(as.numeric(names(branches))>(otree@nnodes-otree@nterm))]
	regs<-paint(otree,subtree=subtrees,branch=branches)
	regs[1]<-regshifts[1]
	regs<-factor(regs)
	regs
	}
