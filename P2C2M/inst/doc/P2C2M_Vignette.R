### R code from vignette source 'P2C2M_Vignette.Rnw'

###################################################
### code chunk number 1: P2C2M_Vignette.Rnw:17-18
###################################################
  options(width=60, size="scriptsize")


###################################################
### code chunk number 2: P2C2M_Vignette.Rnw:43-54
###################################################
require(P2C2M)
require(ggplot2)
require(grid)
data(viz_example_1)
inp = viz_example_1
alpha = 0.05
inData = qnts = df = titles = list()
df$lwr = df$upr = list()
titles = sprintf("gene%02d", c(1:10))
colnames(inp$nReps0) = titles
colnames(inp$nReps10) = titles


###################################################
### code chunk number 3: P2C2M_Vignette.Rnw:64-72
###################################################
inData$nReps10 = stack(as.data.frame(inp$nReps10))
inData$nReps10[,3] = "nReps10"
colnames(inData$nReps10) = c("value", "gene", "n_reps")
qnts$nReps10 = apply(inp$nReps10, 2, quantile, c(alpha, 1-alpha), na.rm=TRUE)
df$lwr$nReps10 = data.frame(lwrQntl=qnts$nReps10[1,], gene=names(qnts$nReps10[1,]),
                            n_reps=rep("nReps10", 10))
df$upr$nReps10 = data.frame(uprQntl=qnts$nReps10[2,], gene=names(qnts$nReps10[2,]),
                            n_reps=rep("nReps10", 10))


###################################################
### code chunk number 4: P2C2M_Vignette.Rnw:82-90
###################################################
inData$nReps0 = stack(as.data.frame(inp$nReps0))
inData$nReps0[,3] = "nReps0"
colnames(inData$nReps0) = c("value", "gene", "n_reps")
qnts$nReps0 = apply(inp$nReps0, 2, quantile, c(alpha, 1-alpha), na.rm=TRUE)
df$lwr$nReps0 = data.frame(lwrQntl=qnts$nReps0[1,], gene=names(qnts$nReps0[1,]),
                           n_reps=rep("nReps0", 10))
df$upr$nReps0 = data.frame(uprQntl=qnts$nReps0[2,], gene=names(qnts$nReps0[2,]),
                           n_reps=rep("nReps0", 10))


###################################################
### code chunk number 5: P2C2M_Vignette.Rnw:99-103
###################################################
inData = rbind(inData$nReps0, inData$nReps10)
dfLwr = rbind(df$lwr$nReps0, df$lwr$nReps10)
dfUpr = rbind(df$upr$nReps0, df$upr$nReps10)
inData$gene = factor(inData$gene, levels = sort(titles))


###################################################
### code chunk number 6: P2C2M_Vignette.Rnw:113-134
###################################################
ggplot(data=inData, aes(x=value)) +
  geom_density() +
  facet_grid(gene ~ n_reps, scales = "free_y") +
  labs(x="Difference values") + 
  ggtitle(expression(atop("Distributions under different N of replicates", 
                     atop(italic("Descriptive Statistic: coal (Liu & Yu 2010)"), "")))) +

  theme_bw() +
  theme(axis.text = element_text(size=5),
        axis.title.y=element_blank(),
        strip.text.y=element_text(angle=0),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        strip.background = element_rect(fill="white"), 
        panel.margin = unit(0.5, "lines"),
        plot.title = element_text(face="bold", size=rel(1.5), vjust=-1)) +
  # Limits on the x-axis improve the visualization
  xlim(-500, 500) + 
  geom_vline(xintercept=0, linetype = "dashed") + 
  geom_vline(aes(xintercept=lwrQntl), dfLwr, color="grey") +
  geom_vline(aes(xintercept=uprQntl), dfUpr, color="grey")


###################################################
### code chunk number 7: P2C2M_Vignette.Rnw:142-163
###################################################
ggplot(data=inData, aes(x=value)) +
  geom_density() +
  facet_grid(gene ~ n_reps, scales = "free_y") +
  labs(x="Difference values") + 
  ggtitle(expression(atop("Distributions under different N of replicates", 
                     atop(italic("Descriptive Statistic: coal (Liu & Yu 2010)"), "")))) +

  theme_bw() +
  theme(axis.text = element_text(size=5),
        axis.title.y=element_blank(),
        strip.text.y=element_text(angle=0),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        strip.background = element_rect(fill="white"), 
        panel.margin = unit(0.5, "lines"),
        plot.title = element_text(face="bold", size=rel(1.5), vjust=-1)) +
  # Limits on the x-axis improve the visualization
  xlim(-500, 500) + 
  geom_vline(xintercept=0, linetype = "dashed") + 
  geom_vline(aes(xintercept=lwrQntl), dfLwr, color="grey") +
  geom_vline(aes(xintercept=uprQntl), dfUpr, color="grey")


###################################################
### code chunk number 8: P2C2M_Vignette.Rnw:184-188
###################################################
require(P2C2M)
require(ggplot2)
data(viz_example_2)
inp = viz_example_2


###################################################
### code chunk number 9: P2C2M_Vignette.Rnw:199-216
###################################################
customfunc = function(inData, simNum){
  handle = inData
  colnames(handle) = c("gtp", "ray", "ndc", "gsi")
  # Convert results into presence/absence matrix
  handle[!grepl("n.s.", handle)] = 1
  handle[grepl("n.s.", handle)] = 0
  # Stack the individual descriptive statistics
  handle = stack(data.frame(handle, stringsAsFactors=FALSE))
  colnames(handle)[1] = "value"
  colnames(handle)[2] = "stat"
  # Add gene identifiers (under the assumption that there are 10 genes)
  handle[,3] = rep(c(1:10), 4)
  colnames(handle)[3] = "gene"
  handle[,4] = simNum
  colnames(handle)[4] = "sim"
return(handle)
}


###################################################
### code chunk number 10: P2C2M_Vignette.Rnw:226-241
###################################################
highL = list()
sims = as.numeric(names(inp$High))
for (i in 1:length(inp$High)) {highL[[i]] = customfunc(inp$High[[i]], sims[i])}
High = do.call("rbind", highL)
High[,ncol(High)+1] = "High_Subst_Rate"
colnames(High)[ncol(High)] = "ratetype"

lowL = list()
sims = as.numeric(names(inp$Low))
for (i in 1:length(inp$Low)) {lowL[[i]] = customfunc(inp$Low[[i]], sims[i])}
Low = do.call("rbind", lowL)
Low[,ncol(Low)+1] = "Low_Subst_Rate"
colnames(Low)[ncol(Low)] = "ratetype"

inData = rbind(High, Low)


###################################################
### code chunk number 11: P2C2M_Vignette.Rnw:250-261
###################################################
ggplot(data=inData, aes(x=sim,y=gene)) + 
  geom_point(aes(colour=value), size = 3) +
  scale_colour_manual(values = c(NA,'black')) + 
  facet_grid(stat~ratetype) + 
  ggtitle(expression(atop("Distribution of False Positive Results",
          atop(italic("Alpha=0.1") , "")))) +
  theme_bw() + 
  theme(axis.text = element_text(size=5),
        strip.background = element_rect(fill="white")) +
  scale_x_discrete(breaks=c(1:5), labels=c(1:5)) +
  scale_y_discrete(breaks=c(10:1), labels=c(10:1))


###################################################
### code chunk number 12: P2C2M_Vignette.Rnw:268-279
###################################################
ggplot(data=inData, aes(x=sim,y=gene)) + 
  geom_point(aes(colour=value), size = 3) +
  scale_colour_manual(values = c(NA,'black')) + 
  facet_grid(stat~ratetype) + 
  ggtitle(expression(atop("Distribution of False Positive Results",
          atop(italic("Alpha=0.1") , "")))) +
  theme_bw() + 
  theme(axis.text = element_text(size=5),
        strip.background = element_rect(fill="white")) +
  scale_x_discrete(breaks=c(1:5), labels=c(1:5)) +
  scale_y_discrete(breaks=c(10:1), labels=c(10:1))


