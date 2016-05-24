obs.v.exp <-
function(mod){
  ptab <- mod$results
  
  ptab$signs <- ifelse(ptab$p_gt>=0.05,0,1) + ifelse(ptab$p_lt>=0.05,0,-1)
  
  exp_cooccur <- ptab$exp_cooccur
  obs_cooccur <- ptab$obs_cooccur
  signs <- ptab$signs
  
  p <- ggplot(ptab, aes(x=exp_cooccur, y=obs_cooccur)) + geom_point(aes(fill=factor(signs,levels=c(-1,0,1))), colour="black",pch=21, size=5)
  p <- p + scale_fill_manual(values = c("#FFCC66","dark gray","light blue"), name = "", labels = c("negative","random","positive"),drop=FALSE) 
  p <- p + theme(plot.title = element_text(vjust=2,size=20, face="bold"),legend.text=element_text(size=18),axis.title = element_text(size = 20),axis.text=element_text(size=18),axis.text.x=element_text(hjust=0,vjust=1)) + xlab("Expected Co-occurrences") + ylab("Observed Co-occurrences")
  p <- p + ggtitle("Observed-Expected Plot") + geom_abline(color="dark gray")
  p
}
