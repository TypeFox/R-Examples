plot.simmr_input <-
function(x,tracers=c(1,2),title='Tracers plot',xlab=colnames(x$mixtures)[tracers[1]],ylab=colnames(x$mixtures)[tracers[2]], sigmas=1, group=1, mix_name = 'Mixtures', colour=TRUE,...) {

# Get mixtures to match current group(s)
curr_rows = which(x$group%in%group)  
curr_mix = x$mixtures[curr_rows,,drop=FALSE]
curr_n_groups = length(group)

# Throw error if too many groups (can only handle max 6 before it runs out of shapes)
#if((length(group)+x$n_sources)>6) stop("Too many groups specified. Total number of groups plus number of sources cannot exceed 6")

# Throw error if plotting only one isotope
if(ncol(curr_mix)==1) stop("This function only works for two or more tracers")

# First get the mean corrected sources and the sd corrected sources
source_means_c = x$source_means + x$correction_means
source_sds_c = sqrt(x$source_sds^2 + x$correction_sds^2)

# Set up data frame for ggplot - have to do it this stupid way because of cran
x2=c(source_means_c[,tracers[1]],curr_mix[,tracers[1]])
y=c(source_means_c[,tracers[2]],curr_mix[,tracers[2]])
x_lower=c(source_means_c[,tracers[1]]-sigmas*source_sds_c[,tracers[1]],curr_mix[,tracers[1]])
x_upper=c(source_means_c[,tracers[1]]+sigmas*source_sds_c[,tracers[1]],curr_mix[,tracers[1]])
y_lower=c(source_means_c[,tracers[2]]-sigmas*source_sds_c[,tracers[2]],curr_mix[,tracers[2]])
y_upper=c(source_means_c[,tracers[2]]+sigmas*source_sds_c[,tracers[2]],curr_mix[,tracers[2]])
if(x$n_groups==1) {
  Source=factor(c(x$source_names,rep(mix_name,nrow(curr_mix))),levels=c(mix_name,x$source_names))
} else {
  Source=factor(c(x$source_names,paste(mix_name,'grp',x$group[curr_rows])),levels=c(paste(mix_name,'grp',unique(x$group[curr_rows])),x$source_names))
}
size=c(rep(0.5,x$n_sources),rep(0.5,nrow(curr_mix)))
df=data.frame(x=x2,y=y,x_lower,y_lower,x_upper,y_upper,Source,size)

if(colour) {
  g=ggplot(data=df, aes(x = x,y = y,colour=Source)) + 
    scale_color_viridis(discrete=TRUE) + 
    theme_bw() +
    labs(x=xlab,y=ylab,title=title) +
    geom_errorbarh(aes(xmax=x_upper,xmin=x_lower,height=0)) +
    geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
    scale_shape_manual(values=1:nlevels(df$Source)) +
    theme(legend.title=element_blank(),legend.key = element_blank()) +
    guides(color=guide_legend(override.aes=list(linetype=c(rep(0,curr_n_groups),rep(1,x$n_sources))))) 
} else {
  g=ggplot(data=df, aes(x = x,y = y,colour=Source)) + 
    theme_bw() +
    labs(x=xlab,y=ylab,title=title) +
    geom_errorbarh(aes(xmax=x_upper,xmin=x_lower,height=0)) +
    geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
    scale_shape_manual(values=1:nlevels(df$Source)) +
    theme(legend.title=element_blank(),legend.key = element_blank()) +
    guides(color=guide_legend(override.aes=list(linetype=c(rep(0,curr_n_groups),rep(1,x$n_sources))))) + 
    scale_colour_grey()
}

print(g)

}
