plots <-
function ()
{
Class<-NULL
loadNamespace("ggplot2")
load("permubiome.RData")
a <- array(, nrow(df_norm))
for (j in 1:nrow(df_norm)){
a[j]<-sum(df_norm[j,3:ncol(df_norm)])
}
for (l in 3:ncol(df_norm)){
for (m in 1:nrow(df_norm)){
df_norm[m,l]<-round((df_norm[m,l]/a[m]), digits=6)
}
}
category <- readline("Type the category you want plotting : ")
if (category == ""){
category <- colnames(df_norm[3])
print(paste("As you declare no categories, the first one of your dataset is plotted!"))
}
else
{
ggplot(df_norm, aes(Class, df_norm[,category]), environment = environment())+geom_boxplot(notch=F, outlier.colour="blue", outlier.shape=1, outlier.size=3, binaxis="y", stackdir="center", dotsize=3)+ggtitle(category)+ylab("Normalized read proportion")+xlab("Classes")+theme(axis.text=element_text(size=12), axis.title=element_text(size=16,face="bold"))+geom_jitter(position=position_jitter(width=0, height=0))
output <- readline("Do you want an output file (yes/no)? : ")
if (substr(output, 1, 1) == "y"){
extension <- readline("What extension do you prefer fo the output plot (ps, pdf, jpeg, tiff, png, bmp )? : ")
ggsave(filename = paste(category, extension, sep = "."), plot = last_plot(), path = NULL, scale = 1, units = c("cm"), dpi = 300, limitsize = TRUE)
} 
else
{
last_plot()
}
}
}
