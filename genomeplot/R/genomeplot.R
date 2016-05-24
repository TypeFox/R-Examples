genomeplot <-
function(data,select){

common_significant <- data #Load file

chrNum <- length(unique(common_significant[,2]))##Determine the chromosome number
common_significant[,4]<--log10(common_significant[,4])
#Position adjustment
for (i in 1:chrNum){
  ndx <- which(common_significant[,2]==i)       #The number of rows in the same chromosomes
  lstcgId <- max(common_significant[ndx,3])      #Maximum position on the same chromosome
  common_significant[ndx,3]<-(common_significant[ndx,3]/lstcgId)+(i-1)
}
#####Find the middle location of each chromosome
bpMidVec <- vector(length=chrNum)
 
for (i in 1:chrNum){
   ndx <- which(common_significant[,2]==i) 
   posSub <- common_significant[ndx,3]
   bpMidVec[i] <- (max(posSub)+min(posSub))/2
   common_significant[ndx,2]<- paste("chr",i,sep="")
}
#Calculate the maximum value and the average value and their related information
ndx <- which(common_significant[,4]==max(common_significant[,4]))
max_mark=common_significant[ndx,]
max_logp=common_significant[ndx,4]
max_scientific=format(max_mark[4],scientific=T)
# mean=colMeans(common_significant[,3:4])
# mean_logp=mean[2]
data <- sort(common_significant[,4])
ndx <- (dim(common_significant)[1])*75/100+1
line_value <- data[ndx]
line_approx <- round(line_value,2)

#The colour scheme
colours_chr1 <- c("#FB545CFF","#29449DFF","#33FF66","#FAB1B8FF","#9A9A9AFF","#E71C30FF","#29449DFF","#238B68FF","#C2C832FF","#AA4E8DFF","#13AAB3FF","#646464FF","#F3424CFF","#5657A7FF","#CD6FA4FF","#C7C7C7FF","#A11021FF","#262873FF","#187C4AFF","#86874DFF","#A2238AFF","#238B68FF")
colours_chr2 <- c("#660900","#661200","#661b00","#662800","#663600","#664300","#664d00","#665500","#606600","#4e6600","#336600","#146600","#00660c","#006642","#006266","#003e66","#001a66","#990066","#130066","#400066","#660056","#660020")
colours_chr_all <- cbind(colours_chr1,colours_chr2)

colours_chr <- colours_chr_all[c(1:chrNum),select]###Select your own colour scheme 
cols <- c("chr1"=colours_chr[1],"chr2"=colours_chr[2],"chr3"=colours_chr[3],"chr4"=colours_chr[4],"chr5"=colours_chr[5],"chr6"=colours_chr[6],"chr7"=colours_chr[7],"chr8"=colours_chr[8],"chr9"=colours_chr[9],"chr10"=colours_chr[10],"chr11"=colours_chr[11],"chr12"=colours_chr[12],"chr13"=colours_chr[13],"chr14"=colours_chr[14],"chr15"=colours_chr[15],"chr16"=colours_chr[16],"chr17"=colours_chr[17],"chr18"=colours_chr[18],"chr19"=colours_chr[19],"chr20"=colours_chr[20],"chr21"=colours_chr[21],"chr22"=colours_chr[22])
cols <- cols[c(1:chrNum)]##Determine the colour number according to the chromosome number
label<-as.character(paste("chr",c(1:chrNum),sep=""))

#Start to plot
#if(datatype=="SNPs"){
p=ggplot(common_significant,aes(x=common_significant[,3], y=common_significant[,4],colour=as.factor(common_significant[,2])))
p <- p + scale_color_manual(values=cols)
p <- p + geom_point(shape=20,size=2)#Set the point size
p <- p + ggtitle('') + xlab('') + ylab(expression(paste('-log'[10],"(P)",sep="")))#Set the transverse and longitudinal axis labels
p <- p + scale_x_continuous(labels=label,limits=c(0,chrNum),breaks=bpMidVec)+scale_y_continuous(limits=c(0,max_logp+1))#Set the transverse and longitudinal axis calibration mark
#p <- p + geom_hline(yintercept=line_value,linetype="dashed", col='black', lwd=1)#Plus a line where the average
p <- p + theme_bw(base_size=15)+ guides(col = guide_legend(ncol=11,byrow=T,reverse = FALSE,override.aes = list(shape=15,size=7,alpha = 1)))
p <- p + theme(panel.grid=element_blank(),text=element_text(),legend.position="bottom",legend.title=element_blank(),legend.justification = "bottom")
p <- p + theme(legend.key=element_rect(colour="black",size=0.5),legend.key.size=unit(0.5,'cm'),legend.key.width=unit(0.5,'cm'),legend.text = element_text(size=10,hjust=3,vjust=3,face='bold'))##Control colour icons and illustrations 
#p <- p + geom_text(aes(x=0,y=line_value,label=line_approx),vjust=0,size=4.5,col="black")
#p <- p + geom_segment(aes(x=max_mark[3],y=max_mark[4]-0.2,xend=max_mark[3],yend=max_mark[4]),col="#5758ABFF",arrow=arrow(length=unit(0.1,"cm"))) #Add the arrow,adjust size
#p <- p + annotate("text",x=max_mark[3],y=max_mark[4]-0.2,label=as.character(max_mark[1]),size=1.5,colour="#5758ABFF")
p
#}



}
