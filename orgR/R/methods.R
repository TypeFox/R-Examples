
getstyle <- function(text_size = 20){
  theme_bw() +
    theme(axis.title.x = element_text(colour="black", size=text_size)) +
    theme(axis.text.x = element_text(size = text_size)) +
    theme(axis.title.y = element_text(colour="black", size=text_size)) +
    theme(axis.text.y = element_text(size = text_size)) +
    theme(legend.position="none") +
    theme(plot.title = element_text(face="bold", size = text_size+2, vjust = 2)) 
}

ggpie <- function(data, category = character(), value = numeric()){
    data$category <- data[, category]
  data$value <- data[, value]
  data$category <- factor(data$category, 
                          levels = data$category[order(data$value, decreasing=TRUE)])
  
  p <- ggplot(data, aes(x = factor(1), fill = factor(category), y = (value)/sum(value),
                        order = (value)/sum(value))) +
    geom_bar(stat = "identity", width = 1) + 
    labs(title = "", x = "", y= "") + 
    getstyle(10) + scale_fill_tableau("colorblind10")+
    coord_polar(theta="y", direction = -1) +
    theme(legend.position="right") +
    theme(axis.ticks=element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), 
          legend.text=element_text(size=14), legend.title=element_text(size=14) )+
    guides(fill = guide_legend(title = category))
  return(p)
}

ggpie2 <- function (dat, by, totals) {
    ggplot(dat, aes_string(x=factor(1), y=totals, fill=by)) +
    geom_bar(stat='identity', color='black') +
    scale_fill_brewer() +
    guides(fill=guide_legend(override.aes=list(colour=NA))) + # removes black borders from legend
    coord_polar(theta='y') +
    theme(axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_text(colour='black'),
          axis.title=element_blank(),
          legend.position="none") +
    scale_y_continuous(breaks=cumsum(dat[[totals]]) - dat[[totals]] / 2, labels=dat[[by]]) 
}
