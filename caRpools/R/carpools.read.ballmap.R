carpools.read.ballmap = function(dataset, fullmatchcolumn=2, namecolumn=1, color="blue", title="title", orderby="name", decreasing=FALSE, labelgenes=NULL, labelcolor="red", extractpattern=expression("^(.+?)_.+"), aggregated=FALSE, normtosgrna=FALSE)
{
  # Setup dataset
  # get information about genes
  if(aggregated)
  { # if already gene data
    aggr = dataset[,namecolumn]
    rownames(dataset) = dataset[,namecolumn]
    normtosgrna=FALSE
  }
  else # extract gene identifier for sgRNAs
  {
    rownames(dataset) = dataset[,namecolumn]
    aggr = sub(extractpattern,"\\1",dataset[,namecolumn],perl=TRUE)
    dataset[,namecolumn] = aggr
  }
  
  # normalize to number of sgRNAs present?
  if(normtosgrna==TRUE)
  {
    # check how many entries for a gene identifier and apply to readcount by dividing
    sgRNA = aggregate.data.frame(dataset[,fullmatchcolumn],list(dataset[,namecolumn]),function(x) as.numeric(length(x[x!=0])) )
    row.names(sgRNA)=sgRNA$Group.1
    

    dataset = aggregatetogenes(dataset, agg.function=sum)
   
    dataset[,fullmatchcolumn] = apply(dataset,1, function(i) return(as.numeric(i[fullmatchcolumn])/as.numeric(sgRNA[sgRNA[,1] == i[1],2])))
    dataset[,fullmatchcolumn] = floor(dataset[,fullmatchcolumn])

  }
  
  # order?
  if(orderby == "readcount")
  {
    dataset = dataset[ order(dataset[,fullmatchcolumn], decreasing=decreasing),]
  }
  else if(orderby == "name")
  {
    dataset = dataset[ order(dataset[,namecolumn], decreasing=decreasing),]
  }
  # calculations
  radius=(sqrt(dataset[,fullmatchcolumn])/2/pi)
  lengthx = ceiling(sqrt(nrow(dataset))*1.5)
  lengthy = ceiling(nrow(dataset)/lengthx)
  
  # Label genes/designs?
  if(!is.null(labelgenes))
  { 
    # add column to dataframe    
    dataset[, "labelgene"] = color
    #dataset[aggr %in% labelgenes,"labelgene"] = labelcolor

    colordf = data.frame(gene = labelgenes,
                         color = labelcolor,
                         stringsAsFactors=FALSE)
    row.colordf = nrow(colordf)
    # add column to dataframe    
    
    
    for(i in 1:row.colordf)
    {
      
      dataset$labelgene = apply(dataset, 1, function(x) 
        
        if(x[1] == colordf[i,"gene"])
        { return(colordf[i,"color"])}
        else if (x[1] != colordf[i,"gene"] && x["labelgene"]==color)
        {return(x["labelgene"])}
        else {return(x["labelgene"])}
      )
      
    }
    
  }
  
  # generate grid
  grid=c()
  for(i in seq(1,lengthx,by=1)){
    for(j in seq(1,lengthy,by=1)){
      grid=rbind(grid,cbind(i,j))
    }
  }
  
# generate legend  
  count=1
  x_val=max(grid[,1])+3
  leg=c()
  text_leg=c()
  text_val=c()
  textgrid=c()
  for(j in seq(min(dataset[,fullmatchcolumn]),max(dataset[,fullmatchcolumn]),length.out=lengthy)){
    grid=rbind(grid,cbind(x_val,count))
    leg=c(leg,sqrt(j)/2/pi)
    text_leg=c(text_leg,j)
    text_val=c(text_val,round(j,digits=1))
    textgrid=rbind(textgrid,cbind(x_val+3,count))
    count=count+1
  }

# create empty plot with increased xlimit to allow legend
plot(1, 1, type = "n", axes = FALSE, ann = FALSE, xlim=c(0,max(grid[,1])+10), ylim=c(0,max(grid[,2])))
par(new=TRUE)
# plot symbols
if(is.null(labelgenes)){
  symbols(grid[,1],grid[,2],circles=c(radius[1:(lengthx*lengthy)],leg),inches=0.06,fg="white",bg=color,xlab="",ylab="",cex=3, add=TRUE)
}
else if(!is.null(labelgenes))
{
  symbols(grid[,1],grid[,2],circles=c(radius[1:(lengthx*lengthy)],leg),inches=0.06,fg="white",bg=dataset$labelgene,xlab="",ylab="",cex=3, add=TRUE)
}
# add text for legend
text(textgrid[,1], textgrid[,2], label=as.integer(text_leg), col="black", cex=0.6)
title(main=title)


}