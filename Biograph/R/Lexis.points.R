Lexis.points <-
function (Bdata,transition,title,cov,group,legend.pos,pdf)
{  	z<- check.par (Bdata) 
  if (missing(cov)) cov=NULL
  if (missing(title)) title <- "Title missing" 
  if (missing(legend.pos)) legend <- "topleft"
  if (missing(group)) group=NULL
 # year <- data.frame(YearTrans(Bdata))
  Bdata2 <- date_b(Bdata=Bdata,selectday=1,format.out="year",covs=NULL)  
  z <- TransitionAB(Bdata=Bdata2,transition=transition,keep=FALSE)
  Bdata2 <- Bdata2[Bdata2$ID%in%z$id,] # select subjects that experienced the transition
  
 date <- z$date
 age <- z$age
 # covariates
# location of selected covariates
  if (!is.null(cov)) { poscov <- which (colnames(Bdata2)%in%cov) 
    # covcat <- sort(unique(Bdata2[,poscov]))
    # colnames(GLHS[poscov])
   }
#group <- Bdata2[,poscov[4]]
#colnames(Bdata2)[which(colnames(Bdata2)==cov)] ="test55"
groupn <- group
if (!is.null(group)) group <- Bdata2[,colnames(Bdata2)==groupn]
if (is.null(cov))
  { covn <- "one"
  	covd <- rep(1,nrow(Bdata2)) }  else
  { covn <- cov
    covd <- Bdata2[,colnames(Bdata2)%in%covn]  }
#namecov1 <- colnames(Bdata2[,covd])
lex <- data.frame(age=age,date=date,covd=covd,gp=group)
#colnames(lex)[3] <- colnames(Bdata2[which (colnames(Bdata2)%in%covn==TRUE)])  
#lex <- subset (lex,!is.na(lex$group))
lex$gp <- factor(lex$gp)
#  unique (levels(lex$gp))

ddc <- round_any(date,5,floor)   # function from Plyr
colours<- c("Red","blue","green","yellow","brown")
ncolours <- length(colours)
ncategories.cov <- length(unique(covd))
if (ncategories.cov > ncolours) 
 { zff <- paste("Error in Lexis.points:", ncolours," colours for ",ncategories.cov," categories of covariate ",cov,sep="" )
  stop(zff)	
 }
 colours <- colours[1:ncategories.cov]
ncolours <- ncategories.cov 
#cov1 <- Bdata2[,poscov[1]]
#cov2 <- Bdata2[,poscov[4]]
# namecov2 <- colnames(Bdata2[poscov[4]])
lex2 <- ggplot(lex,aes(x=date,y=age,colour=covd)) + ggtitle(title) # colour p. 48
lex3 <- lex2+geom_point(aes(colour=covd),size=1.2)+scale_colour_manual(values=colours)+coord_equal() 
lex4 <- lex3+layer(geom="point",stat = "identity", position = "identity", params = list(na.rm = FALSE))+
    theme(legend.direction = "vertical",legend.position = legend.pos,legend.background = element_rect(colour = 'purple', fill = 'pink'))+
    theme(plot.background=element_rect(fill="lightskyblue1",colour=NA),
  panel.background=element_rect("black"),
  axis.text.x=element_text(colour="black"),
  axis.text.y=element_text(colour="black"),
  axis.title.x=element_text(colour="darkgreen",face="bold"),
  axis.title.y=element_text(colour="darkgreen",face="bold",angle=90))   # equal scales p. 136

#lex5 <- lex4 +   theme(axis.text.x=element_text(colour="black",size=6))
if (is.null(groupn)) lex6<- lex4 else 
    { lex6 <- lex4+facet_wrap(~gp,nrow=3)
      lex6 <- lex6+theme(legend.direction = "vertical",legend.position = "right",legend.background = element_rect(colour = 'purple', fill = 'pink'))
    }
print (lex6)
if (pdf)
{ pdf("graph.pdf", width = 8, height = 16)  # portrait
  print (lex6)
  dev.off()  }
  
 return(lex6)
}
