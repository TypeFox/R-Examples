### R code from vignette source 'xkcd-intro.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: packages (eval = FALSE)
###################################################
## install.packages(c("xkcd","splancs","reshape"), dependencies=TRUE)


###################################################
### code chunk number 2: fonts (eval = FALSE)
###################################################
## library(extrafont)
## library(ggplot2)
## if( 'xkcd' %in% fonts()) {
##   p <-  ggplot() + geom_point(aes(x=mpg, y=wt), data=mtcars) + 
##     theme(text = element_text(size = 16, family = "xkcd"))
## } else  {
##   warning("Not xkcd fonts installed!")
##   p <-  ggplot() + geom_point(aes(x=mpg, y=wt), data=mtcars) 
## }
## p


###################################################
### code chunk number 3: xkcd-intro.Rnw:152-157 (eval = FALSE)
###################################################
## p <-  ggplot() + geom_point(aes(x=mpg, y=wt), data=mtcars) 
## ggsave("grnofonts.png",p)
## p <-  ggplot() + geom_point(aes(x=mpg, y=wt), data=mtcars) + 
##   theme(text = element_text(size = 16, family = "xkcd"))
## ggsave("grfonts.png",p)


###################################################
### code chunk number 4: xkcd-intro.Rnw:185-199 (eval = FALSE)
###################################################
## library(extrafont)
## download.file("http://simonsoftware.se/other/xkcd.ttf", 
##               dest="xkcd.ttf", mode="wb")
## system("mkdir ~/.fonts")
## system("cp xkcd.ttf  ~/.fonts")
## font_import(pattern = "[X/x]kcd", prompt=FALSE)
## fonts()
## fonttable()
## if(.Platform$OS.type != "unix") {
##   ## Register fonts for Windows bitmap output
##   loadfonts(device="win") 
## } else {
##   loadfonts() 
## }


###################################################
### code chunk number 5: xkcd-intro.Rnw:208-209 (eval = FALSE)
###################################################
## ggsave("gr1.png", p)


###################################################
### code chunk number 6: embedfonts (eval = FALSE)
###################################################
## ggsave("gr1.pdf", plot=p,  width=12, height=4)
## if(.Platform$OS.type != "unix") {
##    ## Needed for Windows. Make sure you have the correct path
##    Sys.setenv(R_GSCMD = 
##               "C:\\Program Files (x86)\\gs\\gs9.06\\bin\\gswin32c.exe")
## }
## embed_fonts("gr1.pdf")


###################################################
### code chunk number 7: xkcd-intro.Rnw:236-237 (eval = FALSE)
###################################################
## install.packages("xkcd",dependencies = TRUE)


###################################################
### code chunk number 8: xkcd-intro.Rnw:241-244 (eval = FALSE)
###################################################
## help(package="xkcd")
## vignette("xkcd-intro") # It opens the PDF
## browseVignettes(package = "xkcd") # To browse the PDF, R and Rnw


###################################################
### code chunk number 9: library (eval = FALSE)
###################################################
## library(xkcd)


###################################################
### code chunk number 10: axis (eval = FALSE)
###################################################
## xrange <- range(mtcars$mpg)
## yrange <- range(mtcars$wt)
## set.seed(123) # for reproducibility
## p <- ggplot() + geom_point(aes(mpg, wt), data=mtcars) + 
##       xkcdaxis(xrange,yrange)
## p


###################################################
### code chunk number 11: xkcd-intro.Rnw:268-269 (eval = FALSE)
###################################################
## ggsave("graxis.png",p)


###################################################
### code chunk number 12: stickfigure (eval = FALSE)
###################################################
## ratioxy <- diff(xrange)/diff(yrange)
## mapping <- aes(x, y, scale, ratioxy, angleofspine,
##                anglerighthumerus, anglelefthumerus,
##                anglerightradius, angleleftradius,
##                anglerightleg, angleleftleg, angleofneck,
##                linetype=city)
## 
## dataman <- data.frame(x= c(15,30), y=c(3, 4),
##                       scale = c(0.3,0.51) ,
##                       ratioxy = ratioxy,
##                       angleofspine =  -pi/2  ,
##                       anglerighthumerus = c(pi/4, -pi/6),
##                       anglelefthumerus = c(pi/2 + pi/4, pi +pi/6),
##                       anglerightradius = c(pi/3, -pi/3),
##                       angleleftradius = c(pi/3, -pi/3),
##                       anglerightleg = 3*pi/2  - pi / 12,
##                       angleleftleg = 3*pi/2  + pi / 12 ,
##                       angleofneck = runif(1, 3*pi/2-pi/10, 3*pi/2+pi/10),
##                       city=c("Liliput","Brobdingnag"))
## 
## p <- ggplot() + geom_point(aes(mpg, wt, colour=as.character(vs)), data=mtcars) + 
##   xkcdaxis(xrange,yrange) + 
##   xkcdman(mapping, dataman)
## p


###################################################
### code chunk number 13: xkcd-intro.Rnw:322-323 (eval = FALSE)
###################################################
## ggsave("grstickfigure.png",p)


###################################################
### code chunk number 14: xkcd-intro.Rnw:327-328 (eval = FALSE)
###################################################
## p +  facet_grid(.~vs)


###################################################
### code chunk number 15: caritas (eval = FALSE)
###################################################
## volunteers <- data.frame(year=c(2007:2011), 
##                          number=c(56470, 56998, 59686, 61783, 64251))
## xrange <- range(volunteers$year)
## yrange <- range(volunteers$number)
## ratioxy <-  diff(xrange) / diff(yrange)
## 
## datalines <- data.frame(xbegin=c(2008.3,2010.5),ybegin=c(63000,59600), 
##                         xend=c(2008.5,2010.3), yend=c(63400,59000))
## 
## mapping <- aes(x, y, scale, ratioxy, angleofspine,
##                anglerighthumerus, anglelefthumerus,
##                anglerightradius, angleleftradius,
##                anglerightleg, angleleftleg, angleofneck)
## 
## dataman <- data.frame( x= c(2008,2010), y=c(63000, 58850),
##                       scale = 1000 ,
##                       ratioxy = ratioxy,
##                       angleofspine =  -pi/2  ,
##                       anglerighthumerus = c(-pi/6, -pi/6),
##                       anglelefthumerus = c(-pi/2 - pi/6, -pi/2 - pi/6),
##                       anglerightradius = c(pi/5, -pi/5),
##                       angleleftradius = c(pi/5, -pi/5),
##                       angleleftleg = 3*pi/2  + pi / 12 ,
##                       anglerightleg = 3*pi/2  - pi / 12,
##                       angleofneck = runif(1, 3*pi/2-pi/10, 3*pi/2+pi/10))
## 
## p <- ggplot() + geom_smooth(mapping=aes(x=year, y =number), 
##                             data =volunteers, method="loess") +
##   xkcdaxis(xrange,yrange) +
##   ylab("Volunteers at Caritas Spain") +
##   xkcdman(mapping, dataman) +
##   annotate("text", x=2008.7, y = 63700, 
##            label = "We Need\nVolunteers!", family="xkcd" ) +
##   annotate("text", x=2010.5, y = 60000, 
##            label = "Sure\nI can!", family="xkcd" ) +
##   xkcdline(aes(xbegin=xbegin,ybegin=ybegin,xend=xend,yend=yend),
##            datalines, xjitteramount = 0.12) 
## p # Figure 5.a


###################################################
### code chunk number 16: xkcd-intro.Rnw:396-397 (eval = FALSE)
###################################################
## ggsave("grcaritas.png",p)


###################################################
### code chunk number 17: barchart (eval = FALSE)
###################################################
## data <- volunteers
## data$xmin <- data$year - 0.1
## data$xmax <- data$year + 0.1
## data$ymin <- 50000
## data$ymax <- data$number
## xrange <- range(min(data$xmin)-0.1, max(data$xmax) + 0.1)
## yrange <- range(min(data$ymin)+500, max(data$ymax) + 1000)
## 
## mapping <- aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax)
## p <- ggplot() + xkcdrect(mapping,data) + 
##   xkcdaxis(xrange,yrange) +
##   xlab("Year") + ylab("Volunteers at Caritas Spain")
## p # Figure 5.b


###################################################
### code chunk number 18: xkcd-intro.Rnw:419-420 (eval = FALSE)
###################################################
## ggsave("grbar.png",p)


###################################################
### code chunk number 19: help (eval = FALSE)
###################################################
## library(zoo)
## library(xkcd)
## require(splancs) #install.packages("splancs", dependencies = TRUE, repos="http://cran.es.r-project.org/")
## 
## mydatar <- read.table(text="
## 6.202
## 5.965 
## 5.778 
## 5.693 
## 5.639 
## 5.273 
## 4.978 
## 4.833 
## 4.910 
## 4.696 
## 4.574 
## 4.645 
## 4.612
## ")
## 
## mydata1 <- mydatar[dim(mydatar)[1]:1,]
## z <- zooreg(mydata1, end = as.yearqtr("2013-1"), frequency = 4)
## z
## 
## 
## mydata <- data.frame(parados=z)
## mydata$year <- as.numeric(as.Date(as.yearqtr(rownames(mydata))))
## mydata$label <- paste(substr(rownames(mydata),3,4),substr(rownames(mydata),6,7),sep="")
## 
## data <- mydata
## data$xmin <- as.numeric(data$year) -1
## data$xmax <- data$xmin + 90
## data$ymin <- 4.5
## data$ymax <- data$parados
## 
## n <- 3200
## poligono <- mydata[,c("year","parados")]
## names(poligono) <- c("x","y")
## poligono <- rbind(poligono, c(max(poligono$x),4.4))
## poligono <- rbind(poligono, c(min(poligono$x),4.4))
## points <- data.frame(x=runif(n,range(poligono$x)[1],range(poligono$x)[2] ),
##                      y=runif(n,range(poligono$y)[1],range(poligono$y)[2] ))
## kk <- inout(points, poligono)
## points <- points[kk, ]
## points <- rbind(points,poligono)
## 
## x <- points$x
## y <- points$y
## nman <- length(x)
## sizer <-runif(nman, 4, 6)
## nman
## 
## xrange <- c(min(x),max(x))
## yrange <- c(min(y),max(y))
## ratioxy <- diff(xrange)/diff(yrange)
## 
## n <- 2
## set.seed(123)
## twomen <-  xkcdman(mapping= aes(x,  y,
##                  scale,
##                  ratioxy,
##                  angleofspine ,
##                  anglerighthumerus,
##                  anglelefthumerus,
##                  anglerightradius,
##                  angleleftradius,
##                  anglerightleg,
##                  angleleftleg,
##                  angleofneck),
##           data.frame(x=c(15600, 14800) ,
##                      y=c(5.3, 5.7),
##                      scale = 0.2,
##                      ratioxy = ratioxy,
##                      angleofspine = runif(n, - pi/2 - pi/10, -pi/2 + pi/10),
##                      anglerighthumerus = runif(n, -pi/6- pi/10, - pi/6 + pi/10),
##                      anglelefthumerus = runif(n, pi + pi/6 -pi/10, pi + pi/6 + pi/10),
##                      anglerightradius =  runif(n, -pi/4, pi/4),
##                      angleleftradius =  runif(n, pi -pi/4, pi + pi/4),
##                      anglerightleg = runif(n,  3* pi/2 + pi/12 , 3* pi/2  + pi/12 + pi/10),
##                      angleleftleg = runif(n, 3* pi/2  - pi/12 - pi/10, 3* pi/2 - pi/12 ),
##                      angleofneck = runif(n, -pi/2-pi/10, -pi/2 + pi/10)))
## 
## p1 <- ggplot() + geom_text(aes(x,y,label="0"), data.frame(x=x,y=y),family="xkcd",alpha=0.4,size=sizer) +  xkcdaxis(xrange,yrange) +
##    ylab("Unemployed persons (millions)") + xlab("Date") +
##  twomen +
##   annotate("text", x= 15250, y=5.95,label="Help!", family="xkcd",size=7) +
##    xkcdline(aes(xbegin=xbegin,ybegin=ybegin,xend=xend,yend=yend),
##             data=data.frame( xbegin=15600, ybegin=5.42, xend=15250, yend=5.902  )
##             , xjitteramount = 200) + theme(legend.position="none")
## #p1
## p2 <- p1 + scale_x_continuous(breaks=as.numeric(mydata$year),label=mydata$label)
## p2
## 
## ggsave("grhelp.png")
## 


###################################################
### code chunk number 20: homosapiens (eval = FALSE)
###################################################
## library(reshape)
## 
## mydata <- read.table(header=TRUE,sep=",",text="
## year,ministerio,banco,fmi,homo
## 2013,2,1.95,1.96,1.94
## 2014,2.1,1.97,1.93,1.88
## 2015,2.2,2.05,1.90,1.87
## ")
## mydatalong <- melt(mydata, id="year", measure.vars= names(mydata)[-1])
## 
## xrange <- c(2013,2015)
## yrange <- c(1.86,2.21)
## set.seed(123)
## ##p <- ggplot() + geom_smooth(aes(x=year, y=value, group=variable,linetype=variable), data=mydatalong, position = position_jitter(h=0.0001),color="black") + theme(legend.position = "none") + xkcdaxis(xrange,yrange)
## p <- ggplot() + geom_smooth(aes(x=year, y=value, group=variable,linetype=variable), data=mydatalong,color="black") + theme(legend.position = "none") + xkcdaxis(xrange,yrange)
## p2 <- p + ylab("Change in real GDP (%)") + xlab("Economic Projections of several Institutes") + scale_x_continuous(breaks=c(mydata$year),   labels=c(mydata$year))
## datalabel <- data.frame(x=2014.95,
##                         y=t(mydata[mydata$year==2015,c(2,3,4,5)]),
##                         label=c("Ministry of Economy","National Bank","International Monetary Fund","Homo Sapiens Sapiens*"))
## names(datalabel) <- c("x","y","label")
## 
## p3 <- p2 + geom_text(aes(x=x,y=y,label=label), data=datalabel, hjust=1, vjust=1,family="xkcd",size=7) +
##   annotate("text", x=2013.4, y=1.852, label="*Homo Sapiens Sapiens = Doubly Wise Man",family="xkcd",size=3.5)
## ggsave("grhomosapiens.png",p3)


###################################################
### code chunk number 21: sevan (eval = FALSE)
###################################################
## resumen <-
##   structure(list(tonombre = structure(c(1L, 2L, 3L, 11L, 4L, 5L, 
##                    8L, 6L, 7L, 9L, 10L, 14L, 12L, 13L, 15L), .Label = c("Andalucía", 
##                                                                "Aragón", "Asturias", "Canarias", "Cantabria", "C-LaMancha", 
##                                                                "CyLeón", "Cataluña", "Extremadura", "Galicia", "Baleares", 
##                                                                "Madrid", "Murcia", "La Rioja", "Valencia"), class = "factor"), 
##                  persons = c(2743706L, 515772L, 364410L, 399963L, 699410L, 
##                    212737L, 2847377L, 717874L, 894946L, 371502L, 942277L, 119341L, 
##                    2561918L, 493833L, 1661613L), frompersons = c(14266L, 3910L, 
##                                                    3214L, 3283L, 4371L, 1593L, 10912L, 8931L, 9566L, 3231L, 
##                                                    5407L, 940L, 21289L, 3202L, 9939L), topersons = c(10341L, 
##                                                                                          3805L, 2523L, 4039L, 3911L, 1524L, 12826L, 10897L, 7108L, 
##                                                                                          2312L, 4522L, 1066L, 26464L, 3529L, 9187L), llegan = c(0.38, 
##                                                                                                                                        0.74, 0.69, 1.01, 0.56, 0.72, 0.45, 1.52, 0.79, 0.62, 0.48, 
##                                                                                                                                        0.89, 1.03, 0.71, 0.55), sevan = c(0.52, 0.76, 0.88, 0.82, 
##                                                                                                                                                                   0.62, 0.75, 0.38, 1.24, 1.07, 0.87, 0.57, 0.79, 0.83, 0.65, 
##                                                                                                                                                                   0.6)), .Names = c("tonombre", "persons", "frompersons", "topersons", 
##                                                                                                                                                                            "llegan", "sevan"), row.names = c(NA, -15L), class = "data.frame")
## 
## 
## resumenlargo <- melt(resumen[,c("tonombre","llegan","sevan")])
## 
## oo <- order(resumen$llegan)
## nombreordenados <- (resumen$tonombre)[oo]
## nombreordenados
## 
## resumenlargo$tonombre <- factor( resumenlargo$tonombre, levels=nombreordenados, ordered=TRUE)
## 
## 
## set.seed(130613)
## kk <- ggplot() + 
##   geom_bar( aes(y= value, x=tonombre,fill=variable ), data=resumenlargo[resumenlargo$variable=="llegan", ], stat="identity") + 
##   geom_bar(aes(y= (-1)* value, x=tonombre,fill=variable ), data=resumenlargo[resumenlargo$variable=="sevan", ], stat="identity") + 
##  scale_y_continuous(breaks=seq(-1.2,1.5,0.3),labels=abs(seq(-1.2,1.5,0.3))) +
##   ylab("Movilidad de los asalariados (% sobre asalariados residentes)") +
##   coord_flip() +
##   theme_xkcd() + xlab("") +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
## 
## 
## 
## kk2 <- kk +
##   geom_text(aes(x=(tonombre),y=0,label=tonombre), data=resumenlargo,family="xkcd")
## 
## 
## lleganespana <- sum(resumen$topersons)*100/ sum(resumen$persons)
## sevanespana <- sum(resumen$frompersons)*100/ sum(resumen$persons)
## 
## 
## lineaespana1 <-   xkcdline(mapping=aes(xbegin=1-0.5,ybegin=lleganespana,xend=15+0.5, yend=lleganespana, yjitteramount=0.051), data= resumenlargo, linetype=2,mask=FALSE)
## lineaespana2 <- xkcdline(mapping=aes(xbegin=1-0.5,ybegin=-lleganespana,xend=15+0.5, yend=-lleganespana), yjitteramount=0.051, data= resumenlargo,linetype=2,mask=FALSE)
## 
## 
## kk3 <- kk2  +   xkcdline(mapping=aes(xbegin=as.numeric(tonombre)-0.5,ybegin=-1.24,xend=as.numeric(tonombre)-0.5, yend=1.52, xjitteramount=0.151), data= resumenlargo, size=3,color="white") + lineaespana1 + lineaespana2
##  
## 
## kk4 <- kk3 + annotate("text",x=1, y=c(lleganespana,-lleganespana),label="Media de España", hjust=c(-0.11,-0.11), vjust=c(-0.1,0.1),family="xkcd",angle=90)
## 
## 
## kk5 <- kk4 + scale_fill_discrete(name="",
##                          breaks=c("llegan", "sevan"),
##                          labels=c("Llegan", "Se van")) + theme(legend.justification=c(0,0), legend.position=c(0,0))
## 
## 
## 
## 
## xrange <- c(1,15)
## yrange <- c(-1.3,1.6)
## ratioxy <- diff(xrange)/diff(yrange)
## x <- 7
## y <-  1.5
## scale <- 0.35
## mapman  <- aes(x,  y,
##                scale,
##                ratioxy,
##                angleofspine ,
##                anglerighthumerus,
##                anglelefthumerus,
##                anglerightradius,
##                angleleftradius,
##                anglerightleg,
##                angleleftleg,
##                angleofneck)
## n <- 1
## set.seed(130613)
## datamanflip <- data.frame( x= x,
##                        y= y,
##                       scale = scale ,
##                       ratioxy = ratioxy,
##                           angleofspine = runif(n, -pi/2 -pi/2 - pi/10,-pi/2 -pi/2 + pi/10),
##                       ##angleofspine = runif(n, -0 - pi/10,-0 + pi/10),
##                       anglerighthumerus = runif(n, -pi/2-pi/6-pi/10, -pi/2 -pi/6+pi/10),
##                       anglelefthumerus = runif(n, -pi/2-pi/2 - pi/10, -pi/2 -pi/2 + pi/10),
##                       anglerightradius = runif(n, -pi/2-pi/5 - pi/10, -pi/2-pi/5 + pi/10),
##                       angleleftradius = runif(n, -pi/2-pi/5 - pi/10, -pi/2-pi/5 + pi/10),
##                       angleleftleg = runif(n, -pi/2 + 3*pi/2  + pi / 12  -pi/20,-pi/2  +3*pi/2  + pi / 12  +pi/20) ,
##                       anglerightleg =  runif(n, -pi/2 + 3*pi/2  - pi / 12  -pi/20, -pi/2+ 3*pi/2  - pi / 12  +pi/20) ,
##                       angleofneck = runif(n, -pi/2+3*pi/2-pi/10, -pi/2+3*pi/2+pi/10))
## p1 <-  xkcdman(mapman , datamanflip) 
## 
## kk6 <- kk5 + p1
## 
## 
## kk7 <- kk6 + annotate("text", x=9.3, y = 1.3, label="Unos vienen, otros se van",family="xkcd" ) +
##   xkcdline(aes(xbegin=xbegin,xend=xend,yend=yend,ybegin=ybegin), yjitteramount=0.135,data=data.frame(xbegin=9.0, xend=7.2, ybegin=1.2, yend=1.3))
## 
## 
## ggsave(kk7,filename="grsevan.png")
## 


###################################################
### code chunk number 22: motherday (eval = FALSE)
###################################################
## mommy <- read.table(sep=" ",text ="
## 8 100
## 9 0
## 10 0
## 11 0
## 12 0
## 13 0
## 14 100
## 15 100
## 16 500
## 17 420
## 18 75
## 19 50
## 20 100
## 21 40
## 22 0
## ")
## names(mommy) <- c("hour","number")
## data <- mommy
## data$xmin <- data$hour - 0.25
## data$xmax <- data$xmin + 1
## data$ymin <- 0
## data$ymax <- data$number
## xrange <- range(8, 24)
## yrange <- range(min(data$ymin) + 10 , max(data$ymax) + 200)
## ratioxy <- diff(xrange)/diff(yrange)
## timelabel <-  function(text,x,y) {
##     te1 <- annotate("text", x=x, y = y + 65, label=text, size = 6,family ="xkcd")
##   list(te1,
##   xkcdline(aes(xbegin=xbegin, ybegin= ybegin, xend=xend,yend=yend),
##            data.frame(xbegin=x,ybegin= y + 50, xend=x,yend=y), xjitteramount = 0.5))
##   }
## n <- 1800
## set.seed(123)
## x <- runif(n, xrange[1],xrange[2] )
## y <- runif(n, yrange[1],yrange[2] )
## inside <- unlist(lapply(1:n, function(i) any(data$xmin <= x[i] & x[i] < data$xmax &
##                             data$ymin <= y[i] & y[i] < data$ymax)))
## x <- x[inside]
## y <- y[inside]
## nman <- length(x)
## sizer <- round(runif(nman, 1, 10),0)
## angler <- runif(nman, -10,10)
## 
## p <- ggplot() +
##   geom_text(aes(x,y,label="Mummy",angle=angler,hjust=0, vjust=0),
##             family="xkcd",size=sizer,alpha=0.3) +
##   xkcdaxis(xrange,yrange) +
##   annotate("text", x=16, y = 650,
##            label="Happy Mother's day", size = 16,family ="xkcd") +
##   xlab("daily schedule") +
##   ylab("Number of times mothers are called on by their children") +
##   timelabel("Wake up", 9, 125) + timelabel("School", 12.5, 90) +
##   timelabel("Lunch", 15, 130) +
##   timelabel("Homework", 18, 525) +
##   timelabel("Bath", 21, 110) +
##   timelabel("zzz", 23.5, 60)
## 
## 
## p
## ggsave("grmotherday.png",p)


###################################################
### code chunk number 23: PDF (eval = FALSE)
###################################################
## require(tools)
## texi2dvi("xkcd-intro.tex", pdf = TRUE)


