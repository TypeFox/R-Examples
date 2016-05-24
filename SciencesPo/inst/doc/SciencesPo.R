## ----echo=TRUE, message=FALSE, comment=NA--------------------------------
library(SciencesPo)


## Do things 


detach("package:SciencesPo", unload=TRUE)


#You can also use the unloadNamespace command,

unloadNamespace("SciencesPo")


## ----echo=TRUE, message=FALSE--------------------------------------------
help.search('bar.plot')


help.search('twoway', package = 'SciencesPo')

## ----echo=TRUE, message=FALSE--------------------------------------------
vignette(package = "SciencesPo")

## ----echo=TRUE, message=FALSE--------------------------------------------
data(package = "SciencesPo")

## ----echo=FALSE, message=FALSE-------------------------------------------
require(SciencesPo)

set.seed(51)
 w <-sample(4,10, TRUE)
 x <- sample(10, 1000, replace=TRUE, prob=w)
 
skewness(x, type = 1);
kurtosis(x, type = 1);
skewness(x); # Type 2 is the default 
kurtosis(x); # Type 2 is the default 
skewness(x, type = 3);
kurtosis(x, type = 3);


## ----echo=FALSE, message=FALSE-------------------------------------------
aad(pres) 


winsorize(pres)

## ----echo=FALSE, message=FALSE-------------------------------------------
require(SciencesPo)
str(iris)

iris_2 = safe.chars(iris)

str(iris_2)

## ----echo=TRUE, message=FALSE--------------------------------------------
require(SciencesPo)

mylevels <- c('Strongly Disagree', 
              'Disagree', 
              'Neither', 
              'Agree', 
              'Strongly Agree')

myvar <- factor(sample(mylevels[1:5], 10, replace=TRUE))

## ----echo=TRUE, message=FALSE--------------------------------------------
unclass(myvar) # testing the order

## ----echo=TRUE, message=FALSE--------------------------------------------
destring(myvar) 

## ----echo=TRUE, message=FALSE--------------------------------------------
 (x = seq(0, 1, by=.1))
 rounded(x) 

## ----one-way, eval=FALSE, echo=TRUE, message=FALSE, comment=NA-----------
#  CrossTabs(titanic$SURVIVED)

## ----Freq, echo=TRUE, message=FALSE, comment=NA--------------------------
Frequency(titanic, SURVIVED) 

## ----two-way, echo=TRUE, message=FALSE-----------------------------------
crosstable(titanic, SEX, CLASS, SURVIVED) 

## ----politicalDiversity1, echo=TRUE, message=FALSE-----------------------
library("SciencesPo")

# The 1980 presidential election in the US (vote share):

(US1980 <- c("Democratic"=0.410, "Republican"=0.507,
              "Independent"=0.066, "Libertarian"=0.011,
              "Citizens"=0.003, "Others"=0.003));

politicalDiversity(US1980); # ENEP (laakso/taagepera) method 

politicalDiversity(US1980, index= "golosov");


## ----Helsinki-election, echo=TRUE, message=FALSE-------------------------
# Helsinki's 1999

Helsinki <- data.frame(votes = c(68885, 18343, 86448, 21982, 51587,
                                 27227, 8482, 7250, 365, 2734, 1925,
                                 475, 1693, 693, 308, 980, 560, 590, 185),
                       seats.SL=c(5, 1, 6, 1, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0),
                       seats.dH=c(5, 1, 7, 1, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0))

## ----echo=TRUE, message=FALSE, comment=NA--------------------------------
# politicalDiversity(Helsinki$votes); #ENEP Votes

politicalDiversity(Helsinki$seats.SL); #ENP for Saint-Lague

politicalDiversity(Helsinki$seats.dH); #ENP for D'Hondt

## ----highestAverages1, echo=TRUE, message=FALSE, comment=NA--------------
highestAverages(parties=names(Ceara), votes=Ceara,
                seats = 42, method = "dh") 

## ----echo=TRUE, message=FALSE, comment=NA--------------------------------
highestAverages(parties=names(Ceara), votes=Ceara,
                seats = 42, method = "sl") 

## ----echo=TRUE, message=FALSE, comment=NA--------------------------------
highestAverages(parties=names(Ceara), votes=Ceara, 
                seats = 42, method = "hh") 

## ----echo=TRUE, message=FALSE, comment=NA--------------------------------
highestAverages(parties=names(Ceara), votes=Ceara, 
                seats = 42, method = "imperiali") 

## ----echo=TRUE, message=FALSE, comment=NA--------------------------------
highestAverages(parties=names(Ceara), votes=Ceara,
               seats = 42, method = "dh", threshold = 5/100) 

## ----eval=FALSE, echo=TRUE, message=FALSE, comment=NA--------------------
#  largestRemainders(parties=names(Ceara), votes=Ceara,
#                  seats = 42, method = "hare")

## ----eval=FALSE, echo=TRUE, message=FALSE, comment=NA--------------------
#  largestRemainders(parties=names(Ceara), votes=Ceara,
#                  seats = 42, method = "droop")

## ----data-Italy, eval=FALSE, echo=TRUE, message=FALSE--------------------
#  
#  # The 1946 Italian Constituent Assembly election results: parties and unspoilt votes
#  
#  Italy = data.frame(party=c("DC", "PSIUP", "PCI", "UDN", "UQ", "PRI",
#                              "BNL", "PdA", "MIS", "PCd'I", "CDR",
#                             "PSd'Az", "MUI", "PCS", "PDL", "FDPR"),
#                     votes=c(8101004, 4758129, 4356686, 1560638,	1211956,
#                             1003007, 637328, 334748, 171201, 102393,
#                             97690, 78554, 71021, 51088, 40633, 21853))

## ----eval=FALSE, echo=TRUE, message=FALSE, comment=NA--------------------
#  with(Italy, largestRemainders(parties=party, votes=votes,
#                  seats = 556, method = "imperiali.q") )

## ----echo=TRUE, message=FALSE, comment=NA--------------------------------
mytable = highestAverages(parties=names(Ceara), votes=Ceara, 
                seats = 42, method = "dh") 

library(knitr)

kable(mytable, align=c("l","c","c"))

## ----echo=TRUE, message=FALSE, fig.width=4.5, fig.height=4.5, fig.align="center", fig.cap= "2014 Legislative Election in Ceara (M=42)"----

mytable = highestAverages(parties=names(Ceara), votes=Ceara, 
                seats = 42, method = "dh") 

p <- ggplot(mytable, aes(x=reorder(Party, Seats), y=Seats)) + 
  geom_bar(position="dodge", stat = "identity") +
  coord_flip() + labs(x="", y="# Seats")
p + theme_grey() 

## ----eval=TRUE-----------------------------------------------------------
detach("package:SciencesPo")

ggplot(mtcars, aes(mpg, disp,color=factor(carb),size=hp)) + geom_point(alpha=0.7) + labs(title="Bubble Plot") + scale_size_continuous(range = c(3,10))

qplot(1:3, 1:3)

## ----eval=TRUE-----------------------------------------------------------
require(SciencesPo)
qplot(1:3, 1:3)

## ----echo=FALSE, message=FALSE-------------------------------------------
require(SciencesPo)
theme_set(theme_pub(font_size=12)) # default fontsize doesn't work well for online viewing
qplot(1:3, 1:3)

## ----echo=FALSE, message=FALSE-------------------------------------------
require(SciencesPo)
# "Verdana", "serif" and "sans" are also high-readability fonts
theme_set(theme_pub(font_size=12, font_family = "Consolas")) 
qplot(1:3, 1:3)

## ----echo=FALSE, message=FALSE-------------------------------------------
prefs <- theme(axis.text = element_text(size=14, colour=NULL))

qplot(1:3, 1:3) + prefs

## ----echo=FALSE, message=FALSE-------------------------------------------
# Modifying a theme function
themeMod <- theme_gray() +
  theme(text = element_text(family = "Times", colour = "blue", size = 14))

ggplot(mpg, aes(x = cty, y = hwy, colour = factor(cyl))) + 
   geom_point(size = 2.5)

## ----echo=FALSE, message=FALSE-------------------------------------------
# Only change the 'colour' property of theme element 'text'

mytheme1 <- theme_grey() + theme(text = element_text(colour="red"))
mytheme1$text

## ----echo=FALSE, message=FALSE-------------------------------------------
# Replace the 'text' element entirely
mytheme2 <- theme_grey() %+replace% theme(text = element_text(colour="red"))
mytheme2$text

## ----eval=FALSE, message=FALSE-------------------------------------------
#  plot.mpg + background_grid(major = "xy", minor = "none")

## ----eval=FALSE, message=FALSE, fig.width=7, fig.height=5----------------
#  
#  plot.iris <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
#    geom_point() + facet_grid(. ~ Species) + stat_smooth(method = "lm") +
#    background_grid(major = 'y', minor = "none") + # add thin horizontal lines
#    panel_border() # and a border around each panel
#  # plot.mpg and plot.diamonds were defined earlier
#  ggdraw() +
#    draw_plot(plot.iris, 0, .5, 1, .5) +
#    draw_plot(plot.mpg, 0, 0, .5, .5) +
#    draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)

## ----eval=FALSE, echo=FALSE, message=FALSE, fig.width=7, fig.height=5----
#  # Of course, we can also go crazy:
#  ggdraw() +
#    #geom_rect(data = boxes, aes(xmin = x, xmax = x + .15, ymin = y, ymax = y + .15),
#    #          colour = "gray60", fill = "red", alpha=.03) +
#    geom_path(data = spiral, aes(x = x, y = y, colour = t), size = 6, alpha = .4) +
#    draw_plot(plot.mpg, .3, .3, .4, .4) +
#    draw_plot(plot.iris, 0, .7, .7, .35 ) +
#    draw_plot(plot.iris, .45, .0, .6, .3 )

## ----height.matters, fig.width=7, fig.height=5---------------------------
theme_set(theme_pub())

# Generating a ratio winner/opponent measure 
Presidents = transform(Presidents, 
                       height_ratio = winner.height/opponent.height) 

# Avoid missing data
Presidents <- subset(Presidents, !is.na(height_ratio))

fit=lm(winner.vote~height_ratio,data=Presidents)

mylabel=lm2eqn("Presidents","height_ratio","winner.vote")

p1 <- ggplot(Presidents, aes(x=height_ratio, y=winner.vote)) +
      geom_smooth(method=lm, colour="red", fill="gold")+
      geom_point(size = 5, alpha = .7) +
      annotate(geom = 'text', x = 1.1, y = 70, size = 5, label = mylabel, fontface = 'italic') +
      xlim(0.85,1.2) + ylim(25, 70) +
      xlab("Winner/Opponent Height Ratio") + 
      ylab("Relative Support for the Winner")
p1 

geom_foot("Draft Analysis, 2015", color = fade("brown1"))


