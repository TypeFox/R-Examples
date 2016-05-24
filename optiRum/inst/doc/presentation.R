## ----fig.height=15, fig.width=5------------------------------------------
library(ggplot2)
library(optiRum)
 
basicplot<-ggplot(data=iris, 
                  aes(x=Sepal.Width, y=Sepal.Length, colour=Species)) + 
  geom_point()
 
multiplot(basicplot+ggtitle("theme_bw()"), 
          basicplot+theme_minimal()+ggtitle("theme_minimal()"), 
          basicplot+theme_optimum()+ggtitle("theme_optimum()"), layout=matrix(1:3,nrow=3))

