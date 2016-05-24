#####################################################################
#                        Main function                              #
#             Plot TAS diagram (Le Bas et al., 1986)                #
#####################################################################
TASplot <-
function(filename,color="blue",size=3,shape=1) {
 DATA <- AdjRock (filename)


######################################################################
#                          Define plot                               #
######################################################################
## if not already installed, then run:
# install.packages("ggplot2")
#require(ggplot2)


plot.TAS <-
function(DATA) {

 TAS = data.frame(S = DATA$SiO2.adj, NK = DATA$K2O.adj + DATA$Na2O.adj) #define TAS again
 TAS.legend <- data.frame (a=c(74,67,60,54.5,48.5,43,64,64,57.5,52.5,48.86,57.5,53,49,45,43,42.5), b=c(11,4.5,4.0,3.5,3.0,2.0,12.8,10,9,7.5,6,15,11.5,9.3,7.5,5,12.5))
 TAS.names <- c("R","D","A","BA","B","PB","T","TD","TA","BTA","TB","P","TP","PT","Te","Ba","F")
 rownames(TAS.legend) <- TAS.names

with(TAS,
ggplot(data = TAS, aes(x = S, y = NK)) +
     geom_point(aes(x = S, y = NK),color=color,size=size,shape=shape)+
     geom_line(data.frame(x=c(39.2,40,43.2,45,48,50,53.7,55,60,65,77), y=c(0,0.4,2,2.8,4,4.75,6,6.4,8,8.8,10)), mapping=aes(x=x, y=y), colour = "darkred",size=1) + ##Irvine & Baragar##
     geom_line(data=data.frame(x=c(41,41,52.5), y=c(0,7,14)), mapping=aes(x=x, y=y)) +
     geom_line(data=data.frame(x=c(45,45,61), y=c(0,5,13.5)), mapping=aes(x=x, y=y)) +
     geom_line(data=data.frame(x=c(45,52,69), y=c(5,5,8)), mapping=aes(x=x, y=y)) +
     geom_line(data=data.frame(x=c(45,49.4,52,52), y=c(9.4,7.3,5,0)), mapping=aes(x=x, y=y)) +
     geom_line(data=data.frame(x=c(48.4,53,57,57), y=c(11.5,9.3,5.9,0)), mapping=aes(x=x, y=y)) +
     geom_line(data=data.frame(x=c(52.5,57.6,63,63), y=c(14,11.7,7,0)), mapping=aes(x=x, y=y)) +
     geom_line(data=data.frame(x=c(69,69,77.5), y=c(13,8,0)), mapping=aes(x=x, y=y)) +
     geom_line(data=data.frame(x=c(41,45), y=c(3,3)), mapping=aes(x=x, y=y)) +
     geom_line(data=data.frame(x=c(61,64), y=c(13.5,15)), mapping=aes(x=x, y=y)) +
     geom_line(data=data.frame(x=c(52.5,49), y=c(14,15.5)), mapping=aes(x=x, y=y)) +
     geom_text(data=TAS.legend,aes(x=a, y=b, label=c),label=rownames(TAS.legend),size=7, family="Times") +
     xlab(expression(paste(SiO[2],"  (wt.adj %)"))) +
     ylab(expression(paste(Na[2],O," + ",K[2],O, "  (wt.adj %)")))+ 
     geom_text(data = NULL, x = 75, y = 2, label = "Subalkaline", colour = "darkred",size=6, family="Times", fontface="italic", lineheight=.8) +
     geom_text(data = NULL, x = 45, y = 10, label = "Alkaline", colour = "darkred",size=6, fontface="italic",family="Times", lineheight=0.5) +
     theme(text = element_text(size=15,family="Times"))
)
     }

return (plot.TAS(DATA))
}
## End(Not run)
