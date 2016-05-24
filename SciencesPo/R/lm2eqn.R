#' @title Linear model to equation style
#'
#' @description Produces a text equation style to be added in plots.
#'
#' @param .data The data.frame object.
#' @param x The independent variable(s).
#' @param y The dependent variable.
#' @param spaced A logical value indicating if spaces should be added; default is TRUE.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @examples
#' lm2eqn("mtcars","wt","mpg")
#'
#'
#' data(Presidents)
#'
#' Presidents <- transform(Presidents, ratio = winner.height/opponent.height)
#' Presidents <- transform(Presidents, selected = ifelse(winner %in% c("Barack Obama"),1,0))
#'
#' # subsetting election > 1824
#' Presidents = subset(Presidents, election > 1824 & !is.na(ratio))
#'
#' selected=Presidents[Presidents$selected==1,]
#' myeqn=lm2eqn("Presidents","ratio","winner.vote")
#'
#' ggplot(Presidents, aes(x=ratio,y=winner.vote,colour=selected)) +
#'  geom_text(data=selected,aes(label=winner),hjust=-0.1) +
#'  geom_smooth(method=lm, colour="red", fill="gold") +
#'  geom_point(size=5, alpha=.7) +
#' annotate(geom='text',x=1.1,y=64,size=7,label=myeqn,family='Times',fontface='italic') +
#' xlim(.9,1.2) + ylim(38, 65) +
#'  scale_colour_gradientn(guide="none" , colours=c("black","red")) +
#' xlab("Presidential Height Ratio") +
#'  ylab("Relative Support for President")
#'
#' @export
`lm2eqn` <- function(.data, x, y, spaced=TRUE){
  fit=eval(parse(text=paste0("glm(",y,"~",x,", na.action='na.omit', data=",.data,")")))
  intercept=round(stats::coef(fit)[1],1)
  slope=round(stats::coef(fit)[2],1)
  if(spaced) equation=paste0("y = ",slope,"x",ifelse(intercept>=0,' + ',' - '),abs(intercept))
   else equation=paste0("y==",slope,"*x",ifelse(intercept>=0,'+','-'),abs(intercept))
  p=round(summary(fit)$coeff[2,4],3)
  if(p==0) equation=paste(equation,"(p < 0.001)")
  else equation=paste(equation,"(p =",p,")")
  equation
}
NULL
