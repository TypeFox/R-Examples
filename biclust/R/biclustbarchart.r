biclustbarchart<-function(x, Bicres, which=NULL, ...)
{
  TAB<-matrix(0,Bicres@Number,dim(x)[2])
  for (i in 1:Bicres@Number) {
    identper <- Bicres@RowxNumber[, i]
    identq<- Bicres@NumberxCol[i,]
    anz <- sum(identper)
    if (anz == 0) {
      break
      }

    TAB[i,]<-colMeans(x[identper,],na.rm=TRUE)
    TAB[i,!identq]<-NA
    }
  m<-colMeans(x,na.rm=TRUE)
  colnames(TAB) <- names(m)


  b<-Barchart(TAB,m,which=which,...)

  grid.newpage()
  plot(b, more=TRUE)

  pushViewport(viewport(x=0.2,y=0, width=0.8,height=0.05,
                      just=c("left","bottom")))

  grid.text("Population mean:", 0.1, 0.5, just=1)
  grid.segments(x0=0.12, y0=0.5, x1=0.2, y1=0.5,
              gp=gpar(col="darkred"))
  grid.points(0.2, 0.5, pch=16,
            size=unit(0.5, "char"), gp=gpar(col="darkred"))

  grid.text("Segmentwise means:", 0.48, 0.5, just=1)
  grid.rect(0.5,0.75,width=0.1,height=0.25, just=c(0,0.5),
          gp=gpar(fill="#A3E0D8"))
  grid.rect(0.5,0.25,width=0.1,height=0.25, just=c(0,0.5),
          gp=gpar(col="#9E9E9E"))

  grid.text("in bicluster", 0.62, 0.75, just=0)
  grid.text("outside bicluster", 0.62, 0.25, just=0)

}



#### Following Part from Flexclust Package!

### Currently not exported, hence document arguments here:
### x: matrix of cluster centers
### m: vector with location of total population
### labels: text for panel header strips (default is rownames(x))
### REST: see barchart method for kccasimple objects
Barchart <- function(x, m, which=NULL, col=NULL, mcol="darkred",
                     strip.labels=NULL, xlab="", ...)
{
    x <- as.matrix(x)
    m <- as.vector(m)

    if(!is.null(strip.labels))
        rownames(x) <- rep(strip.labels, length=nrow(x))

    if(is.null(col))
        col <- c("#FAC8D1", "#D4D8AE", "#A3E0D8", "#D5D0F6", "#EECEB7", "#B5DFBD", "#B2DAEF","#F1C8EA")

    col <- rep(col, length=nrow(x))

    if(is.null(which))
        which <- seq(1, ncol(x))

    x <- x[,which]
    ## sonst musz man die barplots von unten nach oben lesen
    x <- x[,ncol(x):1]
    m <- rev(m[which])

    x <- as.data.frame(as.table(x))

    panel <- createBarchartPanel(m=m, col=col, mcol=mcol)

    barchart(Var2~Freq|Var1, data=x,
             panel=panel, as.table=TRUE,
             xlab=xlab, ...)
}


createBarchartPanel <- function(m, col, mcol)
{
    KKK <- 1
    KKKplus <- function() KKK <<- KKK+1

    mypanel <- function(x, y, shade=FALSE, diff=NULL, ...)
    {
        if(is.null(diff))
            diff <- c(max(m)/4, 0.5)
        else
            diff <- rep(diff, length=2)

        grey <-  "#9E9E9E"
        COL <- rep("white", length(x))
        MCOL <- rep(grey, length=length(x))
        BCOL <- rep(grey, length=length(x))

        if(length(shade)==1){
            if(shade){
                d1 <- abs(x-m) >= diff[1]
                d2 <- abs((x-m)/m) >= diff[2]
                shade <- d1|d2
            }
            else{
                shade <- rep(TRUE, length(x))
            }
        }
        else{
            if(is.matrix(shade)) shade <- shade[KKK,]
            ### reverse to match reversing in Barchart() above
            shade <- rev(rep(as.logical(shade), length=length(x)))
        }

        COL[shade] <- col[KKK]
        MCOL[shade] <- mcol
        BCOL[shade] <- "black"

        MCOL[is.na(x)] <- NA
        grid.segments(x0=0, y0=1:length(x), x1=m, y1=1:length(x),
                      gp=gpar(col=MCOL),
                      default.units="native")

        panel.barchart(x, y, col=COL, border=BCOL, ...)
        grid.points(m, 1:length(x), pch=16,
                    size=unit(0.5, "char"), gp=gpar(col=MCOL))
        grid.segments(1, 1, 4, 4)
        KKKplus()
    }
    return(mypanel)
}

###**********************************************************


