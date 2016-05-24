segmentLegend <-
function (attribute,
                        colPalette=rainbow(length(attribute)),
                        border=NA,
                        legendName="Legend",
                        bgc='white',
                        temp=FALSE) {


png(filename =ifelse(temp,paste(tempdir(),'/',legendName,'.png',sep=""),paste(legendName,'.png',sep="") ), width=220,
       height=220,units = "px",pointsize = 10, bg="white")
       
 par(mar=c(2.1,3.1,2.1,3.1), bg=bgc)


	niv  <- attribute
  cols <-as.character(substr(colPalette,1,7))

                  x<-rep(1,length(niv))

            pie(  x,
                  labels=niv,
                  clockwise=FALSE,
                  radius=1,
                  col= cols,
                  bg=bgc,
                  init.angle=90)

                graph1 <- dev.cur()
                dev.off(graph1)



    }
