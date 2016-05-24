plot.irates <-
function(x, 
						covar.code=NULL, 
						full.sample = FALSE,
						n.row=1,
						viewport.size = list(w=3.5, h=2.5),
						box.size = list(w=NULL, h=NULL),
						dist = 0.5,
						irates.vbw = NULL,
						arrow.maxlwd = 10,
						display.digits = 2,
						cex=0.9,
						show.values = TRUE,
						mark = NULL,
						mark.col = "red",
						main = "",
						main.dist = 0.4,
						main.gp = gpar(cex=1.2),
						...)
          {
          	if(class(x) != "irates"){stop("Object needs to be of class irates")}
          	object = x
            u <- "inch"
            
            if(is.null(covar.code)) covar.code = object$covar.code
            covar.code = as.character(covar.code)
            
            covar.labels = rep(NA, length(covar.code))
            for(l in 1:length(covar.code)){
  				covar.labels[l] = object$covar.lab[which(object$covar.code == covar.code[l])]
  			}
            
            if(full.sample){
            	covar.code = c(covar.code, object$full.sample.code)
            	covar.labels = c(covar.labels, object$full.sample.lab)
            	}

            n.events = length(object$event.code)
                
            if ( is.null(box.size$w))
              box.size$w <- 1.5*max(convertUnit(stringWidth(c(object$event.lab,object$no.event.lab)),u,valueOnly=TRUE))
            
            if ( is.null(box.size$h))
              box.size$h <- viewport.size$h*3/(n.events*4-1)
                         
             col.vec <- rep("black",n.events)
             col.vec[mark == 1] <- mark.col

             ## dimensions of grid.layout
             n.col <- ceiling(length(covar.code)/n.row)

             ws <- rep(c(viewport.size$w,dist),n.col-1)
             ws <- c(ws,viewport.size$w)

             hs <- rep(c(viewport.size$h,dist),n.row-1)
             hs <- c(hs,viewport.size$h)

             ## check if there exists main
             if( nchar(main) > 0 )
             {
               par <- main.gp
               lh <- convertUnit(unit(1,"lines"),u,valueOnly=TRUE)
               hs <- c(lh,main.dist,hs)
               y.add <- 2
             }
             else
               y.add <- 0
             
             lo <- grid.layout(nrow=2*n.row-1+y.add,ncol=2*n.col-1,widths=unit(ws,u),heights=unit(hs,u))

             grid.newpage()

             pushViewport(viewport(layout=lo))

             if (nchar(main) > 0)
             {
               pushViewport(viewport(layout.pos.row=1, gp=main.gp))
               grid.text(main)
               popViewport()
             }
             

             n <- 1
             for ( k in covar.code)
             {
               posx <- (((n-1) %% n.col)+1)*2-1
               posy <- (ceiling(n/n.col))*2-1+y.add
               pushViewport(viewport(layout.pos.col=posx, layout.pos.row=posy, gp=gpar(cex=cex)))

               n <- n + 1

               Z <- object$ir[k,]
               
               grid.text(covar.labels[which(k == covar.code)],x=0,y=1,just=c("left","top"))

               ## intitial state
               grid.rect(x=unit(box.size$w*1/2,u),y=0.5,width=unit(box.size$w,u),height=unit(box.size$h,u))
               grid.text(object$no.event.lab,x=unit(box.size$w*1/2,u),y=0.5)

               ## competing event states
               y <- seq(viewport.size$h-box.size$h/2,box.size$h/2,by=-(viewport.size$h-box.size$h)/(n.events - 1))
               x <- rep(viewport.size$w - box.size$w*1/2,n.events)
               grid.rect(x=unit(x,u),y=unit(y,u),
                   width=unit(box.size$w,u),height=unit(box.size$h,u),gp=gpar(col=col.vec))
               grid.text(object$event.lab,x=unit(x,u),y=unit(y,u))

               ## arrows
               alpha.max <- max(object$ir)
               Z.lwd <- (Z*arrow.maxlwd)/alpha.max

               offset <- convertUnit(unit(Z.lwd,"points"),u,valueOnly=TRUE)/2
               y <- seq(viewport.size$h-box.size$h/2,box.size$h/2,by=-(viewport.size$h-box.size$h)/(n.events - 1))

               grid.segments(
                 x0=unit(box.size$w+offset,u),
                 x1=unit(viewport.size$w-box.size$w-offset,u),
                 y0=unit(rep(viewport.size$h/2,n.events,n.events),u),
                 y1=unit(y,u),
                 arrow = arrow(),
                 gp = gpar(lwd = Z.lwd, col=col.vec)
               )
               if ( show.values == TRUE )
               {
                 if ( is.null (irates.vbw))
                   irates.vbw <- 1.5*max(convertUnit(stringWidth(round(Z,display.digits)),u,valueOnly=TRUE))
                 
                 ## arrow values
                 y <- seq(3/4*viewport.size$h-box.size$h/4,viewport.size$h/4+box.size$h/4,
                     by=-(viewport.size$h-box.size$h)/2/(n.events - 1))
                 grid.rect(x=unit(viewport.size$w/2,u),y=unit(y,u),width=unit(irates.vbw,u),height=unit(1.5,"lines"),
                     gp=gpar(col="white", fill="white"))

                 for ( i in 1:n.events ) 
                 {
                 #  if(is.na(k))
                 #  {
                 #    if(show.values == 2)
                 #      grid.text(bquote(hat(alpha)[.(i)]==.(round(Z[i],display.digits))),
                 #          x=unit(viewport.size$w/2,u),y=unit(y[i],u))
                 #    else
                 #      grid.text(round(Z[i],display.digits),
                 #          x=unit(viewport.size$w/2,u),y=unit(y[i],u))
                 #  }
                 #  else
                 #  {
                 #    if(show.values == 2)
                 #      grid.text(bquote(hat(alpha)[.(i)]==.(round(Z[1,i],display.digits))),
                 #          x=unit(viewport.size$w/2,u),y=unit(y[i],u))
                 #    else
                       grid.text(round(Z[i],display.digits),
                           x=unit(viewport.size$w/2,u),y=unit(y[i],u))
                 # }
                 }
               }
               popViewport() 
             }
          }

