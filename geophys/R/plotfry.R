plotfry<-function(fry, dis, col=grey(0.7), ann=FALSE, axes=FALSE)
  {
    if(missing(dis)) { dis = diff(range(attr(fry, "X")))/2 }
    if(missing(col)) col = grey(0.7)
    if(missing(ann)) ann=FALSE
    if(missing(axes)) axes=FALSE
    
    mx = attr(fry, "mx")
    my =  attr(fry, "my")

    xlim=c(mx-dis, mx+dis)
    ylim=c(my-dis, my+dis)
    
    plot(attr(fry, "X"), attr(fry, "Y"), type='n' , asp=1, xlim=xlim, ylim=ylim , ann=ann, axes=axes )

    thex = vector()
    they = vector()
    for(i in 1:length(fry))
      {
        flag = fry[[i]]$x>=xlim[1] & fry[[i]]$x<=xlim[2]  & fry[[i]]$y>=ylim[1] & fry[[i]]$y<=ylim[2]
        x = fry[[i]]$x[flag]
        y = fry[[i]]$y[flag]
        points(x, y, pch=".", cex=2, col=col)
        thex = c(thex, x)
        they = c(they, y)

      }

    
    invisible(list(x=thex, y=they, mx=mx, my=my, dis=dis))
    
  }
