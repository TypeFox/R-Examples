`phist` <-
  function(G, Z=1, L=1, col=2, add=FALSE, tit="GEOL 105 Exam 1")
  {
    if(missing(col)) { col = 1 }
    if(missing(Z))  {  Z = G$grades  }
    if(missing(L))  {  L = G$lett  }
    if(missing(add))  {  add=FALSE  }
    if(missing(tit))  {  tit="GEOL  Exam "  }

    ######  add grades to histogram
    
###  gist = histogram structure
###  Z = scores
#### L = letter grades or labels

    u = par("usr")



    ud = unique(G$grades)
    rediv = sort(unique(c(ud-1, ud+1)))
    HH = hist(G$grades,  breaks = rediv , plot = FALSE  )

    G$hist = HH
    gist = G$hist
    
    if(add==FALSE)
      {
        ## plot(gist, col = grey(.8) )
       ## acols = c("#FFDED9", "#E6E0FF", "#F5DFFF", "#FEE0FF", "#D6F6FF")

        acols=c("#FFCBFE","#D7FFF0","#FFF3DC","#CAFFEE","#FFFAD3")

        
      ##  acols =  pastel.colors(5)
       ## gcols = rep(grey(.8), length(A$lett))
        ##gcols[grep("A", A$lett)]= acols[1]
        ##gcols[grep("B", A$lett)]= acols[2]
        ##gcols[grep("C", A$lett)]= acols[3]
        ##gcols[grep("D", A$lett)]= acols[4]
        ##gcols[grep("E", A$lett)]= acols[5]
       
        icol = findInterval(gist$mids, G$divs, all.inside =TRUE)

       gcols = acols[icol]
        
       ##  plot(gist, col =gcols , xlab="grades", main=tit   )
        
     ##    u = par("usr")
        plot(range(gist$breaks),   range(gist$counts) , type='n', axes=FALSE, ann=FALSE)
        axis(1)
        axis(2)
        
        u = par("usr")


        As = grep("A", L)
        Bs = grep("B", L)
        Cs = grep("C", L)
        Ds = grep("D", L)
        Fs = grep("F", L)

        scol  = rep(0, length(L))
        scol[As] = 1
        scol[Bs] = 2
        scol[Cs] = 3
        scol[Ds] = 4
        scol[Fs] = 5

        
        abline(v=G$divs, lwd=2, col=rgb(.5, .5, 1), lty=2)
        
        ddivs = diff(G$divs)
        xgrad = G$divs[1:(length(G$divs)-1)] + ddivs/2
        
        xmin = gist$breaks[1]+(G$divs[2]-gist$breaks[1])/2
        xgrad[1]= xmin
        text(xgrad , rep(u[4], 5), labels=c("F", "D", "C", "B", "A"), pos=1   )
  
        grid(NA, NULL, lwd = 2)


        den = density(Z)
        Y = RPMG::RESCALE(den$y, u[3], u[4], min(den$y), max(den$y))
        lines(den$x, Y)
        

      }

    ## box()
###  abline(v=seq(from=round(min(Z[Z>0])), to=max(Z), by=2))
    axis(1, at=seq(from=10*round(min(Z[Z>0])/10), to=max(Z), by=10), labels=TRUE)
    
    axis(1, at=seq(from=10*round(min(Z[Z>0])/10), to=max(Z), by=2), labels=FALSE)

    J = rep(1, length(Z))

    f1 = findInterval(G$grades,  gist$breaks, all.inside =TRUE)

    x =    gist$mids[f1]

    ###########################  this is fix to get a good looking plot but it is wrong
   ####  fauxgrades = G$grades-.01
   ####  fauxgrades[fauxgrades<0] = 0
fauxgrades = G$grades
    f1 = findInterval(fauxgrades,  sort(G$hist$breaks) , all.inside =TRUE)

    u = sort( unique(f1) )


  ##  x =    gist$mids[f1]

    x = Z

    
    bx = G$hist$breaks[f1]

    fo = order(f1)
### cbind(G$grades[fo], f1[fo],ix[fo], bx[fo], G$hist$breaks[f1[fo]] )

    y = rep(NA, length(G$grades))

    for(i in 1:length(u)){
      w = which(f1==u[i])
      k = length(w)
      sk = seq(from=0.5, to=k-0.5, length=k)
      y[w] = sk
    }


#####  o = order(Z); cbind(Z[o], J[o])
    
  #####   text(x, y, labels=L, col=col, xpd=TRUE, cex=.8, font=2)
  #####   text(x, y, labels=L, col=acols[scol] , xpd=TRUE, cex=.8, font=2)

    RPMG::textrect(x, y, L, textcol = "black", col=acols[scol])


    
    invisible(list(x=x, y=y, L=L)  )
  }

