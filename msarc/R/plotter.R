# In: tree of GO terms
# Out: list of depths, named by GO term
treeDepth <- function(tree,current) {
  out <- list()
  for (n in names(tree)) {
    out[[n]] <- current;
    out <- c(out,treeDepth(tree[[n]],current+1))
  }
  return(out)
}

# radians to degrees:
r2d <- function(r) {
  return(as.integer(r*360/(2.0*pi)))
}

listsum <- function(lst) {
  m <- 0
  for (n in lst) {
    m <- m + n
  }
  return(m)
}

write.background <- function(left,right,radius,theta,conn,large) {
  # write a background arc behind each outer arc in its own color but
  # less intense
  xleft <- radius * cos(left)
  yleft <- radius * sin(left)
  xright <- radius * cos(right)
  yright <- radius * sin(right)
  xinner <- 0.1 * cos(theta)
  yinner <- 0.1 * sin(theta)
  rcolor <- hcl(h=r2d(theta),c=25,l=99)
  color <- col2rgb(rcolor) 
  line <- sprintf("<path d=\"M %f,%f L %f,%f A %d,%d 0 %d,1 %f,%f Z\" style=\"fill: rgb(%d,%d,%d);\" />\n",xinner,yinner,xleft,yleft,radius,radius,large,xright,yright,as.integer(color[1]),as.integer(color[2]),as.integer(color[3]))
  cat(line,file=conn,sep="")
}

write.separator <- function(left,radius,theta,conn,large) {
  # write a background arc behind each outer arc in its own color but
  # less intense
  xleft <- radius * cos(left)
  yleft <- radius * sin(left)
  xinner <- 0.1 * cos(theta)
  yinner <- 0.1 * sin(theta)
  rcolor <- hcl(h=r2d(theta),c=0,l=75)
  color <- col2rgb(rcolor) 
  line <- sprintf("<path d=\"M %f,%f L %f,%f\" style=\"stroke-width: 1; stroke: rgb(%d,%d,%d);\" />\n",xinner,yinner,xleft,yleft,as.integer(color[1]),as.integer(color[2]),as.integer(color[3]))
  cat(line,file=conn,sep="")
}

write.arc <- function(arc,left,right,inner,outer,conn,topLevel,colors) {
  if (right - left > pi) {
    large <- 1
  } else {
    large <- 0
  }
  left <- left + 0.01
  right <- right - 0.01
  if (!is.null(colors[[arc]])) {
    middle <- colors[[arc]]
  } else {
    middle <- left + (right - left) / 2
  }
  rcolor <- hcl(h=r2d(middle),c=150)
  color <- col2rgb(rcolor)
  iouter <- as.integer(outer)
  iinner <- as.integer(inner)
  ilarge <- as.integer(large)
  xil <- inner * cos(left)
  xir <- inner * cos(right)
  xol <- outer * cos(left)
  xor <- outer * cos(right)
  yil <- inner * sin(left)
  yir <- inner * sin(right)
  yol <- outer * sin(left)
  yor <- outer * sin(right)
  if (topLevel) {
    write.background(left - 0.01,right + 0.01,iouter,middle,conn,ilarge)
#    write.separator(left - 0.01,iouter,middle,conn,ilarge)
  }
  line <- sprintf("<path d=\"M %f,%f A %d,%d 0 %d,1 %f,%f L %f,%f A %d,%d 0 %d,0 %f,%f Z\" style=\"stroke-width: 5; stroke: rgb(0,0,0); fill: rgb(%d,%d,%d);\" />\n",
                  xol,yol,iouter,iouter,ilarge,xor,yor,xir,yir,iinner,iinner,ilarge,xil,yil,as.integer(color[1]),as.integer(color[2]),as.integer(color[3]))
  cat(line,file=conn,sep="")
  return(list(rcolor,middle))
}

draw.header <- function(radius,padding,conn) {
    middle = radius + padding
    size = middle * 2
    cat("<?xml version=\"1.0\" standalone=\"no\"?>\n",file=conn,sep="")
    cat("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n",file=conn,sep="")
    cat(sprintf("<svg width=\"%spx\" height=\"%spx\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n",size,size),file=conn,sep="")
    cat(sprintf("<g transform=\"translate(%d,%d)\"><g transform=\"rotate(-90)\">\n",middle,middle),file=conn,sep="")
}

draw.footer <- function(conn) {
  cat("</g></g></svg>\n",file=conn)
}

index2letters <- function(index) {
  radix <- 26
  asciiOffset <- 97
  rem <- index %% radix
  quot <- (index %/% radix)
  c <- rawToChar(as.raw(rem+asciiOffset))
  if (quot == 0) {
    s = c
  } else {
    s = paste(index2letters(quot-1),c,sep="",collapse="")
  }
  return(s)
}

write.label <- function(left,right,radius,depth,id,descr,conn,index) {
  if (depth == 1) {
    label <- descr
  } else {
    label <- index2letters(index-1)
  }
  theta <- left + (right - left) / 2
  offset <- theta + pi/2
  if (offset > pi && offset < 2*pi) {
    offset <- offset + pi
    radius <- radius + 15
  }
  offset_deg <- r2d(offset)
  offset_deg <- as.integer(offset_deg %% 360)
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  line <- sprintf("<text id=\"%s\" x=\"%f\" y=\"%f\" text-anchor=\"middle\" style=\"font-size: 25pt;\" transform=\"rotate(%f,%f,%f)\"> %s </text>\n",
                  id,x,y,offset_deg,x,y,label)
  cat(line,file=conn,sep="")
  if (descr != "other") {
    message(sprintf("%s\t%s\t%s",index2letters(index-1),id,descr))
  }
}

draw.arc.set <- function(tree,gotbl,depths,parent,others,mdepth,inner,outer,left,right,conn,index,topLevel,colors) {
  if (length(tree) == 0) {
    return(list())
  }
  arcs <- names(tree)
  arcs <- sort(arcs)
  tot <- sum(gotbl[arcs,]$count) + others[[parent]]
  cleft <- left
  depth <- depths[[arcs[[1]]]]
  bw <- (outer - inner) / mdepth
  arc_inner <- outer - bw * depth
  arc_outer <- arc_inner + bw / 3
  angles <- data.frame()
  for (arc in arcs) {
    # draw arc
    frac <- gotbl[arc,]$count / tot
    width <- (right - left) * frac
    cright <- cleft + width
    rv <- write.arc(arc,cleft,cright,arc_inner,arc_outer,conn,topLevel,colors)
    color <- rv[[1]]
    mid_angle <- rv[[2]]
    write.label(cleft,cright,arc_outer+10,depth,arc,gotbl[arc,]$description,conn,index)
    # draw kids
    kangles <- draw.arc.set(tree[[arc]],gotbl,depths,arc,others,mdepth,inner,outer,cleft,cright,conn,index+1,F,colors)
    if (is.null(dim(kangles))) {
      index <- index + 1
    } else {
      index <- index + 1 + dim(kangles)[1]
    }
    aangles <- data.frame(arc,cleft,cright,color,mid_angle,stringsAsFactors=F)
    angles <- rbind(angles,kangles,aangles)
    cleft <- cright
  } 
  if (topLevel) {
    # draw "other"
    frac <- others[[parent]] / tot
    width <- (right - left) * frac
    cright <- cleft + width
    rv <- write.arc('other',cleft,cright,arc_inner,arc_outer,conn,topLevel,colors)
    write.label(cleft,cright,arc_outer+10,depth,arc,"other",conn,index)
#    aangles <- data.frame("other",cleft,cright,color,stringsAsFactors=F)
#    angles <- rbind(angles,aangles)
#    cleft <- cright
  }
  return(angles)
}

# In: a. Tree of GO terms
#     b. Counts for each GO term
#     c. Radius of inner boundary
#     d. Radius of outer boundary
# Out: left, right angles for each GO term, plus side effect of writing to
#      output file.
draw.arcs <- function(tree,gotbl,others,inner,outer,conn,parent,colors) {
  depths <- treeDepth(tree,1)
  mdepth = listmax(depths)
  left = 0
  right = 2*pi
  angles <- draw.arc.set(tree,gotbl,depths,parent,others,mdepth,inner,outer,left,right,conn,1,T,colors)
  colnames(angles) = c('id','left','right','color','middle')
  return(angles)
}

countOthers <- function(parent,kids,gosets) {
  kset <- list()
  for (k in kids) {
    kset <- c(kset,gosets[[k]])
  }
  notOther <- gosets[[parent]] %in% kset
  other <- gosets[[parent]][!notOther]
  return(length(other))
}

getOthers <- function(tree,gosets,parent="all") {
  kids <- names(tree)
  oc <- countOthers(parent,kids,gosets)
  counts <- list()
  counts[[parent]] = oc
  for (k in kids) {
    counts <- c(counts,getOthers(tree[[k]],gosets,k))
  }
  return(counts)
}

uniprotInTree <- function(tree,gosets) {
  kids <- names(tree)
  ids <- vector()
  for (kid in names(tree)) {
    ids <- c(ids,gosets[[kid]])
  }
  return(unique(ids))
}
  
getMissing <- function(tree,gosets,known) {
  inTree <- uniprotInTree(tree,gosets)
  missFlag <- !(known %in% inTree)
  return(known[missFlag])
}

getBestGoCat <- function(gos,uniprots,chosenFew) {
  final <- list()
  for (id in chosenFew) {
    final[[id]] <- vector()
  }
  final[["other"]] = vector()
  for (up in uniprots) {
    small = 100000
    target = "other" 
    for (id in chosenFew) {
      gset <- gos[[id]]
      if (up %in% gset && length(gset) < small) {
        small <- length(gset)
        target <- id
      }
    }
    final[[target]] = c(final[[target]],up)
  }
  return(final)
}

getSmallestGoCat <- function(gos,u2g) {
  final <- list()
  ups <- names(u2g)
  for (u in ups) {
    small <- 100000
    winner <- "other"
    cans <- u2g[[u]]
    for (c in cans) {
      if (length(gos[[c]]) < small) {
        small <- length(gos[[c]])
        winner <- c
      }
    }
    final[[u]] <- winner
  }
  return(final)
}

mergeByName <- function(l1,l2) {
  n1 = names(l1)
  n2 = names(l2)
  for (n in n2) {
    if (n %in% n1) {
      l1[[n]] <- c(l1[[n]],l2[[n]])
    } else {
      l1[[n]] <- l2[[n]]
    }
  }
  return(l1)
}

getLowestGoCats <- function(tree,gos,uniprots) {
  cats <- names(tree)
  found <- list()
  for (cat in cats) {
    found <- mergeByName(found,getLowestGoCats(tree[[cat]],gos,uniprots))
  }
  fnames = names(found)
  for (cat in cats) {
    us = uniprots[uniprots %in% gos[[cat]]]
    for (u in us) {
      if (!(u %in% fnames)) {
        found[[u]] = cat
      }
    }
  }
  return(found)
}

invertCats <- function(u2g) {
  cats <- list()
  for (u in names(u2g)) {
    c <- u2g[[u]]
    cats[[c]] = c(cats[[c]],u)
  }
  return(cats)
}

findAngle <- function(kids,angles) {
  maxR <- -1
  for (x in kids) {
    r <- angles[x,]$right
    if (r > maxR) {
      maxR <- r
    }
  }
  return(maxR)
}

draw.tag <- function(conn,theta,dist,symbol) {
  if (nchar(symbol) > 6) {
    symbol <- substr(symbol,1,6)
  }
  if (theta > pi && theta < 2*pi) {
    offset = theta + pi
    anchor = "end"
  } else {
    offset = theta
    anchor = "start"
  }
  x <- dist * cos(theta)
  y <- dist * sin(theta)
  td <- r2d(offset)
  line <- sprintf("<text x=\"%f\" y=\"%f\" text-anchor=\"%s\" style=\"font-size: 20pt;\" transform=\"rotate(%f,%f,%f)\">%s</text>\n",
                  x,y,anchor,td,x,y,symbol)
  cat(line,file=conn,sep="")
}

draw.scores <- function(tree,inner,outer,scale,cats,raw,angles,conn) {
  for (id in names(tree)) {
    draw.scores(tree[[id]],inner,outer,scale,cats,raw,angles,conn)
    ups <- cats[[id]]
    kids <- names(tree[[id]])
    if (length(kids) > 0) {
      langle <- findAngle(kids,angles)
    } else {
      langle <- angles[id,]$left
    }
    rangle <- angles[id,]$right
    count <- length(cats[[id]])
    segment <- (rangle - langle) / count
    j = 0
    rcolor <- angles[id,]$color
    color <- col2rgb(rcolor)
    for (u in cats[[id]]) {
      rd <- raw[raw$uniprot==u,]
      dist <- sqrt(as.numeric(rd$score)) * scale
      if (is.na(dist) || dist < 1) {
        dist <- 1
      }
      l <- langle + segment * j
      r <- langle + segment * (j+1)
      m <- (l+r)/2
      xi <- inner * cos(m)
      yi <- inner * sin(m)
      xo <- (inner + dist) * cos(m)
      yo <- (inner + dist) * sin(m)
      line <- sprintf("<path d=\"M %f,%f L %f,%f\" style=\"stroke-width:3; stroke: rgb(%d,%d,%d);\" />\n",
                      xi,yi,xo,yo,color[1],color[2],color[3])
      cat(line,file=conn,sep="")
      if (dist > (outer - inner) * 0.05) {
        draw.tag(conn,r,inner+dist+20,rd$symbol)
      }
      j <- j+1
    }
  }
}

draw.missing <- function(missing,inner,outer,scale,raw,left,right,conn) {
  count <- length(missing)
  segment <- (right - left) / count
  middle <- left + (right - left) / 2
  rcolor <- hcl(h=r2d(middle),c=150)
  color <- col2rgb(rcolor)
  j <- 0
  for (u in missing) {
    rd <- raw[raw$uniprot==u,]
    dist <- sqrt(as.numeric(rd$score)) * scale
    if (is.na(dist) || dist < 1) {
      dist <- 1
    }
    l <- left + segment * j
    r <- left + segment * (j+1)
    m <- (l+r)/2
    xi <- inner * cos(m)
    yi <- inner * sin(m)
    xo <- (inner + dist) * cos(m)
    yo <- (inner + dist) * sin(m)
    line <- sprintf("<path d=\"M %f,%f L %f,%f\" style=\"stroke-width:3; stroke: rgb(%d,%d,%d);\" />\n",
                      xi,yi,xo,yo,color[1],color[2],color[3])
    cat(line,file=conn,sep="")
    if (dist > (outer - inner) * 0.05) {
      draw.tag(conn,r,inner+dist+20,rd$symbol)
    }
    j <- j+1
  }
}

msarc.copyColors <- function(target,source) {
  target$palette <- source$palette
  return(target)
}
  
msarc.plotSVG <- function(msarc,file="msarc.svg") {
  if (requireNamespace("GO.db",quietly=TRUE)) {
    msarc$tree <- makeHierarchy(msarc$candidates,msarc$counts)
    tree <- msarc$tree
    conf <- msarc$conf
    gotbl <- msarc$gotbl
    radius <- conf$radius
    conn <- file(file,open="w")
    gosets <- msarc$go2uniAll
    rawdata <- msarc$data
    few <- msarc$candidates
    colors <- msarc$palette
    padding=10 
    inner = radius * 0.7
    outer = radius
    draw.header(radius,padding,conn)
    others <- getOthers(tree,gosets)
    angles <- draw.arcs(tree,gotbl,others,inner,outer,conn,"all",colors)
    rownames(angles) <- angles$id
    palette <- as.list(angles$middle)
    names(palette) <- angles$id
    msarc$palette <- palette
    sc_inner = radius * 0.1
    sc_outer = radius * 0.6
    lowCats <- getLowestGoCats(tree,gosets,rawdata$uniprot)
    bestCats <- getSmallestGoCat(gosets,lowCats)
    cat2u <- invertCats(bestCats)
    scores <- sqrt(as.numeric(rawdata$score))
    mscore <- max(scores)
    scale <- (sc_outer - sc_inner) / mscore
    draw.scores(tree,sc_inner,sc_outer,scale,cat2u,rawdata,angles,conn)
    missing <- getMissing(tree,gosets,rawdata$uniprot)
    missing <- sort(missing)
    left <- angles$right[length(angles$right)]
    right <- 2 * pi
    draw.missing(missing,sc_inner,sc_outer,scale,rawdata,left,right,conn)
    draw.footer(conn)
    close(conn)
    return(msarc)
  } else {
    stop("Failed to load required package 'GO.db'.")
  }
}

draw.blind.scores <- function(inner,outer,scale,raw,conn,color) {
  lines <- dim(raw)[1]
  perm <- sample(1:lines)
  dist <- (sqrt(as.numeric(raw$score)) * scale)[perm]
  dist[is.na(dist) | dist < 1] = 1
  sym <- raw$symbol[perm]
  for (i in 1:lines) {
    m <- (2*pi) * (i / lines)
    xi <- inner * cos(m)
    yi <- inner * sin(m)
    xo <- (inner + dist[i]) * cos(m)
    yo <- (inner + dist[i]) * sin(m)
    line <- sprintf("<path d=\"M %f,%f L %f,%f\" style=\"stroke-width:3; stroke: rgb(%d,%d,%d);\" />\n",
                      xi,yi,xo,yo,color[1],color[2],color[3])
    cat(line,file=conn,sep="")
    if (dist[i] > (outer - inner) * 0.1) {
      draw.tag(conn,m,inner+dist[i]+20,sym[i])
    }
  }
}

msarc.plotCircle <- function(msarc,file="nocat.svg",col=c(0,255,0)) {
  conf <- msarc$conf
  radius <- conf$radius
  conn <- file(file,open="w")
  rawdata <- msarc$data
  padding=10 
  draw.header(radius,padding,conn)
  sc_inner <- radius * 0.1
  sc_outer <- radius * 0.9
  scores <- sqrt(as.numeric(rawdata$score))
  mscore <- max(scores)
  scale <- (sc_outer - sc_inner) / mscore
  draw.blind.scores(sc_inner,sc_outer,scale,rawdata,conn,col)
  draw.footer(conn)
  close(conn)
  return(msarc)
}

msarc.tagCloud <- function(msarc,pal=NA,...) {
  if (is.na(pal)) {
    pal = brewer.pal(10,"Spectral")
  }
  wordcloud(msarc$data$symbol,sqrt(as.numeric(msarc$data$score)),colors=pal,...)
}
