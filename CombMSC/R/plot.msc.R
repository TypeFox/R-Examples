`plot.msc` <-
function(x, sumfun.Index = "def", fancy = FALSE, grayscale = FALSE, draw.Grid = TRUE, ...){
  # 3-D plotting function.  Used in plot.msc only when the "fancy" option is TRUE and
  # the number of MSCs is 3, which allows this particular type of display
  showResults <- function(fram, resp, maintitle = "", obj){
    require(rgl)
    x = unique(as.vector(fram[,1]))
    y = unique(as.vector(fram[,2]))
    ht <- x %o% y
    for(i in 1:length(x)){
    for(j in 1:length(y)){
    ifelse(x[i] + y[j] > 1, ht[i,j] <-0, ht[i,j] <- resp[which(fram[,1]==x[i] & fram[,2]==y[j])])
    }}
    ht <- ht * (100/max(resp))
    rgl.open()
    rgl.clear("all")
    rgl.bg(sphere = TRUE, color = c("black", "white"), lit = FALSE, size=.5, alpha=0.25,bg=lines)
    rgl.light()
    rgl.bbox()
    rgl.surface(100*x,100*y,ht,color="blue",alpha=.95,shininess=128,tick=FALSE)
    ind = which.min(resp)
    rgl.spheres(100*fram[,1][ind],(100/max(resp))*resp[ind],100*fram[,2][ind],color="red",radius=2,
    specular="red",texmipmap=T, texminfilter="linear.mipmap.linear")
    axes3d(labels = c(names(obj$msc.List)[1], names(obj$msc.List)[2], "Rank"), tick=FALSE, ntick=0, color="green")
    title3d(xlab = names(obj$msc.List)[1], zlab = names(obj$msc.List)[2], main = maintitle, color = "green")
    rgl.viewpoint(theta = 225, phi = 30)
  }
  if(fancy && length(x$msc.List)==3){
  if(sumfun.Index == "all"){
      for(i in 1:length(x$summary.Functions))showResults(fram= x$Sum.Stats,
      resp = x$Sum.Stats[,names(x$summary.Function)[i]],
    maintitle = names(x$summary.Functions)[i], obj = x)}
    else {if(sumfun.Index=="def")sumfun.Index = 1
           showResults(fram= x$Sum.Stats, resp = x$Sum.Stats[,names(x$summary.Function)[sumfun.Index]],
           maintitle = names(x$summary.Functions)[sumfun.Index], obj = x)}}
  else{
    if(sumfun.Index == "all"){
      for(i in 1:length(x$summary.Functions)){
         win.graph()
         print(plot(x, sumfun.Index=i))

      }
    }
    else{
      if(sumfun.Index == "def"){
        if(length(x$summary.Functions) < 6) plot(x, 1)
        else{
          print(plot(x, sumfun.Index=1, grayscale = grayscale, draw.Grid = draw.Grid),
          split = c(1, 1, 3, 2), more=TRUE)
          print(plot(x, sumfun.Index=2, grayscale = grayscale, draw.Grid = draw.Grid),
          split = c(2, 1, 3, 2), more=TRUE)
          print(plot(x, sumfun.Index=3, grayscale = grayscale, draw.Grid = draw.Grid),
          split = c(3, 1, 3, 2), more=TRUE)
          print(plot(x, sumfun.Index=4, grayscale = grayscale, draw.Grid = draw.Grid),
          split = c(1, 2, 3, 2), more=TRUE)
          print(plot(x, sumfun.Index=5, grayscale = grayscale, draw.Grid = draw.Grid),
          split = c(2, 2, 3, 2), more=TRUE)
          print(plot(x, sumfun.Index=6, grayscale = grayscale, draw.Grid = draw.Grid),
          split = c(3, 2, 3, 2), more=FALSE)
        }
      }
      else{
        require(lattice)
        obj = x
        sumfun = x$summary.Functions[sumfun.Index]
        fram <- x$Sum.Stats
        num.Mscs <- length(x$msc.List)
        if(num.Mscs==1) stop ("Something's wrong... you should have more than one MSC!")
        if(num.Mscs==2) {out <- xyplot(fram[,names(sumfun)] ~ fram[,1], main = paste(
           names(sumfun)), xlab = names(x$msc.List)[1],
           ylab = names(sumfun))
           }
        if(num.Mscs==3) {
          mypanel <- function(x,y,z, ...){
            panel.fill(col="black")
            if(draw.Grid)panel.grid(col.line = "blue", h = 4, v = 4)
            panel.levelplot(x, y, z, background = list(alpha = 1, col = "blue"), ...)
          }
          if(!grayscale){
          out <- levelplot(fram[,names(sumfun)] ~ fram[[1]]
          + fram[[2]], xlab = paste(names(x$msc.List)[1], "Weight"),
          ylab= paste(names(x$msc.List)[2], "Weight"), zlab = names(sumfun), auto.key=TRUE
          , main = paste(names(sumfun)), bg="black", panel=mypanel,
          axis = function(side, ...) {
             axis.default(side, ...)
             if (side == "bottom")
                 panel.text(0, 0, lab = names(obj$msc.List)[[3]], adj = c(1, 1))})
          }

          else out <- levelplot(fram[,names(sumfun)] ~ fram[[1]]
          + fram[[2]], xlab = paste(names(x$msc.List)[1], "Weight"),
          ylab= paste(names(x$msc.List)[2], "Weight"), zlab = names(sumfun), auto.key=TRUE
          , main = paste(names(sumfun)), col.regions = gray(1:100/100),
          axis = function(side, ...) {
             axis.default(side, ...)
             if (side == "bottom")
                 panel.text(0, 0, lab = names(obj$msc.List)[[3]], adj = c(1, 1))})


        }
        if(num.Mscs==4) {
          ints = (1/9)*cbind(0:8, 1:9)
          dim(ints) <- c(9,2)
          Name <- shingle(fram[[3]], intervals=ints)
          out <- (levelplot(fram[[names(sumfun)]] ~ fram[[1]] * fram[[2]]
          |  Name,  xlab = names(x$msc.List)[1],
          ylab= names(x$msc.List)[2], zlab= names(sumfun), auto.key=TRUE,
          strip=strip.custom(var.name=names(fram)[3], shingle.intervals=ints,
          strip.levels=c(FALSE, FALSE)), region=TRUE, main = paste(names(sumfun))))
          out
          }

        if(num.Mscs>4) {print("No default method available."); return(NULL)}
        out
      }
    }
  }
}

