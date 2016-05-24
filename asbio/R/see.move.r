see.move<-function(){
local({
    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if(have_ttk) {
        tkbutton <- ttkbutton
        tkframe <- ttkframe
        tklabel <- ttklabel
    }
    
    tclServiceMode(FALSE) 
    top <- tktoplevel()
    tktitle(top) <- "Demonstration of Least Squares Regression -- Move points"
 

    buttons <- tkframe(top)
    tkpack(buttons, side="bottom", fill="x", pady="2m")
    dismiss <- tkbutton(buttons, text="Exit",
                        command=function()tkdestroy(top))
    tkpack(dismiss, side="right", expand=TRUE)
    
    canvas <- tkcanvas(top, relief="raised", width=550, height=350)
    tkpack(canvas, side="top", fill="x")

    plotFont <- "Helvetica 9"
    plotFont2 <- "Helvetica 11"
    

    tkcreate(canvas, "polygon", 100, 50, 100, 250, 400, 250, 400, 50, width=0, fill = "white")
    tkcreate(canvas, "line", 100, 250, 400, 250, width=2)
    tkcreate(canvas, "line", 100, 250, 100, 50, width=2)
    tkcreate(canvas, "text", 275, 20, text="Moving points in simple linear regression",
             font=plotFont2)
    tkcreate(canvas, "text", 250, 290, text="x",
             font=c("Helvetica", "11","italic"))
    tkcreate(canvas, "text", 30, 150, text="y",
             font=c("Helvetica", "11","italic"))
             
     tkcreate(canvas, "polygon", 430, 60, 430, 180, 530, 180, 530, 60, width=2, fill = "white")

    # X tickmarks & labels; from x = 100 to 400 
        for (i in 0:10) {
        x <- 100 + i * 30
        tkcreate(canvas, "line", x, 250, x, 245, width=2)
        tkcreate(canvas, "text", x, 254,
                 text=10*i, anchor="n", font=plotFont)
    }
    # Y tickmarks & labels; from y = 50 to 250
    for (i in 0:5) {
        y <- 250 - i * 40
        tkcreate(canvas, "line", 100, y, 105, y, width=2)
        tkcreate(canvas, "text", 96, y,
                 text=formatC(50*i,format="f",digits=1),
                 anchor="e", font=plotFont)
    }

    # The (original) data
    points <- matrix(c(12, 56,
                       20, 94,
                       33, 98,
                       32, 120,
                       61, 180,
                       75, 160,
                       98, 223), ncol=2, byrow=TRUE)

    ## `self-drawing' point object
    point.items <- apply(points, 1, function(row) {
        x <- 100 + 3 * row[1]
        y <- 250 - 4/5 * row[2]
        item <- tkcreate(canvas, "oval", x - 5, y - 5, x + 5, y + 5,
                         width=1, outline="black",
                         fill="SkyBlue2")
        tkaddtag(canvas, "point", "withtag", item)
        item
    })

    plotDown <- function(x, y) {
      ## Arguments:
      ## x, y -	The coordinates of the mouse press.
      x <- as.numeric(x)
      y <- as.numeric(y)
      tkdtag(canvas, "selected")
      tkaddtag(canvas, "selected", "withtag", "current")
      tkitemraise(canvas,"current")
      lastX <<- x
      lastY <<- y
    }

    plotMove <- function(x, y) {
        ## Arguments:
        ## x, y -	The coordinates of the mouse.
        x <- as.numeric(x)
        y <- as.numeric(y)
        tkmove(canvas, "selected", x - lastX, y - lastY)
        lastX <<- x
        lastY <<- y
    }

    step1 <- function(){lapply(point.items,
                         function(item)
                         as.double(tkcoords(canvas,item)))}
    
    plotLine <- function(){
        coords <-step1()
        x <- sapply(coords, function(z) (z[1]+z[3])/2)
        y <- sapply(coords, function(z) (z[2]+z[4])/2)
        lm.out <- lm(y~x)
        x0 <- range(x)
        y0 <- predict(lm.out, data.frame(x=x0))
        tkcreate(canvas, "line", x0[1], y0[1], x0[2], y0[2], width=3)
        
    }

    line <- plotLine()

    slope.func<-function(){
    coords <-step1()
        x <- sapply(coords, function(z) (z[1]+z[3])/2)
        x <- (x-100)/3
        y <- sapply(coords, function(z) (z[2]+z[4])/2)
        y <- -1*(y-250)*5/4
        lm.out<-lm(y~x)
        slope<-lm.out$coefficients[2]
        tkcreate(canvas, "text", 500, 80, text=formatC(slope,format="f",digits=2),
             font=plotFont)      
        }
    
    slope<-slope.func()    
    
    pre.func0 <- function(){
     tkcreate(canvas, "text", 460, 80, text= "Slope = ",
             font=plotFont)}
             
    pre0 <-pre.func0()
    
    yint.func<-function(){
    coords <-step1()
        x <- sapply(coords, function(z) (z[1]+z[3])/2)
        x <- (x-100)/3
        y <- sapply(coords, function(z) (z[2]+z[4])/2)
        y <- -1*(y-250)*5/4
        lm.out<-lm(y~x)
        yint<-lm.out$coefficients[1]
        tkcreate(canvas, "text", 505, 105, text=formatC(yint,format="f",digits=2),
             font=plotFont)      
        }
    
    yint<-yint.func()    
    
    pre.func01 <- function(){
     tkcreate(canvas, "text", 460, 105, text= "Y int. = ",
             font=plotFont)}
             
    pre01 <-pre.func01()
    
    
    r.func <-function(){
        coords <-step1()
        x <- sapply(coords, function(z) (z[1]+z[3])/2)
        y <- sapply(coords, function(z) (z[2]+z[4])/2)
        r<-round(cor(y,x),2)
        tkcreate(canvas, "text", 485, 140, text=formatC(-1*r,format="f",digits=2),
             font=plotFont)
               }
    
    r <-r.func()
     
    pre.func1 <- function(){
     tkcreate(canvas, "text", 460, 140, text= "r = ",
             font=plotFont)}
     
    pre1 <-pre.func1() 
        
    r2.func <-function(){
        coords <-step1()
        x <- sapply(coords, function(z) (z[1]+z[3])/2)
        y <- sapply(coords, function(z) (z[2]+z[4])/2)
        r<-round(cor(y,x),2)
        tkcreate(canvas, "text", 485, 165, text=formatC(r*r,format="f",digits=2),
             font=plotFont)
              }
    
     r2 <-r2.func()
    
     pre.func2 <- function(){
     tkcreate(canvas, "text", 460, 165, text= 'r\u00b2 = ',
             font=plotFont)}
     
    pre2 <-pre.func2() 
    
    lastX <- 0
    lastY <- 0

    tkitembind(canvas, "point", "<Any-Enter>",
               function() tkitemconfigure(canvas, "current",
                                          fill="red"))
    tkitembind(canvas, "point", "<Any-Leave>",
               function() tkitemconfigure(canvas, "current",
                                          fill="SkyBlue2"))
    tkitembind(canvas, "point", "<1>", plotDown)
    tkitembind(canvas, "point", "<ButtonRelease-1>",
               function(x){
                   tkdtag(canvas, "selected")
                   tkdelete(canvas, "withtag", line)
                   tkdelete(canvas, "withtag", slope)
                   tkdelete(canvas, "withtag", yint)
                   tkdelete(canvas, "withtag", r)
                   tkdelete(canvas, "withtag", r2)
                   tkdelete(canvas, "withtag", pre0)
                   tkdelete(canvas, "withtag", pre01)
                   tkdelete(canvas, "withtag", pre1)
                   tkdelete(canvas, "withtag", pre2)
                   
                   line <<- plotLine()
                   slope <<- slope.func()
                   yint <<- yint.func()
                   pre0 <<- pre.func0()
                   pre01 <<- pre.func01()
                   r <<- r.func()
                   r2 <<- r2.func()
                   pre1 <<- pre.func1()
                   pre2 <<- pre.func2()
                      })
    tkbind(canvas, "<B1-Motion>", plotMove)
    tclServiceMode(TRUE)
})
}


    