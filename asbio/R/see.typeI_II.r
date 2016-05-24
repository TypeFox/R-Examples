see.typeI_II<-function(){

local({
    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if(have_ttk) {
        tkbutton <- ttkbutton
        tkframe <- ttkframe
        tklabel <- ttklabel
    }
tclServiceMode(FALSE)
   top <- tktoplevel() 
   tkwm.geometry(top, "+500+4")
   tktitle(top) <- "Demonstration of type I and type II error"
    canvas <- tkcanvas(top, relief="raised", width=550, height=550)
    tkpack(canvas, side="top", fill="x")

    plotFont <- "Helvetica 9"
    plotFont2 <- "Calibri 11"
    
#------------------------ Message boxes ----------------------#

I<-function(){
 tkmessageBox(message="This constitutes a type I error, i.e. the null hypothesis is rejected, but it is actually true.  Type I error will occur with a probability equivelent to the significance level (alpha) of the hypothesis test.  Type I error is generally considered to be more serious than type II error, and as a result it is emphasized more in hypothesis testing.\n
 
The significance level is often set at 0.05, based on an early reccomendation by R.A. Fisher.  Fisher later argued against using set significance levels")}

II<-function(){
 tkmessageBox(message="This constitutes a correct decision, i.e. the null hypothesis is rejected, and it is false.  The scenario is referred to as power.  Power will occur with a probability equivelent to 1 - the probability of type II error.  A conventional level of power used by working statisticians is 0.8.")}
    
III<-function(){
 tkmessageBox(message="This constitutes a correct decision, i.e. we fail to reject the null hypothesis, and it is true.  The scenario will occur with a probability equivelent to 1 - the probability of type I error.")}
 
IV<-function(){
 tkmessageBox(message="This constitutes a type II error, i.e. we fail to reject the null hypothesis, but it is actually false.  Just as a level of 0.05 is often used for an acceptable level of type I error, a level of 0.2 is often used as a cutoff for type II error.")}
 

    tkcreate(canvas, "polygon", 100, 100, 100, 500, 500, 500, 500, 100, width=0, fill = "white")
    tkcreate(canvas, "line", 100, 100, 500, 100, width=2)
    tkcreate(canvas, "line", 100, 500, 500, 500, width=2)
    tkcreate(canvas, "line", 500, 100, 500, 500, width=2)
    tkcreate(canvas, "line", 100, 500, 100, 100, width=2)
    tkcreate(canvas, "line", 300, 100, 300, 500, width=2)
    tkcreate(canvas, "line", 100, 300, 500, 300, width=2)
    tkcreate(canvas, "text", 300, 20, text="Type I and type II error",
             font=plotFont2)
    tkcreate(canvas, "text", 50, 200, text="Reject H\u2080",
             font=plotFont2)
    tkcreate(canvas, "text", 50, 400, text="FTR H\u2080",
             font=plotFont2)
    tkcreate(canvas, "text", 200, 80, text="H\u2080 True",
             font=plotFont2)
    tkcreate(canvas, "text", 400, 80, text="H\u2080 False",
             font=plotFont2)

#------------------------ alpha ----------------------#
    tkcreate(canvas, "text", 200, 150, text="Type I error\n          \u03b1",
             font=plotFont2)
             
            item1 <- tkcreate(canvas, "polygon", 160, 235, 160, 265, 240, 265, 240, 235,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "point1", "withtag", item1)
         
            tkcreate(canvas, "text", 200, 250, text="More info",
                     font=plotFont2)
                     
#------------------------ 1-beta ----------------------#    
    tkcreate(canvas, "text", 400, 160, text="Correct decision\n         Power\n           1 - \u03b2",
             font=plotFont2)
             
             item2 <- tkcreate(canvas, "polygon", 360, 235, 360, 265, 440, 265, 440, 235,
                           width=1, outline="black",
                           fill="SkyBlue2")        
            tkaddtag(canvas, "point2", "withtag", item2)
            tkcreate(canvas, "text", 400, 250, text="More info",
                     font=plotFont2)

#------------------------ 1-alpha ----------------------#                 
    tkcreate(canvas, "text", 200, 350, text="Correct decision\n            1 - \u03b1",
             font=plotFont2)
    
            item3 <- tkcreate(canvas, "polygon", 160, 435, 160, 465, 240, 465, 240, 435,
                           width=1, outline="black",
                           fill="SkyBlue2")        
            tkaddtag(canvas, "point3", "withtag", item3)
            tkcreate(canvas, "text", 200, 450, text="More info",
                     font=plotFont2)

#------------------------ beta ----------------------#    
    tkcreate(canvas, "text", 400, 350, text="Type II error\n          \u03b2",
             font=plotFont2)
             
            item4 <- tkcreate(canvas, "polygon", 360, 435, 360, 465, 440, 465, 440, 435,
                           width=1, outline="black",
                           fill="SkyBlue2")        
            tkaddtag(canvas, "point4", "withtag", item4)
            tkcreate(canvas, "text", 400, 450, text="More info",
                     font=plotFont2)         

#------------------------ Drive widgets ----------------------#             

        tkitembind(canvas, "point1", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "point1", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "point1", "<1>", substitute(I()))


        tkitembind(canvas, "point2", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "point2", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
                                             
        tkitembind(canvas, "point2", "<1>", substitute(II()))
        
        tkitembind(canvas, "point3", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "point3", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "point3", "<1>", substitute(III()))
        
        tkitembind(canvas, "point4", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "point4", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "point4", "<1>", substitute(IV()))
        

tkpack(tkbutton(top, text = "Exit", command = function()tkdestroy(top)),side="right")
tclServiceMode(TRUE)
})
}


    