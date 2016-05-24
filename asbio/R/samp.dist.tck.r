samp.dist.method.tck<-function(){

tclRequire("BWidget")
local({
    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if(have_ttk) {
        tkbutton <- ttkbutton
        tkcheckbutton <- ttkcheckbutton
        tkentry <- ttkentry
        tkframe <- ttkframe
        tklabel <- ttklabel
        tkradiobutton <- ttkradiobutton
    }

tclServiceMode(FALSE)
dialog.ci <- function(){
        tt <- tktoplevel()
        tkwm.title(tt,"Sampling distributions")
        Par1<- tclVar("snapshot")
        Stats<-tclVar("mean") 
        done <- tclVar(0)
    
reset <- function()
        {
        tclvalue(Par)<-"sanpshot"
        tclvalue(Stats)<-"mean"
        } 

distrib<-function(Par,Stats){
        if(Par=="animation") samp.dist.tck(Stats)
        if(Par=="snapshot"){
            if(Stats == "mean"|Stats == "median"|Stats == "trimmed mean"| Stats == "Winsorized mean"| Stats == "sd"| Stats == "var"| Stats == "Huber estimator" | Stats == "H-L estimator" | Stats == "IQR" | Stats == "MAD") samp.dist.snap.tck1(Stats)
            if(Stats != "mean"  & Stats != "median" & Stats != "trimmed mean" & Stats != "Winsorized mean" & Stats != "sd" & Stats != "var" & Stats != "Huber estimator" & Stats != "H-L estimator" & Stats != "IQR" & Stats != "MAD") samp.dist.snap.tck2(Stats)
        }
        }
  
build <- function()
        {
            Par <-tclvalue(Par1)
            Stats <-tclvalue(Stats)
            distrib(Par, Stats)
        }
               
        
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",
                               command=function()tclvalue(done)<-0)

        
        alt.rbuts <- tkframe(tt)
        tkpack(tklabel(alt.rbuts, text="Depiction method"))
        for ( i in c("animation", "snapshot")){
            tmp <- tkradiobutton(alt.rbuts, text=i, variable=Par1, value=i)
            tkpack(tmp,anchor="w")
        }

        tkgrid(tklabel(tt, text="Statistic"),columnspan=2)
        stats<- c("custom", "mean","median", "trimmed mean", "Winsorized mean", "Huber estimator", "H-L estimator", "sd", "var", "IQR", "MAD","(n-1)S^2/sigma^2", "F*", "t* (1 sample)", "t* (2 sample)", "Pearson correlation", "covariance")
            comboBox <- tkwidget(tt,"ComboBox", editable=FALSE, values=stats, textvariable = Stats, width = 20)
        tkgrid(comboBox,columnspan=2)
        tkgrid(tklabel(tt,text=""),columnspan=2)
        tkgrid(alt.rbuts,columnspan=2)
        tkgrid(tklabel(tt,text=""),columnspan=2)
        tkgrid(submit.but,reset.but)
            
tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)
tkwait.variable(done)
if(tclvalue(done)=="2") stop("aborted")
tkdestroy(tt) 
cmd <- build()
eval.parent(cmd)
tclServiceMode(TRUE)
    }                            
Par1 <-tclVar("snapshot")
Stats <- tclVar("mean")
dialog.ci()

})
}

 



samp.dist.tck<-function(statc = "mean"){

local({
    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if(have_ttk) {
        tkbutton <- ttkbutton
        tkcheckbutton <- ttkcheckbutton
        tkentry <- ttkentry
        tkframe <- ttkframe
        tklabel <- ttklabel
        tkradiobutton <- ttkradiobutton
    }
    tclServiceMode(FALSE)
    dialog.sd <- function(){
        tt <- tktoplevel()
        tkwm.title(tt,"Sampling distributions")
        biv.parent.entry <- tkentry(tt, textvariable=Biv.parent)
        parent.entry1 <- tkentry(tt, textvariable=Parent1)
        parent.entry2 <- tkentry(tt, textvariable=Parent2)
        s.size1.entry<-tkentry(tt, textvariable=S.size1)
        s.size2.entry<-tkentry(tt, textvariable=S.size2)
        stat.entry1<-tkentry(tt, textvariable=Stat1)
        stat.entry2<-tkentry(tt, textvariable=Stat2)
        stat.entry3<-tkentry(tt, textvariable=Stat3)
        stat.entry4<-tkentry(tt, textvariable=Stat4)
        func.entry<-tkentry(tt, textvariable=Func)
        R.entry<-tkentry(tt, textvariable=Rep)
        x.entry<-tkentry(tt, textvariable=Xlab)
        int.entry<-tkentry(tt, textvariable=Int)
        col.entry<-tkentry(tt, textvariable=Col)
        
           
    done <- tclVar(0)
    show.SE<-tclVar(1)
        
        reset <- function(){
            tclvalue(Biv.parent)<-"NULL"
            tclvalue(Parent1)<-"rexp(100)"
            tclvalue(Parent2)<-"NULL"
            tclvalue(S.size1)<-""
            tclvalue(S.size2)<-"NULL"
            tclvalue(Stat1)<-"mean"
            tclvalue(Stat2)<-"NULL"
            tclvalue(Stat3)<-"NULL"
            tclvalue(Stat4)<-"NULL"
            tclvalue(Func)<-"NULL"
            tclvalue(Rep)<-""
            tclvalue(Xlab)<-"expression(Mean)"
            tclvalue(Int)<-"0.1"
            tclvalue(show.SE)<-"0"
            tclvalue(Col)<-"rainbow"
        }

        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)

        build <- function()
        {
            biv.parent <- parse(text=tclvalue(Biv.parent))[[1]]
            parent <-parse(text=tclvalue(Parent1))[[1]]
            parent2 <-parse(text=tclvalue(Parent2))[[1]]
            s.size <-tclvalue(S.size1)
            s.size2 <-parse(text=tclvalue(S.size2))[[1]]
            stat <-parse(text=tclvalue(Stat1))[[1]]
            stat2 <-parse(text=tclvalue(Stat2))[[1]]
            stat3 <-parse(text=tclvalue(Stat3))[[1]]
            stat4 <-parse(text=tclvalue(Stat4))[[1]]
            R <-tclvalue(Rep)
            x<-parse(text=tclvalue(Xlab))[[1]]
            func<-parse(text=tclvalue(Func))[[1]]
            se <- as.logical(tclObj(show.SE))
            interval<-tclvalue(Int)
            col.anim <-tclvalue(Col)
            substitute(samp.dist(biv.parent=biv.parent,parent=parent,parent2=parent2,s.size=as.numeric(s.size),s.size2=s.size2,stat=stat,stat2=stat2,stat3=stat3,stat4=stat4, func=func,R=as.numeric(R),xlab=x,interval=as.numeric(interval),show.SE=se,col.anim=col.anim))
        }
       
        se.cbut <- tkcheckbutton(tt, text="Show SE", variable=show.SE)
        tkgrid(tklabel(tt,text="Sampling distribution animations"),columnspan=4)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="Bivariate parent"), biv.parent.entry)
        tkgrid(tklabel(tt,text="Parent 1",font=c("Helvetica","9","bold")), parent.entry1,tklabel(tt,text="Parent 2"), parent.entry2)
        tkgrid(tklabel(tt,text="Sample size 1",font=c("Helvetica","9","bold")), s.size1.entry,tklabel(tt,text="Sample size 2"), s.size2.entry)
        tkgrid(tklabel(tt,text="Stat",font=c("Helvetica","9","bold")), stat.entry1,tklabel(tt,text="Stat 2"), stat.entry2)
        tkgrid(tklabel(tt,text="Stat 3"), stat.entry3, tklabel(tt,text="Stat 4"), stat.entry4)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="Iterations"), R.entry)
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="Function"), func.entry)
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="Anim. color"), col.entry)
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="Anim. interval"), int.entry)
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="X-axis label"), x.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text=""),se.cbut,sticky="w")
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text=""),submit.but, reset.but)
        tkgrid(tklabel(tt,text=""))
       
        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")

        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
    tclServiceMode(TRUE)
    }                            
if(statc == "custom"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("NULL")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("NULL")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("NULL")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("NULL")
      Ne<-tclVar("NULL")
      Xlab<- tclVar("NULL")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}
if(statc == "mean"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("mean")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(bar(x))")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}
if(statc == "median"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("median")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(Median)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}
if(statc == "trimmed mean"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("NULL")
      tr.mean <- function(x)mean(x,trim = 0.2)
      Stat1<-tclVar("tr.mean")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(Trimmed.mean)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}      
if(statc == "Winsorized mean"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("NULL")
      win.mean <- function(x)mean(win(x))
      Stat1<-tclVar("win.mean")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(Winsorized.mean)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}      
if(statc == "Huber estimator"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("huber.mu")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(Huber.estimator)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}      
if(statc == "H-L estimator"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("HL.mean")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(H-L.estimator)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}  
if(statc == "sd"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("sd")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(S)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}   
if(statc == "var"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("var")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(S^2)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}   
if(statc == "MAD"){

      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("mad")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(MAD))")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}           
if(statc == "IQR"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("10")
      S.size2<-tclVar("10")
      Stat1<-tclVar("IQR")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(IQR)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}    
if(statc == "(n-1)S^2/sigma^2"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("8")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("var")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      xlabel <- expression(paste("(n - 1)", s^2, "/", sigma^2))
      Func<-tclVar("function(s.dist, s.dist2, s.size, s.size2, sigma.sq = 1)((s.size - 1) * s.dist)/sigma.sq")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("xlabel")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}    
if(statc == "F*"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("expression(rnorm(s.size2))")
      S.size1<-tclVar("10")
      S.size2<-tclVar("8")
      Stat1<-tclVar("var")
      Stat2<-tclVar("var")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("function(s.dist,s.dist2,s.size, s.size2) s.dist/s.dist2")
      Rep<-tclVar("2000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(F.star)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}  
if(statc == "t* (1 sample)"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("8")
      S.size2<-tclVar("NULL")
      Stat1<-tclVar("mean")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("var")
      Stat4<-tclVar("NULL")
      Func<-tclVar("function(s.dist, s.dist3, s.size, s.size2)s.dist/(sqrt(s.dist3/s.size))")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(t.star)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}    
if(statc == "t* (2 sample)"){      
      Biv.parent<-tclVar("NULL")
      Parent1 <- tclVar("expression(rnorm(s.size))")
      Parent2 <- tclVar("expression(rnorm(s.size))")
      SS <- tclVar("10")
      SS2<-tclVar("12")
      Stat1<-tclVar("mean")
      Stat2<-tclVar("mean")
      Stat3<-tclVar("var")
      Stat4<-tclVar("var")
      Func<-tclVar("function(s.dist1, s.dist2, s.dist3, s.dist4, s.size = 6, s.size2 = s.size2)({MSE<-(((s.size - 1) * s.dist3) + ((s.size2 - 1) * s.dist4))/((s.size + s.size2) - 2);(s.dist1 - s.dist2)/(sqrt(MSE) * sqrt((1/s.size) + (1/s.size2)))})")
      Rep<-tclVar("10000")
      Xlab<- tclVar("expression(t.star)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}    
if(statc == "Pearson correlation"){
  
      Biv.parent<-tclVar("expression(rmvnorm(s.size, c(0,0), sigma = matrix(nrow=2,ncol=2, data =c(1,0,0,1))))")
      Parent1<-tclVar("NULL")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("8")
      S.size2<-tclVar("8")
      Stat1<-tclVar("NULL")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("cor")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(r)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")} 
if(statc == "covariance"){
      
      Biv.parent<-tclVar("expression(rmvnorm(s.size, c(0,0), sigma = matrix(nrow=2,ncol=2, data =c(1,0,0,1))))")
      Parent1<-tclVar("NULL")
      Parent2<-tclVar("NULL")
      S.size1<-tclVar("8")
      S.size2<-tclVar("8")
      Stat1<-tclVar("NULL")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("cor")
      Rep<-tclVar("1000")
      Ne<-tclVar("seq(1,30)")
      Xlab<- tclVar("expression(Covariance)")
      Col<- tclVar("gray")
      Int<-tclVar("0.01")}      
dialog.sd()
})
}
