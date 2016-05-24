require(rneos) || stop("rneos is not installed")
require(tcltk) || stop("tcltk support is absent")
require(tcltk2) || stop("tcltk2 support is absent")
setClass("hieroClass", S3methods = TRUE)
hieroEnv <- new.env(hash= FALSE, parent = baseenv())
assign("wfile", ".hie", envir = hieroEnv)
attr(hieroEnv$wfile, "class") <- "hieroClass"
assign("wID", ".neo", envir = hieroEnv) 
attr(hieroEnv$wID, "class") <- "hieroClass"
assign("nameOld", "", envir = hieroEnv) 
attr(hieroEnv$nameOld, "class") <- "hieroClass"
assign("jobID_old", "", envir = hieroEnv) 
attr(hieroEnv$jobID_old, "class") <- "hieroClass"
assign("crows", 1, envir = hieroEnv) 
attr(hieroEnv$crows, "class") <- "hieroClass"
assign("nLevs", 1, envir = hieroEnv) 
attr(hieroEnv$nLevs, "class") <- "hieroClass"
assign("rbValue1", "", envir = hieroEnv) 
attr(hieroEnv$rbValue1, "class") <- "hieroClass"
assign("rbValue2", "", envir = hieroEnv) 
attr(hieroEnv$rbValue2, "class") <- "hieroClass"
assign("err", 0, envir = hieroEnv) 
attr(hieroEnv$err, "class") <- "hieroClass"
assign("newProbName", "name", envir = hieroEnv) 
attr(hieroEnv$newProbName, "class") <- "hieroClass"
assign("jobID", "name", envir = hieroEnv) 
attr(hieroEnv$jobID, "class") <- "hieroClass"
assign("FresultsTemp", "name", envir = hieroEnv) 
attr(hieroEnv$FresultsTemp, "class") <- "hieroClass"
assign("jobSource", "name", envir = hieroEnv) 
attr(hieroEnv$jobSource, "class") <- "hieroClass"
assign("res.n1", 0, envir = hieroEnv) 
attr(hieroEnv$res.n1, "class") <- "hieroClass"
assign("res.n2", 0, envir = hieroEnv) 
attr(hieroEnv$res.n2, "class") <- "hieroClass"
assign("res.n3", 0, envir = hieroEnv) 
attr(hieroEnv$res.n3, "class") <- "hieroClass"
assign("res.n4", 0, envir = hieroEnv) 
attr(hieroEnv$res.n4, "class") <- "hieroClass"
assign("res.n5", 0, envir = hieroEnv) 
attr(hieroEnv$res.n5, "class") <- "hieroClass"
assign("res.m1", 0, envir = hieroEnv) 
attr(hieroEnv$res.m1, "class") <- "hieroClass"
assign("res.m2", 0, envir = hieroEnv) 
attr(hieroEnv$res.m2, "class") <- "hieroClass"
assign("res.m3", 0, envir = hieroEnv) 
attr(hieroEnv$res.m3, "class") <- "hieroClass"
assign("res.m4", 0, envir = hieroEnv) 
attr(hieroEnv$res.m4, "class") <- "hieroClass"
assign("res.m5", 0, envir = hieroEnv) 
attr(hieroEnv$res.m5, "class") <- "hieroClass"
assign("res.obj", 0, envir = hieroEnv) 
attr(hieroEnv$res.obj, "class") <- "hieroClass"
assign("res.con", 0, envir = hieroEnv) 
attr(hieroEnv$res.con, "class") <- "hieroClass"
assign("res.groups", 0, envir = hieroEnv) 
attr(hieroEnv$res.groups, "class") <- "hieroClass"
assign("res.alpha", 0, envir = hieroEnv) 
attr(hieroEnv$res.alpha, "class") <- "hieroClass"
assign("res.Delta", 0, envir = hieroEnv) 
attr(hieroEnv$res.Delta, "class") <- "hieroClass"
assign("res.power", 0, envir = hieroEnv) 
attr(hieroEnv$res.power, "class") <- "hieroClass"

delta2<-function(size=size,power=power){
    K1<-function(size=size,power=power){
        z<-qnorm(1-size)
        k<-z-qnorm(1-power)
        k^2}
        {k<-K1(size=size,power=power)
        k<-sqrt(k)
        z<-qnorm(1-size/2)
        delta<-5
        while (abs(delta)>0.01){
            A<-1-pnorm(z-k)+pnorm(-z-k)-power
            B<-dnorm(z-k)-dnorm(-z-k)
            delta<- -A/B
            k<-k+delta}
        k^2}
}
calcPower<-function(alpha=hieroEnv$res.alpha, Delta=hieroEnv$res.Delta, ncpar=hieroEnv$res.con)
    {1-pchisq(qchisq(1-alpha,df=1),df=1,ncp=Delta^2/ncpar)
}    
HierO <- function(){
loadJob <- function() {
  file <- tclvalue(tkgetOpenFile(initialfile=tclvalue(tclfile.tail(hieroEnv$wfile)),initialdir=tclvalue(tclfile.dir(hieroEnv$wfile)),filetypes="{{HierO files} {.hie}} {{All files} *}"))
  if (!length(file)) return()
  chn <- tclopen(file, "r")
  tkdelete(formulas_txt,"0.0","end")
  tkinsert(formulas_txt, "0.0", tclvalue(tclread(chn)))
  tclclose(chn)
  code <- tclvalue(tkget(formulas_txt,"0.0","end"))
  e <- try(parse(text=code))
  if (inherits(e, "try-error")) {
    tkmessageBox(message="Syntax error",icon="error")
    return()
   }
   print(eval(e))
   tkdelete(formulas_txt,"0.0","end")
   tk2notetab.select(nb, "Define")
   hieroEnv$wfile <- file  
}
loadJobID <- function()  {
    nc2 <- CreateNeosComm()
    nameOld <- tclvalue(tkgetOpenFile(initialfile=tclvalue(tclfile.tail(hieroEnv$wID)),initialdir=tclvalue(tclfile.dir(hieroEnv$wID)),filetypes="{{NEOS job ID files} {.neo}} {{All files} *}"))
    if (!length(nameOld)) return()
    load(nameOld)
    jobID_old_sav@nc <- nc2
    jobID_old <- jobID_old_sav
    if(regexpr("/",nameOld)[1] > 0){
        nameOld2 <- substr(nameOld, start=regexpr("/",nameOld)[1]+1, stop=100)
        if(regexpr("/",nameOld2)[1] > 0){
            nameOld2 <- substr(nameOld2, start=regexpr("/",nameOld2)[1]+1, stop=100)
            if(regexpr("/",nameOld2)[1] > 0){
                nameOld2 <- substr(nameOld2, start=regexpr("/",nameOld2)[1]+1, stop=100)
                if(regexpr("/",nameOld2)[1] > 0){
                    nameOld2 <- substr(nameOld2, start=regexpr("/",nameOld2)[1]+1, stop=100)
                    if(regexpr("/",nameOld2)[1] > 0){
                        nameOld2 <- substr(nameOld2, start=regexpr("/",nameOld2)[1]+1, stop=100)
                        if(regexpr("/",nameOld2)[1] > 0){
                            nameOld2 <- substr(nameOld2, start=regexpr("/",nameOld2)[1]+1, stop=100)
                            if(regexpr("/",nameOld2)[1] > 0){
                                nameOld2 <- substr(nameOld2, start=regexpr("/",nameOld2)[1]+1, stop=100)}}}}}}}
    tkdelete(oldJob_txt,"0.0","end")
    tkinsert(oldJob_txt,"end", paste("loaded \n ", nameOld2))  
    print(eval(jobID_old))
    tk2notetab.select(nb, "Results")
    hieroEnv$wID <- nameOld
    hieroEnv$nameOld <- nameOld
    hieroEnv$jobID_old <- jobID_old
}
saveJob <- function() {
    file <- tclvalue(tkgetSaveFile(initialfile=tclvalue(tclfile.tail(hieroEnv$wfile)),initialdir=tclvalue(tclfile.dir(hieroEnv$wfile)),filetypes="{{HierO files} {.hie}} {{All files} *}"))
    if (!length(file)) return()
    chn <- tclopen(file, "w")
    sink(file)
    cat("hieroEnv$nLevs <-", hieroEnv$nLevs, "\n", sep=" ")
    cat("tkconfigure(levs, text=tclVar(hieroEnv$nLevs)) \n")
    rbVal_sav <- as.character(tclvalue(hieroEnv$rbValue1))
    cat("hieroEnv$rbValue1 <- tclVar(", shQuote(rbVal_sav), ") \n",sep="")
    cat("tkconfigure(rbs1,variable=hieroEnv$rbValue1, value=", shQuote("obj.s1"), ")\n", sep="")
    cat("tkconfigure(rbs2,variable=hieroEnv$rbValue1, value=", shQuote("obj.s2"), ")\n", sep="")
    rbVal2_sav <- as.character(tclvalue(hieroEnv$rbValue2))
    cat("hieroEnv$rbValue2 <- tclVar(", shQuote(rbVal2_sav), ") \n",sep="")
    cat("tkconfigure(rb1,variable=hieroEnv$rbValue2, value=", shQuote("obj.ct"), ")\n",sep="")
    cat("tkconfigure(rb2,variable=hieroEnv$rbValue2, value=", shQuote("obj.ce"), ")\n",sep="")
    cat("tkconfigure(rb3,variable=hieroEnv$rbValue2, value=", shQuote("obj.vt"), ")\n",sep="")
    cat("tkconfigure(rb4,variable=hieroEnv$rbValue2, value=", shQuote("obj.ve"), ")\n",sep="")
    SliderValue_sav <- as.numeric(tkget(levsT))
    size_sav <- tkget(size)
    Delta_sav <- tkget(Delta)
    Limpow_sav <- tkget(Limpow)
    LimVC_sav <- tkget(LimVC)    
    cat("if (as.character(tclvalue(hieroEnv$rbValue1)) != ", shQuote("obj.s2"), "){oneSample()}\n", sep="")
    cat("if (as.character(tclvalue(hieroEnv$rbValue1)) == ", shQuote("obj.s2"), "){twoSample()}\n", sep="")   
    cat("if (as.character(tclvalue(hieroEnv$rbValue2)) == ", shQuote("obj.ct"), " || as.character(tclvalue(hieroEnv$rbValue2)) ==", shQuote("obj.ce"), "){minCosts()}\n", sep="")
    cat("if (as.character(tclvalue(hieroEnv$rbValue2)) == ", shQuote("obj.vt"), " || as.character(tclvalue(hieroEnv$rbValue2)) ==", shQuote("obj.ve"), "){minVar()}\n", sep="")     
    cat("SliderValue <- tclVar(", shQuote(SliderValue_sav), ")\n",sep="")
    cat("tkconfigure(levsT, from=1, to=", hieroEnv$nLevs, ", variable=SliderValue)\n",sep="") 
    cat("grLevelConstr()\n")
    cat("tkconfigure(size, text= tclVar(", shQuote(size_sav),"))\n",sep="")
    cat("tkconfigure(Delta, text= tclVar(", shQuote(Delta_sav),"))\n",sep="")
    cat("tkconfigure(Limpow, text= tclVar(", shQuote(Limpow_sav),"))\n",sep="")
    cat("tkconfigure(LimVC, text= tclVar(", shQuote(LimVC_sav),"))\n",sep="")
    cat("hieroEnv$crows <- ", hieroEnv$crows, "\n", sep="")
    cat("if (hieroEnv$crows > 2){tkgrid(LL_3, oper1_3, funct_3, oper2_3, UL_3)} \n")
    cat("if (hieroEnv$crows > 3){tkgrid(LL_4, oper1_4, funct_4, oper2_4, UL_4)} \n")                                        
    cat("if (hieroEnv$crows > 4){tkgrid(LL_5, oper1_5, funct_5, oper2_5, UL_5)} \n")  
    cat("if (hieroEnv$crows > 5){tkgrid(LL_6, oper1_6, funct_6, oper2_6, UL_6)} \n")
    cat("if (hieroEnv$crows > 6){tkgrid(LL_7, oper1_7, funct_7, oper2_7, UL_7)} \n")            
    cat("if (hieroEnv$crows > 7){tkgrid(LL_8, oper1_8, funct_8, oper2_8, UL_8)} \n")              
    cat("if (hieroEnv$crows > 8){tkgrid(LL_9, oper1_9, funct_9, oper2_9, UL_9)} \n")
    cat("if (hieroEnv$crows > 9){tkgrid(LL_10, oper1_10, funct_10, oper2_10, UL_10)} \n")
    LL_1_sav <- tkget(LL_1)
    cat("tkconfigure(LL_1,textvariable=tclVar(", shQuote(LL_1_sav), ")) \n", sep="")
    LL_2_sav <- tkget(LL_2)
    cat("tkconfigure(LL_2,textvariable=tclVar(", shQuote(LL_2_sav), ")) \n", sep="")
    LL_3_sav <- tkget(LL_3)
    cat("tkconfigure(LL_3,textvariable=tclVar(", shQuote(LL_3_sav), ")) \n", sep="")
    LL_4_sav <- tkget(LL_4)
    cat("tkconfigure(LL_4,textvariable=tclVar(", shQuote(LL_4_sav), ")) \n", sep="")
    LL_5_sav <- tkget(LL_5)
    cat("tkconfigure(LL_5,textvariable=tclVar(", shQuote(LL_5_sav), ")) \n", sep="")
    LL_6_sav <- tkget(LL_6)
    cat("tkconfigure(LL_6,textvariable=tclVar(", shQuote(LL_6_sav), ")) \n", sep="")
    LL_7_sav <- tkget(LL_7)
    cat("tkconfigure(LL_7,textvariable=tclVar(", shQuote(LL_7_sav), ")) \n", sep="")
    LL_8_sav <- tkget(LL_8)
    cat("tkconfigure(LL_8,textvariable=tclVar(", shQuote(LL_8_sav), ")) \n", sep="")
    LL_9_sav <- tkget(LL_9)
    cat("tkconfigure(LL_9,textvariable=tclVar(", shQuote(LL_9_sav), ")) \n", sep="")
    LL_10_sav <- tkget(LL_10)
    cat("tkconfigure(LL_10,textvariable=tclVar(", shQuote(LL_10_sav), ")) \n", sep="")
    oper1_1_sav <- tkget(oper1_1)
    cat("tkconfigure(oper1_1,textvariable=tclVar(", shQuote(oper1_1_sav), ")) \n", sep="")
    oper1_2_sav <- tkget(oper1_2)
    cat("tkconfigure(oper1_2,textvariable=tclVar(", shQuote(oper1_2_sav), ")) \n", sep="")
    oper1_3_sav <- tkget(oper1_3)
    cat("tkconfigure(oper1_3,textvariable=tclVar(", shQuote(oper1_3_sav), ")) \n", sep="")
    oper1_4_sav <- tkget(oper1_4)
    cat("tkconfigure(oper1_4,textvariable=tclVar(", shQuote(oper1_4_sav), ")) \n", sep="")
    oper1_5_sav <- tkget(oper1_5)
    cat("tkconfigure(oper1_5,textvariable=tclVar(", shQuote(oper1_5_sav), ")) \n", sep="")
    oper1_6_sav <- tkget(oper1_6)
    cat("tkconfigure(oper1_6,textvariable=tclVar(", shQuote(oper1_6_sav), ")) \n", sep="")
    oper1_7_sav <- tkget(oper1_7)
    cat("tkconfigure(oper1_7,textvariable=tclVar(", shQuote(oper1_7_sav), ")) \n", sep="")
    oper1_8_sav <- tkget(oper1_8)
    cat("tkconfigure(oper1_8,textvariable=tclVar(", shQuote(oper1_8_sav), ")) \n", sep="")
    oper1_9_sav <- tkget(oper1_9)
    cat("tkconfigure(oper1_9,textvariable=tclVar(", shQuote(oper1_9_sav), ")) \n", sep="")
    oper1_10_sav <- tkget(oper1_10)
    cat("tkconfigure(oper1_10,textvariable=tclVar(", shQuote(oper1_10_sav), ")) \n", sep="")
    funct_1_sav <- tkget(funct_1)
    cat("tkconfigure(funct_1,textvariable=tclVar(", shQuote(funct_1_sav), ")) \n", sep="")
    funct_2_sav <- tkget(funct_2)
    cat("tkconfigure(funct_2,textvariable=tclVar(", shQuote(funct_2_sav), ")) \n", sep="")
    funct_3_sav <- tkget(funct_3)
    cat("tkconfigure(funct_3,textvariable=tclVar(", shQuote(funct_3_sav), ")) \n", sep="")
    funct_4_sav <- tkget(funct_4)
    cat("tkconfigure(funct_4,textvariable=tclVar(", shQuote(funct_4_sav), ")) \n", sep="")
    funct_5_sav <- tkget(funct_5)
    cat("tkconfigure(funct_5,textvariable=tclVar(", shQuote(funct_5_sav), ")) \n", sep="")
    funct_6_sav <- tkget(funct_6)
    cat("tkconfigure(funct_6,textvariable=tclVar(", shQuote(funct_6_sav), ")) \n", sep="")
    funct_7_sav <- tkget(funct_7)
    cat("tkconfigure(funct_7,textvariable=tclVar(", shQuote(funct_7_sav), ")) \n", sep="")
    funct_8_sav <- tkget(funct_8)
    cat("tkconfigure(funct_8,textvariable=tclVar(", shQuote(funct_8_sav), ")) \n", sep="")
    funct_9_sav <- tkget(funct_9)
    cat("tkconfigure(funct_9,textvariable=tclVar(", shQuote(funct_9_sav), ")) \n", sep="")
    funct_10_sav <- tkget(funct_10)
    cat("tkconfigure(funct_10,textvariable=tclVar(", shQuote(funct_10_sav), ")) \n", sep="")
    oper2_1_sav <- tkget(oper2_1)
    cat("tkconfigure(oper2_1,textvariable=tclVar(", shQuote(oper2_1_sav), ")) \n", sep="")
    oper2_2_sav <- tkget(oper2_2)
    cat("tkconfigure(oper2_2,textvariable=tclVar(", shQuote(oper2_2_sav), ")) \n", sep="")
    oper2_3_sav <- tkget(oper2_3)
    cat("tkconfigure(oper2_3,textvariable=tclVar(", shQuote(oper2_3_sav), ")) \n", sep="")
    oper2_4_sav <- tkget(oper2_4)
    cat("tkconfigure(oper2_4,textvariable=tclVar(", shQuote(oper2_4_sav), ")) \n", sep="")
    oper2_5_sav <- tkget(oper2_5)
    cat("tkconfigure(oper2_5,textvariable=tclVar(", shQuote(oper2_5_sav), ")) \n", sep="")
    oper2_6_sav <- tkget(oper2_6)
    cat("tkconfigure(oper2_6,textvariable=tclVar(", shQuote(oper2_6_sav), ")) \n", sep="")
    oper2_7_sav <- tkget(oper2_7)
    cat("tkconfigure(oper2_7,textvariable=tclVar(", shQuote(oper2_7_sav), ")) \n", sep="")
    oper2_8_sav <- tkget(oper2_8)
    cat("tkconfigure(oper2_8,textvariable=tclVar(", shQuote(oper2_8_sav), ")) \n", sep="")
    oper2_9_sav <- tkget(oper2_9)
    cat("tkconfigure(oper2_9,textvariable=tclVar(", shQuote(oper2_9_sav), ")) \n", sep="")
    oper2_10_sav <- tkget(oper2_10)
    cat("tkconfigure(oper2_10,textvariable=tclVar(", shQuote(oper2_10_sav), ")) \n", sep="")
    UL_1_sav <- tkget(UL_1)
    cat("tkconfigure(UL_1,textvariable=tclVar(", shQuote(UL_1_sav), ")) \n", sep="")
    UL_2_sav <- tkget(UL_2)
    cat("tkconfigure(UL_2,textvariable=tclVar(", shQuote(UL_2_sav), ")) \n", sep="")
    UL_3_sav <- tkget(UL_3)
    cat("tkconfigure(UL_3,textvariable=tclVar(", shQuote(UL_3_sav), ")) \n", sep="")
    UL_4_sav <- tkget(UL_4)
    cat("tkconfigure(UL_4,textvariable=tclVar(", shQuote(UL_4_sav), ")) \n", sep="")
    UL_5_sav <- tkget(UL_5)
    cat("tkconfigure(UL_5,textvariable=tclVar(", shQuote(UL_5_sav), ")) \n", sep="")
    UL_6_sav <- tkget(UL_6)
    cat("tkconfigure(UL_6,textvariable=tclVar(", shQuote(UL_6_sav), ")) \n", sep="")
    UL_7_sav <- tkget(UL_7)
    cat("tkconfigure(UL_7,textvariable=tclVar(", shQuote(UL_7_sav), ")) \n", sep="")
    UL_8_sav <- tkget(UL_8)
    cat("tkconfigure(UL_8,textvariable=tclVar(", shQuote(UL_8_sav), ")) \n", sep="")
    UL_9_sav <- tkget(UL_9)
    cat("tkconfigure(UL_9,textvariable=tclVar(", shQuote(UL_9_sav), ")) \n", sep="")
    UL_10_sav <- tkget(UL_10)
    cat("tkconfigure(UL_10,textvariable=tclVar(", shQuote(UL_10_sav), ")) \n", sep="")
    vars.1.sav <- tkget(vars.1)
    cat("tkconfigure(vars.1,textvariable=tclVar(", shQuote(vars.1.sav), ")) \n", sep="")
    vars.2.sav <- tkget(vars.2)
    cat("tkconfigure(vars.2,textvariable=tclVar(", shQuote(vars.2.sav), ")) \n", sep="")
    vars.3.sav <- tkget(vars.3)
    cat("tkconfigure(vars.3,textvariable=tclVar(", shQuote(vars.3.sav), ")) \n", sep="")
    vars.4.sav <- tkget(vars.4)
    cat("tkconfigure(vars.4,textvariable=tclVar(", shQuote(vars.4.sav), ")) \n", sep="")
    vars.5.sav <- tkget(vars.5)
    cat("tkconfigure(vars.5,textvariable=tclVar(", shQuote(vars.5.sav), ")) \n", sep="")
    costs.1.sav <- tkget(costs.1)
    cat("tkconfigure(costs.1,textvariable=tclVar(", shQuote(costs.1.sav), ")) \n", sep="")
    costs.2.sav <- tkget(costs.2)
    cat("tkconfigure(costs.2,textvariable=tclVar(", shQuote(costs.2.sav), ")) \n", sep="")
    costs.3.sav <- tkget(costs.3)
    cat("tkconfigure(costs.3,textvariable=tclVar(", shQuote(costs.3.sav), ")) \n", sep="")
    costs.4.sav <- tkget(costs.4)
    cat("tkconfigure(costs.4,textvariable=tclVar(", shQuote(costs.4.sav), ")) \n", sep="")
    costs.5.sav <- tkget(costs.5)
    cat("tkconfigure(costs.5,textvariable=tclVar(", shQuote(costs.5.sav), ")) \n", sep="")
    costs.12.sav <- tkget(costs.12)
    cat("tkconfigure(costs.12,textvariable=tclVar(", shQuote(costs.12.sav), ")) \n", sep="")
    costs.22.sav <- tkget(costs.22)
    cat("tkconfigure(costs.22,textvariable=tclVar(", shQuote(costs.22.sav), ")) \n", sep="")
    costs.32.sav <- tkget(costs.32)
    cat("tkconfigure(costs.32,textvariable=tclVar(", shQuote(costs.32.sav), ")) \n", sep="")
    costs.42.sav <- tkget(costs.42)
    cat("tkconfigure(costs.42,textvariable=tclVar(", shQuote(costs.42.sav), ")) \n", sep="")
    costs.52.sav <- tkget(costs.52)
    cat("tkconfigure(costs.52,textvariable=tclVar(", shQuote(costs.52.sav), ")) \n", sep="")
    problemName_sav <- tkget(problemName)
    optSolver_sav <- tkget(optSolver)
    tolerance_sav <- tkget(tolerance)
    cat("tkconfigure(problemName,textvariable=tclVar(", shQuote(problemName_sav), ")) \n", sep="")
    cat("tkconfigure(optSolver,textvariable=tclVar(", shQuote(optSolver_sav), ")) \n", sep="")
    cat("tkconfigure(tolerance,textvariable=tclVar(", shQuote(tolerance_sav), "))", sep="")
    sink() 
    hieroEnv$wfile <- file 
    tclclose(chn)
}    
saveJobID <- function() {
    nameOld <- tclvalue(tkgetSaveFile(initialfile=tclvalue(tclfile.tail(hieroEnv$wID)),initialdir=tclvalue(tclfile.dir(hieroEnv$wID)),filetypes="{{NEOS job ID files} {.neo}} {{All files} *}"))
    if (!length(nameOld)) return()
    jobID_old_sav <- hieroEnv$jobID
    save(jobID_old_sav, ascii=TRUE, compress = FALSE, file=nameOld)
    hieroEnv$wID <- nameOld
    hieroEnv$nameOld <- nameOld
    tkinsert(model_txt,"end", paste("Job ID saved to ", nameOld))
    tkinsert(model_txt,"end", paste("\n-------------------------------- \n"))   
    print(eval(jobID_old_sav))
}    
helpAbout <- function() {
ttTemp <- tktoplevel(bg = "white")
tktitle(ttTemp) <- "About HierO"
tempTxt <- tktext(ttTemp, bg="white",font=fontText)
tkgrid(tempTxt)
tkinsert(tempTxt,"end", paste("HierO is a tool for optimizing sample size for hierarchical data.\n",
                                "1. Define your setting by filling all of the required fields at the Define tab and press continue. \n",
                                "2. Check at the Final problem tab that the setting is correct \n",
                                "  and send the optimization problem to NEOS server by pressing the button. \n",
                                "3. Get info and solution of the current or previously send job at the Results tab. \n \n",
                                "You can save or load a job at any stage from the file menu. \n",
                                "You can save or load a job id of a job that was send to NEOS and request for the results." ))   
}
helpDef <- function() {
ttTemp <- tktoplevel(bg = "white")
tktitle(ttTemp) <- "Help on Define tab"
yscrt <- tkscrollbar(ttTemp, repeatinterval=5,command=function(...)tkyview(tempTxt,...))
tempTxt <- tktext(ttTemp, bg="white",font=fontText, yscrollcommand=function(...)tkset(yscrt,...),wrap="none",width=80, height=30)
tkgrid(tempTxt,yscrt)
tkgrid.configure(yscrt,rowspan=4,sticky="nsw")
tkinsert(tempTxt,"end", paste("NUMBER OF LEVELS \n",
                                "   -Give the number of hierarchical levels of the data design \n \n",
                                "NUMBER OF GROUPS \n",
                                "   -Select the number of groups (samples) \n \n",
                                "OBJECTIVE \n",
                                "   -Choose the objective for the optimization process \n \n",
                                "CONSTANTS \n",
                                "   -Set the values for required constants \n",
                                "    - Level of grouping; at which level the randomization is done \n",
                                "    - Effect size; minimum difference you want to detect \n",
                                "    - Size; limit for type I error \n",
                                "    - Power; limit for (1 - type II) error \n",
                                "    - Limit for total costs; set the upper limit for the costs of the study \n",
                                "    - Limit for the variance of the estimate; set the upper limit for the variation of the estimate \n \n",
                                "VARIANCE COMPONENTS AND COSTS \n",
                                "   -Give the (positive) estimates of the variance components and costs for each level (and group) \n \n",
                                "DEFINE CONSTRAINTS \n",
                                "   -Set the lower and upper limit for the number of units at each level (and group) as well as\n",
                                "     other possible constraints. \n",
                                "   -Note that each level (for both groups) has to have positive lower limit"))     
}
helpFb <- function() {
ttTemp <- tktoplevel(bg = "white")
tktitle(ttTemp) <- "Help on Final problem tab"
tempTxt <- tktext(ttTemp, bg="white",font=fontText)
tkgrid(tempTxt)
tkinsert(tempTxt,"end", paste("The Final problem tab is divided into 3 sections; left, right and bottom. \n \n",
                                "LEFT SECTION \n",
                                "   -Here is shown the optimization problem defined at the Define tab. \n \n",
                                "RIGHT SECTION \n",
                                "   -Here is shown the optimization problem defined at the Define tab in GAMS form. \n",
                                "   -The code can be modified manually. \n \n",
                                "BOTTOM SECTION \n",
                                "   -Name of the optimization problem; change the name of your problem by \n",
                                "     pressing the Set button.\n", 
                                "   -Select solver; You can choose one of three solvers at the NEOS server to solve your problem. \n",
                                "     BARON and LINDOGlobal are Global optimizers and AlphaECP should also find solution \n", 
                                "     quite close to global optimum. \n", 
                                "   -Relative gap for global optimum; you can set how close to the global optimum the \n",
                                "     iterations search for the solution. Confirm by pressing Set button. \n",
                                "   -Ping NEOS; Check if NEOS server is online. \n",
                                "   -Send optimization problem to NEOS; Sent the optimization problem \n",
                                "     as shown at the RIGHT SECTION to NEOS server \n"))    
}
helpRes <- function() {
ttTemp <- tktoplevel(bg = "white")
tktitle(ttTemp) <- "Help on Results tab" 
tempTxt <- tktext(ttTemp, bg="white",font=fontText)
tkgrid(tempTxt)
tkinsert(tempTxt,"end", paste("The Results tab is divided into 3 sections; CURRENT JOB, SAVED JOB and OUTPUT. \n \n",
                                "CURRENT JOB \n",
                                "   -Current job means the last optimization problem you have send to NEOS since starting HierO \n",
                                "   -You can save the job id returned by NEOS to require results later \n",
                                "   -Get job info from NEOS; NEOS returns info about the optimization problem \n",
                                "   -Get job status from NEOS; NEOS returns the status of the optimization problem \n",
                                "   -Get results from NEOS; Require solution for the optimization problem \n \n",
                                "SAVED JOB \n",
                                "   -You can require for the results of a previously send job with saved job id file \n",
                                "   -Load old job id; load a saved job id earlier returned by NEOS \n",
                                "   -Get old job info from NEOS; NEOS returns info about the job that refers to the loaded old job id \n",
                                "   -Get old job results from NEOS; Require solution for the optimization problem referred \n",
                                "    by the loaded old job id\n \n",
                                "TEXT WINDOW \n",
                                "   -The most important parts of the NEOS output will be shown here"))    
}
helpLevs <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on Levels"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=50, height=4)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Number of Levels: \n Set number of hierarchical levels in data structure" ))   
}
helpGroups <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on Groups"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=50, height=4)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Number of groups: \n Select number of groups or samples" ))   
}
helpObjective <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on Objective"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=70, height=25)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Objective: Choose the objective of the optimization problem \n\n",
                                    "Minimize costs for a given power when testing: \n",
                                    "   Test the difference of location parameter from some constant (one group) or\n",
                                    "   test the difference of two location parameters (two groups). \n",
                                    "   Total costs of the study are minimized for limited power.\n\n",
                                    "Minimize costs for a limited length of confidence interval when estimating: \n",
                                    "   Estimate location parameter (one group) or \n",
                                    "   difference of two location parameters (two groups). \n",
                                    "   Total costs of the study are minimized for limited variance of the estimate.\n\n",
                                    "Maximize power for a limited costs when testing: \n",
                                    "   Test the difference of location parameter from some constant (one group) or\n",
                                    "   test the difference of two location parameters (two groups).\n",
                                    "   Power is maximized for limited costs.\n\n", 
                                    "Minimize length of a confidence interval for limited costs when estimating: \n",
                                    "   Estimate location parameter (one group) or \n",
                                    "   difference of two location parameters (two groups). \n",
                                    "   Length of a confidence interval is minimized for limited costs."))   
}
helpAlpha <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on Size"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=50, height=4)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Size (alpha): \n Limit for 'type-1' error" ))   
}
helpRnd <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on level of grouping"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=50, height=4)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Level of grouping: \n At wich level of hierarchy the randomization is done. \n",
                                    " "))   
}
helpDelta <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on Effect size"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=50, height=4)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Effect size (Delta): \n Minimum difference of sample means you want to detect." ))   
}
helpPower <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on Power"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=50, height=4)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Power (1-beta): \n Lower limit for 1 - 'type-2' error." ))   
}
helpVC <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on Limit for variance of the estimate and Limit for total costs"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=70, height=10)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Depending on the objective: \n\n",
                                    "1. Limit for variance of the estimate: \n    upper limit for the variance of the estimate\n\n",
                                    "2. Limit for total costs: upper limit for total costs of the study" ))   
}
helpVarCo <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on Variance components and costs"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=50, height=10)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Variance component: \n",
                                    "   Set variance component for each level of hierarchy. \n\n",
                                    "Costs for a single unit at group 1: \n",
                                    "   Cost for single unit at each level for the first group. \n\n",
                                    "Costs for a single unit at group 2: \n Cost for single unit at each level for the second group." ))   
}
helpConstr <- function() {
    ttTemp <- tktoplevel(bg = "white")
    tktitle(ttTemp) <- "Help on Constraints"
    tempTxt <- tktext(ttTemp, bg="white",font=fontText,wrap="none",width=60, height=22)
    tkgrid(tempTxt)
    tkinsert(tempTxt,"end", paste("Define constraints for variables. \n\n",
                                    "Lower limit: \n",
                                    "  Lower limit for constraint function, Each variable should have set \n",
                                    "  lower limit > 0 \n",
                                    "Operator 1:\n",
                                    "   Select one of two operators '>=' or '='\n",
                                    "Constraint function:\n",
                                    "   One sample: \n",
                                    "     n1 = lowest level unit \n",
                                    "     n2 = second lowest level unit \n",
                                    "      ... \n",
                                    "   Two samples: \n",
                                    "     n1, n2,... refers to sample 1 and \n",
                                    "     m1, m2,... refers to sample 2 \n",
                                    "   Standard mathematical operators in R can be used \n",
                                    "     (for example (n1+m1)/n3*n4) \n",
                                    "Operator 2:\n",
                                    "   Leave empty or select one of two operators '>=' or '='\n",
                                    "Upper limit: \n",
                                    "  Upper limit for constraint function"))   
}
tt <- tktoplevel(bg = "white")
fontHeading <- tkfont.create(family="serif",size=14,weight="bold")
fontText <- tkfont.create(family="serif",size=11)
fontFormula <- tkfont.create(family="courier",size=12)
fontNotification <- tkfont.create(family="serif",size=11,slant="italic")   
topMenu <- tkmenu(tt)
tkconfigure(tt, menu=topMenu)
tktitle(tt) <- "HierO - Sample size optimization for hierarchical data"
fileMenu1 <- tkmenu(topMenu, tearoff=FALSE)
fileMenu2 <- tkmenu(topMenu, tearoff=FALSE)
fileMenu3 <- tkmenu(topMenu, tearoff=FALSE)
tkadd(topMenu, "cascade", label="File", menu=fileMenu1)
tkadd(topMenu, "cascade", label="Job ID", menu=fileMenu2)
tkadd(topMenu, "cascade", label="Help", menu=fileMenu3)
tkadd(fileMenu1, "command", label="Load problem", command=loadJob)
tkadd(fileMenu2, "command", label="Load job ID", command=loadJobID)
tkadd(fileMenu1, "command", label="Save problem", command=saveJob)
tkadd(fileMenu2, "command", label="Save job ID", command=saveJobID)
tkadd(fileMenu3, "command", label="About", command=helpAbout)
tkadd(fileMenu3, "command", label="Help on 'Define'", command=helpDef)
tkadd(fileMenu3, "command", label="Help on 'Final problem'", command=helpFb)
tkadd(fileMenu3, "command", label="Help on 'Results'", command=helpRes)
nb <- tk2notebook(tt, tabs = c("Define", "Final problem", "Results"),width=1162,height=806)
tkpack(nb, fill = "both", expand = 1)
tb1 <- tk2notetab(nb, "Define")
tb2 <- tk2notetab(nb, "Final problem")
tb3 <- tk2notetab(nb, "Results")
tkconfigure(nb,cursor="hand1") 
frameWhole1  <- tkframe(tb1,relief="groove",borderwidth=5, width=1160, height=804)
frameUpper  <- tkframe(frameWhole1, relief="ridge",borderwidth=5)
frameUpperW <- tkframe(frameUpper,relief="groove",borderwidth=2)
frameUpperC <- tkframe(frameUpper,relief="groove",borderwidth=2)
frameUpperE <- tkframe(frameUpper,relief="groove",borderwidth=2)    
frame2  <- tkframe(frameWhole1, relief="ridge",borderwidth=5)  
frame2W <- tkframe(frame2,relief="groove",borderwidth=2)   
frame2CE <- tkframe(frame2,relief="groove",borderwidth=2)  
frame3  <- tkframe(frameWhole1, relief="ridge",borderwidth=5)
frame3B  <- tkframe(frameWhole1, relief="ridge",borderwidth=0)
frame3W <- tkframe(frame3,relief="groove",borderwidth=2)
frame3W1 <- tkframe(frame3W,relief="groove",borderwidth=0)
frame3W2 <- tkframe(frame3W,relief="groove",borderwidth=0)
frame3C <- tkframe(frame3,relief="groove",borderwidth=0)
frame3E <- tkframe(frame3,relief="groove",borderwidth=0)
tkgrid(frameWhole1, padx=5, pady=5, sticky="wens")
tkgrid(frameUpper, padx=10, pady=5, sticky="we", columnspan=2)
tkgrid(frameUpperW, frameUpperC, frameUpperE, ipadx=4, ipady=4, padx=10, pady=5, sticky="ns")
tkgrid(frame2, padx=10, pady=5, sticky="we", columnspan=2)
tkgrid(frame2W, frame2CE, ipadx=4, ipady=4, padx=10, pady=5, sticky="ns")
tkgrid(frame3, frame3B, padx=10, pady=5, sticky="wesn")
tkgrid(frame3W, frame3C, frame3E, ipadx=4, ipady=4, padx=10, pady=5, sticky="wesn")
tkgrid(frame3W1,frame3W2, sticky="wesn")
frameWhole3  <- tkframe(tb2,relief="groove",borderwidth=5)
frameWhole3U  <- tkframe(frameWhole3,relief="ridge",borderwidth=0)
frameWhole3W  <- tkframe(frameWhole3U,relief="groove",borderwidth=1)
frameWhole3C  <- tkframe(frameWhole3U,relief="groove",borderwidth=1)
frameWhole3D  <- tkframe(frameWhole3,relief="ridge",borderwidth=0)
frameWhole3D1  <- tkframe(frameWhole3D,relief="groove",borderwidth=0)
frameWhole3D2  <- tkframe(frameWhole3D,relief="groove",borderwidth=0)
frameWhole3D3  <- tkframe(frameWhole3D,relief="groove",borderwidth=0)
frameWhole3D4  <- tkframe(frameWhole3D,relief="groove",borderwidth=0)
frameWhole3D5  <- tkframe(frameWhole3D,relief="groove",borderwidth=0)
tkgrid(frameWhole3, padx=5, pady=5)
tkgrid(frameWhole3U, padx=1, pady=1, sticky="we")
tkgrid(frameWhole3W, frameWhole3C, padx=1, pady=1, sticky="ns")
tkgrid(frameWhole3D, sticky="we")
tkgrid(frameWhole3D1, frameWhole3D2, frameWhole3D3, frameWhole3D4, frameWhole3D5, padx=1, pady=1)
frameWhole4  <- tkframe(tb3,relief="groove",borderwidth=5)
frameWhole4W <- tkframe(frameWhole4,relief="groove",borderwidth=2)
frameWhole4E <- tkframe(frameWhole4,relief="groove",borderwidth=2)
frameWhole40W <- tkframe(frameWhole4W,relief="groove",borderwidth=0)
frameWhole41W <- tkframe(frameWhole4W,relief="ridge",borderwidth=4)
frameWhole42W <- tkframe(frameWhole4W,relief="ridge",borderwidth=4)
tkgrid(frameWhole4, padx=5, pady=5, sticky="nswe")
tkgrid(frameWhole4W, frameWhole4E, padx=5, pady=5, sticky="ns")
tkgrid(frameWhole40W, ipady=5, ipadx=0, padx=0, pady=c(10,5))
tkgrid(frameWhole41W, ipady=10, ipadx=5, padx=5, pady=10)
tkgrid(frameWhole42W, ipady=5, ipadx=5, padx=5, pady=5)
levelsLabel <- tklabel(frameUpperW,text="Number of levels",font=fontHeading)
helpButLevs <- tkbutton(frameUpperW,text="?", relief="solid",borderwidth=1,borderwidth=1, command=function() helpLevs())
tk2tip(helpButLevs, "Number of hierarchical levels \n in data structure") 
levs <- tk2spinbox(frameUpperW, from = 1, to = 5, increment = 1, width=5, text=tclVar("0"),justify='right',command=function() changeLevs())
hieroEnv$nLevs <- as.integer(tkget(levs))
labelText <- tclVar("Select number of levels! ")
label1 <- tklabel(frameUpperW,text=paste(tclvalue(labelText)),fg="red",font=fontNotification)
level_guide11_value <- tclVar("")
level_guide11 <- tklabel(frameUpperW,text=paste(tclvalue(level_guide11_value)),fg="darkgreen",font=fontNotification)
level_guide1mid_v <- tclVar("")
level_guide1mid <- tklabel(frameUpperW,text=paste(tclvalue(level_guide1mid_v)),fg="darkgreen",font=fontNotification)
level_guide1K_value <- tclVar("")
level_guide1K <- tklabel(frameUpperW,text=paste(tclvalue(level_guide1K_value)),fg="darkgreen",font=fontNotification)
tkgrid(levelsLabel,helpButLevs,padx=10) 
tkgrid(levs, padx=5, pady=5,columnspan=2)          
tkgrid(label1,padx=5, pady=5,columnspan=2)
tkgrid(level_guide11,padx=5, pady=5,columnspan=2)
tkgrid(level_guide1mid,padx=5, pady=1,columnspan=2)
tkgrid(level_guide1K,padx=5, pady=5,columnspan=2)
tkconfigure(label1,textvariable=labelText)
tkconfigure(level_guide11,textvariable=level_guide11_value)
tkconfigure(level_guide1mid,textvariable=level_guide1mid_v)
tkconfigure(level_guide1K,textvariable=level_guide1K_value)
samplesLabel <- tklabel(frameUpperC,text="Number of groups",font=fontHeading)
helpButGr <- tkbutton(frameUpperC,text="?", relief="solid",borderwidth=1, command=function() helpGroups())
tk2tip(helpButGr, "Select number of groups") 
rbs1 <- tkradiobutton(frameUpperC,padx=5)
rbs2 <- tkradiobutton(frameUpperC)
hieroEnv$rbValue1 <- tclVar("obj.s1")
onesample <- tklabel(frameUpperC,text="One group",font=fontText)
twosample <- tklabel(frameUpperC,text="Two groups",font=fontText)
level_guide21_value <- tclVar("")
level_guide21 <- tklabel(frameUpperC,text=paste(tclvalue(level_guide21_value)),fg="darkgreen",font=fontNotification)
level_guide2mid_v <- tclVar("")
level_guide2mid <- tklabel(frameUpperC,text=paste(tclvalue(level_guide2mid_v)),fg="darkgreen",font=fontNotification)
level_guide2K_value <- tclVar("")
level_guide2K <- tklabel(frameUpperC,text=paste(tclvalue(level_guide2K_value)),fg="darkgreen",font=fontNotification)
tkconfigure(rbs1,variable=hieroEnv$rbValue1,value="obj.s1", command=function() oneSample())
tkconfigure(rbs2,variable=hieroEnv$rbValue1,value="obj.s2", command=function() twoSample())
tkgrid(samplesLabel, helpButGr, padx=30, columnspan=2)
tkgrid(onesample,rbs1, padx=2, pady=5)
tkgrid(twosample, rbs2, padx=2, pady=5)
tkgrid(level_guide21,padx=5, pady=5,columnspan=3)
tkgrid(level_guide2mid,padx=5, pady=1,columnspan=3)
tkgrid(level_guide2K,padx=5, pady=1,columnspan=3)
tkconfigure(level_guide21,textvariable=level_guide21_value)
tkconfigure(level_guide2mid,textvariable=level_guide2mid_v)
tkconfigure(level_guide2K,textvariable=level_guide2K_value) 
objLabel <- tklabel(frameUpperE,text="Objective",font=fontHeading)
helpButObj <- tkbutton(frameUpperE,text="?", relief="solid",borderwidth=1, command=function() helpObjective())
tk2tip(helpButObj, "Choose the objective of the optimization problem \n\n Estimate location parameter (one group)/\n difference of two location parameters (two groups) or \n\n
test the difference of location parameter from some constant (one group)/\ntest the difference of two location parameters (two groups)") 
rb1 <- tkradiobutton(frameUpperE,padx=3)
rb2 <- tkradiobutton(frameUpperE)
rb3 <- tkradiobutton(frameUpperE,padx=3)
rb4 <- tkradiobutton(frameUpperE)
hieroEnv$rbValue2 <- tclVar("obj")
ObjCt <- tklabel(frameUpperE,text="Minimize costs for a given power \n when testing",font=fontText)
ObjCe <- tklabel(frameUpperE,text="Minimize costs for a limited length of \n confidence interval when estimating",font=fontText)
ObjVt <- tklabel(frameUpperE,text="Maximize power for a limited costs \n when testing ",font=fontText)
ObjVe <- tklabel(frameUpperE,text="Minimize length of a confidence interval for \n limited costs when estimating ",font=fontText)
tkconfigure(rb1,variable=hieroEnv$rbValue2,value="obj.ct", command=function() minCosts())
tkconfigure(rb2,variable=hieroEnv$rbValue2,value="obj.ce", command=function() minCosts())
tkconfigure(rb3,variable=hieroEnv$rbValue2,value="obj.vt", command=function() minVar())
tkconfigure(rb4,variable=hieroEnv$rbValue2,value="obj.ve", command=function() minVar())
tkgrid(objLabel, helpButObj)
tkgrid(ObjCt,rb1, padx=10, pady=3)
tkgrid(ObjCe,rb2, padx=10, pady=3)
tkgrid(ObjVt,rb3, padx=10, pady=3)
tkgrid(ObjVe,rb4, padx=10, pady=3)
sizeLabel <- tklabel(frame2W,text="Size (alpha)",font=fontText)
helpButAlpha <- tkbutton(frame2W,text="?", relief="solid",borderwidth=1, command=function() helpAlpha())
tk2tip(helpButAlpha, "Limit for 'type-1' error") 
size <- tk2spinbox(frame2W, from = 0.01, to = 0.99, increment = 0.01, width=7, text=tclVar("0.05"),justify='right')  
DeltaLabel <- tklabel(frame2W,text="Effect size (Delta)",font=fontText)                       
helpButRnd <- tkbutton(frame2W,text="?", relief="solid",borderwidth=1, command=function() helpRnd())
tk2tip(helpButRnd, "At wich level of hierarchy the randomization is done.") 
Delta <- tk2spinbox(frame2W, from = 0, to = 100, increment = 0.05, width=7,justify='right') 
helpButDelta <- tkbutton(frame2W,text="?", relief="solid",borderwidth=1, command=function() helpDelta())
tk2tip(helpButDelta, "Minimum difference you want to detect") 
SliderValue <- tclVar(hieroEnv$nLevs)
levsGrLabel <- tklabel(frame2W,text="Level of grouping",font=fontText)
levsT <- tkscale(frame2W, from=1, to=hieroEnv$nLevs, showvalue=T, variable=SliderValue, resolution=1, orient="horizontal")
helpButPow <- tkbutton(frame2W,text="?", relief="solid",borderwidth=1, command=function() helpPower())
tk2tip(helpButPow, "Limit for 1 - 'type-2' error")
powLabel <- tklabel(frame2W,text="Power (1-beta)",font=fontText)
Limpow <- tk2spinbox(frame2W, from = 0.01, to = 0.99, increment = 0.01, width=7, text=tclVar("0.80"),justify='right')  
LimVCLabel <- tklabel(frame2W,text="",font=fontText)
helpButVC <- tkbutton(frame2W,text="?", relief="solid",borderwidth=1, command=function() helpVC())
tk2tip(helpButVC, "Depending on the objective \n the upper limit for total costs \n or the variance of the estimate")
LimVC <- tkentry(frame2W,width="10",justify='right')
tkgrid(levsGrLabel, levsT, helpButRnd, padx=26, pady=5)
tkgrid(DeltaLabel, Delta, helpButDelta, padx=10, pady=5)
tkgrid(sizeLabel, size, helpButAlpha, padx=10, pady=5)
tkgrid(powLabel, Limpow, helpButPow, padx=10, pady=5)
tkgrid(LimVCLabel, LimVC, helpButVC, padx=10, pady=5)
tkconfigure(levsGrLabel, fg="grey")
tkconfigure(levsT, fg="grey")
tkconfigure(Delta,fg="grey", bg="grey")
tkconfigure(DeltaLabel, fg="grey")
tkconfigure(sizeLabel, fg="grey")
tkconfigure(size, fg="grey", bg="grey")
tkconfigure(powLabel,fg="grey")
tkconfigure(Limpow,fg="grey", bg="grey")
tkconfigure(LimVCLabel,fg="grey")
tkconfigure(LimVC,fg="grey",bg = "grey")
varcostLabel <- tklabel(frame2CE,text="   Define variance components and single unit cost at each level   ",font=fontHeading)
helpButVC <- tkbutton(frame2CE,text="?", relief="solid",borderwidth=1, command=function() helpVarCo())
tk2tip(helpButVC, "Set variance component for each level of hierarchy and \n cost for single unit at each level") 
varscostsHl <- tklabel(frame2CE,text="Level",font=fontText)
varscostsHv <- tklabel(frame2CE,text="Variance \n component",font=fontText)
varscostsHc1 <- tklabel(frame2CE,text="Cost for single unit \n at group 1",font=fontText)
varscostsHc2 <- tklabel(frame2CE,text="Cost for single unit \n at group 2",font=fontText)
varscosts1 <- tklabel(frame2CE,text="Level 1 (lowest)",fg="grey",font=fontText)
vars.1 <- tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.1 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.12 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
varscosts2 <- tklabel(frame2CE,text="Level 2",fg="grey",font=fontText)
vars.2 <- tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.2 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.22 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
varscosts3 <- tklabel(frame2CE,text="Level 3",fg="grey",font=fontText)
vars.3 <- tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.3 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.32 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
varscosts4 <- tklabel(frame2CE,text="Level 4",fg="grey",font=fontText)
vars.4 <- tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.4 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.42 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
varscosts5 <- tklabel(frame2CE,text="Level 5",fg="grey",font=fontText)
vars.5 <- tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.5 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
costs.52 <-tkentry(frame2CE,width="7",textvariable=tclVar("-1"),justify='right',fg="grey",bg = "grey")
tkgrid(varcostLabel, helpButVC, padx=6, columnspan=4)
tkgrid(varscostsHl,varscostsHv,varscostsHc1,varscostsHc2,padx=4, pady=5)
tkgrid(varscosts1, vars.1, costs.1, costs.12, padx=1, pady=1)
tkgrid(varscosts2, vars.2, costs.2, costs.22, padx=1, pady=1)
tkgrid(varscosts3, vars.3, costs.3, costs.32, padx=1, pady=1)
tkgrid(varscosts4, vars.4, costs.4, costs.42, padx=1, pady=1)
tkgrid(varscosts5, vars.5, costs.5, costs.52, padx=1, pady=1)        
costrLabel <- tklabel(frame3W1,text="Define constraints",font=fontHeading)
helpButCon <- tkbutton(frame3W1,text="?", relief="solid",borderwidth=1, command=function() helpConstr())
tk2tip(helpButCon, "Define constraints for variables. \n Each variable should have lower limit > 0. \n (Default upper limit for each variable is 100) \n\n Constraint function: \n\n One sample: \n n1 = lowest level unit \n n2 = second lowest level unit \n ... \n
Two samples: \n n1, n2,... refers to sample 1 and \n m1, m2,... refers to sample 2 \n
Standard mathematical operators in R can be used \n (for example (n1+m1)/n3*n4) ") 
tkgrid(costrLabel,helpButCon,columnspan=5)
Operators <- c("","<=","=")     
head1 <- tklabel(frame3W1,text="Lower limit ",font=fontText) 
head2 <- tklabel(frame3W1,text="Operator 1 ",font=fontText) 
head3 <- tklabel(frame3W1,text="Constraint function ",font=fontText) 
head4 <- tklabel(frame3W1,text="Operator 2 ",font=fontText) 
head5 <- tklabel(frame3W1,text="Upper Limit ",font=fontText)
tclRequire("BWidget")
LL_1 <- tkentry(frame3W1,width="12",justify='right')
oper1_1 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_1 <- tkentry(frame3W1,width="30",justify='right')
oper2_1 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_1 <- tkentry(frame3W1,width="12",justify='right')
LL_2 <- tkentry(frame3W1,width="12",justify='right')
oper1_2 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_2 <- tkentry(frame3W1,width="30",justify='right')
oper2_2 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_2 <- tkentry(frame3W1,width="12",justify='right') 
LL_3 <- tkentry(frame3W1,width="12",justify='right')
oper1_3 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_3 <- tkentry(frame3W1,width="30",justify='right')
oper2_3 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_3 <- tkentry(frame3W1,width="12",justify='right') 
LL_4 <- tkentry(frame3W1,width="12",justify='right')
oper1_4 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_4 <- tkentry(frame3W1,width="30",justify='right')
oper2_4 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_4 <- tkentry(frame3W1,width="12",justify='right') 
LL_5 <- tkentry(frame3W1,width="12",justify='right')
oper1_5 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_5 <- tkentry(frame3W1,width="30",justify='right')
oper2_5 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_5 <- tkentry(frame3W1,width="12",justify='right') 
LL_6 <- tkentry(frame3W1,width="12",justify='right')
oper1_6 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_6 <- tkentry(frame3W1,width="30",justify='right')
oper2_6 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_6 <- tkentry(frame3W1,width="12",justify='right') 
LL_7 <- tkentry(frame3W1,width="12",justify='right')
oper1_7 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_7 <- tkentry(frame3W1,width="30",justify='right')
oper2_7 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_7 <- tkentry(frame3W1,width="12",justify='right') 
LL_8 <- tkentry(frame3W1,width="12",justify='right')
oper1_8 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_8 <- tkentry(frame3W1,width="30",justify='right')
oper2_8 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_8 <- tkentry(frame3W1,width="12",justify='right') 
LL_9 <- tkentry(frame3W1,width="12",justify='right')
oper1_9 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_9 <- tkentry(frame3W1,width="30",justify='right')
oper2_9 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_9 <- tkentry(frame3W1,width="12",justify='right') 
LL_10 <- tkentry(frame3W1,width="12",justify='right')
oper1_10 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
funct_10 <- tkentry(frame3W1,width="30",justify='right')
oper2_10 <- tkwidget(frame3W1,"ComboBox",values=Operators,width="8")
UL_10 <- tkentry(frame3W1,width="12",justify='right')   
delete.but <- tkbutton(frame3W2,text="Delete",command=function() deleteRowC())
add.but <- tkbutton(frame3W2,text="Add",command=function() addRowC()) 
e.constraints <- tklabel(frame3E,text="Extra constraints",fg="grey", font=fontHeading)
c.extra.l1t <- tclVar("")
c.extra.l2t <- tclVar("")
c.extra.l3t <- tclVar("")
c.extra.l4t <- tclVar("")
c.extra.l1 <- tklabel(frame3E,text=paste(tclvalue(c.extra.l1t)),fg="darkgreen",font=fontNotification)
c.extra.l2 <- tklabel(frame3E,text=paste(tclvalue(c.extra.l2t)),fg="darkgreen",font=fontNotification)
c.extra.l3 <- tklabel(frame3E,text=paste(tclvalue(c.extra.l3t)),fg="darkgreen",font=fontNotification)
c.extra.l4 <- tklabel(frame3E,text=paste(tclvalue(c.extra.l4t)),fg="darkgreen",font=fontNotification)
tkgrid(head1, head2, head3,head4,head5)
tkgrid(LL_1, oper1_1, funct_1, oper2_1, UL_1, padx=4)
tkgrid(add.but, ipadx=8, ipady=4, padx=4, pady=4, sticky="swe")
tkgrid(delete.but, ipadx=8, ipady=4, padx=4, pady=4, sticky="swe")
tkgrid(e.constraints)
tkgrid(c.extra.l1)
tkgrid(c.extra.l2)
tkgrid(c.extra.l3)
tkgrid(c.extra.l4)
problemLoad <- tkbutton(frame3B,font=fontHeading,bg = "white", height=2, text=" Load \n problem ",command=loadJob,font=fontHeading)
problemSave <- tkbutton(frame3B,font=fontHeading,bg = "white", height=2, text=" Save \n problem ",command=saveJob,font=fontHeading)
def1 <- tkbutton(frame3B,font=fontHeading,bg = "green", height=2, text=" Continue ",command=function()DefGams(),font=fontHeading)
tkgrid(problemLoad, padx=10, pady=15, sticky="we")
tkgrid(problemSave, padx=10, pady=15, sticky="we")
tkgrid(def1, padx=10, pady=15, sticky="we")
HeadProb <- tklabel(frameWhole3W,text="Problem",font=fontHeading)
f.funct.lt <- tclVar("")
f.funct.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.funct.lt)),fg="black",font=fontText)
f.funct.ct <- tclVar("")
f.funct.c <- tklabel(frameWhole3W,text=paste(tclvalue(f.funct.ct)),fg="black",font=fontText)
f.funct.ot <- tclVar("")
f.funct.o <- tklabel(frameWhole3W,text=paste(tclvalue(f.funct.ot)),fg="black",font=fontText)
tkgrid(HeadProb, columnspan=4)
tkgrid(tklabel(frameWhole3W, text="  ----------------------------------------------------------  "), columnspan=4)
tkgrid(f.funct.l, columnspan=4)
tkgrid(f.funct.o, columnspan=4)
tkgrid(f.funct.c, columnspan=4)
HeadConst <- tklabel(frameWhole3W,text="Constants",font=fontHeading)
f.1.lt <- tclVar("") 
f.1.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.1.lt)),fg="black",font=fontText)
f.2.lt <- tclVar("") 
f.2.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.2.lt)),fg="black",font=fontText)
f.3.lt <- tclVar("") 
f.3.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.3.lt)),fg="black",font=fontText)
f.4.lt <- tclVar("") 
f.4.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.4.lt)),fg="black",font=fontText)
tkgrid(tklabel(frameWhole3W, text="  ----------------------------------------------------------  "), columnspan=4)
tkgrid(HeadConst, columnspan=4) 
tkgrid(tklabel(frameWhole3W, text="  ----------------------------------------------------------  "), columnspan=4)
tkgrid(f.1.l, columnspan=4)
tkgrid(f.2.l, columnspan=4)
tkgrid(f.3.l, columnspan=4)
tkgrid(f.4.l, columnspan=4)
tkgrid(tklabel(frameWhole3W, text="  ----------------------------------------------------------  "), columnspan=4)
f.lev.l <- tklabel(frameWhole3W,fg="black",font=fontText)
f.var.l <- tklabel(frameWhole3W,fg="black",font=fontText)
f.c1.lt <- tclVar("")
f.c1.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.c1.lt)),fg="black",font=fontText)
f.c2.lt <- tclVar("")
f.c2.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.c2.lt)),fg="black",font=fontText)
tkgrid(f.lev.l,f.var.l, f.c1.l, f.c2.l)
f.lev1.l <- tklabel(frameWhole3W,fg="black",font=fontText)
f.var1.lt <- tclVar("")
f.var1.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.var1.lt)),fg="black",font=fontText)
f.cost1.lt <- tclVar("")
f.cost1.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost1.lt)),fg="black",font=fontText)
f.cost12.lt <- tclVar("")
f.cost12.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost12.lt)),fg="black",font=fontText)
f.lev2.l <- tklabel(frameWhole3W,fg="black",font=fontText)
f.var2.lt <- tclVar("")
f.var2.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.var2.lt)),fg="black",font=fontText)
f.cost2.lt <- tclVar("")
f.cost2.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost2.lt)),fg="black",font=fontText)
f.cost22.lt <- tclVar("")
f.cost22.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost22.lt)),fg="black",font=fontText)
f.lev3.l <- tklabel(frameWhole3W,fg="black",font=fontText)
f.var3.lt <- tclVar("")
f.var3.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.var3.lt)),fg="black",font=fontText)
f.cost3.lt <- tclVar("")
f.cost3.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost3.lt)),fg="black",font=fontText)
f.cost32.lt <- tclVar("")
f.cost32.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost32.lt)),fg="black",font=fontText)
f.lev4.l <- tklabel(frameWhole3W,fg="black",font=fontText)
f.var4.lt <- tclVar("")
f.var4.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.var4.lt)),fg="black",font=fontText)
f.cost4.lt <- tclVar("")
f.cost4.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost4.lt)),fg="black",font=fontText)
f.cost42.lt <- tclVar("")
f.cost42.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost42.lt)),fg="black",font=fontText)
f.lev5.l <- tklabel(frameWhole3W,fg="black",font=fontText)
f.var5.lt <- tclVar("")
f.var5.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.var5.lt)),fg="black",font=fontText)
f.cost5.lt <- tclVar("")
f.cost5.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost5.lt)),fg="black",font=fontText)
f.cost52.lt <- tclVar("")
f.cost52.l <- tklabel(frameWhole3W,text=paste(tclvalue(f.cost52.lt)),fg="black",font=fontText)
tkgrid(f.lev1.l, f.var1.l, f.cost1.l, f.cost12.l)
tkgrid(f.lev2.l, f.var2.l, f.cost2.l, f.cost22.l)
tkgrid(f.lev3.l, f.var3.l, f.cost3.l, f.cost32.l)
tkgrid(f.lev4.l, f.var4.l, f.cost4.l, f.cost42.l)
tkgrid(f.lev5.l, f.var5.l, f.cost5.l, f.cost52.l)
gamsLabel <- tklabel(frameWhole3C, text="Optimization problem in GAMS form",font=fontText)
yscr3 <- tkscrollbar(frameWhole3C, repeatinterval=5,command=function(...)tkyview(formulas_txt,...))
formulas_txt <- tktext(frameWhole3C, bg="white",font=fontFormula, yscrollcommand=function(...)tkset(yscr3,...),wrap="none",width=80, height=38)
tkgrid(gamsLabel)
tkgrid(formulas_txt,yscr3)
tkgrid.configure(yscr3,rowspan=4,sticky="nsw")
probNameLabel <- tklabel(frameWhole3D1, text="Name of optimization problem",font=fontText)
problemName <- tkentry(frameWhole3D1,width="25",justify='right')
tkconfigure(problemName,textvariable=tclVar("opt_problem1"))
nameAcc <- tkbutton(frameWhole3D1,text="Set",command=function() DefGams(),font=fontText)
optSolverChoices <- c("BARON","AlphaECP","LINDOGlobal")
optSolverLabel <- tklabel(frameWhole3D2, text="Select solver",font=fontText)
optSolver <- tkwidget(frameWhole3D2,"ComboBox",editable=FALSE,values=optSolverChoices, text="BARON", width=20,justify='right')
tkconfigure(optSolver,values=optSolverChoices)
toleranceLabel <- tklabel(frameWhole3D3,text="Relative gap for global optimum",font=fontText)
tolerance <- tk2spinbox(frameWhole3D3, from = 0.0001, to = 0.50, increment = 0.00001, width=10, text=tclVar("0.0001"),justify='right',font=fontText)
tolAcc <- tkbutton(frameWhole3D3,text="Set",command=function() DefGams(),font=fontText)
problemSave2 <- tkbutton(frameWhole3D4,font=fontHeading,bg = "white", height=2, text=" Save \n problem",command=saveJob)
checkNeos <- tkbutton(frameWhole3D4,font=fontHeading,bg = "white", height=2, text=" Ping \n NEOS", command=function() pingNeos())
neosStatus <- tclVar("                     ")
neosStatusL <- tklabel(frameWhole3D4,text=paste(tclvalue(neosStatus)),fg="darkgreen",font=fontNotification)
toNeosBut <- tkbutton(frameWhole3D5, bg = "green", text="Send optimization \n problem to NEOS",font=fontHeading, command=function() sentToNEOS())
tkgrid(probNameLabel, padx=5, pady=5, columnspan=2)
tkgrid(problemName,nameAcc, padx=5, pady=5)
tkgrid(optSolverLabel, padx=15, pady=5)
tkgrid(optSolver, padx=15, pady=5)
tkgrid(toleranceLabel, padx=15, pady=5, columnspan=2)
tkgrid(tolerance, tolAcc, padx=0, pady=5)
tkgrid(neosStatusL,checkNeos, problemSave2,padx=5, pady=5, sticky="e") #rowspan=4, 
tkgrid(toNeosBut, ipadx=8, ipady=4, padx=4, pady=4,sticky="ens",rowspan=2)
tkconfigure(neosStatusL,textvariable=neosStatus)

neosLabel <- tklabel(frameWhole40W, text="Communication with NEOS",font=fontHeading)
tkgrid(neosLabel, sticky="wn")
neosNew <- tklabel(frameWhole41W, text="Current job",font=fontHeading,justify='center')
storeJobID <- tkbutton(frameWhole41W,text="Save job ID",command=function() saveJobID())
jobInfo <- tkbutton(frameWhole41W, text="Get job info \n from NEOS", command=function() getJobInfo())
jobStatus <- tkbutton(frameWhole41W, text="Get job status \n from NEOS", command=function() getJobStatus())
jobResults <- tkbutton(frameWhole41W, text="Get results \n from NEOS", command=function() getResults())
tkgrid(tklabel(frameWhole41W, text="  -----------------------------------------  "))
tkgrid(neosNew, sticky="wn", columnspan=2, sticky="we")
tkgrid(tklabel(frameWhole41W, text="  -----------------------------------------  "))
tkgrid(storeJobID, ipadx=5, ipady=5, padx=5, pady=5, sticky="swe")
tkgrid(jobInfo, ipadx=5, ipady=5, padx=5, pady=5, sticky="swe")
tkgrid(jobStatus, ipadx=5, ipady=5, padx=5, pady=5, sticky="swe")
tkgrid(jobResults, ipadx=5, ipady=5, padx=5, pady=5, sticky="swe")
neosOldLabel <- tklabel(frameWhole42W, text="Saved jobs",font=fontHeading,justify='center')
getOldJobID <- tkbutton(frameWhole42W,text="Load old job ID",command=function() loadJobID())
oldJob_txt <- tktext(frameWhole42W, bg="white",fg="darkgreen", font=fontText,wrap="none",width=8, height=3)
jobInfo_old <- tkbutton(frameWhole42W, text="Get old job \n info from NEOS", command=function() getJobInfo_old())
jobResults_old <- tkbutton(frameWhole42W, text="Get old job \n results from NEOS", command=function() getJobResults_old())
tkgrid(tklabel(frameWhole42W, text="  -----------------------------------------  "))
tkgrid(neosOldLabel, sticky="wn", columnspan=2, sticky="we")
tkgrid(tklabel(frameWhole42W, text="  -----------------------------------------  "))
tkgrid(getOldJobID, ipadx=5, ipady=5, padx=5, pady=5, sticky="s", sticky="we")
tkgrid(oldJob_txt, ipadx=5, ipady=5, padx=5, pady=5, sticky="s", sticky="we")
tkgrid(jobInfo_old, ipadx=5, ipady=5, padx=5, pady=5, sticky="s", sticky="we")
tkgrid(jobResults_old, ipadx=5, ipady=5, padx=5, pady=5, sticky="s", sticky="we")
glob.l <- tklabel(frameWhole4E, text="Optimization problem will be solved as Global Optimization problem with GAMS modelling language")
yscr <- tkscrollbar(frameWhole4E, repeatinterval=5,command=function(...)tkyview(model_txt,...))
model_txt <- tktext(frameWhole4E, bg="white",font=fontFormula, yscrollcommand=function(...)tkset(yscr,...),wrap="none",width=84, height=41)
tkgrid(glob.l)
tkgrid(model_txt,yscr)
tkgrid.configure(yscr,rowspan=4,sticky="nsw")

changeLevs <- function(){
    hieroEnv$nLevs <- as.integer(tkget(levs))  
    tkconfigure(varscosts1,fg="black")
    tkconfigure(vars.1,fg="black",bg = "white")
    tkconfigure(costs.1,fg="black",bg = "white")
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2") tkconfigure(costs.12,fg="black",bg = "white")
        tclvalue(labelText) <- paste("You have selected ", hieroEnv$nLevs," Levels", sep="")
        tkconfigure(label1,fg="darkgreen")
    if (hieroEnv$nLevs > 1) tclvalue(level_guide11_value) <- paste("Highest level unit = n",hieroEnv$nLevs, sep="")
    else tclvalue(level_guide11_value) <- paste("Lowest level unit = n1", sep="")
    if (hieroEnv$nLevs > 2) tclvalue(level_guide1mid_v) <- paste(". . .")
    else  tclvalue(level_guide1mid_v) <- paste("")
    if (hieroEnv$nLevs > 1) tclvalue(level_guide1K_value) <- paste("Lowest level unit = n1", sep="")
    else tclvalue(level_guide1K_value) <- paste("")
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
        tclvalue(level_guide21_value) <- paste("")
        tclvalue(level_guide2mid_v) <- paste("")
        tclvalue(level_guide2K_value) <- paste("")}
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
        if (hieroEnv$nLevs > 1) tclvalue(level_guide21_value) <- paste("Highest level unit for second group = m",hieroEnv$nLevs, sep="")
        else tclvalue(level_guide21_value) <- paste("Lowest level unit for second group = m1", sep="")
        if (hieroEnv$nLevs > 2) tclvalue(level_guide2mid_v) <- paste(". . .")
        else  tclvalue(level_guide2mid_v) <- paste("")
        if (hieroEnv$nLevs > 1) tclvalue(level_guide2K_value) <- paste("Lowest level unit for the second group = m1", sep="")
        else tclvalue(level_guide2K_value) <- paste("")} 
    tkconfigure(size,fg="black", bg="white")
    tkconfigure(sizeLabel,fg="black")
    SliderValue <- tclVar(hieroEnv$nLevs)
    tkconfigure(levsT, from=1, to=hieroEnv$nLevs, variable=SliderValue)
    changeGrids()
    grLevelConstr()
    }
oneSample <- function(){
    tkconfigure(levsGrLabel, fg="grey")
    tkconfigure(levsT, fg="grey")
    tkconfigure(varscostsHc2,fg="grey")
    tkconfigure(costs.12,fg="grey",bg="grey")
    tkconfigure(costs.22,fg="grey",bg="grey")
    tkconfigure(costs.32,fg="grey",bg="grey")
    tkconfigure(costs.42,fg="grey",bg="grey")
    tkconfigure(costs.52,fg="grey",bg="grey")
    tclvalue(level_guide21_value) <- paste("")
    tclvalue(level_guide2mid_v) <- paste("") 
    tclvalue(level_guide2K_value) <- paste("")
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){
        tkconfigure(Delta,fg="grey", bg="grey")
        tkconfigure(DeltaLabel,fg="grey")}
    changeLevs()
    }
twoSample <- function(){
    tkconfigure(levsGrLabel, fg="black")
    tkconfigure(levsT, fg="black")
    tkconfigure(varscostsHc2,fg="black")
    tclvalue(level_guide21_value) <- paste("Lowest level unit for second group = m1", sep="")
    if (hieroEnv$nLevs >= 3) tclvalue(level_guide2mid_v) <- paste(". . .")
    if (hieroEnv$nLevs >= 2) tclvalue(level_guide2K_value) <- paste("Highest level unit for second group = m",hieroEnv$nLevs, sep="")
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){
        tkconfigure(Delta,fg="grey", bg="grey")
        tkconfigure(DeltaLabel,fg="grey")}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt"){    
        tkconfigure(Delta,fg="black", bg="white")
        tkconfigure(DeltaLabel,fg="black")}
    tkconfigure(levsGrLabel, fg="black")
    tkconfigure(levsT, fg="black")
    changeLevs()
    }      
changeGrids <- function(){
    if (hieroEnv$nLevs < 6 ){
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
            tkgrid.forget(LL_6, oper1_6, funct_6, oper2_6, UL_6)
            tkconfigure(LL_6,textvariable=tclVar(""))
            tkconfigure(oper1_6,textvariable=tclVar(""))
            tkconfigure(funct_6,textvariable=tclVar(""))
            tkconfigure(oper2_6,textvariable=tclVar(""))
            tkconfigure(UL_6,textvariable=tclVar(""))
            tkgrid.forget(LL_7, oper1_7, funct_7, oper2_7, UL_7)
            tkconfigure(LL_7,textvariable=tclVar(""))
            tkconfigure(oper1_7,textvariable=tclVar(""))
            tkconfigure(funct_7,textvariable=tclVar(""))
            tkconfigure(oper2_7,textvariable=tclVar(""))
            tkconfigure(UL_7,textvariable=tclVar(""))
            tkgrid.forget(LL_8, oper1_8, funct_8, oper2_8, UL_8)
            tkconfigure(LL_8,textvariable=tclVar(""))
            tkconfigure(oper1_8,textvariable=tclVar(""))
            tkconfigure(funct_8,textvariable=tclVar(""))
            tkconfigure(oper2_8,textvariable=tclVar(""))
            tkconfigure(UL_8,textvariable=tclVar(""))
            tkgrid.forget(LL_9, oper1_9, funct_9, oper2_9, UL_9)
            tkconfigure(LL_9,textvariable=tclVar(""))
            tkconfigure(oper1_9,textvariable=tclVar(""))
            tkconfigure(funct_9,textvariable=tclVar(""))
            tkconfigure(oper2_9,textvariable=tclVar(""))
            tkconfigure(UL_9,textvariable=tclVar(""))
            tkgrid.forget(LL_10, oper1_10, funct_10, oper2_10, UL_10)
            tkconfigure(LL_10,textvariable=tclVar(""))
            tkconfigure(oper1_10,textvariable=tclVar(""))
            tkconfigure(funct_10,textvariable=tclVar(""))
            tkconfigure(oper2_10,textvariable=tclVar(""))
            tkconfigure(UL_10,textvariable=tclVar(""))}}
    if (hieroEnv$nLevs < 5 ){
        tkconfigure(varscosts5,fg="grey")
        tkconfigure(vars.5,fg="grey",bg = "grey")
        tkconfigure(costs.5,fg="grey",bg = "grey")
        tkconfigure(costs.52,fg="grey",bg = "grey")
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
            tkgrid.forget(LL_5, oper1_5, funct_5, oper2_5, UL_5)
            tkconfigure(LL_5,textvariable=tclVar(""))
            tkconfigure(oper1_5,textvariable=tclVar(""))
            tkconfigure(funct_5,textvariable=tclVar(""))
            tkconfigure(oper2_5,textvariable=tclVar(""))
            tkconfigure(UL_5,textvariable=tclVar(""))}
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkgrid.forget(LL_10, oper1_10, funct_10, oper2_10, UL_10)
            tkconfigure(LL_10,textvariable=tclVar(""))
            tkconfigure(oper1_10,textvariable=tclVar(""))
            tkconfigure(funct_10,textvariable=tclVar(""))
            tkconfigure(oper2_10,textvariable=tclVar(""))
            tkconfigure(UL_10,textvariable=tclVar(""))
            tkgrid.forget(LL_9, oper1_9, funct_9, oper2_9, UL_9)
            tkconfigure(LL_9,textvariable=tclVar(""))
            tkconfigure(oper1_9,textvariable=tclVar(""))
            tkconfigure(funct_9,textvariable=tclVar(""))
            tkconfigure(oper2_9,textvariable=tclVar(""))
            tkconfigure(UL_9,textvariable=tclVar(""))}}
    if (hieroEnv$nLevs < 4){
        tkconfigure(varscosts4,fg="grey")           
        tkconfigure(vars.4,fg="grey",bg = "grey")
        tkconfigure(costs.4,fg="grey",bg = "grey")
        tkconfigure(costs.42,fg="grey",bg = "grey")
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
            tkgrid.forget(LL_4, oper1_4, funct_4, oper2_4, UL_4)
            tkconfigure(LL_4,textvariable=tclVar(""))
            tkconfigure(oper1_4,textvariable=tclVar(""))
            tkconfigure(funct_4,textvariable=tclVar(""))
            tkconfigure(oper2_4,textvariable=tclVar(""))
            tkconfigure(UL_4,textvariable=tclVar(""))}
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkgrid.forget(LL_8, oper1_8, funct_8, oper2_8, UL_8)
            tkconfigure(LL_8,textvariable=tclVar(""))
            tkconfigure(oper1_8,textvariable=tclVar(""))
            tkconfigure(funct_8,textvariable=tclVar(""))
            tkconfigure(oper2_8,textvariable=tclVar(""))
            tkconfigure(UL_8,textvariable=tclVar(""))
            tkgrid.forget(LL_7, oper1_7, funct_7, oper2_7, UL_7)
            tkconfigure(LL_7,textvariable=tclVar(""))
            tkconfigure(oper1_7,textvariable=tclVar(""))
            tkconfigure(funct_7,textvariable=tclVar(""))
            tkconfigure(oper2_7,textvariable=tclVar(""))
            tkconfigure(UL_7,textvariable=tclVar(""))}}
    if (hieroEnv$nLevs < 3){
        tkconfigure(varscosts3,fg="grey")          
        tkconfigure(vars.3,fg="grey",bg = "grey")
        tkconfigure(costs.3,fg="grey",bg = "grey")
        tkconfigure(costs.32,fg="grey",bg = "grey")
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
            tkgrid.forget(LL_3, oper1_3, funct_3, oper2_3, UL_3)
            tkconfigure(LL_3,textvariable=tclVar(""))
            tkconfigure(oper1_3,textvariable=tclVar(""))
            tkconfigure(funct_3,textvariable=tclVar(""))
            tkconfigure(oper2_3,textvariable=tclVar(""))
            tkconfigure(UL_3,textvariable=tclVar(""))}
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkgrid.forget(LL_6, oper1_6, funct_6, oper2_6, UL_6)
            tkconfigure(LL_6,textvariable=tclVar(""))
            tkconfigure(oper1_6,textvariable=tclVar(""))
            tkconfigure(funct_6,textvariable=tclVar(""))
            tkconfigure(oper2_6,textvariable=tclVar(""))
            tkconfigure(UL_6,textvariable=tclVar(""))
            tkgrid.forget(LL_5, oper1_5, funct_5, oper2_5, UL_5)
            tkconfigure(LL_5,textvariable=tclVar(""))
            tkconfigure(oper1_5,textvariable=tclVar(""))
            tkconfigure(funct_5,textvariable=tclVar(""))
            tkconfigure(oper2_5,textvariable=tclVar(""))
            tkconfigure(UL_5,textvariable=tclVar(""))}}
    if (hieroEnv$nLevs < 2){
        tkconfigure(varscosts2,fg="grey")          
        tkconfigure(vars.2,fg="grey",bg = "grey")
        tkconfigure(costs.2,fg="grey",bg = "grey")
        tkconfigure(costs.22,fg="grey",bg = "grey")
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
            tkgrid.forget(LL_2, oper1_2, funct_2, oper2_2, UL_2)
            tkconfigure(LL_2,textvariable=tclVar(""))
            tkconfigure(oper1_2,textvariable=tclVar(""))
            tkconfigure(funct_2,textvariable=tclVar(""))
            tkconfigure(oper2_2,textvariable=tclVar(""))
            tkconfigure(UL_2,textvariable=tclVar(""))}
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkgrid.forget(LL_4, oper1_4, funct_4, oper2_4, UL_4)
            tkconfigure(LL_4,textvariable=tclVar(""))
            tkconfigure(oper1_4,textvariable=tclVar(""))
            tkconfigure(funct_4,textvariable=tclVar(""))
            tkconfigure(oper2_4,textvariable=tclVar(""))
            tkconfigure(UL_4,textvariable=tclVar(""))
            tkgrid.forget(LL_3, oper1_3, funct_3, oper2_3, UL_3)
            tkconfigure(LL_3,textvariable=tclVar(""))
            tkconfigure(oper1_3,textvariable=tclVar(""))
            tkconfigure(funct_3,textvariable=tclVar(""))
            tkconfigure(oper2_3,textvariable=tclVar(""))
            tkconfigure(UL_3,textvariable=tclVar(""))}
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
        tkconfigure(costs.12,fg="grey",bg = "grey")
        tkconfigure(costs.22,fg="grey",bg = "grey")
        tkconfigure(costs.32,fg="grey",bg = "grey")
        tkconfigure(costs.42,fg="grey",bg = "grey")
        tkconfigure(costs.52,fg="grey",bg = "grey")}         
    }
    if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
        tkconfigure(costs.12,fg="grey",bg = "grey")
        tkconfigure(costs.22,fg="grey",bg = "grey")
        tkconfigure(costs.32,fg="grey",bg = "grey")
        tkconfigure(costs.42,fg="grey",bg = "grey")
        tkconfigure(costs.52,fg="grey",bg = "grey")}
    if (hieroEnv$nLevs > 0){
        tkconfigure(varscosts1,fg="black") 
        tkconfigure(vars.1,fg="black",bg = "white")
        tkconfigure(costs.1,fg="black",bg = "white")
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
        tkgrid(LL_1, oper1_1, funct_1, oper2_1, UL_1)
        tkconfigure(LL_1,textvariable=tclVar("1"))
        tkconfigure(oper1_1,textvariable=tclVar("<="))
        tkconfigure(funct_1,textvariable=tclVar("n1"))}
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkconfigure(costs.12,fg="black",bg = "white")
            tkgrid(LL_2, oper1_2, funct_2, oper2_2, UL_2)
            tkconfigure(LL_2,textvariable=tclVar("1"))
            tkconfigure(oper1_2,textvariable=tclVar("<="))
            tkconfigure(funct_2,textvariable=tclVar("m1"))}}
    if (hieroEnv$nLevs > 1){
        tkconfigure(varscosts2,fg="black") 
        tkconfigure(vars.2,fg="black",bg = "white")
        tkconfigure(costs.2,fg="black",bg = "white")
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
            tkgrid(LL_2, oper1_2, funct_2, oper2_2, UL_2)
            tkconfigure(LL_2,textvariable=tclVar("1"))
            tkconfigure(oper1_2,textvariable=tclVar("<="))
            tkconfigure(funct_2,textvariable=tclVar("n2"))}
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkconfigure(costs.12,fg="black",bg = "white")
            tkconfigure(costs.22,fg="black",bg = "white")
            tkgrid(LL_3, oper1_3, funct_3, oper2_3, UL_3)
            tkconfigure(LL_3,textvariable=tclVar("1"))
            tkconfigure(oper1_3,textvariable=tclVar("<="))
            tkconfigure(funct_3,textvariable=tclVar("n2"))
            tkgrid(LL_4, oper1_4, funct_4, oper2_4, UL_4)
            tkconfigure(LL_4,textvariable=tclVar("1"))
            tkconfigure(oper1_4,textvariable=tclVar("<="))
            tkconfigure(funct_4,textvariable=tclVar("m2"))}}
    if (hieroEnv$nLevs > 2){
        tkconfigure(varscosts3,fg="black") 
        tkconfigure(vars.3,fg="black",bg = "white")
        tkconfigure(costs.3,fg="black",bg = "white")
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
            tkgrid(LL_3, oper1_3, funct_3, oper2_3, UL_3)
            tkconfigure(LL_3,textvariable=tclVar("1"))
            tkconfigure(oper1_3,textvariable=tclVar("<="))
            tkconfigure(funct_3,textvariable=tclVar("n3"))}
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkconfigure(costs.12,fg="black",bg = "white")
            tkconfigure(costs.22,fg="black",bg = "white")
            tkconfigure(costs.32,fg="black",bg = "white")
            tkgrid(LL_5, oper1_5, funct_5, oper2_5, UL_5)
            tkconfigure(LL_5,textvariable=tclVar("1"))
            tkconfigure(oper1_5,textvariable=tclVar("<="))
            tkconfigure(funct_5,textvariable=tclVar("n3"))
            tkgrid(LL_6, oper1_6, funct_6, oper2_6, UL_6)
            tkconfigure(LL_6,textvariable=tclVar("1"))
            tkconfigure(oper1_6,textvariable=tclVar("<="))
            tkconfigure(funct_6,textvariable=tclVar("m3"))}}  
    if ( hieroEnv$nLevs > 3){
        tkconfigure(varscosts4,fg="black") 
        tkconfigure(vars.4,fg="black",bg = "white")
        tkconfigure(costs.4,fg="black",bg = "white")
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
            tkgrid(LL_4, oper1_4, funct_4, oper2_4, UL_4)
            tkconfigure(LL_4,textvariable=tclVar("1"))
            tkconfigure(oper1_4,textvariable=tclVar("<="))
            tkconfigure(funct_4,textvariable=tclVar("n4"))}                                        
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkconfigure(costs.12,fg="black",bg = "white")
            tkconfigure(costs.22,fg="black",bg = "white")
            tkconfigure(costs.32,fg="black",bg = "white")        
            tkconfigure(costs.42,fg="black",bg = "white") 
            tkgrid(LL_7, oper1_7, funct_7, oper2_7, UL_7)
            tkconfigure(LL_7,textvariable=tclVar("1"))
            tkconfigure(oper1_7,textvariable=tclVar("<="))
            tkconfigure(funct_7,textvariable=tclVar("n4"))
            tkgrid(LL_8, oper1_8, funct_8, oper2_8, UL_8)
            tkconfigure(LL_8,textvariable=tclVar("1"))
            tkconfigure(oper1_8,textvariable=tclVar("<="))
            tkconfigure(funct_8,textvariable=tclVar("m4"))}}  
    if  (hieroEnv$nLevs > 4){
        tkconfigure(varscosts5,fg="black") 
        tkconfigure(vars.5,fg="black",bg = "white")
        tkconfigure(costs.5,fg="black",bg = "white")
        if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
            tkgrid(LL_5, oper1_5, funct_5, oper2_5, UL_5)
            tkconfigure(LL_5,textvariable=tclVar("1"))
            tkconfigure(oper1_5,textvariable=tclVar("<="))
            tkconfigure(funct_5,textvariable=tclVar("n5"))} 
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkconfigure(costs.12,fg="black",bg = "white")
            tkconfigure(costs.22,fg="black",bg = "white")
            tkconfigure(costs.32,fg="black",bg = "white")        
            tkconfigure(costs.42,fg="black",bg = "white")        
            tkconfigure(costs.52,fg="black",bg = "white")
            tkgrid(LL_9, oper1_9, funct_9, oper2_9, UL_9)
            tkconfigure(LL_9,textvariable=tclVar("1"))
            tkconfigure(oper1_9,textvariable=tclVar("<="))
            tkconfigure(funct_9,textvariable=tclVar("n5"))
            tkgrid(LL_10, oper1_10, funct_10, oper2_10, UL_10)
            tkconfigure(LL_10,textvariable=tclVar("1"))
            tkconfigure(oper1_10,textvariable=tclVar("<="))
            tkconfigure(funct_10,textvariable=tclVar("m5"))}}
    if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2"){
        hieroEnv$crows <- hieroEnv$nLevs}
    else if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
        hieroEnv$crows <- 2*hieroEnv$nLevs}            
}
addRowC <- function(){ 
    if (hieroEnv$crows == 0){tkgrid(LL_1, oper1_1, funct_1, oper2_1, UL_1)
                            tkconfigure(LL_1,textvariable=tclVar("1"))
                            tkconfigure(oper1_1,textvariable=tclVar("<="))}
    if (hieroEnv$crows == 1){tkgrid(LL_2, oper1_2, funct_2, oper2_2, UL_2)
                            tkconfigure(LL_2,textvariable=tclVar("1"))
                            tkconfigure(oper1_2,textvariable=tclVar("<="))}
    if (hieroEnv$crows == 2){tkgrid(LL_3, oper1_3, funct_3, oper2_3, UL_3)
                            tkconfigure(LL_3,textvariable=tclVar("1"))
                            tkconfigure(oper1_3,textvariable=tclVar("<="))}
    if (hieroEnv$crows == 3){tkgrid(LL_4, oper1_4, funct_4, oper2_4, UL_4)
                            tkconfigure(LL_4,textvariable=tclVar("1"))
                            tkconfigure(oper1_4,textvariable=tclVar("<="))}                                         
    if (hieroEnv$crows == 4){tkgrid(LL_5, oper1_5, funct_5, oper2_5, UL_5)
                            tkconfigure(LL_5,textvariable=tclVar("1"))
                            tkconfigure(oper1_5,textvariable=tclVar("<="))}  
    if (hieroEnv$crows == 5){tkgrid(LL_6, oper1_6, funct_6, oper2_6, UL_6)
                            tkconfigure(LL_6,textvariable=tclVar("1"))
                            tkconfigure(oper1_6,textvariable=tclVar("<="))}
    if (hieroEnv$crows == 6){tkgrid(LL_7, oper1_7, funct_7, oper2_7, UL_7)
                            tkconfigure(LL_7,textvariable=tclVar("1"))
                            tkconfigure(oper1_7,textvariable=tclVar("<="))}            
    if (hieroEnv$crows == 7){tkgrid(LL_8, oper1_8, funct_8, oper2_8, UL_8)
                            tkconfigure(LL_8,textvariable=tclVar("1"))
                            tkconfigure(oper1_8,textvariable=tclVar("<="))}              
    if (hieroEnv$crows == 8){tkgrid(LL_9, oper1_9, funct_9, oper2_9, UL_9)
                            tkconfigure(LL_9,textvariable=tclVar("1"))
                            tkconfigure(oper1_9,textvariable=tclVar("<="))}
    if (hieroEnv$crows == 9){tkgrid(LL_10, oper1_10, funct_10, oper2_10, UL_10)
                            tkconfigure(LL_10,textvariable=tclVar("1"))
                            tkconfigure(oper1_10,textvariable=tclVar("<="))}
    hieroEnv$crows <- hieroEnv$crows+1
    }
deleteRowC <- function(){
    if (hieroEnv$crows == 1){tkgrid.forget(LL_1, oper1_1, funct_1, oper2_1, UL_1)
                    tkconfigure(LL_1,textvariable=tclVar(""))
                    tkconfigure(oper1_1,textvariable=tclVar(""))
                    tkconfigure(funct_1,textvariable=tclVar(""))
                    tkconfigure(oper2_1,textvariable=tclVar(""))
                    tkconfigure(UL_1,textvariable=tclVar(""))}
    if (hieroEnv$crows== 2){tkgrid.forget(LL_2, oper1_2, funct_2, oper2_2, UL_2)
                    tkconfigure(LL_2,textvariable=tclVar(""))
                    tkconfigure(oper1_2,textvariable=tclVar(""))
                    tkconfigure(funct_2,textvariable=tclVar(""))
                    tkconfigure(oper2_2,textvariable=tclVar(""))
                    tkconfigure(UL_2,textvariable=tclVar(""))}
    if (hieroEnv$crows == 3){tkgrid.forget(LL_3, oper1_3, funct_3, oper2_3, UL_3)
                    tkconfigure(LL_3,textvariable=tclVar(""))
                    tkconfigure(oper1_3,textvariable=tclVar(""))
                    tkconfigure(funct_3,textvariable=tclVar(""))
                    tkconfigure(oper2_3,textvariable=tclVar(""))
                    tkconfigure(UL_3,textvariable=tclVar(""))}
    if (hieroEnv$crows == 4){tkgrid.forget(LL_4, oper1_4, funct_4, oper2_4, UL_4)
                    tkconfigure(LL_4,textvariable=tclVar(""))
                    tkconfigure(oper1_4,textvariable=tclVar(""))
                    tkconfigure(funct_4,textvariable=tclVar(""))
                    tkconfigure(oper2_4,textvariable=tclVar(""))
                    tkconfigure(UL_4,textvariable=tclVar(""))}                                         
    if (hieroEnv$crows == 5){tkgrid.forget(LL_5, oper1_5, funct_5, oper2_5, UL_5)
                    tkconfigure(LL_5,textvariable=tclVar(""))
                    tkconfigure(oper1_5,textvariable=tclVar(""))
                    tkconfigure(funct_5,textvariable=tclVar(""))
                    tkconfigure(oper2_5,textvariable=tclVar(""))
                    tkconfigure(UL_5,textvariable=tclVar(""))}  
    if (hieroEnv$crows == 6){tkgrid.forget(LL_6, oper1_6, funct_6, oper2_6, UL_6)
                    tkconfigure(LL_6,textvariable=tclVar(""))
                    tkconfigure(oper1_6,textvariable=tclVar(""))
                    tkconfigure(funct_6,textvariable=tclVar(""))
                    tkconfigure(oper2_6,textvariable=tclVar(""))
                    tkconfigure(UL_6,textvariable=tclVar(""))}
    if (hieroEnv$crows == 7){tkgrid.forget(LL_7, oper1_7, funct_7, oper2_7, UL_7)
                    tkconfigure(LL_7,textvariable=tclVar(""))
                    tkconfigure(oper1_7,textvariable=tclVar(""))
                    tkconfigure(funct_7,textvariable=tclVar(""))
                    tkconfigure(oper2_7,textvariable=tclVar(""))
                    tkconfigure(UL_7,textvariable=tclVar(""))}            
    if (hieroEnv$crows == 8){tkgrid.forget(LL_8, oper1_8, funct_8, oper2_8, UL_8)
                    tkconfigure(LL_8,textvariable=tclVar(""))
                    tkconfigure(oper1_8,textvariable=tclVar(""))
                    tkconfigure(funct_8,textvariable=tclVar(""))
                    tkconfigure(oper2_8,textvariable=tclVar(""))
                    tkconfigure(UL_8,textvariable=tclVar(""))}              
    if (hieroEnv$crows == 9){tkgrid.forget(LL_9, oper1_9, funct_9, oper2_9, UL_9)
                    tkconfigure(LL_9,textvariable=tclVar(""))
                    tkconfigure(oper1_9,textvariable=tclVar(""))
                    tkconfigure(funct_9,textvariable=tclVar(""))
                    tkconfigure(oper2_9,textvariable=tclVar(""))
                    tkconfigure(UL_9,textvariable=tclVar(""))}
    if (hieroEnv$crows == 10){tkgrid.forget(LL_10, oper1_10, funct_10, oper2_10, UL_10)
                    tkconfigure(LL_10,textvariable=tclVar(""))
                    tkconfigure(oper1_10,textvariable=tclVar(""))
                    tkconfigure(funct_10,textvariable=tclVar(""))
                    tkconfigure(oper2_10,textvariable=tclVar(""))
                    tkconfigure(UL_10,textvariable=tclVar(""))}
    hieroEnv$crows <- hieroEnv$crows-1
    }
minCosts <- function(){
    tkconfigure(LimVCLabel,fg="grey")
    tkconfigure(LimVC,fg="grey",bg = "grey")
    tkconfigure(powLabel,fg="grey")
    tkconfigure(Limpow,fg="grey",bg = "grey")
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct"){
        tkconfigure(powLabel,fg="black")
        tkconfigure(Limpow,fg="black",bg = "white")
        tkconfigure(Delta,fg="black", bg="white")
        tkconfigure(DeltaLabel,fg="black")}
    else if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce"){
        tkconfigure(LimVCLabel,fg="black",text="Limit for variance \n of the estimate")
        tkconfigure(LimVC,fg="black",bg = "white",text=tclVar(""))
        tkconfigure(Delta,fg="grey", bg="grey")
        tkconfigure(DeltaLabel,fg="grey")}        
}
minVar <- function(){
    tkconfigure(powLabel,fg="grey")
    tkconfigure(Limpow,fg="grey",bg = "grey")
    tkconfigure(LimVCLabel,fg="black",text="Limit for total costs")
    tkconfigure(LimVC,fg="black",bg = "white",text=tclVar(""))
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt"){
        tkconfigure(Delta,fg="black", bg="white")
        tkconfigure(DeltaLabel,fg="black")}
    else if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){
        tkconfigure(Delta,fg="grey", bg="grey")
        tkconfigure(DeltaLabel,fg="grey")}
}
grLevelConstr <- function(){
    tkconfigure(e.constraints,fg="grey") 
    tclvalue(c.extra.l1t) <- paste("")
    tclvalue(c.extra.l2t) <- paste("")
    tclvalue(c.extra.l3t) <- paste("")
    tclvalue(c.extra.l4t) <- paste("")            
    if (hieroEnv$nLevs > as.numeric(tkget(levsT)) && as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
        tkconfigure(e.constraints,fg="black") 
        if (hieroEnv$nLevs == 5){
            tclvalue(c.extra.l1t) <- paste("n5 = m5")
            tclvalue(c.extra.l2t) <- paste("")
            tclvalue(c.extra.l3t) <- paste("")
            tclvalue(c.extra.l4t) <- paste("")
            if (as.numeric(tkget(levsT)) < 4){
                tclvalue(c.extra.l2t) <- paste("n4 = m4")}
            if (as.numeric(tkget(levsT)) < 3){
                tclvalue(c.extra.l3t) <- paste("n3 = m3")}
            if (as.numeric(tkget(levsT)) < 2){
                tclvalue(c.extra.l4t) <- paste("n2 = m2")}}
        if (hieroEnv$nLevs == 4){
            tclvalue(c.extra.l1t) <- paste("n4 = m4")
            tclvalue(c.extra.l2t) <- paste("")
            tclvalue(c.extra.l3t) <- paste("")
            if (as.numeric(tkget(levsT)) < 3){
                tclvalue(c.extra.l2t) <- paste("n3 = m3")}
            if (as.numeric(tkget(levsT)) < 2){
                tclvalue(c.extra.l3t) <- paste("n2 = m2")}
                tclvalue(c.extra.l4t) <- paste("")}
        if (hieroEnv$nLevs == 3){
                tclvalue(c.extra.l1t) <- paste("n3 = m3")
                tclvalue(c.extra.l2t) <- paste("")
                if (as.numeric(tkget(levsT)) < 2){
                    tclvalue(c.extra.l2t) <- paste("n2 = m2")}  
                    tclvalue(c.extra.l4t) <- paste("")
                    tclvalue(c.extra.l3t) <- paste("")}
        if (hieroEnv$nLevs == 2){
            tclvalue(c.extra.l1t) <- paste("n2 = m2")
            tclvalue(c.extra.l4t) <- paste("")
            tclvalue(c.extra.l3t) <- paste("")
            tclvalue(c.extra.l2t) <- paste("")}}           
    tkconfigure(c.extra.l1,text=paste(tclvalue(c.extra.l1t)))
    tkconfigure(c.extra.l2,text=paste(tclvalue(c.extra.l2t)))
    tkconfigure(c.extra.l3,text=paste(tclvalue(c.extra.l3t)))
    tkconfigure(c.extra.l4,text=paste(tclvalue(c.extra.l4t)))
    }
tkbind(levsT, "<ButtonRelease-1>", grLevelConstr)
DefGams <- function(){
    hieroEnv$err <- as.numeric(0)
    tkdelete(formulas_txt,"0.0","end")
    errorSamples()
    errorTrLevel()
    errorObj()
    errorVariance()
    errorCosts()
    errorConstr()
    errorNumeric()
    hieroEnv$newProbName <- tkget(problemName)
    if (hieroEnv$err == 0){
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
        tclvalue(f.funct.lt) <- paste("You have constructed a ", hieroEnv$nLevs,"-level, \n one group optimization problem.", sep=" ")}
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){                           
        tclvalue(f.funct.lt) <- paste("You have constructed a ", hieroEnv$nLevs,"-level, \n two groups optimization problem.", sep="")}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt"){
        tclvalue(f.funct.ot) <- paste("Objective is to optimize sample size so that \n the statistical power is maximized when \n testing difference of means and \n total costs are restricted to ", tkget(LimVC), ".", sep="")}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){
        tclvalue(f.funct.ot) <- paste("Objective is to optimize sample size so that \n the length of confidence interval is minimized \n when estimating the mean and \n total costs are restricted to ", tkget(LimVC), ".", sep="")}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct"){
        tclvalue(f.funct.ot) <- paste("Objective is to optimize sample size so that \n the total costs are minimized when \n testing  difference of means and \n statistical power is limited to (minumum) ", tkget(Limpow), ".", sep="")}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce"){
        tclvalue(f.funct.ot) <- paste("Objective is to optimize sample size so that \n the total costs are minimized when \n estimating the mean and variance of the \n estimate is limited to (maximum) ", tkget(LimVC), ".", sep="")}
    tclvalue(f.funct.ct) <- paste("Constraints are set as \n",
                                    tkget(LL_1), tkget(oper1_1), tkget(funct_1), tkget(oper2_1), tkget(UL_1), "\n",
                                    tkget(LL_2), tkget(oper1_2), tkget(funct_2), tkget(oper2_2), tkget(UL_2), "\n", 
                                    tkget(LL_3), tkget(oper1_3), tkget(funct_3), tkget(oper2_3), tkget(UL_3), "\n",                                                                        
                                    tkget(LL_4), tkget(oper1_4), tkget(funct_4), tkget(oper2_4), tkget(UL_4), "\n", 
                                    tkget(LL_5), tkget(oper1_5), tkget(funct_5), tkget(oper2_5), tkget(UL_5), "\n", 
                                    tkget(LL_6), tkget(oper1_6), tkget(funct_6), tkget(oper2_6), tkget(UL_6), "\n", 
                                    tkget(LL_7), tkget(oper1_7), tkget(funct_7), tkget(oper2_7), tkget(UL_7), "\n", 
                                    tkget(LL_8), tkget(oper1_8), tkget(funct_8), tkget(oper2_8), tkget(UL_8), "\n", 
                                    tkget(LL_9), tkget(oper1_9), tkget(funct_9), tkget(oper2_9), tkget(UL_9), "\n", 
                                    tkget(LL_10), tkget(oper1_10), tkget(funct_10), tkget(oper2_10), tkget(UL_10), sep=" ")   
    tkconfigure(f.funct.l,text=paste(tclvalue(f.funct.lt)))    
    tkconfigure(f.funct.c,text=paste(tclvalue(f.funct.ct)))
    tkconfigure(f.funct.o,text=paste(tclvalue(f.funct.ot)))                                
    tclvalue(f.1.lt) <- paste("Size (alpha) =", tkget(size), sep=" ")
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct"){
        tclvalue(f.2.lt) <- paste("Power (1-beta) =", tkget(Limpow), sep=" ")}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce"){
        tclvalue(f.2.lt) <- paste("Limit for total variance V<=", tkget(LimVC), sep=" ")}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){
        tclvalue(f.2.lt) <- paste("Limit for total costs C<=", tkget(LimVC), sep=" ")}    
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt"){
        tclvalue(f.3.lt) <- paste("Effect size =", tkget(Delta), sep=" ")}
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
        if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){ 
         tclvalue(f.3.lt) <- paste("Level of grouping =", tkget(levsT), sep=" ")
         tclvalue(f.4.lt) <- paste("")}
        if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt"){ 
         tclvalue(f.4.lt) <- paste("Level of grouping =", tkget(levsT), sep=" ")}}
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
         if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){ 
            tclvalue(f.3.lt) <- paste("")
            tclvalue(f.4.lt) <- paste("")}
        if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt"){    
         tclvalue(f.4.lt) <- paste("")}}        
    tkconfigure(f.1.l,text=paste(tclvalue(f.1.lt)))
    tkconfigure(f.2.l,text=paste(tclvalue(f.2.lt)))        
    tkconfigure(f.3.l,text=paste(tclvalue(f.3.lt))) 
    tkconfigure(f.4.l,text=paste(tclvalue(f.4.lt)))    
    tkconfigure(f.lev.l,textvariable=tclVar("  Level  "))
    tkconfigure(f.lev1.l,textvariable=tclVar("Level 1"))
    if (hieroEnv$nLevs < 2) tkconfigure(f.lev2.l,textvariable=tclVar(""))
    else tkconfigure(f.lev2.l,textvariable=tclVar("Level 2"))
    if (hieroEnv$nLevs < 3) tkconfigure(f.lev3.l,textvariable=tclVar(""))     
    else tkconfigure(f.lev3.l,textvariable=tclVar("Level 3")) 
    if (hieroEnv$nLevs < 4) tkconfigure(f.lev4.l,textvariable=tclVar(""))     
    else tkconfigure(f.lev4.l,textvariable=tclVar("Level 4")) 
    if (hieroEnv$nLevs < 5) tkconfigure(f.lev5.l,textvariable=tclVar(""))             
    else tkconfigure(f.lev5.l,textvariable=tclVar("Level 5")) 
    tclvalue(f.var1.lt) <- paste("v1 =", tkget(vars.1), sep=" ")
    tclvalue(f.var2.lt) <- paste("v2 =", tkget(vars.2), sep=" ")
    tclvalue(f.var3.lt) <- paste("v3 =", tkget(vars.3), sep=" ")
    tclvalue(f.var4.lt) <- paste("v4 =", tkget(vars.4), sep=" ")
    tclvalue(f.var5.lt) <- paste("v5 =", tkget(vars.5), sep=" ")
    if (hieroEnv$nLevs < 2) tclvalue(f.var2.lt) <- paste("")
    if (hieroEnv$nLevs < 3) tclvalue(f.var3.lt) <- paste("")
    if (hieroEnv$nLevs < 4) tclvalue(f.var4.lt) <- paste("") 
    if (hieroEnv$nLevs < 5) tclvalue(f.var5.lt) <- paste("")  
    tkconfigure(f.var.l,textvariable=tclVar("  Variance  "))
    tkconfigure(f.var1.l,text=paste(tclvalue(f.var1.lt)))    
    tkconfigure(f.var2.l,text=paste(tclvalue(f.var2.lt)))      
    tkconfigure(f.var3.l,text=paste(tclvalue(f.var3.lt))) 
    tkconfigure(f.var4.l,text=paste(tclvalue(f.var4.lt)))              
    tkconfigure(f.var5.l,text=paste(tclvalue(f.var5.lt)))         
    tclvalue(f.cost1.lt) <- paste("c11 =", tkget(costs.1), sep=" ")
    tclvalue(f.cost2.lt) <- paste("c21 =", tkget(costs.2), sep=" ")
    tclvalue(f.cost3.lt) <- paste("c31 =", tkget(costs.3), sep=" ")
    tclvalue(f.cost4.lt) <- paste("c41 =", tkget(costs.4), sep=" ")
    tclvalue(f.cost5.lt) <- paste("c51 =", tkget(costs.5), sep=" ")
    if (hieroEnv$nLevs < 2) tclvalue(f.cost2.lt) <- paste("")
    if (hieroEnv$nLevs < 3) tclvalue(f.cost3.lt) <- paste("")
    if (hieroEnv$nLevs < 4) tclvalue(f.cost4.lt) <- paste("") 
    if (hieroEnv$nLevs < 5) tclvalue(f.cost5.lt) <- paste("") 
    tkconfigure(f.cost1.l,text=paste(tclvalue(f.cost1.lt)))
    tkconfigure(f.cost2.l,text=paste(tclvalue(f.cost2.lt)))    
    tkconfigure(f.cost3.l,text=paste(tclvalue(f.cost3.lt)))    
    tkconfigure(f.cost4.l,text=paste(tclvalue(f.cost4.lt)))    
    tkconfigure(f.cost5.l,text=paste(tclvalue(f.cost5.lt)))
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
        tclvalue(f.c1.lt) <- paste("  Costs  ")
        tclvalue(f.c2.lt) <- paste("")      
        tclvalue(f.cost12.lt) <- paste("")
        tclvalue(f.cost22.lt) <- paste("")
        tclvalue(f.cost32.lt) <- paste("")
        tclvalue(f.cost42.lt) <- paste("")
        tclvalue(f.cost52.lt) <- paste("")                        
        } 
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
        tclvalue(f.c1.lt) <- paste("  Costs for \n group 1  ")  
        tclvalue(f.c2.lt) <- paste("  Costs for \n group 2  ")
        tclvalue(f.cost12.lt) <- paste("c12 =", tkget(costs.12), sep=" ")        
        tclvalue(f.cost22.lt) <- paste("c22 =", tkget(costs.22), sep=" ")
        tclvalue(f.cost32.lt) <- paste("c32 =", tkget(costs.32), sep=" ") 
        tclvalue(f.cost42.lt) <- paste("c42 =", tkget(costs.42), sep=" ") 
        tclvalue(f.cost52.lt) <- paste("c52 =", tkget(costs.52), sep=" ")
        }     

    if (hieroEnv$nLevs < 2) tclvalue(f.cost22.lt) <- paste("")
    if (hieroEnv$nLevs < 3) tclvalue(f.cost32.lt) <- paste("")
    if (hieroEnv$nLevs < 4) tclvalue(f.cost42.lt) <- paste("") 
    if (hieroEnv$nLevs < 5) tclvalue(f.cost52.lt) <- paste("")
    
    tkconfigure(f.c1.l,text=paste(tclvalue(f.c1.lt)))
    tkconfigure(f.c2.l,text=paste(tclvalue(f.c2.lt))) 
    tkconfigure(f.cost12.l,text=paste(tclvalue(f.cost12.lt)))
    tkconfigure(f.cost22.l,text=paste(tclvalue(f.cost22.lt)))
    tkconfigure(f.cost32.l,text=paste(tclvalue(f.cost32.lt)))
    tkconfigure(f.cost42.l,text=paste(tclvalue(f.cost42.lt)))
    tkconfigure(f.cost52.l,text=paste(tclvalue(f.cost52.lt)))
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce"){
        tkinsert(formulas_txt, "end", paste("VARIABLE COSTS;"))}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){
        tkinsert(formulas_txt, "end", paste("VARIABLE VARIANCE;"))}    
    tkinsert(formulas_txt, "end", "\n \n")
    tkinsert(formulas_txt, "end", paste("INTEGER VARIABLES n1"))
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
        tkinsert(formulas_txt, "end", ", m1")}
    if (hieroEnv$nLevs >= 2){tkinsert(formulas_txt, "end",  ", n2")
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", ", m2")}}
    if (hieroEnv$nLevs >= 3){tkinsert(formulas_txt, "end", ", n3")
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", ", m3")}}
    if (hieroEnv$nLevs >= 4){tkinsert(formulas_txt, "end", ", n4")
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", ", m4")}}
    if (hieroEnv$nLevs >= 5){tkinsert(formulas_txt, "end", ", n5")
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", ", m5")}}
    tkinsert(formulas_txt, "end", ";\n \n")
    n.cons <- 1
    for (i in 1:hieroEnv$crows){
        if(nchar(paste(tkget(get(paste("funct_", i, sep = ""))))) > 3 ){
            if(as.numeric(tkget(get(paste("LL_", i, sep = ""))))>=0){
                n.cons <- n.cons+1}
                if(length(as.numeric(tkget(get(paste("UL_", i, sep = "")))))>=1){
                    n.cons <- n.cons+1}}}
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            if (hieroEnv$nLevs > as.numeric(tkget(levsT))){
            n.cons <- n.cons+(hieroEnv$nLevs-as.numeric(tkget(levsT)))}}        
    tkinsert(formulas_txt, "end", "EQUATIONS CON")
    for(i in 1:n.cons){
        tkinsert(formulas_txt, "end", ", CON")
        tkinsert(formulas_txt, "end", i)}
    tkinsert(formulas_txt, "end", ";\n \n")
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce"){
        tkinsert(formulas_txt, "end", paste("CON..     COSTS =E= "))
        transCosts()
        tkinsert(formulas_txt, "end", paste("; \n \n","CON1..     ", sep=""))
        transVariation()
        if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce"){
            tkinsert(formulas_txt, "end", paste(" =L= ", tkget(LimVC), ";\n", sep=""))}
        if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct"){
            size.f <- as.numeric(tkget(size))
            power.f <- as.numeric(tkget(Limpow))
            Delta.f <- as.numeric(tkget(Delta))
            sigma2 = Delta.f^2/delta2(size.f,power.f)
            tkinsert(formulas_txt, "end", paste(" =L= ", sigma2, ";\n", sep=""))}}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){
        tkinsert(formulas_txt, "end", paste("CON..     VARIANCE =E= "))
        transVariation()
        tkinsert(formulas_txt, "end", paste("; \n \n","CON1..     ", sep=""))
        transCosts()
        tkinsert(formulas_txt, "end", paste(" =L= ", tkget(LimVC), "; \n", sep=""))}     
    cons_tempo=2
    if (paste("funct_", hieroEnv$crows, sep = "") == ""){
        errorConstr()}
    else     
        for (i in 1:hieroEnv$crows){
            if(nchar(paste(tkget(get(paste("funct_", i, sep = ""))))) > 3 ){
                if(as.numeric(tkget(get(paste("LL_", i, sep=""))))>=0){
                    tkinsert(formulas_txt, "end", paste("CON", cons_tempo, "..      ", tkget(get(paste("LL_", i, sep=""))), sep=""))
                    if (Operators[as.numeric(tclvalue(tcl(get(paste("oper1_",i, sep="")),"getvalue")))+1] == "<="){tkinsert(formulas_txt, "end", " =L= ")}
                    if (Operators[as.numeric(tclvalue(tcl(get(paste("oper1_",i, sep="")),"getvalue")))+1] == "="){tkinsert(formulas_txt, "end", " =E= ")}
                        tkinsert(formulas_txt, "end", paste(tkget(get(paste("funct_",i, sep=""))), ";\n", sep=""))
                        cons_tempo=cons_tempo+1}
            if(length(as.numeric(tkget(get(paste("UL_", i, sep="")))))>=1){
                tkinsert(formulas_txt, "end", paste("CON", cons_tempo, "..      ", tkget(get(paste("funct_", i, sep=""))),sep=""))
                    if (Operators[as.numeric(tclvalue(tcl(get(paste("oper2_",i, sep="")),"getvalue")))+1] == "<="){tkinsert(formulas_txt, "end", " =L= ")}
                    if (Operators[as.numeric(tclvalue(tcl(get(paste("oper2_",i, sep="")),"getvalue")))+1] == "="){tkinsert(formulas_txt, "end", " =E= ")}     
                        tkinsert(formulas_txt, "end", paste(tkget(get(paste("UL_",i, sep=""))),";\n", sep=""))
                        cons_tempo=cons_tempo+1}}}
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
        if (hieroEnv$nLevs > as.numeric(tkget(levsT))){
            if (hieroEnv$nLevs >= 5){
                tkinsert(formulas_txt, "end", paste("CON", cons_tempo, "..      n5/m5 =E= 1;\n", sep=""))
                cons_tempo=cons_tempo+1}
            if (hieroEnv$nLevs >= 4 && as.numeric(tkget(levsT)) < 4){
                tkinsert(formulas_txt, "end", paste("CON", cons_tempo, "..      n4/m4 =E= 1;\n", sep=""))
                cons_tempo=cons_tempo+1}        
            if (hieroEnv$nLevs >= 3 && as.numeric(tkget(levsT)) < 3){
                tkinsert(formulas_txt, "end", paste("CON", cons_tempo, "..      n3/m3 =E= 1;\n", sep=""))
                cons_tempo=cons_tempo+1} 
            if (hieroEnv$nLevs >= 2 && as.numeric(tkget(levsT)) < 2){
                tkinsert(formulas_txt, "end", paste("CON", cons_tempo, "..      n2/m2 =E= 1;\n", sep=""))
                cons_tempo=cons_tempo+1}}}          
    if(hieroEnv$crows >= 1 && length(paste(tkget(funct_1))) !=0 ){
        if(nchar(paste(tkget(funct_1))) < 4 ){
            if(as.numeric(tkget(LL_1))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_1), ".lo = ",tkget(LL_1), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_1)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_1), ".up = ",tkget(UL_1), "; ", sep=""))}}}
    if(hieroEnv$crows >= 2){
        if(nchar(paste(tkget(funct_2))) < 4 ){
            if(as.numeric(tkget(LL_2))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_2), ".lo = ",tkget(LL_2), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_2)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_2), ".up = ",tkget(UL_2), "; ", sep=""))}}}
    if(hieroEnv$crows >= 3){  
        if(nchar(paste(tkget(funct_3))) < 4 ){
            if(as.numeric(tkget(LL_3))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_3), ".lo = ",tkget(LL_3), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_3)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_3), ".up = ",tkget(UL_3), "; ", sep=""))}}}
    if(hieroEnv$crows >= 4){  
        if(nchar(paste(tkget(funct_4))) < 4 ){
            if(as.numeric(tkget(LL_4))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_4), ".lo = ",tkget(LL_4), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_4)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_4), ".up = ",tkget(UL_4), "; ", sep=""))}}}
    if(hieroEnv$crows >= 5){      
        if(nchar(paste(tkget(funct_5))) < 4 ){
            if(as.numeric(tkget(LL_5))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_5), ".lo = ",tkget(LL_5), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_5)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_5), ".up = ",tkget(UL_5), "; ", sep=""))}}}
    if(hieroEnv$crows >= 6){                 
        if(nchar(paste(tkget(funct_6))) < 4 ){
            if(as.numeric(tkget(LL_6))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_6), ".lo = ",tkget(LL_6), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_6)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_6), ".up = ",tkget(UL_6), "; ", sep=""))}}}
    if(hieroEnv$crows >= 7){
        if(nchar(paste(tkget(funct_7))) < 4 ){
            if(as.numeric(tkget(LL_7))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_7), ".lo = ",tkget(LL_7), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_7)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_7), ".up = ",tkget(UL_7), "; ", sep=""))}}}
    if(hieroEnv$crows >= 8){
        if(nchar(paste(tkget(funct_8))) < 4 ){
            if(as.numeric(tkget(LL_8))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_8), ".lo = ",tkget(LL_8), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_8)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_8), ".up = ",tkget(UL_8), "; ", sep=""))}}}
    if(hieroEnv$crows >= 9){
        if(nchar(paste(tkget(funct_9))) < 4 ){
            if(as.numeric(tkget(LL_9))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_9), ".lo = ",tkget(LL_9), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_9)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_9), ".up = ",tkget(UL_9), "; ", sep=""))}}}
    if(hieroEnv$crows >= 10){
        if(nchar(paste(tkget(funct_10))) < 4 ){
            if(as.numeric(tkget(LL_10))>=0){ 
                tkinsert(formulas_txt, "end", paste("\n", tkget(funct_10), ".lo = ",tkget(LL_10), "; ", sep=""))} 
            if(length(as.numeric(tkget(UL_10)))>=1){
                tkinsert(formulas_txt, "end", paste(tkget(funct_10), ".up = ",tkget(UL_10), "; ", sep=""))}}}                       
    tkinsert(formulas_txt, "end", "\n \n")
    tkinsert(formulas_txt, "end",paste("option optcr =", tkget(tolerance),";\n\n", sep=" "))
    tkinsert(formulas_txt, "end", paste("MODEL", tkget(problemName), "/ALL/ ;", sep=" "))
    tkinsert(formulas_txt, "end", "\n \n")
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce"){
        tkinsert(formulas_txt, "end", paste("SOLVE", tkget(problemName), "USING MINLP MINIMIZING COSTS;", sep=" "))}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve"){    
        tkinsert(formulas_txt, "end", paste("SOLVE", tkget(problemName), "USING MINLP MINIMIZING VARIANCE;", sep=" "))}
    tkinsert(formulas_txt, "end", "\n \n")
    tkinsert(formulas_txt, "end", paste("display n1.l"))
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
        tkinsert(formulas_txt, "end", ", m1.l")}
    if (hieroEnv$nLevs >= 2){tkinsert(formulas_txt, "end",  ", n2.l")
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", ", m2.l")}}
    if (hieroEnv$nLevs >= 3){tkinsert(formulas_txt, "end", ", n3.l")
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", ", m3.l")}}
    if (hieroEnv$nLevs >= 4){tkinsert(formulas_txt, "end", ", n4.l")
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", ", m4.l")}}
    if (hieroEnv$nLevs >= 5){tkinsert(formulas_txt, "end", ", n5.l")
        if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", ", m5.l")}}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce"){
        tkinsert(formulas_txt, "end", ", COSTS.l")}
    else tkinsert(formulas_txt, "end", ", VARIANCE.l")            
        tkinsert(formulas_txt, "end", paste(", CON1.l, \n \t '....END...., ", as.character(tclvalue(hieroEnv$rbValue2)), ", size.l=", as.numeric(tkget(size)),"_1_, ",sep=""))
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct"){
        tkinsert(formulas_txt, "end", paste("power.l=", as.numeric(tkget(Limpow)),"_3_, ",sep=""))}
    if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt"){
        tkinsert(formulas_txt, "end", paste("Delta.l=", as.numeric(tkget(Delta)),"_2_, ",sep=""))}
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
        tkinsert(formulas_txt, "end", paste("groups.l=1"))} 
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
        tkinsert(formulas_txt, "end", paste("groups.l=2"))}        
    tkinsert(formulas_txt, "end", paste("';\n \n"))
    if (hieroEnv$err == 0) tk2notetab.select(nb, "Final problem")
    if (hieroEnv$err >= 1) tk2notetab.select(nb, "Define")      
    hieroEnv$wfile <- tkget(problemName)}
}       
transCosts <- function(){
if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
    if (hieroEnv$nLevs == 1){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1",sep=""))}
    if (hieroEnv$nLevs == 2){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1*n2+", tkget(costs.2), "*n2", sep=""))}
    if (hieroEnv$nLevs == 3){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1*n2*n3+", tkget(costs.2), "*n2*n3+", tkget(costs.3), "*n3", sep=""))}
    if (hieroEnv$nLevs == 4){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1*n2*n3*n4+", tkget(costs.2), "*n2*n3*n4 \n                  +", tkget(costs.3), "*n3*n4+", tkget(costs.4), "*n4", sep=""))}
    if (hieroEnv$nLevs == 5){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1*n2*n3*n4*n5+", tkget(costs.2), "*n2*n3*n4*n5 \n                  +", tkget(costs.3), "*n3*n4*n5+", tkget(costs.4), "*n4*n5+", tkget(costs.5), "*n5",sep=""))}}
if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
    if (hieroEnv$nLevs == 1){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1", sep=""))
        tkinsert(formulas_txt, "end", paste("+",tkget(costs.12), "*m1",sep=""))}
    if (hieroEnv$nLevs == 2){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1*n2+", tkget(costs.2), "*n2", sep=""))
        tkinsert(formulas_txt, "end", paste("\n                   +", tkget(costs.12), "*m1*m2+", tkget(costs.22), "*m2", sep=""))}
    if (hieroEnv$nLevs == 3){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1*n2*n3+", tkget(costs.2), "*n2*n3+", tkget(costs.3), "*n3", sep=""))
        tkinsert(formulas_txt, "end", paste("\n                   +",tkget(costs.12), "*m1*m2*m3+", tkget(costs.22), "*m2*m3+", tkget(costs.32), "*m3", sep=""))}
    if (hieroEnv$nLevs == 4){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1*n2*n3*n4+", tkget(costs.2), 
                    "*n2*n3*n4 \n                   +", tkget(costs.3), "*n3*n4+", tkget(costs.4), "*n4", sep=""))
        tkinsert(formulas_txt, "end", paste("\n                   +",tkget(costs.12), "*m1*m2*m3*m4+", tkget(costs.22), 
                    "*m2*m3*m4 \n                   +", tkget(costs.32), "*m3*m4+", tkget(costs.42), "*m4", sep=""))}
    if (hieroEnv$nLevs == 5){
        tkinsert(formulas_txt, "end", paste(tkget(costs.1), "*n1*n2*n3*n4*n5+", tkget(costs.2), "*n2*n3*n4*n5 \n                  +",
            tkget(costs.3), "*n3*n4*n5+", tkget(costs.4), "*n4*n5+", tkget(costs.5), "*n5",sep=""))
        tkinsert(formulas_txt, "end", paste("\n                   +",tkget(costs.12), "*m1*m2*m3*m4*m5+", tkget(costs.22), "*m2*m3*m4*m5 \n                  +",
            tkget(costs.32), "*m3*m4*m5+", tkget(costs.42), "*m4*m5+", tkget(costs.52), "*m5",sep=""))}}
}
transVariation <- function(){
    if(hieroEnv$nLevs == 1){
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "/n1", sep=""))}
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "*(1/n1 +1/m1)", sep=""))}}        
    if(hieroEnv$nLevs == 2){
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "/(n1*n2)+", tkget(vars.2),"/n2", sep=""))}
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "*(1/(n1*n2)+1/(m1*m2)) \n                  ", sep=""))
            if (as.numeric(tkget(levsT))>=2){        
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.2),"*(1/n2+1/m2)", sep=""))}}}       
    if(hieroEnv$nLevs == 3){
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "/(n1*n2*n3)+", tkget(vars.2),"/(n2*n3)+", tkget(vars.3), "/n3", sep=""))}
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "*(1/(n1*n2*n3)+1/(m1*m2*m3)) \n                  ", sep=""))
            if (as.numeric(tkget(levsT))>=2){
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.2),"*(1/(n2*n3)+1/(m2*m3)) \n                  ", sep=""))}
            if (as.numeric(tkget(levsT))>=3){
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.3), "*(1/n3+1/m3)", sep=""))}}}
    if(hieroEnv$nLevs == 4){
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "/(n1*n2*n3*n4)+", tkget(vars.2),"/(n2*n3*n4) \n                  +", 
            tkget(vars.3), "/(n3*n4)+", tkget(vars.4), "/n4 ", sep=""))}
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "*(1/(n1*n2*n3*n4)+1/(m1*m2*m3*m4))\n                  ", sep=""))
            if (as.numeric(tkget(levsT))>=2){
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.2),"*(1/(n2*n3*n4)+1/(m2*m3*m4))\n                  ", sep=""))}
            if (as.numeric(tkget(levsT))>=3){
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.3), "*(1/(n3*n4)+1/(m3*m4))\n                  ",sep=""))}
            if (as.numeric(tkget(levsT))>=4){
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.4), "*(1/n4+1/m4) ", sep=""))}}}
    if(hieroEnv$nLevs == 5){
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s1"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "/(n1*n2*n3*n4*n5)+", tkget(vars.2),"/(n2*n3*n4*n5) \n                  +",
            tkget(vars.3), "/(n3*n4*n5)+", tkget(vars.4), "/(n4*n5)+", tkget(vars.5), "/n5 ", sep=""))}
        if(as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2"){
            tkinsert(formulas_txt, "end", paste(tkget(vars.1), "*(1/(n1*n2*n3*n4*n5)+1/(m1*m2*m3*m4*m5))\n                  ",sep=""))
            if (as.numeric(tkget(levsT))>=2){    
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.2),"*(1/(n2*n3*n4*n5)+1/(m2*m3*m4*m5)) \n                  ",sep=""))}
            if (as.numeric(tkget(levsT))>=3){        
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.3), "*(1/(n3*n4*n5)+1/(m3*m4*m5))\n                  ",sep=""))}
            if (as.numeric(tkget(levsT))>=4){           
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.4), "*(1/(n4*n5)+1/(m4*m5))\n                  ",sep=""))}
            if (as.numeric(tkget(levsT))>=5){           
                tkinsert(formulas_txt, "end", paste("+",tkget(vars.5), "*(1/n5+1/m5) ", sep=""))}}}
}
pingNeos <- function(){
    neosStatus_value <- Nping()
    neosStatus_v <- sub("\n", "", neosStatus_value@ans)
    neosStatus_v <- sub("is", "\n is", neosStatus_v)
    tclvalue(neosStatus) <- paste(neosStatus_v)
    print(eval(neosStatus_value))  
    }          
sentToNEOS <- function(){
    hieroEnv$err <- 0
    errorName()
    if (hieroEnv$err == 0){
        tk2notetab.select(nb, "Results")
        nc <- CreateNeosComm()
        SolverChoice <- optSolverChoices[as.numeric(tclvalue(tcl(optSolver,"getvalue")))+1]
        if (SolverChoice == "AlphaECP"){
            solver.template <- NgetSolverTemplate(category = "minco", solvername = SolverChoice, inputMethod = "GAMS")}
        else solver.template <- NgetSolverTemplate(category = "go", solvername = SolverChoice, inputMethod = "GAMS")    
        modct <- paste(paste(tclvalue(tkget(formulas_txt,"0.0","end")), collapse = "\n"), "\n")
        argslist <- list(model = modct, options = "", wantlog = "", comments = "")
        xmls2 <- CreateXmlString(neosxml = solver.template, cdatalist = argslist)
        hieroEnv$jobID <- NsubmitJob(xmlstring = xmls2, user = "rneos", interface = "", id = 0)
        tkinsert(model_txt,"end", paste("Job was sent to NEOS. \n-------------------------------- \n"))
        return(hieroEnv$jobID)}
    }
getJobInfo <- function(){
    Jinfo <- NgetJobInfo(hieroEnv$jobID, convert = TRUE)
    tkinsert(model_txt,"end", paste("Job info for ", hieroEnv$newProbName, ": "))  
    tkinsert(model_txt,"end", paste(Jinfo@ans))
    tkinsert(model_txt,"end", paste(" \n-------------------------------- \n"))  
    print(eval(Jinfo))
    }
getJobStatus <- function(){
    Jstatus <- NgetJobStatus(hieroEnv$jobID)
    tkinsert(model_txt,"end", paste("Job status for ", hieroEnv$newProbName, ": "))  
    tkinsert(model_txt,"end", paste(Jstatus@ans))
    tkinsert(model_txt,"end", paste("\n-------------------------------- \n"))  
    print(eval(Jstatus))
    }
getResults <- function(){
    Fresults <- NgetFinalResults(hieroEnv$jobID)
    hieroEnv$FresultsTemp <- gsub("[[:space:]]","", Fresults@ans)
    print(eval(Fresults))
    hieroEnv$jobSource <- "new"
    errorCosts()
    organizeResults()
    }
getJobInfo_old <- function(){
    Jstatus_old <- NgetJobStatus(hieroEnv$jobID_old, convert = TRUE)
    tkinsert(model_txt,"end", paste("Status of ", hieroEnv$nameOld, ": ", sep=""))  
    tkinsert(model_txt,"end", paste(Jstatus_old@ans))
    tkinsert(model_txt,"end", paste("\n-------------------------------- \n"))   
    print(eval(Jstatus_old))
    }
getJobResults_old <- function(){
    Fresults_old <- NgetFinalResults(hieroEnv$jobID_old)
    hieroEnv$FresultsTemp <- gsub("[[:space:]]","", Fresults_old@ans)
    print(eval(Fresults_old))
    hieroEnv$jobSource <- "old"
    errorCosts()
    organizeResults()
}    
organizeResults <- function(){
    if(regexpr("divisionbyzero",hieroEnv$FresultsTemp)[1] > 0){
            tkinsert(model_txt,"end", paste("\n \nError found. \nPlease check your model constraints and R console for more information!"))    
            tkinsert(model_txt,"end", paste("\n--------------------------------\n"))}    
    if(regexpr("divisionbyzero",hieroEnv$FresultsTemp)[1] <= 0){
    if (regexpr("m1.L",hieroEnv$FresultsTemp)[1] > 0)
    num.groups <- 2 else num.groups <- 1
    if (regexpr("n1.L",hieroEnv$FresultsTemp)[1] > 0){
        num.levels <- 1
        if (regexpr("n2.L",hieroEnv$FresultsTemp)[1] > 0){
            num.levels <- 2
            if (regexpr("n3.L",hieroEnv$FresultsTemp)[1] > 0){
                num.levels <- 3
                if (regexpr("n4.L",hieroEnv$FresultsTemp)[1] > 0){
                    num.levels <- 4
                    if (regexpr("n5.L",hieroEnv$FresultsTemp)[1] > 0){
                        num.levels <- 5}}}}}
    if(regexpr("obj.ct",hieroEnv$FresultsTemp)[1] > 0){
        objective <- "obj.ct"}
    if(regexpr("obj.ce",hieroEnv$FresultsTemp)[1] > 0){
        objective <- "obj.ce"}
    if(regexpr("obj.vt",hieroEnv$FresultsTemp)[1] > 0){
        objective <- "obj.vt"}
    if(regexpr("obj.ve",hieroEnv$FresultsTemp)[1] > 0){
        objective <- "obj.ve"}                 
    if(num.levels == 1){
        if (num.groups == 1){
            if(objective == "obj.ct" || objective == "obj.ce"){
                hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
            if(objective == "obj.vt" || objective == "obj.ve"){    
                hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}
        if(num.groups == 2){
            hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m1.L",hieroEnv$FresultsTemp)[1]-9))
        if(objective == "obj.ct" || objective == "obj.ce"){
            hieroEnv$res.m1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
        if(objective == "obj.vt" || objective == "obj.ve"){    
            hieroEnv$res.m1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}}
    if(num.levels == 2){
        if (num.groups == 1){
            hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n2.L",hieroEnv$FresultsTemp)[1]-9))
            if(objective == "obj.ct" || objective == "obj.ce"){
                hieroEnv$res.n2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
            if(objective == "obj.vt" || objective == "obj.ve"){    
                hieroEnv$res.n2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}
        if(num.groups == 2){
            hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m1.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.m1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n2.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m2.L",hieroEnv$FresultsTemp)[1]-9)) 
            if(objective == "obj.ct" || objective == "obj.ce"){
                hieroEnv$res.m2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
            if(objective == "obj.vt" || objective == "obj.ve"){    
                hieroEnv$res.m2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}}
    if(num.levels == 3){
        if (num.groups == 1){
            hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n2.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n3.L",hieroEnv$FresultsTemp)[1]-9)) 
            if(objective == "obj.ct" || objective == "obj.ce"){
                hieroEnv$res.n3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
            if(objective == "obj.vt" || objective == "obj.ve"){    
                hieroEnv$res.n3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}
        if(num.groups == 2){
            hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m1.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.m1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n2.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m2.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.m2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n3.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m3.L",hieroEnv$FresultsTemp)[1]-9)) 
            if(objective == "obj.ct" || objective == "obj.ce"){
                hieroEnv$res.m3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
            if(objective == "obj.vt" || objective == "obj.ve"){    
                hieroEnv$res.m3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}}
    if(num.levels == 4){
        if (num.groups == 1){
            hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n2.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n3.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.n3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n4.L",hieroEnv$FresultsTemp)[1]-9)) 
            if(objective == "obj.ct" || objective == "obj.ce"){
                hieroEnv$res.n4 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n4.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
            if(objective == "obj.vt" || objective == "obj.ve"){    
                hieroEnv$res.n4 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n4.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}
        if(num.groups == 2){
            hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m1.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.m1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n2.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m2.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.m2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n3.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m3.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.m3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n4.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n4 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n4.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m4.L",hieroEnv$FresultsTemp)[1]-9)) 
            if(objective == "obj.ct" || objective == "obj.ce"){
                hieroEnv$res.m4 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m4.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
            if(objective == "obj.vt" || objective == "obj.ve"){    
                hieroEnv$res.m4 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m4.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}}
    if(num.levels == 5){
        if (num.groups == 1){
            hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n2.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n3.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.n3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n4.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.n4 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n4.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n5.L",hieroEnv$FresultsTemp)[1]-9)) 
            if(objective == "obj.ct" || objective == "obj.ce"){
                hieroEnv$res.n5 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n5.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
            if(objective == "obj.vt" || objective == "obj.ve"){    
                hieroEnv$res.n5 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n5.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}
        if(num.groups == 2){
            hieroEnv$res.n1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m1.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.m1 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m1.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n2.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m2.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.m2 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m2.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n3.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m3.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.m3 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m3.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n4.L",hieroEnv$FresultsTemp)[1]-9))
            hieroEnv$res.n4 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n4.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m4.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.m4 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m4.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("n5.L",hieroEnv$FresultsTemp)[1]-9)) 
            hieroEnv$res.n5 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("n5.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("m5.L",hieroEnv$FresultsTemp)[1]-9)) 
            if(objective == "obj.ct" || objective == "obj.ce"){
                hieroEnv$res.m5 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m5.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLECOSTS.L",hieroEnv$FresultsTemp)[1]-1))}
            if(objective == "obj.vt" || objective == "obj.ve"){    
                hieroEnv$res.m5 <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("m5.L",hieroEnv$FresultsTemp)[1]+5, stop=regexpr("VARIABLEVARIANCE.L",hieroEnv$FresultsTemp)[1]-1))}}}
    if(regexpr("COSTS.L",hieroEnv$FresultsTemp)[1] > 0){
        hieroEnv$res.obj <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("COSTS.L",hieroEnv$FresultsTemp)[1]+8, stop=regexpr("EQUATIONCON1.L",hieroEnv$FresultsTemp)[1]))}
    if(regexpr("VARIANCE.L",hieroEnv$FresultsTemp)[1] > 0){ 
        hieroEnv$res.obj <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("VARIANCE.L",hieroEnv$FresultsTemp)[1]+11, stop=regexpr("EQUATIONCON1.L",hieroEnv$FresultsTemp)[1]))}
    hieroEnv$res.con <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("CON1.L",hieroEnv$FresultsTemp)[1]+7, stop=regexpr("....END....",hieroEnv$FresultsTemp)[1]-1))
    hieroEnv$res.groups <- as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("groups.l",hieroEnv$FresultsTemp)[1]+9, stop=regexpr("groups.l",hieroEnv$FresultsTemp)[1]+10))
    hieroEnv$res.alpha <-as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("size.l",hieroEnv$FresultsTemp)[1]+7, stop=regexpr("_1_",hieroEnv$FresultsTemp)[1]-1))
    if(regexpr("Delta.l",hieroEnv$FresultsTemp)[1] > 0){
        hieroEnv$res.Delta <-as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("Delta.l",hieroEnv$FresultsTemp)[1]+8, stop=regexpr("_2_",hieroEnv$FresultsTemp)[1]-1))}
    if(regexpr("power.l",hieroEnv$FresultsTemp)[1] > 0){    
        hieroEnv$res.power <-as.numeric(substr(hieroEnv$FresultsTemp, start=regexpr("power.l",hieroEnv$FresultsTemp)[1]+8, stop=regexpr("_3_",hieroEnv$FresultsTemp)[1]-1))}
    if(hieroEnv$jobSource == "old"){   
        tkinsert(model_txt,"end", paste("Results for ", hieroEnv$nameOld, ": \n \n"))} 
    if(hieroEnv$jobSource == "new"){
        tkinsert(model_txt,"end", paste("Results for ", hieroEnv$newProbName, ": \n"))}
    if (num.levels == 5){ 
        tkinsert(model_txt,"end", paste("\n n5 = ", hieroEnv$res.n5))
        if (num.groups == 2){
            tkinsert(model_txt,"end", paste("\t\t m5 = ", hieroEnv$res.m5))}}   
    if (num.levels >= 4){        
        tkinsert(model_txt,"end", paste("\n n4 = ", hieroEnv$res.n4))
        if (num.groups == 2){
            tkinsert(model_txt,"end", paste("\t\t m4 = ", hieroEnv$res.m4))}}
    if (num.levels >= 3){
        tkinsert(model_txt,"end", paste("\n n3 = ", hieroEnv$res.n3))
        if (num.groups == 2){
            tkinsert(model_txt,"end", paste("\t\t m3 = ", hieroEnv$res.m3))}}
    if (num.levels >= 2){
        tkinsert(model_txt,"end", paste("\n n2 = ", hieroEnv$res.n2))
        if (num.groups == 2){
            tkinsert(model_txt,"end", paste("\t\t m2 = ", hieroEnv$res.m2))}}
    tkinsert(model_txt,"end", paste("\n n1 = ", hieroEnv$res.n1))    
    if (num.groups == 2){
        tkinsert(model_txt,"end", paste("\t\t m1 = ", hieroEnv$res.m1))}
    tkinsert(model_txt,"end", paste("\n\n Total number of units:")) 
    if (num.levels == 1){
        tkinsert(model_txt,"end", paste("\n N = ", hieroEnv$res.n1))
        if (num.groups == 2){
            tkinsert(model_txt,"end", paste("\t\t M = ", hieroEnv$res.m1))}}
    if (num.levels == 2){
        tkinsert(model_txt,"end", paste("\n N = ", hieroEnv$res.n1*hieroEnv$res.n2))
        if (num.groups == 2){
            tkinsert(model_txt,"end", paste("\t\t M = ", hieroEnv$res.m1*hieroEnv$res.m2))}}            
    if (num.levels == 3){
        tkinsert(model_txt,"end", paste("\n N = ", hieroEnv$res.n1*hieroEnv$res.n2*hieroEnv$res.n3))
        if (num.groups == 2){
            tkinsert(model_txt,"end", paste("\t\t M = ", hieroEnv$res.m1*hieroEnv$res.m2*hieroEnv$res.m3))}} 
    if (num.levels == 4){
        tkinsert(model_txt,"end", paste("\n N = ", hieroEnv$res.n1*hieroEnv$res.n2*hieroEnv$res.n3*hieroEnv$res.n4))
        if (num.groups == 2){
            tkinsert(model_txt,"end", paste("\t\t M = ", hieroEnv$res.m1*hieroEnv$res.m2*hieroEnv$res.m3*hieroEnv$res.m4))}}
    if (num.levels == 5){
        tkinsert(model_txt,"end", paste("\n N = ", hieroEnv$res.n1*hieroEnv$res.n2*hieroEnv$res.n3*hieroEnv$res.n4*hieroEnv$res.n5))
        if (num.groups == 2){
            tkinsert(model_txt,"end", paste("\t\t M = ", hieroEnv$res.m1*hieroEnv$res.m2*hieroEnv$res.m3*hieroEnv$res.m4*hieroEnv$res.m5))}}                                          
    if(objective == "obj.ct"){
        tkinsert(model_txt,"end", paste("\n\n Total costs \t\t=\t", hieroEnv$res.obj))
        tkinsert(model_txt,"end", paste("\n Power \t\t=\t", round(calcPower(ncpar=hieroEnv$res.con),digits=5)))}
    if(objective == "obj.ce"){
        tkinsert(model_txt,"end", paste("\n\n Total costs \t\t=\t", hieroEnv$res.obj))
        tkinsert(model_txt,"end", paste("\n Variance of the estimate \t=\t", hieroEnv$res.con))}
    if(objective == "obj.vt"){
        tkinsert(model_txt,"end", paste("\n\n Power \t\t=\t", round(calcPower(ncpar=hieroEnv$res.obj),digits=3)))
        tkinsert(model_txt,"end", paste("\n Total costs \t\t=\t", hieroEnv$res.con))}
    if(objective == "obj.ve"){
        tkinsert(model_txt,"end", paste("\n\n Variance of the estimate \t=\t", hieroEnv$res.obj))
        tkinsert(model_txt,"end", paste("\n Total costs \t\t=\t", hieroEnv$res.con))}
    tkinsert(model_txt,"end", paste("\n \n Check R console for more information!"))    
    tkinsert(model_txt,"end", paste("\n--------------------------------\n"))}
} 
errorSamples <- function(){
    if (as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s1" && as.character(tclvalue(hieroEnv$rbValue1)) != "obj.s2" && hieroEnv$err == 0){
        tkmessageBox(message="Please select number \n of samples!",icon="warning")
        hieroEnv$err <- hieroEnv$err+1}
}
errorTrLevel <- function(){
    if (length(as.numeric(tkget(levsT))) == 0 && as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2" && hieroEnv$err == 0){
        tkmessageBox(message="Please give \n level of grouping!",icon="warning")
        hieroEnv$err <- hieroEnv$err+1}
}
errorNumeric <- function(){
    if (hieroEnv$err == 0){
        if(is.na(as.numeric(hieroEnv$nLevs))){
            tkmessageBox(message="Please check number of levels!",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}
        if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct"){
            if(is.na(as.numeric(tkget(Delta)))){ 
                tkmessageBox(message="Please check effect size!",icon="warning")
                hieroEnv$err <- hieroEnv$err+1}
            else if(as.numeric(tkget(Delta)) == 0){
                tkmessageBox(message="Please check effect size!",icon="warning")
                hieroEnv$err <- hieroEnv$err+1}}  
        if(is.na(as.numeric(tkget(size))) || as.numeric(tkget(size)) > 0.9999 || as.numeric(tkget(size)) == 0){
            tkmessageBox(message="Please check Size (alpha)!",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}    
        if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct" && (as.numeric(tkget(Limpow)) == 0 || as.numeric(tkget(Limpow)) > 0.9999 || is.na(as.numeric(tkget(Limpow))))){
            tkmessageBox(message="Please check Power (1-beta)!",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}
        if (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce" && (as.numeric(tkget(LimVC)) == 0 || is.na(as.numeric(tkget(LimVC))))){
            tkmessageBox(message="Please check Limit for variance of the estimate",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}
        if ((as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve") && (as.numeric(tkget(LimVC)) == 0 || is.na(as.numeric(tkget(LimVC))))){
            tkmessageBox(message="Please check Limit for total costs",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}}
}
errorVariance <- function(){
    for (i in 1:hieroEnv$nLevs){
    if ((as.numeric(tkget(get(paste("vars.",i,sep="")))) == -1 || is.na(as.numeric(tkget(get(paste("vars.",i,sep="")))))) && hieroEnv$err == 0){ 
        hieroEnv$err <- hieroEnv$err+1
        tkmessageBox(message="Please check variance estimates for each level!",icon="warning")}}
}
errorCosts <- function(){
    for (i in 1:hieroEnv$nLevs){
    if ((as.numeric(tkget(get(paste("costs.",i,sep="")))) == -1 || is.na(as.numeric(tkget(get(paste("costs.",i,sep="")))))) && hieroEnv$err == 0){ 
        hieroEnv$err <- hieroEnv$err+1
        tkmessageBox(message="Please check costs for the levels of the first group!",icon="warning")}
    if (as.character(tclvalue(hieroEnv$rbValue1)) == "obj.s2" && (as.numeric(tkget(get(paste("costs.",i,"2",sep="")))) == -1 || is.na(as.numeric(tkget(get(paste("costs.",i, "2",sep=""))))))&& hieroEnv$err == 0){        
         hieroEnv$err <- hieroEnv$err+1
        tkmessageBox(message="Please check costs for the levels of the second group!",icon="warning")}}
}
errorObj <- function(){
    if (hieroEnv$err == 0){
        if (as.character(tclvalue(hieroEnv$rbValue2)) != "obj.ct" && as.character(tclvalue(hieroEnv$rbValue2)) != "obj.ce" && as.character(tclvalue(hieroEnv$rbValue2)) != "obj.vt" && as.character(tclvalue(hieroEnv$rbValue2)) != "obj.ve"){
            tkmessageBox(message="Please set your objective!",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}
        if (length(as.numeric(tkget(LimVC))) == 0 && (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.vt" || as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ve")){
            tkmessageBox(message="Please give limit for \n total costs!",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}
        if (length(as.numeric(tkget(Limpow))) == 0 && (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ct")){
            tkmessageBox(message="Please give limit \n for power!",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}
        if (length(as.numeric(tkget(LimVC))) == 0 && (as.character(tclvalue(hieroEnv$rbValue2)) == "obj.ce")){
            tkmessageBox(message="Please give limit for \n total variance!",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}}
}
errorConstr <- function(){
    if (length(paste(tkget(get(paste("funct_", hieroEnv$crows, sep = ""))))) == 0 && hieroEnv$err == 0){
        tkmessageBox(message="Please delete empty constraints!",icon="warning")
            hieroEnv$err <- hieroEnv$err+1}
}
errorName <- function(){
    if (length(as.numeric(tkget(problemName))) == 0){
        tkmessageBox(message="Please give name to the optimization problem!",icon="warning")
        hieroEnv$err <- hieroEnv$err+1}} 
} 
HierO()
