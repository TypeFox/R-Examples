
toLatex.yuima <- function (object, ...) 
{
    mod <- NULL
    if (class(object) == "yuima.model") 
	mod <- object
    if (class(object) == "yuima.carma") 
  mod <- object
	if (class(object) == "yuima.cogarch") 
	  mod <- object
    if (class(object) == "yuima") 
	mod <- object@model
    #if(class(mod) =="yuima.carma" && length(mod@info@lin.par)==0 )
      if((class(mod) =="yuima.carma") || (class(mod) =="yuima.cogarch")  )
      { 
#         yuima.warn("")
        
        
        mysymb <- c("*", "alpha", "beta", "gamma", "delta", "rho", 
                    "theta","sigma","mu", "sqrt")
        #     myrepl <- c(" \\cdot ", "\\alpha ", "\\beta ", "\\gamma ", 
        #   			"\\delta ", "\\rho ", "\\theta ", "\\sqrt ")
        myrepl <- c(" \\cdot ", "\\alpha ", "\\beta ", "\\gamma ", 
                    "\\delta ", "\\rho ", "\\theta ","\\sigma","\\mu", "\\sqrt ")
        ns <- length(mysymb)
        
        
        n.eq <- mod@equation.number
        info <- mod@info
        noise.var<-"W"
        # We construc the system that describes the CARMA(p,q) process
        
        if (!length(mod@jump.variable)==0){noise.var <- mod@jump.variable}
        dr <- paste("\\left\\{\\begin{array}{l} \n")
        main.con <- info@ma.par
        if(class(mod)=="yuima.carma"){
          if(length(info@loc.par)==0 && !length(info@scale.par)==0){
            main.con<-paste(info@scale.par,"* \\ ", info@ma.par)
          }
        
          if(!length(info@loc.par)==0 && length(info@scale.par)==0){
            main.con<-paste(info@loc.par,"+ \\ ", info@ma.par)
          }
          
          if(!length(info@loc.par)==0 && !length(info@scale.par)==0){
            main.con<-paste(info@loc.par,"+ \\ ",info@scale.par,"* \\ ", info@ma.par)
          }
        }else{
          if(class(mod)=="yuima.cogarch"){
            main.con<-paste(info@loc.par,"+ \\ ", info@ma.par)  
          }
        }
        if((class(mod) =="yuima.carma")){
          dr <- paste(dr, info@Carma.var,
                      "\\left(", sprintf("%s", mod@time.variable),"\\right) = ",main.con, "'" , 
                     info@Latent.var,"\\left(", sprintf("%s", mod@time.variable),"\\right) \\\\ \n")
        }else{
          if((class(mod) =="yuima.cogarch")){
            dr <- paste(dr, sprintf("d%s", info@Cogarch.var),
                        "\\left(", sprintf("%s", mod@time.variable),"\\right) = \\ sqrt{",info@V.var, 
                        "\\left(", sprintf("%s", mod@time.variable),"\\right)} \\ ",
                        sprintf("d%s", noise.var),"\\left(", sprintf("%s", mod@time.variable),"\\right) \\\\ \n")
            dr <- paste(dr, info@V.var,
                        "\\left(", sprintf("%s", mod@time.variable),"\\right) = ",main.con, "'" , 
                        info@Latent.var,"\\left(", sprintf("%s", mod@time.variable),"\\right) \\\\ \n")
                        
          }
        }
        
        if((class(mod) =="yuima.carma")){
          noise.latent.var <- noise.var
        }else{
          if((class(mod) =="yuima.cogarch")){
            noise.latent.var <- paste0("\\left[",noise.var,",",noise.var,"\\right]^{q}")
          }
            
        }
        dr <- paste(dr, sprintf("d%s", info@Latent.var),
                       "\\left(", sprintf("%s", mod@time.variable),"\\right)",
                       "=","A",info@Latent.var,
                        "\\left(", sprintf("%s", mod@time.variable),"\\right)",
                        sprintf("d%s", mod@time.variable),
                       "+ e",sprintf("d%s", noise.latent.var),"\\left(",
                       mod@time.variable, "\\right) \\\\ \n")
        
        dr<- paste(dr, "\\end{array}\\right.")
#11/12
        for (i in 1:ns) {
          dr <- gsub(mysymb[i], myrepl[i], dr, fixed = "TRUE")
        }
        
        body <- paste("%%% Copy and paste the following output in your LaTeX file")
        body <- c(body, paste("$$"))
        body <- c(body, dr)
        body <- c(body, paste("$$"))
        # Vector Latent Variable.
        
        body <- c(body, paste("$$"))
        if(class(mod)=="yuima.carma"){
          latent.lab0<-paste(info@Latent.var,0:(info@p-1),sep="_")
        }else{
          if(class(mod)=="yuima.cogarch"){
              latent.lab0<-paste(info@Latent.var,1:info@q,sep="_")
          }
        }

        if(length(latent.lab0)==1){latent.lab<-latent.lab0}
        if(length(latent.lab0)==2){
          latent.lab0[1]<-paste(latent.lab0[1],"(",mod@time.variable,")",",\\ ",sep="")
          latent.lab0[2]<-paste(latent.lab0[2],"(",mod@time.variable,")",sep="")
          latent.lab<-latent.lab0
        }
        if(length(latent.lab0)>2){
          latent.lab<-paste(latent.lab0[1],"(",mod@time.variable,")",
                            ",\\ ","\\ldots \\ ",
                            ",\\ ",tail(latent.lab0,n=1),
                            "(",mod@time.variable,")")
        }
        latent.lab<-paste(latent.lab,collapse="") 
        X<-paste(info@Latent.var,"(",mod@time.variable,")",
                 "=\\left[",latent.lab,
                 "\\right]'")
        
        for (i in 1:ns) {
          X <- gsub(mysymb[i], myrepl[i], X, fixed = "TRUE")
        }
        
        body <- c(body, X)
        body <- c(body, paste("$$"))
        # Vector Moving Average Coefficient.
        body <- c(body, paste("$$"))
        
        #b.nozeros <-c(0:info@q)
        
      #  ma.lab0<-paste(paste(info@ma.par,0:(info@q),sep="_"),collapse=", \\ ")
        if(class(mod)=="yuima.carma"){
          ma.lab0<-paste(info@ma.par,0:(info@q),sep="_")
        }else{
          if(class(mod)=="yuima.cogarch"){
            ma.lab0<-paste(info@ma.par,1:(info@p),sep="_")
          }
        }
        #if(length(ma.lab0)==1){ma.lab1<-ma.lab0}
        if(class(mod)=="yuima.carma"){
          if(info@q>=0 && info@q<=1){
            ma.lab1<-paste(ma.lab0,collapse=", \\ ")}
        #if(length(ma.lab0)==2){
#         if(info@q==1){
#           ma.lab0[1]<-paste(ma.lab0[1],",\\ ",sep="")
#       #    ma.lab0[2]<-paste(ma.lab0[2],"(",mod@time.variable,")",sep="")
#           ma.lab1<-ma.lab0
#         }
        #if(length(ma.lab0)>2){
          if(info@q>1){
            ma.lab1<-paste(ma.lab0[1],
                              ",\\ ","\\ldots",
                              " \\ , \\ ",tail(ma.lab0,n=1))
          }
        }else{
          if(class(mod)=="yuima.cogarch"){
            if(info@p>=0 && info@p<=2){
              ma.lab1<-paste(ma.lab0,collapse=", \\ ")
            }
            if(info@p>2){
              ma.lab1<-paste(ma.lab0[1],
                             ",\\ ","\\ldots",
                             " \\ , \\ ",tail(ma.lab0,n=1))
            }
          }  
        }
        if(class(mod)=="yuima.carma"){  
          numb.zero<-(info@p-(info@q+1))
        }else{
          if(class(mod)=="yuima.cogarch"){
            numb.zero<-(info@q-info@p)
          }
        }
        if (numb.zero==0){ma.lab <- ma.lab1}
        if (numb.zero>0&&numb.zero<=2){
          zeros<- 0*c(1:numb.zero)
          zero.el <- paste(zeros, collapse=", \\ ")
          ma.lab <- paste(ma.lab1," ,\\ ", zero.el)
        }
         if (numb.zero>2 ){
           ma.lab <- paste(ma.lab1," ,\\ 0, \\ \\ldots \\ , \\ 0")
         }         
        Vector.ma <- paste(info@ma.par,"=","\\left[",ma.lab,"\\right]'")
        
        for (i in 1:ns) {
          Vector.ma <- gsub(mysymb[i], myrepl[i], Vector.ma, fixed = "TRUE")
        }
        
        body <- c(body, Vector.ma)
        body <- c(body, paste("$$"))
        
        # e vector
        body <- c(body, paste("$$"))
        
        noise.coef<-mod@diffusion
        vect.e0 <- substr(tail(noise.coef,n=1), 13, nchar(tail(noise.coef,n=1)) -2)
        vect.e1 <- substr(vect.e0, 2, nchar(vect.e0) -1)
        dummy.e0<-strsplit(vect.e1,split="+",fixed=TRUE)
        
        
        if (!length(mod@jump.variable)==0){
          noise.coef <- mod@jump.coeff
          if(class(mod)=="yuima.carma"){
            vect.e0 <- substr(tail(noise.coef,n=1), 18, nchar(tail(noise.coef,n=1)) -2)
          }else{
            vect.e0 <- substr(tail(noise.coef,n=1), 18, nchar(tail(noise.coef,n=1)) -2)
          }
          #vect.e0 <- substr(tail(noise.coef,n=1), 2, nchar(tail(noise.coef,n=1)) -1)
        } else{ 
          if(length(info@lin.par) != 0){
                
            if (info@lin.par != info@ma.par){
              dummy.e0<-c(dummy.e0[[1]][1], paste(dummy.e0[[1]][(2:length(dummy.e0[[1]]))],
                                                 "(",mod@time.variable,")",sep=""))
              dummy.e0<-gsub(info@Latent.var,paste0(info@Latent.var,"_",collapse=""),dummy.e0)
              dummy.e0<-gsub(info@lin.par,paste0(info@lin.par,"_",collapse=""),dummy.e0)  
              if(length(dummy.e0)>3){
                vect.e0<-paste(dummy.e0[1],"+",dummy.e0[2],
                               "+ \\ldots +",tail(dummy.e0,n=1))
                vect.e0<-paste("(",vect.e0,")")
              } else{vect.e0<-paste("(",paste(dummy.e0,collapse="+"),")")}
            } 
  #           else{
  #             yuima.warm("The case of loc.par and scale.par will be implemented as soon as possible")
  #             return(NULL)
  #           }
          }  
        }
        if(class(mod)=="yuima.carma"){
          if (info@p==1){vect.e <- vect.e0}
          if (info@p==2){vect.e <- paste("0, \\ ",vect.e0)}
          if (info@p==3){vect.e <- paste("0, \\ 0, \\ ",vect.e0)}
          if (info@p>3){vect.e <- paste("0, \\ \\ldots \\ , \\ 0, \\  ",vect.e0)}
        }else{
          if(class(mod)=="yuima.cogarch"){
            if (info@q==1){vect.e <- vect.e0}
            if (info@q==2){vect.e <- paste("0, \\ ",vect.e0)}
            if (info@q==3){vect.e <- paste("0, \\ 0, \\ ",vect.e0)}
            if (info@q>3){vect.e <- paste("0, \\ \\ldots \\ , \\ 0, \\  ",vect.e0)}
          }
        }
        coeff.e<- paste("e","=","\\left[",  vect.e , "\\right]'")
        
        for (i in 1:ns) {
          coeff.e <- gsub(mysymb[i], myrepl[i], coeff.e, fixed = "TRUE")
        }
        
        
        
        
        body <- c(body, coeff.e)
        body <- c(body, paste("$$"))
        # Matrix A        
        body <- c(body, paste("$$"))

        if(class(mod)=="yuima.cogarch"){
          Up.A<-NULL
        }
        
        if(class(mod)=="yuima.carma"){
          if(info@p==1){
            cent.col<-"c"
            last.A<-paste(paste(paste("",info@ar.par,sep=" -"),info@p:1,sep="_"),collapse=" &")
          }
        }else{
          if(class(mod)=="yuima.cogarch"){ 
            if(info@q==1){
              cent.col<-"c"
              last.A<-paste(paste(paste("",info@ar.par,sep=" -"),info@q:1,sep="_"),collapse=" &")
            }
          }
        }
       
        if(class(mod)=="yuima.carma"){
            if(info@p==2){
              cent.col<-"cc"
              Up.A <-" 0 & 1 \\\\ \n"
              last.A<-paste(paste(paste("",info@ar.par,sep=" -"),info@p:1,sep="_"),collapse=" &")
            }
            
            if(info@p==3){
              cent.col<-"ccc"
              Up.A <-" 0 & 1 & 0 \\\\ \n 0 & 0 & 1 \\\\ \n"
              last.A<-paste(paste(paste("",info@ar.par,sep=" -"),info@p:1,sep="_"),collapse=" &")
              
            }
            
            if(info@p>3){
              cent.col<-"cccc"
              Up.A <-" 0 & 1 & \\ldots & 0 \\\\ \n \\vdots & \\vdots & \\ddots & \\vdots \\\\ \n 0 & 0 & \\ldots & 1 \\\\ \n"
              dummy.ar<-paste(paste("",info@ar.par,sep=" -"),info@p:1,sep="_")
              last.A <- paste(dummy.ar[1]," & ", dummy.ar[2]," & \\ldots &", tail(dummy.ar,n=1) )
            
            }
        }else{
          if(class(mod)=="yuima.cogarch"){ 
              if(info@q==2){
                cent.col<-"cc"
                Up.A <-" 0 & 1 \\\\ \n"
                last.A<-paste(paste(paste("",info@ar.par,sep=" -"),info@q:1,sep="_"),collapse=" &")
              }
              
              if(info@q==3){
                cent.col<-"ccc"
                Up.A <-" 0 & 1 & 0 \\\\ \n 0 & 0 & 1 \\\\ \n"
                last.A<-paste(paste(paste("",info@ar.par,sep=" -"),info@q:1,sep="_"),collapse=" &")
                
              }
              
              if(info@q>3){
                cent.col<-"cccc"
                Up.A <-" 0 & 1 & \\ldots & 0 \\\\ \n \\vdots & \\vdots & \\ddots & \\vdots \\\\ \n 0 & 0 & \\ldots & 1 \\\\ \n"
                dummy.ar<-paste(paste("",info@ar.par,sep=" -"),info@q:1,sep="_")
                last.A <- paste(dummy.ar[1]," & ", dummy.ar[2]," & \\ldots &", tail(dummy.ar,n=1) )
                
              }          
          }
        }

        matrix.A <-paste(Up.A ,last.A," \\\\ \n",sep="")
        
        array.start<-paste0("\\begin{array}{",cent.col,"}\n",collapse="")
        MATR.A<-paste("A ","=","\\left[",array.start,  matrix.A,  "\\end{array}\\right]" )
        
        for (i in 1:ns) {
          MATR.A <- gsub(mysymb[i], myrepl[i], MATR.A, fixed = "TRUE")
        }
        
        
        body <- c(body, MATR.A)
        body <- c(body, paste("$$"))
        body <- structure(body, class = "Latex")
        
        return(body)
        
    } else{ 
    n.eq <- mod@equation.number
    dr <- paste("\\left(\\begin{array}{c}\n")
    for (i in 1:n.eq) {
            dr <- paste(dr, substr(mod@drift[i], 2, nchar(mod@drift[i]) -1), "\\\\ \n")
        #      dr <- paste(dr, substr(mod@drift[i], 3, nchar(mod@drift[i]) - 2), "\\\\ \n")
    }
    #
    dr <- paste(dr, "\\end{array}\\right)", sprintf("d%s", mod@time.variable))
    n.n <- mod@noise.number
    df <- paste(sprintf("\\left[\\begin{array}{%s}\n",paste(rep("c",n.n),sep="",collapse="")))
    for (i in 1:n.eq) {
        #df <- paste(df, paste(mod@diffusion[[i]], collapse = "&"))
                df <- paste(df, paste(substr(mod@diffusion[[i]], 2, nchar(mod@diffusion[[i]]) - 1)  , collapse = "&"))
        df <- paste(df, "\\\\ \n")
    }
    df <- paste(df, "\\end{array}\\right]")
# We consider the Jump 6/11
    if (length(mod@jump.coeff)>=1){
      dj<-paste("\\left(\\begin{array}{c}\n")
      for (i in 1:n.eq) {
        dj <- paste(dj, substr(mod@jump.coeff[i], 2, nchar(mod@jump.coeff[i]) - 1), "\\\\ \n")
        #dj <- paste(dj, mod@jump.coeff[i], "\\\\ \n")
        
      }
      dj <- paste(dj, "\\end{array}\\right)", sprintf("d %s", mod@jump.variable))
    }
    
    wn <- paste("\\left(\\begin{array}{c}\n")
    wn <- paste(wn, paste(sprintf("dW%s", 1:n.n), sep = "", collapse = "\\\\ "))
    wn <- paste(wn, "\n \\end{array}\\right)")
    st <- paste("\\left(\\begin{array}{c}\n")
    st <- paste(st, paste(sprintf("d%s", mod@solve.variable), 
						  sep = "", collapse = "\\\\ "))
    st <- paste(st, "\n \\end{array}\\right)")
    mysymb <- c("*", "alpha", "beta", "gamma", "delta", "rho", 
				"theta","sigma","mu", "sqrt")
#     myrepl <- c(" \\cdot ", "\\alpha ", "\\beta ", "\\gamma ", 
# 				"\\delta ", "\\rho ", "\\theta ", "\\sqrt ")
    myrepl <- c(" \\cdot ", "\\alpha ", "\\beta ", "\\gamma ", 
                "\\delta ", "\\rho ", "\\theta ","\\sigma","\\mu", "\\sqrt ")
    ns <- length(mysymb)
    for (i in 1:ns) {
        dr <- gsub(mysymb[i], myrepl[i], dr, fixed = "TRUE")
        df <- gsub(mysymb[i], myrepl[i], df, fixed = "TRUE")
# for Jump         
        if (length(mod@jump.coeff)>=1){
          dj <- gsub(mysymb[i], myrepl[i], dj, fixed = "TRUE")
        }
    }
    body <- paste("%%% Copy and paste the following output in your LaTeX file")
    body <- c(body, paste("$$"))
    body <- c(body, paste(st))
    body <- c(body, paste(" = "))
    body <- c(body, paste(dr))
    body <- c(body, paste(" +"))
    body <- c(body, paste(df))
    body <- c(body, paste("'"))
    body <- c(body, paste(wn))
    # For Jump 6/11
    if (length(mod@jump.coeff)>=1){
      body <- c(body, paste(" +"))
      body <- c(body, paste(dj))
    }
    
    body <- c(body, paste("$$"))

    body <- c(body, paste("$$"))
# For Initial Conditions     
    numb.solve.var <- length(mod@solve.variable)
    bodyaus <-c( paste("\\left(\\begin{array}{c}\n"))
    for (i in 1:numb.solve.var) {
      bodyaus <-paste(bodyaus, paste(paste(mod@solve.variable[i],"(0)",sep=""),substr(mod@xinit[i], 2, nchar(mod@xinit[i]) - 1),sep="="), "\\\\ \n")
      # paste(bodyaus, paste(paste(mod@solve.variable[i],"(0)",sep=""),substr(mod@xinit[i], 3, nchar(mod@xinit[i]) - 2),sep="="), "\\\\ \n")
      #     paste(bodyaus, paste(paste(mod@solve.variable[i],"(0)",sep=""),substr(mod@xinit[i], 2, nchar(mod@xinit[i]) - 1),sep="="), "\\\\ \n")
    }
    
    bodyaus <- paste(bodyaus, "\\end{array}\\right)")
    for (i in 1:ns) {
      bodyaus <- gsub(mysymb[i], myrepl[i], bodyaus, fixed = "TRUE")
    }
    
    body<-c(body,paste(bodyaus))
    
#     body <- c(body, paste(sprintf("%s(0)=%f,\\quad", mod@solve.variable, 
# 								  mod@xinit)))
    body <- c(body, paste("$$"))
    structure(body, class = "Latex")
    }
}



toLatex.yuima.model <- toLatex.yuima 

toLatex.yuima.carma <- toLatex.yuima

toLatex.yuima.cogarch <- toLatex.yuima

