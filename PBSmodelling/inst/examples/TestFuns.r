#TestFuns-------------------------------2013-07-03
# GUI menu to test various PBSmodelling functions.
#-------------------------------------------RH/ACB

local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

TestFuns <- function(funs=ls(pos=grep("package:PBSmodelling",search()))){
	setPBSoptions("testWins", list(
		createVector = "vector",
		focusWin = c("master","slave","drudge","employee"),
		testWidgets = c("widWin","testW"),
		chooseWinVal = "choisir",
		closeWin = c("window","vector","testW","widWin") ) )
	#is1 <- ifelse(length(funs)==1,TRUE,FALSE);  # --- NOT USED ?
	resetGraph(); getWinVal(scope="L"); print(funs)
	testWins <- getPBSoptions("testWins")
	if (any(funs=="testWidgets")) closeWin(setdiff(sort(unique(unlist(testWins))),"window"))
	else closeWin(setdiff(sort(unique(unlist(testWins[setdiff(names(testWins),"createVector")]))),"window"))
	#if (!any(funs=="createVector")) closeWin("vector");
	#if (!any(funs=="focusWin")) closeWin(c("master","slave","drudge","employee"));
	#if (!any(funs=="chooseWinVal")) closeWin(c("choisir"));
	if (any(funs=="closeWin")){
		#closeWin(c("window","vector","testW","widWin"));
		setWinVal(list(eN = 0, wtxt = "No examples chosen"), winName = "runE") }
	if (any(funs=="addArrows")) {
		tt=seq(from=-5,to=5,by=0.01)
		plot(sin(tt), cos(tt)*(1-sin(tt)), type="l")
		addArrows(0.2,0.5,0.8,0.5)
		addArrows(0.8,0.95,0.95,0.55, col="#FF0066")
	}
	if (any(funs=="addLabel")) {
		addLabel(0.75,seq(from=0.9,to=0.1,by=-0.10),c('a','b','c'), col="#0033AA")
	}
	if (any(funs=="addLegend")) {
		n <- sample(1:length(colors()),20); clrs <- colors()[n]
		addLegend(.2,.9,fill=clrs,leg=clrs,cex=1.5)
	}
	if (any(funs=="calcFib")) {
		addLabel(.3,.9,paste(paste(1:20,calcFib(20,20),sep=" : "),collapse="\n"),cex=1.5,adj=c(0,1))
	}
	if (any(funs=="calcGM")) {
		ex1 <- "calcGM(c(0,1,100))"; x1 <- eval(parse(text=ex1))
		ex2 <- "calcGM(c(0,1,100),offset=0.01,exzero=FALSE)"; x2 <- eval(parse(text=ex2))
		addLabel(.5,.8,paste(c(ex1,x1,ex2,x2),collapse="\n\n"),cex=1.2)
	}
	if (any(funs=="calcMin")) {
		Ufun <- function(P) {
			Linf <- P[1]; K <- P[2]; t0 <- P[3]; obs <- afile$len
			pred <- Linf * (1 - exp(-K*(afile$age-t0)))
			n <- length(obs); ssq <- sum((obs-pred)^2 )
			return(n*log(ssq)) }
		afile <- data.frame(age=1:16,len=c(7.36,14.3,21.8,27.6,31.5,35.3,39,41.1,43.8,45.1,47.4,48.9,
			50.1,51.7,51.7,54.1))
		pvec <- data.frame(val=c(70,0.5,0),min=c(40,0.01,-2),max=c(100,2,2),active=c(TRUE,TRUE,TRUE),
			row.names=c("Linf","K","t0"),stringsAsFactors=FALSE)
		alist <- calcMin(pvec=pvec,func=Ufun,method="nlm",steptol=1e-4,repN=10)
		print(alist[-1]); P <- alist$Pend
		expandGraph()
		xnew <- seq(afile$age[1],afile$age[nrow(afile)],len=100)
		ynew <- P[1] * (1 - exp(-P[2]*(xnew-P[3])) )
		plot(afile); lines(xnew,ynew,col="red",lwd=2)
		addLabel(.05,.88,paste(paste(c("Linf","K","t0"),round(P,c(2,4,4)),sep=" = "),collapse="\n"),adj=0,cex=0.9)
	}
	if (any(funs=="chooseWinVal")) {
		dfnam <-
		c("airquality","attitude","ChickWeight","faithful","freeny",
		"iris","LifeCycleSavings","longley","morley","Orange",
		"quakes","randu","rock","stackloss","swiss","trees")
		wlist <- c(
		"window name=\"choisir\" title=\"Test chooseWinVal\"",
		"label text=\"Press <ENTER> in the green entry box
		\\nto choose a file, then press <GO>\" sticky=W pady=5",
		"grid 1 3 sticky=W",
		"label text=File: sticky=W",
		"entry name=fnam mode=character width=23 value=\"\" 
		func=.win.chFile entrybg=darkolivegreen1 pady=5",
		"button text=GO bg=green sticky=W func=.win.chTest",
		"")
		chFile <- function(ch=dfnam,fn="fnam") {chooseWinVal(ch,fn,winname="choisir")}
		tput(chFile)
		chTest <- function() {
			getWinVal(winName="choisir",scope="L")
			if (fnam!="" && any(fnam==dfnam)) {
				file <- get(fnam)
				pairs(file,gap=0) }
			else {
				resetGraph()
				addLabel(.5,.5,"Press <ENTER> in the green entry box\nto choose a file, then press <GO>",
					col="red",cex=1.5)}}
		tput(chTest)
		createWin(wlist,astext=TRUE); chTest()
	}
	if (any(funs=="cleanProj")) {
		cleanProj(prefix="TestFuns", suffix=c("_junk.aaa","_junk.bbb"), files=c("TestFuns_JUNK1", "TestFuns_JUNK2"))
	}
	if (any(funs=="createVector")) {
		createVector(c(m=2,n=3,phi=0,k=1000),vectorLabels=c("x cycles","y cycles", "y phase", "points"),
			#func=".win.tget",action="drawLiss",windowname="vector")
			func="drawLiss",windowname="vector")
	}
	if (any(funs=="drawBars")) {
		plot(0:10,0:10,type="n")
		drawBars(x=1:9,y=9:1,col="deepskyblue4",lwd=3)
	}
	if (any(funs=="expandGraph")) {
		expandGraph(mfrow=c(2,1))
		tt=seq(from=-10, to=10, by=0.05)
		plot(tt,sin(tt), xlab="this is the x label",  ylab="this is the y label", 
			main="main title", sub="sometimes there is a \"sub\" title")
		plot(cos(tt),sin(tt*2), xlab="cos(t)", ylab="sin(2 t)", main="main title", 
			sub="sometimes there is a \"sub\" title")
	}
	if (any(funs=="findPat")) {
		junk<-findPat(c("[aeoiy]", "^[0-9]"), c("hello", "WRLD", "11b"))
		addLabel(.5,.7,"Find [aeoiy] or [1-9] in",cex=1.5)
		addLabel(.5,.5,paste("c(",paste(c("hello", "WRLD", "11b"),collapse=", "),")"),cex=1.5)
		addLabel(.5,.3,paste(junk,collapse=" "),cex=1.5,col="blue")
	}
	if (any(funs=="focusWin")) {
		createWin(c("window name=master title=Master onclose=.win.closeSDE","grid 4 1 padx=25", 
			"radio name=it value=0 text=master sticky=W function=grab", 
			"radio name=it value=1 text=slave sticky=W function=grab", 
			"radio name=it value=2 text=drudge sticky=W function=grab",
			"radio name=it value=3 text=employee sticky=W function=grab"),astext=TRUE)
		createWin(c("window name=slave title=Slave", "button text=\"Yes master\" font=\"times italic bold\" fg=darkgreen padx=25 pady=10 function=grab"), astext=TRUE)
		createWin(c("window name=drudge title=Drudge", "button text=\"Who cares\" font=\"times italic bold\" fg=blue padx=25 pady=10 function=grab"), astext=TRUE)
		createWin(c("window name=employee title=Employee", "button text=\"I`m on a break\" font=\"times italic bold\" fg=red padx=25 pady=10 function=grab"), astext=TRUE)
	}
	if (any(funs=="genMatrix")) {
		plotBubbles(genMatrix(20,6))
	}
	if (any(funs=="getChoice")) {
		getChoice(c("Fame","Fortune","Health","Beauty","Lunch"),"What do you want?",qcolor="red",gui=TRUE)
		cat("\n====================\nYou chose: ")
	}
	if (any(funs=="getYes")) {
		if(getYes("Click Yes or No:"))
			showAlert("YES :)")
		else
			showAlert("NO :(")
	}
	if (any(funs=="GT0")) {
		plotGT0 <- function(eps=1,x1=-2,x2=10,n=1000,col="black") {
			x <- seq(x1,x2,len=n); y <- GT0(x,eps)
			lines(x,y,col=col); invisible(list(x=x,y=y)) }
		testGT0 <- function(eps=c(7,5,3,1,.1),x1=-2,x2=10,n=1000) {
			x <- seq(x1,x2,len=n); y <- x
			plot(x,y,type="l")
			mycol <- c("red","blue","green","brown","violet","orange","pink")
			for (i in 1:length(eps)) 
				plotGT0(eps=eps[i],x1=x1,x2=x2,n=n,col=mycol[i])
			invisible() }
		testGT0()
	}
	if (any(funs=="openFile")) {
		sys=Sys.info()[["sysname"]]
		if (sys=="Windows") prog="notepad.exe"
		else if (sys=="Darwin" && .Platform$GUI=="AQUA") prog="TextEdit"
		else prog="gedit"
		setPBSext("html", paste(prog,"%f") )
		openFile("TestFuns_openFile.html")
		clearPBSext("html")
	}
	if (any(funs=="openProjFiles")) {
		openProjFiles("vonB", c(".r", "data.txt", "Win.txt"))
	}
	if (any(funs=="pad0")) {
		x <- pad0(x=123,n=10,f=0:7)
		addLabel(.5,.5,paste(x,collapse="\n"),cex=1.5)
	}
	if (any(funs=="pickCol")) {
		junk<-pickCol()
		addLabel(.5,.5,junk,cex=4,col=junk)
	}
	if (any(funs=="plotACF")) {
		plotACF(trees,lwd=2,lags=30)
	}
	if (any(funs=="plotAsp")) {
		x <- seq(0,10,0.1); y <- sin(x)
		par(mfrow=2:1)
		plotAsp(x,y,asp=1,xlim=c(0,10),ylim=c(-2,2), main="sin(x)")
		plotAsp(x,y^2,asp=1,xlim=c(0,10),ylim=c(-2,2), main="sin^2(x)")
	}
	if (any(funs=="plotBubbles")) {
		plotBubbles(genMatrix(20,6))
	}
	if (any(funs=="plotCsum")) {
		x <- rgamma(n=1000,shape=2)
		plotCsum(x)
	}
	if (any(funs=="plotDens")) {
		z <- data.frame(y1=rnorm(50,sd=2),y2=rnorm(50,sd=1),y3=rnorm(50,sd=.5))
		plotDens(z,lwd=3)
	}
	if (any(funs=="plotFriedEggs")) {
		x=rnorm(5000,10,3); y=-x+rnorm(5000,1,4); z=x+rnorm(5000,1,3)
		A=data.frame(x=x,y=y,z=z)
		plotFriedEggs(A,eggs=TRUE,rings=FALSE);  cat("Here are the eggs...\n")
		oldAsk=par("ask"); par(ask=TRUE)
		plotFriedEggs(A,eggs=FALSE,rings=TRUE);  cat("Here are the rings...\n")
		plotFriedEggs(A,eggs=FALSE,rings=FALSE); cat("And the pepper alone. END.\n")
		par(ask=oldAsk)
	}
	if (any(funs=="plotTrace")) {
		z <- data.frame(x=1:50,y1=rnorm(50,sd=3),y2=rnorm(50,sd=1),y3=rnorm(50,sd=.25))
		plotTrace(z,lwd=3)
	}
	if (any(funs=="promptOpenFile")) {
		scan(promptOpenFile(),what=character(),sep="\n")
	}
	if (any(funs=="promptSaveFile")) {
		promptSaveFile("trees")
	}
	if (any(funs=="readList")) {
		junk <- readList("TestFuns_Data-readList.txt"); print(junk)
	}
	if (any(funs=="resetGraph")) {
		addLabel(.5,.5,"Extent of new plot space",col="darkgreen",cex=1.5); box()
	}
	if (any(funs=="restorePar")) {
		pvec <- data.frame(val=c(1,100,10000),min=c(0,0,0),max=c(5,500,50000),active=c(TRUE,TRUE,TRUE))
		S <- c(.5,.5,.5); P <- restorePar(S,pvec); print(cbind(pvec,S,P))
	}
	if (any(funs=="runExamples")) {
	runExamples()
	}
	if (any(funs=="scalePar")) {
		pvec <- data.frame(val=c(1,100,10000),min=c(0,0,0),max=c(5,500,50000),active=c(TRUE,TRUE,TRUE))
		S <- scalePar(pvec); print(cbind(pvec,S))
	}
	if (any(funs=="setPBSext")) {
		setPBSext("goo", "notepad.exe %f"); openFile("TestFuns_setPBSext.goo")
	}
	if (any(funs=="show0")) {
		frame(); par(cex=2)
		addLabel(0.25,0.75,show0(15.2,4)); addLabel(0.25,0.7,show0(15.1,4)); addLabel(0.25,0.65,show0(15,4))
		addLabel(0.25,0.55,show0(15.2,4,TRUE)); addLabel(0.25,0.5,show0(15.1,4,TRUE)); addLabel(0.25,0.45,show0(15,4,TRUE))
	}
	if (any(funs=="showArgs")) {
		showArgs()
	}
	if (any(funs=="testAlpha")) {
		testAlpha()
	}
	if (any(funs=="testCol")) {
		#testCol(c("plum","tomato","olive","peach","honeydew"))
		testCol()
	}
	if (any(funs=="testLty")) {
		testLty(newframe = TRUE)
	}
	if (any(funs=="testLwd")) {
		testLwd(lwd=1:20, col=c("black","blue"), newframe=TRUE)
	}
	if (any(funs=="testPch")) {
		testPch(pch=1:100, ncol=10, grid=TRUE, newframe=TRUE)
	}
	if (any(funs=="testWidgets")) {
		testWidgets()
	}
	if (any(funs=="unpackList")) {
		x<-list(a=view(swiss),b=view(trees)); unpackList(x)
		print("Object a:"); print(a); print(""); print("Object b:"); print(b)
	}
	if (any(funs=="view")) {
		print("View(swiss)"); print(view(swiss))
	}
	if (any(funs=="viewCode")) {
		viewCode("PBSmodelling")
	}
	if (any(funs=="writeList")) {
		writeList(list(trees=view(trees),swiss=view(swiss)),fname="wList.txt",format="P")
		openFile("wList.txt")
	}
} # end TestFuns

#Auxilliary functions-----------------------
drawLiss <- function() {
	getWinVal(scope="L");
	ti <- 2*pi*(0:k)/k;;  x <- sin(2*pi*m*ti);  y <- sin(2*pi*(n*ti+phi));
	plot(x,y,type="l");
	invisible(NULL); }
grab <- function() {
	it <- getWinVal(scope="L")$it;
	if (is.null(it) || it==0){
		focusWin("master"); setWinVal(list(it=0),winName="master"); }
	else focusWin(switch(it,"slave","drudge","employee"));
	invisible()}
# Close the slave/drudge/employee windows controlled by master
closeSDE <- function(){ closeWin(getPBSoptions("testWins")[["focusWin"]][-1]) }
# Close all the windows associated with TestFuns
closeALL <- function(){ closeWin(sort(unique(unlist(getPBSoptions("testWins"))))) }

#require(PBSmodelling)
createWin("TestFunsWin.txt")

}) # end local scope
