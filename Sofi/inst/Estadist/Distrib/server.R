pkg <- c("VGAM")
new.pkg <- pkg[!(pkg %in% installed.packages())]
if (length(new.pkg)) {
  install.packages(new.pkg)
}
#Distribuciones de variables aleatorias (test para UAA)

library(shiny)
library(VGAM)
load(system.file("Estadist/Distrib/samplingApp.RData", package="Sofi"), envir=.GlobalEnv)
options(shiny.deprecation.messages=FALSE)

shinyServer(function(input, output, session){

	output$distName <- renderUI({
		if(input$disttype=="Discreta"){
			radioButtons("dist","Distribución:",selected="bern",
				list("Bernoulli"="bern","Binomial"="bin","Uniforme discreta"="dunif","Geométrica"="geom","Hipergeométrica"="hgeom","Binomial negativa"="nbin","Poisson"="poi") # 
			)
		} else if(input$disttype=="Continua"){
			radioButtons("dist","Distribución:",selected="beta",
				list("Beta"="beta","Cauchy"="cauchy","Chi cuadrado"="chisq","Exponencial"="exp","F"="F","Gamma"="gam","De Laplace"="lap", # Continua
					"Logística"="logi","Log-Normal"="lognorm","Normal"="norm","De Pareto"="pareto","t de Student"="t","Uniforme continua"="unif","Weibull"="weib")
			)
		}
	})
		
	dat <- reactive({
		dist <- switch(input$dist,
			bern=rbern, bin=rbinom2, dunif=drunif, geom=rgeom2, hgeom=rhyper2, nbin=rnbinom2, poi=rpois2, # discrete
			beta=rbeta2, cauchy=rcauchy2, chisq=rchisq2, exp=rexp2, F=rf2, gam=rgamma2, lap=rlaplace2, # Continua
			logi=rlogis2, lognorm=rlnorm, norm=rnorm, pareto=rpareto2, t=rt2, unif=runif, weib=rweibull2
			)

		def.args <- switch(input$dist,
			# discrete
			bern=c(input$bern.prob),
			bin=c(input$binom.size,input$binom.prob),
			dunif=c(input$drunif.min,input$drunif.max,input$drunif.step),
			geom=c(input$geom.prob),
			hgeom=c(input$hyper.M,input$hyper.N,input$hyper.K),
			nbin=c(input$nbin.size,input$nbin.prob),
			poi=c(input$poi.lambda),
			# Continua
			beta=c(input$beta.shape1,input$beta.shape2),
			cauchy=c(input$cau.location,input$cau.scale),
			chisq=c(input$chisq.df),
			exp=c(input$exp.rate),
			F=c(input$F.df1,input$F.df2),
			gam=c(input$gam.shape,input$gam.rate),
			lap=c(input$lap.location,input$lap.scale),
			logi=c(input$logi.location,input$logi.scale),
			lognorm=c(input$meanlog,input$sdlog),
			norm=c(input$mean,input$sd),
			pareto=c(input$pareto.location,input$pareto.shape),
			t=c(input$t.df),
			unif=c(input$min,input$max),
			weib=c(input$weib.shape,input$weib.scale)
			)

		f <- formals(dist)
		f <- f[names(f)!="nn" & names(f)!="n"]
		if(any(input$dist==c("dunif","hgeom"))){ len <- min(length(f),4-1); f <- f[1:len]
		} else { len <- min(length(f),3-1); f <- f[1:len] }
		argList <- list(n=input$n)
		for(i in 1:len) argList[[names(f)[i]]] <- def.args[i]
		return(list(do.call(dist,argList),names(f)))
	})

	output$dist1 <- renderUI({
		input$dist
		isolate({
			if(length(input$dist)){
				lab <- switch(input$dist,
					bern="Probabilidad:", bin="Tamaño:", dunif="Discreto secuencia mínima:", geom="Probabilidad:", hgeom="M:", nbin="Número de éxitos:",	poi="Media y varianza:", # discrete
					beta="Alfa:", cauchy="Ubicación:", chisq="Grados de libertad:", exp="Proporción", F="Grados de libertad del numerador:", gam="Forma:", lap="Ubicación:",
					logi="Ubicación:", lognorm="Media (log):", norm="Media:", pareto="Ubicación:",	t="Grados de libertad:", unif="Mínimo:", weib="Forma:"
					)
				ini <- switch(input$dist,
					bern=0.5, bin=10, dunif=0, geom=0.5, hgeom=10, nbin=10, poi=10,	# discrete
					beta=2, cauchy=0, chisq=1, exp=1, F=1, gam=1, lap=0, logi=0, lognorm=0, norm=0, pareto=1,	t=15, unif=0, weib=1 # Continua
					)
        if(lab=="Probabilidad:"){
          sliderInput(dat()[[2]][1],lab,
                      min = 0, max = 1, value = ini, step = 0.025)
        } else {numericInput(dat()[[2]][1],lab,ini)}
				
			}
		})
	})
	
	output$dist2 <- renderUI({
		input$dist
		isolate({
			if(length(input$dist)){
				lab <- switch(input$dist,
					bin="Probabilidad:",	dunif="Discreto secuencia máxima:",	hgeom="N:",	nbin="Probabilidad:", # discrete
					beta="Beta:", cauchy="Escala:", F="Grados de libertad del denominador:", gam="Proporción", lap="Escala:", # Continua
					logi="Escala:", lognorm="Desviación estándar (log)", norm="Desviación estándar:", pareto="Forma:", unif="Máximo:", weib="Escala:"
					)
				ini <- switch(input$dist,
					bin=0.5, dunif=100, hgeom=20, nbin=0.5,
					beta=2, cauchy=1, F=15, gam=1, lap=1, logi=1, lognorm=1, norm=1, pareto=3, unif=1, weib=1
					)
				if(any(input$dist==c("bin","dunif","hgeom","nbin","cauchy","lap","logi","pareto","weib",
									"beta","F","gam","lognorm","norm","unif"))){
				  if(lab=="Probabilidad:"){
				    sliderInput(dat()[[2]][2],lab,
				                min = 0, max = 1, value = ini, step = 0.025)
				  } else {numericInput(dat()[[2]][2],lab,ini)}
				} #numericInput(dat()[[2]][2],lab,ini)
			}
		})
	})
	
	output$dist3 <- renderUI({
		input$dist
		isolate({
			if(length(input$dist)){
				lab <- switch(input$dist,
					dunif="Tamaño del paso:",	hgeom="K:")
				ini <- switch(input$dist,
					dunif=1, hgeom=5)
				if(any(input$dist==c("dunif","hgeom"))) numericInput(dat()[[2]][3],lab,ini)
			}
		})
	})
	
	output$sampDens <- renderUI({
		if(input$disttype=="Continua") checkboxInput("density","Curva de densidad de la muestra",FALSE)
	})
	
	output$BW <- renderUI({
		if(length(input$density)){
			if(input$density) numericInput("bw","bandwidth:",1)
		}
	})
	
	doPlot <- function(margins){
		if(length(input$dist)){
			d <- dat()[[1]]
			dist <- input$dist
			n <- input$n
			expr <- get(paste("expr",dist,sep="."))
			par(mar=margins)
			if(input$disttype=="Discreta"){
				barplot(as.numeric(table(d))/input$n,names.arg=names(table(d)),main=expr,xlab="Observaciones",ylab="Densidad",col="orange",cex.main=1.5,cex.axis=1.3,cex.lab=1.3)
			}
			if(input$disttype=="Continua"){
				hist(d,main=expr,xlab="Observaciones",ylab="Densidad",col="orange",cex.main=1.5,cex.axis=1.3,cex.lab=1.3,prob=T)
				if(length(input$density)) if(input$density & length(input$bw)) lines(density(d,adjust=input$bw),lwd=2)
			}
		}
	}
	
	output$plot <- renderPlot({
		doPlot(margins=c(4,4,10,1))
	},
	height=function(){ w <- session$clientData$output_plot_width; round((0.75*w)) }, width="auto"
	)
	
	output$dlCurPlot <- downloadHandler(
		filename = 'curPlot.pdf',
		content = function(file){
			pdf(file = file, width=11, height=8.5)
			doPlot(margins=c(6,6,10,2))
			dev.off()
		}
	)

	output$dldat <- downloadHandler(
		filename = function() { paste(input$dist, '.csv', sep='') },
		content = function(file) {
			write.csv(data.frame(x=dat()[[1]]), file)
		}
	)
	
	output$summary <- renderPrint({
		summary(dat()[[1]])
	})
	
	output$table <- renderTable({
		data.frame(x=dat()[[1]])
	})
	
	output$pageviews <-	renderText({
		if (!file.exists("pageviews.Rdata")) pageviews <- 0 else load(file="pageviews.Rdata")
		pageviews <- pageviews + 1
		save(pageviews,file="pageviews.Rdata")
		paste("Visits:",pageviews)
	})
	
	observe({
	  if (input$sal == 1) stopApp()
	})

})
