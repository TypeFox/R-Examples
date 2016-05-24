library(wppExplorer)
library(reshape2)
library(googleVis)
library(plyr)
library(ggplot2)

shinyServer(function(input, output, session) {
  indicatorData <- reactive({
    wppExplorer:::lookupByIndicator(input$indicator, input$indsexmult, input$indsex, input$selagesmult, input$selages)
  })
  
	indicator.fun <- reactive({
		wppExplorer:::ind.fun(as.integer(input$indicator))
	})
	
   indicatorDataLow <- reactive({
    wppExplorer:::getUncertainty(input$indicator, input$uncertainty, 'low', input$indsexmult, input$indsex, input$selagesmult, input$selages)
  })
  
  indicatorDataHigh <- reactive({
    wppExplorer:::getUncertainty(input$indicator, input$uncertainty, 'high', input$indsexmult, input$indsex, input$selagesmult, input$selages)
  })
  
  data <- reactive({
    wpp.by.year(indicatorData(), input$year)
  })
  
  rangeForAllYears <- reactive({
    range(indicatorData()$value, na.rm=TRUE)
  })
  
  pyramid.data <- reactive({
  	wppExplorer:::get.pyramid.data(input$year, input$seltcountries)
  })
  
  pyramid.data.low <- reactive({
  	wppExplorer:::get.pyramid.data(input$year, input$seltcountries, input$uncertainty, bound='low')
  })
  
   pyramid.data.high <- reactive({
  	wppExplorer:::get.pyramid.data(input$year, input$seltcountries, input$uncertainty, bound='high')
  })
  
  age.profile.mortM <- reactive({
  	wppExplorer:::get.pyramid.data(input$year, input$seltcountries, indicators=c(M='mxM'), load.pred=FALSE)
  })
  
  age.profile.mortF <- reactive({
  	wppExplorer:::get.pyramid.data(input$year, input$seltcountries, indicators=c(F='mxF'), load.pred=FALSE)
  })
  
  age.profile.popM <- reactive({
  	wppExplorer:::get.pyramid.data(input$year, input$seltcountries, indicators=c(M='popM'))
  })
  
  age.profile.popF <- reactive({
  	wppExplorer:::get.pyramid.data(input$year, input$seltcountries, indicators=c(F='popF'))
  })
  
  age.profile.fert <- reactive({
  	wppExplorer:::get.age.profile.fert(input$year, input$seltcountries)
  })
  
  age.profile.pfert <- reactive({
  	wppExplorer:::get.age.profile.pfert(input$year, input$seltcountries)
  })

  
  data.env <- function() wppExplorer:::wpp.data.env
    
  year.range <- reactiveValues(min=1955, max=2100)
    
  output$yearUI <- renderUI({
  	animationOptions(interval = 3000)
    sliderInput('year', h5('Year:'), sep="",
    			animate=TRUE,
                min=isolate(year.range$min), max=isolate(year.range$max), value = 2015, step=5)
  })
  
  observe({
  	# update the year slider if the indicator changes and its data have a different time range
  	if(is.null(input$year)) return(NULL)
  	data <- indicatorData()
  	if(nrow(data)==0) return(NULL)
  	yearRange <- range(data$Year)
  	if(year.range$min == yearRange[1]) return(NULL)
    isolate(year.range$min <- yearRange[1])
    value <- max(yearRange[1], isolate(input$year))	
  	updateSliderInput(session, "year", min=yearRange[1], value=value) 	
  })
  
  output$indicatorDesc <- renderText({
  	wppExplorer:::ind.definition(as.integer(input$indicator))
  })
  
  output$uncertaintyNote <- renderText({
  	if(wppExplorer:::ind.is.low.high(as.integer(input$indicator)) && wppExplorer:::get.wpp.year()>2010) return("")
    "No uncertainty available for this indicator."
  })
  
  year.output <- reactive(paste('Year:', input$year))
  output$mapyear <- renderText(year.output())
  #output$year1 <- renderText(year.output())
  output$year2 <- renderText(year.output())
  output$year3 <- renderText(year.output())
  
  get.country.charcodes <- function()
  	return(data.env()$iso3166[data.env()$iso3166$is.country, 'charcode'])
  	
  output$map <- reactive({
    if (is.null(input$year))
      return(NULL)
    df <- data()
    if (nrow(df) == 0)
      return(NULL)
    #df <- cbind(df, hover=rep('xxx', nrow(df)))
    country.codes <- get.country.charcodes()
    df <- df[df$charcode %in% country.codes,]
    options <- NULL
    if (input$normalizeMapAndCountryPlot) {
    	inddata <- indicatorData()
    	inddata <- inddata[inddata$charcode %in% country.codes, 'value']
    	options <- list( # fixed color scale
           colorAxis = list(
             minValue = min(inddata),
             maxValue = max(inddata)
           )
        )
	}
    list(data = df, options=options)
  })
  
  output$mapgvis <- renderGvis({
  	if (is.null(input$year))
      return(NULL)
    df <- data()
    if (nrow(df) == 0)
      return(NULL)
	#browser()
	col <- c('0x0000CC', '0x00CCFF', '0x33FF66', '0xFFFF66', '0xFF9900', '0xFF3300')
	gvisGeoChart(df, locationvar="charcode", #numvar="value", #hovervar="value", 
				colorvar="value", chartid="map",
				options=list(height=500, width=900, 
				dataMode='regions'#,
				#colors=paste('[', paste(col, collapse=', '), ']'
				))
  })
  
  output$countryPlot <- renderPlot({
    if (is.null(input$map_selection))
      return(NULL)

	data <- indicatorData()
    data.l <- indicatorDataLow()
    data.h <- indicatorDataHigh()
    df <- wpp.by.country(data, input$map_selection)
    low <- wpp.by.country(data.l, input$map_selection)
    high <- wpp.by.country(data.h, input$map_selection)
    idx.col.val <- grep('\\.', colnames(low))
    ylim <- range(c(df$value, low[,idx.col.val], high[,idx.col.val]), na.rm=TRUE)
    if (input$normalizeMapAndCountryPlot) {
    	idx.col.val.all <- grep('\\.', colnames(data.l))
    	ylim <- range(c(data$value, data.l[,idx.col.val.all], data.h[,idx.col.val.all]), na.rm=TRUE)
    }
    plot(df, type='n', ylim=ylim)
    title(main = paste(data.env()$iso3166$name[data.env()$iso3166$charcode==input$map_selection]))
    lines(df$Year, df$value, type='l')
    if(!is.null(low)) {
    	ipres <- which.max(df$Year[df$Year<min(low$Year)])
    	for(i in 3:1) {
			idx <- grep(paste0('\\.',i), colnames(low))
			if(length(idx)==0) next
    		lines(c(df$Year[ipres], low$Year), c(df[ipres, 'value'], low[,idx]), lty=i+1)
    		lines(c(df$Year[ipres], high$Year),c(df[ipres, 'value'], high[,idx]), lty=i+1)
    	}
    }
    abline(v=input$year, col=3, lty=3)
  })
	
 # output$table <- renderTable({
 	# iso <- data.env()$iso3166
 	# if(!input$includeAggr1) iso <- iso[iso$is.country,]
 	# data <- merge(data(), iso[,c('charcode', 'name')], by='charcode')#
	# data[,c('charcode', 'name', 'value')]
 # }, include.rownames = FALSE, include.colnames = FALSE)
  
  output$stable <- DT::renderDataTable({ #renderGvis({
  	year.data <- data()
  	if(nrow(year.data)==0) invalidateLater(1000, session)
  	iso <- data.env()$iso3166
 	if(!input$includeAggr2) iso <- iso[iso$is.country,]
  	year.data <- merge(year.data, iso[,c('charcode', 'name')], by='charcode')
  	low <- indicatorDataLow()
  	data <- cbind(year.data[,c('charcode', 'name', 'value')], rank=rank(year.data$value)) # add rank column
  	if(!is.null(low)) { # add intervals
  		data.l <- wpp.by.year(low, input$year)
  		 if(nrow(data.l) > 0) {		
    		data.h <- wpp.by.year(indicatorDataHigh(), input$year)
    		for(i in 1:3) {
    			colnames(data.l) <- sub(paste0('value.',i), paste('low', wppExplorer:::.get.pi.name.for.label(i)), colnames(data.l))
    			colnames(data.h) <- sub(paste0('value.',i), paste('high', wppExplorer:::.get.pi.name.for.label(i)), colnames(data.h))
    		}
    		ncoldata <- ncol(data)
  			data <- merge(data, data.l, by='charcode')
  			data <- merge(data, data.h, by='charcode')
  			# rearrange, so that columns corresponding to (low, high) pairs is always beside one another
  			if (ncol(data.l) > 2) {
  				l <- ncol(data.l) - 1
  				col.idx <- matrix(1:(2*l), nrow=l)
  				col.idx <- as.vector(t(col.idx))
  				data <- data[,c(1:ncoldata, col.idx+ncoldata)]
  			}
  		}
  	}
  	colnames(data)[1] <- 'code'
	data
  	})
  	
  output$ghist <- renderGvis({
  	data <- data()
  	if(is.null(data) || nrow(data)<=0) return(NULL)
  	# filter out aggregations
  	data <- merge(data[,c('charcode', 'value')], data.env()$iso3166[data.env()$iso3166$is.country, c('charcode', 'name'), drop=FALSE], by='charcode')
  	data <- data[,c('name', 'value')]
  	xlim <- if(input$fiXscaleHist) rangeForAllYears() else range(data$value, na.rm=TRUE)
  	#browser()
  	options <- list(title=paste(input$year), legend="{ position: 'none' }", colors="['green']", 
  						height="500px", width="900px", histogram=paste0("{bucketSize: ", diff(xlim)/30, "}"))
  						
  	options$hAxis <- paste0("{maxAlternation: 1,  minValue:", xlim[1], ", maxValue:", xlim[2], "}")
  	gvisHistogram(data, options=options)
  })
  
  output$hist <- renderPlot({
  	data <- data()
  	if(is.null(data) || nrow(data)<=0) return(NULL)
  	# filter out aggregations
  	data <- merge(data[,c('charcode', 'value')], data.env()$iso3166[data.env()$iso3166$is.country, 'charcode', drop=FALSE], by='charcode')$value
  	xlim <- if(input$fiXscaleHist) rangeForAllYears() else range(data, na.rm=TRUE)
  	binw <- diff(xlim)/20
    qplot(data()$value, binwidth=binw, xlim=c(xlim[1]-binw/2, xlim[2]+binw/2), xlab='value')
  })
        
  output$ageselection <- renderUI({
  	if(indicator.fun() %in% c('fertage', 'pfertage')){
  		ages <- paste(seq(15, by=5, length=7), seq(19, by=5, length=7), sep='-')
  	} else {
  		if(indicator.fun()=='mortagesex')
  			ages <- c(0,1,seq(5, by=5, length=19))
  		else ages <- paste(seq(0, by=5, length=20), seq(4, by=5, length=20), sep='-')
  		ages <- c(ages, '100+')
  	}
  	if (wppExplorer:::ind.no.age.sum(as.integer(input$indicator))) { # no multiple choices allowed
  		multiple <- FALSE
		name <- 'selages'
		selected<-NULL
  	} else { 		
  		multiple <- TRUE
  		name <- 'selagesmult'
  		selected <- ages[1]
  	}
  	selectInput(name, 'Age:', ages, multiple=multiple, selected=selected, selectize = FALSE)
	})
	
	output$sexselection <- renderUI({
		choices<-if(indicator.fun() %in% c('fertage', 'pfertage')) c(Female="F") else c(Female="F", Male="M")
		if(wppExplorer:::ind.no.age.sum(as.integer(input$indicator))){
  			multiple <- FALSE
  			selected <- NULL
  			name <- 'indsex'
  		} else {
  			multiple <- TRUE
  			selected <- choices
  			name <- 'indsexmult'
  		}
  		selectInput(name, 'Sex:', choices=choices, selected=selected, selectize = FALSE, multiple=multiple)
	})
	
  output$cselection <- renderUI({
  	o <- order(data.env()$iso3166[,'name'])
  	codes <- as.character(data.env()$iso3166[o,'charcode'])
  	names <- paste(codes, data.env()$iso3166[o,'name'])
  	countries <- structure(
  		codes,
  		names = names
  	)
  	do.call('selectInput', list('seltcountries', 'Select countries/areas:', countries, multiple=TRUE, selectize = FALSE,
  					selected=countries[1] #names[1]
  					))
	})
	
  cast.profile.data <- function(data) {
    vrange <- range(data$value, na.rm=TRUE)
    hrange <- if(is.element('15-19', data$age)) c(0, length(unique(data$age))) else range(data$age)
    #browser()
    data <- dcast(data, age.num + age ~ charcode, mean)
    data$age.num <- NULL
    list(casted=data, hrange=hrange, vrange=vrange)
  	
  }	
	
  filter.trend.data <- function(data, countries, cast=TRUE){
  	data <- wpp.by.countries(data, countries)
    if(is.null(data) || nrow(data) <= 0) return(NULL)
   	hrange <- range(data$Year, na.rm=TRUE)
    vrange <- range(data[,grep('value', colnames(data))], na.rm=TRUE)
    if(cast)
    	data <- dcast(data, Year ~ charcode, mean)
    list(casted=data, hrange=hrange, vrange=vrange)
  }
  
  get.trends <- reactive({
  	if(is.null(input$seltcountries)) return(NULL)
  	filter.trend.data(indicatorData(), input$seltcountries)
  })
  

  get.age.profiles <- reactive({
  	if(is.null(input$seltcountries)) return(NULL)
  	filter.age.profiles(indicatorData(), input$seltcountries, input$year)
  })
  
  get.trends.nocast <- reactive({
  	if(is.null(input$seltcountries)) return(NULL)
  	filter.trend.data(indicatorData(), input$seltcountries, cast=FALSE)
  })
  
  get.trends.low <- reactive({
  	if(is.null(input$seltcountries)) return(NULL)
  	filter.trend.data(indicatorDataLow(), input$seltcountries, cast=FALSE)
  })
  
  get.trends.high <- reactive({
  	if(is.null(input$seltcountries)) return(NULL)
  	filter.trend.data(indicatorDataHigh(), input$seltcountries, cast=FALSE)
  })
  
  output$trends <- reactive({
	data <- get.trends()
	if(is.null(data)) return(data)
    list(data = wppExplorer:::preserveStructure(data$casted),
         options = list(
           hAxis = list(viewWindowMode = 'explicit', viewWindow = list(
             min = data$hrange[1], max = data$hrange[2]
           ), format="####"),
           vAxis = list(viewWindowMode = 'explicit', viewWindow = list(
             min = data$vrange[1], max = data$vrange[2]
           ), logScale=input$median.logscale)
         )
    )
  })
  
  show.age.profile <- function(sex, fun, logscale=FALSE, year) {
  	if(fun=='mortagesex')
		data <- do.call(paste0('age.profile.mort', sex), list())
	else {
		if(fun %in% c('tpop', 'tpopF', 'tpopM', 'popagesex'))
			data <- do.call(paste0('age.profile.pop', sex), list())
		else {
			if(fun %in% c('fert', 'fertage') && sex=='F') {
				data <- age.profile.fert()
			} else {
				if(fun == 'pfertage' && sex=='F') {
					data <- age.profile.pfert()
				} else {
  					return(list(data=data.frame(age=c(0,0), v=c(0,0)),
  						#data.frame(age=seq(0,100, by=5), value=rep(0, 21)), 
  						options=list(title=paste(c(F='Female', M='Male')[sex], ': No age profiles for this indicator.'),
  							legend= list(position="none"),
  							hAxis = list(viewWindow = list(min=-1, max=1)),
  							vAxis = list(viewWindow = list(min=-1, max=1))
  						)))
  				}
			}
		}
	}
	if(is.null(data)) return(NULL)
	data <- cast.profile.data(data)
    list(data = wppExplorer:::preserveStructure(data$casted),
         options = list(
           hAxis = list(slantedText=fun!='mortagesex',
           				viewWindowMode = 'explicit', viewWindow = list(
           				min = data$hrange[1], max = data$hrange[2])),
           #), format="####"),
           vAxis = list(viewWindowMode = 'explicit', viewWindow = list(
             min = data$vrange[1], max = data$vrange[2]
           ), logScale=logscale),
           legend = list(position="right"),
           title = paste(year, c(F='Female', M='Male')[sex])
         )
    )
  }
  
  output$age.profileM <- reactive({
  	show.age.profile('M', indicator.fun(), input$aprofile.logscale, input$year)
  })
  output$age.profileF <- reactive({
	show.age.profile('F', indicator.fun(), input$aprofile.logscale, input$year)
  })
  reset.trend.data <- function() ggplot.data$trends <- NULL
  reset.pyramid.data <- function() ggplot.data$pyramid <- NULL
  
  ggplot.data <- reactiveValues(trends=NULL, pyramid=NULL)
  trends.ranges <- reactiveValues(year=NULL, value=NULL)
  
  output$probtrends <- renderPlot({
  	reset.trend.data()
  	data <- get.trends.nocast()
  	if(is.null(data)) return(data)
  	data <- data$casted
  	low <- NULL
  	if (!input$trend.logscale) low <- get.trends.low()
  	if(!is.null(low)) {
  		high <- get.trends.high()
  		colnames(low$casted) <- sub('value', 'low', colnames(low$casted))
  		colnames(high$casted) <- sub('value', 'high', colnames(high$casted))
  		low.high <- merge(low$casted, high$casted, by=c('charcode', 'Year'))
  		#browser()
  		#colnames(low.high)[3:4] <- c('low', 'high')
  		min.year <- min(low.high$Year, na.rm=TRUE)
  		data <- merge(data, low.high, by=c('charcode', 'Year'), all=TRUE)
  		idx <- which(data$Year == min.year-5)
  		for (col in grep('low|high', colnames(data)))
  			data[idx,col] <- data$value[idx]
  	}
  	if(!is.null(trends.ranges$year)) {
  		data.zoom <- subset(data, Year >= trends.ranges$year[1] & Year <= trends.ranges$year[2] & 
  								value >= trends.ranges$value[1] & value <= trends.ranges$value[2])
  		if(nrow(data.zoom) == 0) {
  			trends.ranges$year <- trends.ranges$value <- NULL
  		} else data <- data.zoom
  	}
  	g <- ggplot(data, aes(x=Year,y=value,colour=charcode, fill=charcode)) + geom_line() + theme(legend.title=element_blank())
  	if (input$trend.logscale) g <- g + coord_trans(y="log2")
  	isolate(ggplot.data$trends <- data)
  	
  	if(!is.null(low)) {
  		line.data <- NULL
  		for(i in 3:1) {
  			idx <- grep(paste0('\\.',i), colnames(data))
  			if(length(idx)==0) next
  			g <- g + geom_ribbon(aes_string(ymin=colnames(data)[idx][1], ymax=colnames(data)[idx][2], 
  											fill="charcode", colour="charcode", linetype=NA), alpha=c(0.3, 0.2, 0.1)[i])
  			if(!is.null(line.data)) colnames(line.data) <- c('charcode', 'Year', 'low', 'high', 'variant')
  			line.data <- rbind(line.data, setNames(cbind(data[,c(1,2,idx)], wppExplorer:::.get.pi.name.for.label(i)), colnames(line.data)))
  			tmp1 <- tmp2 <- data
  			tmp1[,'value'] <- data[,colnames(data)[idx][1]]
  			tmp2[,'value'] <- data[,colnames(data)[idx][2]]
  			isolate(ggplot.data$trends <- rbind(ggplot.data$trends, tmp1, tmp2))
  		}
  		colnames(line.data) <- c('charcode', 'Year', 'low', 'high', 'variant')
		
  		g <- g + geom_line(data=line.data, aes(y=low, linetype=variant, colour=charcode))
  		g <- g + geom_line(data=line.data, aes(y=high, linetype=variant, colour=charcode))
  		g <- g + scale_linetype_manual(values=c("80%"=2, '1/2child'=4, "95%"=3), na.value=0)
  	}
  	isolate(ggplot.data$trends <- merge(ggplot.data$trends, data.env()$iso3166[,c('charcode', 'name')], by='charcode', sort=FALSE))
  	g
  })
  
  observeEvent(input$probtrends_zoom, {
  	trends.ranges$year <- c(round(input$probtrends_zoom$xmin,0), round(input$probtrends_zoom$xmax,0))
  	trends.ranges$value <- c(input$probtrends_zoom$ymin, input$probtrends_zoom$ymax)
  })

  observeEvent(input$probtrends_zoom_reset, {
  	trends.ranges$year <- NULL
  	trends.ranges$value <- NULL
  })
  
  output$probtrends_selected <- renderText({
  	if(is.null(input$probtrends_values)) return(" ")
  	selected <- nearPoints(ggplot.data$trends, input$probtrends_values, maxpoints = 1, threshold=10)
	if(nrow(selected) == 0) return(" ")
	paste0("Year = ", selected$Year, ', Value = ', round(selected$value,3), ", Country: ", selected$name)
	})
		
  .is.pyramid.indicator <- function() {
  	 if(!indicator.fun() %in% c('tpop', 'tpopF', 'tpopM', 'popagesex')) {
  		df <- data.frame(x=0, y=0, lab='No pyramid data for this indicator.')
  		g <- ggplot(df, aes(x=x, y=y, label=lab)) + geom_text() + scale_y_continuous(name='') + scale_x_continuous(name='')
  		print(g)
  		return(FALSE)
  	}
  	TRUE
  }
  .get.prop.data <- function(data, tpop) {
	tpop <- wppExplorer::wpp.by.countries(wppExplorer::wpp.by.year(tpop, input$year), input$seltcountries)
  	colnames(tpop)[2] <- 'tpop' 
  	data <- merge(data, tpop, by='charcode')  		
  	data <- ddply(data, 'charcode', mutate, value = value/tpop)
  	data$tpop <- NULL
  	data
  }
  
  .get.pyramid.data <- reactive({
  #.get.pyramid.data <- function(proportion=FALSE) {
  	proportion <- input$proppyramids
  	data <- pyramid.data()
  	if(proportion) {
  		tpop <- wppExplorer::wpp.indicator('tpop')
  		data <- .get.prop.data(data, tpop)
  	}
	low <- pyramid.data.low()
  	if(!is.null(low) && nrow(low)>0) {
  		high <- pyramid.data.high()
  		if(proportion) {
			#browser()
			which.pi <- wppExplorer:::.get.pi.name(as.integer(input$uncertainty))
			which.pi <- if('half.child' %in% which.pi) 'half.child' else NULL # currently only half child for pyramid available
			if(!is.null(which.pi)) {
  				tpop <- wppExplorer::wpp.indicator('tpop.ci', which.pi=which.pi, bound='low')
  				low <- .get.prop.data(low, tpop)
  				lowval <- low$value
  				tpop <- wppExplorer::wpp.indicator('tpop.ci', which.pi=which.pi, bound='high')
  				high <- .get.prop.data(high, tpop)
  				low$value <- pmin(low$value, high$value, na.rm=TRUE)
  				high$value <- pmax(high$value, lowval, na.rm=TRUE)
  			} else low <- NULL
  		}
  		if(!is.null(low)) {
  			low.high <- merge(low, high, by=c('charcode', 'age', 'age.num', 'sex'), sort=FALSE)
  			colnames(low.high)[5:6] <- c('low', 'high')
  			data <- merge(data, low.high, all=TRUE, sort=FALSE)
  		}
	}
  	data <- data[order(data$age.num),]
  	data
  })
  
  .print.pyramid <- function(data) {
  	data.range <- range(abs(data$value), na.rm=TRUE)
  	g <- ggplot(data, aes(y=value, x=reorder(age, age.num), group=charcode, colour=charcode)) + 
  			geom_line(data=subset(data, sex=='F')) + geom_line(data=subset(data, sex=='M')) + 
  			scale_x_discrete(name="") + scale_y_continuous(labels=function(x)abs(x)) + 
  			coord_flip() + ggtitle(input$year) + theme(legend.title=element_blank())
  	g <- g + geom_text(data=NULL, y=-data.range[2]/2, x=20, label="Male", colour='black')
  	g <- g + geom_text(data=NULL, y=data.range[2]/2, x=20, label="Female", colour='black')
  	g <- g + geom_hline(yintercept = 0)
  	if(is.element('low', colnames(data))) {
  		g <- g + geom_ribbon(data=subset(data, sex=='F'), aes(ymin=low, ymax=high, linetype=NA), alpha=0.3)
  		g <- g + geom_ribbon(data=subset(data, sex=='M'), aes(ymin=high, ymax=low, linetype=NA), alpha=0.3)
  		line.data <- cbind(data, variant=wppExplorer:::.get.pi.name.for.label(3)) # only half-child variant available 
  		g <- g + geom_line(data=subset(line.data, sex=='F'), aes(y=low, linetype=variant, colour=charcode, group=charcode)) # female low
  		g <- g + geom_line(data=subset(line.data, sex=='F'), aes(y=high, linetype=variant, colour=charcode, group=charcode)) # female high
  		g <- g + geom_line(data=subset(line.data, sex=='M'), aes(y=low, linetype=variant, colour=charcode, group=charcode)) # male low
  		g <- g + geom_line(data=subset(line.data, sex=='M'), aes(y=high, linetype=variant, colour=charcode, group=charcode)) # male high
        g <- g + scale_linetype_manual(values=c("80%"=2, '1/2child'=4, "95%"=3), na.value=0)
        tmp1 <- tmp2 <- data
		tmp1[,'value'] <- data[,'low']
  		tmp2[,'value'] <- data[,'high']
  		isolate(ggplot.data$pyramid <- rbind(ggplot.data$pyramid, tmp1, tmp2))
  	}
  	g
  }
  
  pyramid.ranges <- reactiveValues(age=NULL, value=NULL)
  output$pyramids <- renderPlot({
  	if(!.is.pyramid.indicator()) return()
  	reset.pyramid.data()
	#data <- .get.pyramid.data(proportion=input$proppyramids)
	data <- .get.pyramid.data()
	male.idx <- which(data$sex=='M')
	for(col in c('value', 'low', 'high'))
		if(!is.null(data[[col]])) data[[col]][male.idx] <- -data[[col]][male.idx]
	if(!is.null(pyramid.ranges$age)) {
		#print(input$pyramid_zoom)
  		data.zoom <- subset(data, age.num >= pyramid.ranges$age[1] & age.num <= pyramid.ranges$age[2] & 
  								value >= pyramid.ranges$value[1] & value <= pyramid.ranges$value[2])
  		if(nrow(data.zoom) == 0) {
  			pyramid.ranges$age <- pyramid.ranges$value <- NULL
  		} else data <- data.zoom
  	}
	isolate(ggplot.data$pyramid <- data)
	g <- .print.pyramid(data)
	isolate(ggplot.data$pyramid <- merge(ggplot.data$pyramid, data.env()$iso3166[,c('charcode', 'name')], by='charcode', sort=FALSE))
	g
  })
  
  observeEvent(input$pyramid_zoom, {
  	pyramid.ranges$age <- c(round(input$pyramid_zoom$ymin,0), round(input$pyramid_zoom$ymax,0))
  	pyramid.ranges$value <- c(input$pyramid_zoom$xmin, input$pyramid_zoom$xmax)
  })

  observeEvent(input$pyramid_zoom_reset, {
  	pyramid.ranges$age <- NULL
  	pyramid.ranges$value <- NULL
  })
  
  output$pyramid_selected <- renderText({
  	if(is.null(input$pyramid_values)) return(" ")
  	selected <- nearPoints(ggplot.data$pyramid, input$pyramid_values, xvar='value', yvar='age.num', maxpoints = 1, threshold=10)
	if(nrow(selected) == 0) return(" ")
	paste0(if(selected$value < 0) "Male " else "Female ", selected$age, ', Value = ', round(abs(selected$value),3), ", Country: ", selected$name)
	})
  
  # .get.digits <- reactive({
  	# print(wppExplorer:::ind.digits(as.integer(input$indicator)))
  	# wppExplorer:::ind.digits(as.integer(input$indicator))
  # })
  
  # format_num <- function(col, digits) {
  	# format <- paste0("%.", digits, 'f')
	# if (is.numeric(col)) sprintf(format, col)
	# else col
# }

  output$trendstable <- renderTable({
	data <- get.trends()
	if(is.null(data)) return(data)
	df <- as.data.frame(data$casted[,-1])
	if(ncol(df) > 1) {
		if(wppExplorer:::ind.sum.in.table(as.integer(input$indicator))) {
			df <- cbind(df, rowSums(df))
			colnames(df)[ncol(df)] <- 'Sum'
		}
	} else colnames(df) <- input$seltcountries # one country selected
	
	df <- t(df)
	#df <- t(as.data.frame(lapply(df, format_num, digits=wppExplorer:::ind.digits(as.integer(input$indicator)))))
	# df <- t(data$casted[,-1]) # remove year column
	#browser()
	colnames(df) <- as.integer(data$casted[,'Year'])
	df
	}, include.rownames = TRUE)

  output$trendstabletitle <- renderText({
 	wppExplorer:::get.indicator.title(input$indicator, input$indsexmult, input$indsex, input$selagesmult, input$selages) 	
   })
   
   all.data <- reactive({
   		if(is.null(wppExplorer:::wpp.data.env$mchart.data)) {
   			inds <- unique(c(input$indicator, 1,2,0,4))[1:4]
   		} else inds <- input$indicator
		wppExplorer:::lookupByIndicator.mchart(inds, input$indsexmult, input$indsex, input$selagesmult, input$selages)
	})
	
  output$graphgvis <- renderGvis({
  	# Take a dependency on input$AddIndicator button
  	input$AddIndicator

    df <- isolate(all.data())
	#browser()
	gvisMotionChart(df,
                      idvar="name", 
                      timevar="Year",
                      #xvar="xvalue", yvar="yvalue",
                      colorvar="UN Areas", 
                      #sizevar="zvalue",
                      options=list(width=700, height=600))
  })
  output$AddIndicatorText <- renderText({"\nAdd indicator from the left panel\nto chart axes:"})
})