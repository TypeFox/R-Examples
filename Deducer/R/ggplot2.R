## TODO: Add comment
## 
## Author: Ian
################################################################################


#ggExtras is no longer compatible with ggplot2 0.9


#
#
#ggExtraGuiElements <-  function(){
#	
#	eval(parse(text="require('ggExtra')"))
#	cat("\nnote: the ggExtra GUI interface is experimental\n")
#	
#	
#	PlottingElement <- J("org.rosuda.deducer.plots.PlottingElement");
#	PlotController <- J("org.rosuda.deducer.plots.PlotController");
#	Theme <- J("org.rosuda.deducer.plots.Theme")
#	Layer <- J("org.rosuda.deducer.plots.Layer")
#	Geom <- J("org.rosuda.deducer.plots.Geom")
#	Stat <- J("org.rosuda.deducer.plots.Stat")
#	Scale <- J("org.rosuda.deducer.plots.Scale")
#	Position <- J("org.rosuda.deducer.plots.Position")
#	Aes <- J("org.rosuda.deducer.plots.Aes")
#	
#	Param <- J("org.rosuda.deducer.widgets.param.Param")
#	ParamNumeric <- J("org.rosuda.deducer.widgets.param.ParamNumeric")
#	ParamCharacter <- J("org.rosuda.deducer.widgets.param.ParamCharacter")
#	
#	Color <- J("java.awt.Color")
#	
#	#initialize plot elements
#	PlotController$init()
#	
#	#######################################	#######################################
#										#themes
#	#######################################	#######################################
#	
#	#a paramerter controlling base text size
#	pn <- new(ParamNumeric,"base_size")
#	pn$setViewType(Param$VIEW_ENTER)
#	pn$setLowerBound(0)
#	pn$setValue(14)
#	pn$setValue(14)
#	pn$setDefaultValue(14)	
#	
#	#add theme bb
#	pe <- PlottingElement$createElement("theme","bb");
#	theme <- new(Theme)
#	theme$setName("theme_bb")
#	theme$params$add(pn$clone())
#	pe$setModel(theme)
#	PlotController$addTheme(pe)
#	
#	#add theme bw2
#	pe <- PlottingElement$createElement("theme","bw2");
#	theme <- new(Theme)
#	theme$setName("theme_bw2")
#	theme$params$add(pn$clone())
#	pe$setModel(theme)
#	PlotController$addTheme(pe)
#	
#	#add theme dark
#	pe <- PlottingElement$createElement("theme","dark");
#	theme <- new(Theme)
#	theme$setName("theme_dark")
#	theme$params$add(pn$clone())
#	pe$setModel(theme)
#	PlotController$addTheme(pe)
#	
#	#add theme flashy
#	pe <- PlottingElement$createElement("theme","flashy");
#	theme <- new(Theme)
#	theme$setName("theme_flashy")
#	theme$params$add(pn$clone())
#	pe$setModel(theme)
#	PlotController$addTheme(pe)
#	
#	#add theme gray2
#	pe <- PlottingElement$createElement("theme","gray2");
#	theme <- new(Theme)
#	theme$setName("theme_gray2")
#	theme$params$add(pn$clone())
#	pe$setModel(theme)
#	PlotController$addTheme(pe)
#	
#	#add theme minimal
#	pe <- PlottingElement$createElement("theme","minimal");
#	theme <- new(Theme)
#	theme$setName("theme_minimal")
#	theme$params$add(pn$clone())
#	pe$setModel(theme)
#	PlotController$addTheme(pe)
#	
#	#add theme talk
#	pe <- PlottingElement$createElement("theme","talk");
#	theme <- new(Theme)
#	theme$setName("theme_talk")
#	theme$params$add(pn$clone())
#	pe$setModel(theme)
#	PlotController$addTheme(pe)
#	
#	
#	#######################################	#######################################
#										#geoms
#	#######################################	#######################################
#	
#	
#	
#	#######################################
#	#add geom_ellipse
#	#######################################
#	geom <- new(Geom)
#	geom$name <- "ellipse"
#	geom$defaultStat <- "identity"
#	geom$defaultPosition <- "identity"
#	
#	aes <- Aes$makeAes("x")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("y")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("size")
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("angle")
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("ar")
#	aes$dataType <- Aes$DATA_NUMERIC
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("colour")
#	aes$value <- .jcast(Color$black)
#	aes$defaultValue <- .jcast(Color$black)
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("fill")
#	aes$value <- .jcast(Color$black)
#	aes$defaultValue <- .jcast(Color$black)
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("alpha")
#	geom$aess$add(aes)
#	
#	geomLayer <- new(Layer)
#	geomLayer$isGeom = TRUE
#	geomLayer$name <- "geom_ellipse"
#	geomLayer$geom <- geom
#	geomLayer$stat <- Stat$makeIdentity()
#	geomLayer$pos <- Position$makePosition("identity")
#	
#	pe <- PlottingElement$createElement("geom","ellipse");
#	pe$setModel(geomLayer)
#	PlotController$addGeom(pe)
#	
#	#######################################
#	#add geom_field
#	#######################################
#	geom <- new(Geom)
#	geom$name <- "field"
#	geom$defaultStat <- "identity"
#	geom$defaultPosition <- "identity"
#	
#	aes <- Aes$makeAes("x")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("y")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("angle")
#	aes$required <- TRUE
#	aes$defaultUseVariable <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("length")
#	aes$required <- TRUE
#	aes$defaultUseVariable <- TRUE
#	aes$dataType <- Aes$DATA_NUMERIC
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("size")
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("colour")
#	aes$value <- .jcast(Color$black)
#	aes$defaultValue <- .jcast(Color$black)
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("linetype")
#	geom$aess$add(aes)
#	
#	geomLayer <- new(Layer)
#	geomLayer$isGeom = TRUE
#	geomLayer$name <- "geom_field"
#	geomLayer$geom <- geom
#	geomLayer$stat <- Stat$makeIdentity()
#	geomLayer$pos <- Position$makePosition("identity")
#	
#	pe <- PlottingElement$createElement("geom","field");
#	pe$setModel(geomLayer)
#	PlotController$addGeom(pe)
#
#	#######################################
#	#add geom_fielduv
#	#######################################
#	geom <- new(Geom)
#	geom$name <- "fielduv"
#	geom$defaultStat <- "identity"
#	geom$defaultPosition <- "identity"
#	
#	aes <- Aes$makeAes("x")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("y")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("abscissa")
#	aes$required <- TRUE
#	aes$defaultUseVariable <- TRUE
#	aes$dataType <- Aes$DATA_NUMERIC
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("ordinate")
#	aes$required <- TRUE
#	aes$defaultUseVariable <- TRUE
#	aes$dataType <- Aes$DATA_NUMERIC
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("colour")
#	aes$value <- .jcast(Color$black)
#	aes$defaultValue <- .jcast(Color$black)
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("size")
#	aes$value <- .jcast(new(J("java.lang.Double"),1))
#	aes$defaultValue <- .jcast(new(J("java.lang.Double"),1))
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("linetype")
#	geom$aess$add(aes)
#	
#	geomLayer <- new(Layer)
#	geomLayer$isGeom = TRUE
#	geomLayer$name <- "geom_fielduv"
#	geomLayer$geom <- geom
#	geomLayer$stat <- Stat$makeIdentity()
#	geomLayer$pos <- Position$makePosition("identity")
#	
#	pe <- PlottingElement$createElement("geom","fielduv");
#	pe$setModel(geomLayer)
#	PlotController$addGeom(pe)
#	
#	#######################################
#	#add geom_ngon
#	#######################################
#	geom <- new(Geom)
#	geom$name <- "ngon"
#	geom$defaultStat <- "identity"
#	geom$defaultPosition <- "identity"
#	
#	aes <- Aes$makeAes("x")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("y")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("sides")
#	aes$required <- FALSE
#	aes$dataType <- Aes$DATA_NUMERIC
#	aes$value <- .jcast(new(J("java.lang.Double"),5))
#	aes$defaultValue <- .jcast(new(J("java.lang.Double"),5))
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("size")
#	aes$value <- .jcast(new(J("java.lang.Double"),1))
#	aes$defaultValue <- .jcast(new(J("java.lang.Double"),1))
#	geom$aess$add(aes)
#		
#	
#	aes <- Aes$makeAes("angle")
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("ar")
#	aes$required <- FALSE
#	aes$dataType <- Aes$DATA_NUMERIC
#	aes$value <- .jcast(new(J("java.lang.Double"),1))
#	aes$defaultValue <- .jcast(new(J("java.lang.Double"),1))	
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("colour")
#	aes$value <- .jcast(Color$gray)
#	aes$defaultValue <- .jcast(Color$gray)
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("fill")
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("alpha")
#	geom$aess$add(aes)	
#	
#	aes <- Aes$makeAes("lex")
#	aes$required <- FALSE
#	aes$dataType <- Aes$DATA_NUMERIC
#	aes$value <- .jcast(new(J("java.lang.Double"),1))
#	aes$defaultValue <- .jcast(new(J("java.lang.Double"),1))	
#	geom$aess$add(aes)
#	
#	geomLayer <- new(Layer)
#	geomLayer$isGeom = TRUE
#	geomLayer$name <- "geom_ngon"
#	geomLayer$geom <- geom
#	geomLayer$stat <- Stat$makeIdentity()
#	geomLayer$pos <- Position$makePosition("identity")
#	
#	pe <- PlottingElement$createElement("geom","ngon");
#	pe$setModel(geomLayer)
#	PlotController$addGeom(pe)
#	
#	#######################################
#	#add geom_point2
#	#######################################
#	geom <- new(Geom)
#	geom$name <- "point2"
#	geom$defaultStat <- "identity"
#	geom$defaultPosition <- "identity"
#	
#	aes <- Aes$makeAes("x")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("y")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("shape")
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("size")
#	aes$value <- .jcast(new(J("java.lang.Double"),2))
#	aes$defaultValue <- .jcast(new(J("java.lang.Double"),2))
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("colour")
#	aes$value <- .jcast(Color$black)
#	aes$defaultValue <- .jcast(Color$black)
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("fill")
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("alpha")
#	geom$aess$add(aes)	
#	
#	aes <- Aes$makeAes("lex")
#	aes$required <- FALSE
#	aes$dataType <- Aes$DATA_NUMERIC
#	aes$value <- .jcast(new(J("java.lang.Double"),1))
#	aes$defaultValue <- .jcast(new(J("java.lang.Double"),1))	
#	geom$aess$add(aes)
#	
#	pc <- new(ParamCharacter,"size.unit")
#	pc$setOptions(c("native","cm","mm"))
#	pc$setViewType(Param$VIEW_EDITABLE_COMBO)
#	geom$params$add(pc)
#	
#	geomLayer <- new(Layer)
#	geomLayer$isGeom = TRUE
#	geomLayer$name <- "geom_point2"
#	geomLayer$geom <- geom
#	geomLayer$stat <- Stat$makeIdentity()
#	geomLayer$pos <- Position$makePosition("identity")
#	
#	pe <- PlottingElement$createElement("geom","point2");
#	pe$setModel(geomLayer)
#	PlotController$addGeom(pe)
#	
#	
#	#######################################
#	#add geom_star
#	#######################################
#	geom <- new(Geom)
#	geom$name <- "star"
#	geom$defaultStat <- "identity"
#	geom$defaultPosition <- "identity"
#	
#	aes <- Aes$makeAes("x")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("y")
#	aes$required <- TRUE
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("sides")
#	aes$required <- FALSE
#	aes$dataType <- Aes$DATA_NUMERIC
#	aes$value <- .jcast(new(J("java.lang.Double"),5))
#	aes$defaultValue <- .jcast(new(J("java.lang.Double"),5))
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("size")
#	aes$value <- .jcast(new(J("java.lang.Double"),1))
#	aes$defaultValue <- .jcast(new(J("java.lang.Double"),1))
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("angle")
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("colour")
#	aes$value <- .jcast(Color$black)
#	aes$defaultValue <- .jcast(Color$black)
#	geom$aess$add(aes)
#	
#	aes <- Aes$makeAes("alpha")
#	geom$aess$add(aes)	
#	
#	geomLayer <- new(Layer)
#	geomLayer$isGeom = TRUE
#	geomLayer$name <- "geom_star"
#	geomLayer$geom <- geom
#	geomLayer$stat <- Stat$makeIdentity()
#	geomLayer$pos <- Position$makePosition("identity")
#	
#	pe <- PlottingElement$createElement("geom","star");
#	pe$setModel(geomLayer)
#	PlotController$addGeom(pe)
#	
#	#######################################	#######################################
#										#scales
#	#######################################	#######################################
#	
#	#scale_length
#	scale <- Scale$makeContinuous("length")
#	scale$setName("scale_length")
#	pe <- PlottingElement$createElement("scale","length");
#	pe$setModel(scale)
#	PlotController$addScale(pe)
#	
#	#scale_angle
#	scale <- Scale$makeContinuous("angle")
#	scale$setName("scale_angle")
#	pe <- PlottingElement$createElement("scale","angle");
#	pe$setModel(scale)
#	PlotController$addScale(pe)
#	
#	#scale_ar
#	scale <- Scale$makeContinuous("ar")
#	scale$setName("scale_ar")
#	pe <- PlottingElement$createElement("scale","ar");
#	pe$setModel(scale)
#	PlotController$addScale(pe)
#	
#	#scale_ar
#	scale <- Scale$makeDiscrete("ar")
#	pe <- PlottingElement$createElement("scale","ar");
#	pe$setModel(scale)
#	PlotController$addScale(pe)
#	
#	#scale_lex
#	scale <- Scale$makeContinuous("lex")
#	scale$setName("scale_lex")
#	pe <- PlottingElement$createElement("scale","lex");
#	pe$setModel(scale)
#	PlotController$addScale(pe)
#	
#	#scale_lex
#	scale <- Scale$makeDiscrete("lex")
#	pe <- PlottingElement$createElement("scale","lex");
#	pe$setModel(scale)
#	PlotController$addScale(pe)
#	
#	#scale_sides
#	scale <- Scale$makeContinuous("sides")
#	pe <- PlottingElement$createElement("scale","sides");
#	pe$setModel(scale)
#	PlotController$addScale(pe)
#	
#	#scale_sides
#	scale <- Scale$makeDiscrete("sides")
#	scale$setName("scale_sides")
#	pe <- PlottingElement$createElement("scale","sides");
#	pe$setModel(scale)
#	PlotController$addScale(pe)
#	
#	
#}
#
#
#
#
#

