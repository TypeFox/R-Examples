# =========================================================================
# R package: rplotengine.R
# =========================================================================

rplotengine = function (args_file = "mygraph.arg") {

# -----------------------------------------------------------------------
# Importamos paquetes externos
# Buscamos en directorios especificados en la variable de entorno R_LIBS_USER
# -----------------------------------------------------------------------
#library(xtable)
#require(xtable)
#do.call(library, list("xtable")) 

# -------------------------------------------------------------------------
# Leemos archivos de funciones (file.path para ser multiplataforma)
# -------------------------------------------------------------------------
#source(...)

# ------------------------------------------------------------------------
# Definimos los parametros del script (argumentos para un grafico)
# ------------------------------------------------------------------------

# Variables de R con los valores de los parametros
R_VARS = c( "title", "subtitle", "x_axis_title", "y_axis_title",
    "col_x_values", "col_y_values", "series_names",
    "x_factor", "y_factor", "total_serie", "total_value",
    "x_min", "x_max", "y_min", "y_max", "y_log",
    "show_titles", "show_grid",
    "show_confint", "confint_as_percentage",
    "text_size_title", "text_size_subtitle",
    "text_size_axis_ticks", "text_size_axis_titles", "text_size_legend",
    "pos_legend", "graph_type", "width_factor", "height_factor", "latex_digits", "verbose",
    "data_filename", "graph_filename", "graph_fileext_seq" )

# Para evitar los siguientes mensajes al chequear el paquete:
# rplotengine: no visible binding for global variable 'xxx'
title = ""
subtitle = ""
title = ""
subtitle = ""
x_axis_title = ""
y_axis_title = ""
col_x_values = ""
col_y_values = ""
series_names = ""
x_factor = ""
y_factor = ""
total_serie = ""
total_value = ""
x_min = ""
x_max = ""
y_min = ""
y_max = ""
y_log = ""
show_titles = ""
show_grid = ""
#show_hotspots = ""
show_confint = ""
confint_as_percentage = ""
text_size_title = ""
text_size_subtitle = ""
text_size_axis_ticks = ""
text_size_axis_titles = ""
text_size_legend = ""
pos_legend = ""
graph_type = ""
width_factor = ""
height_factor = ""
latex_digits = ""
verbose = ""
data_filename = ""
graph_filename = ""
graph_fileext_seq = ""

# ------------------------------------------------------------------------
# Leemos parametros del grafico desde archivo de entrada
# ------------------------------------------------------------------------

NUM_ARGS = length(R_VARS)

# Obtenemos lista de argumentos desde linea de ordenes
# First read in the arguments listed at the command line
args  = commandArgs()
nargs = length(args)

# Nos aseguramos de que existe el archivo de datos de entrada
if (!file.exists(args_file)) {
   # Puede que sea cuando se chequean los ejemplos del paquete.
   args_file = system.file(args_file, package = "rplotengine")
   if (length(args_file) == 0) {
      # Cuando se chequea el paquete el directorio actual es .../rplotengine.Rcheck (print(getwd()),
	  # y existe en .../rplotengine.Rcheck/rplotengine (copiado desde /rplotengine/inst).
      if (verbose == "1") {
         writeLines (paste("Graph parameters file '", args_file, "' not found.", sep=""))  # stop
      }
      return (FALSE)
   }
} else {
   args_file = paste(getwd(), "/", args_file, sep="")
}
writeLines (paste("Graph parameters file ...", sep=""))
writeLines (paste("   File: ", args_file, sep=""))

args_lines = readLines (args_file)
nlines = length (args_lines)

# Asignamos parametros
var = 1
for (l in 1:nlines) {
   args_lines[l] = trim(args_lines[l])
   if (args_lines[l] == "") next
   if (substr (args_lines[l],1,1) == '#') next
   arg = unlist(strsplit(args_lines[l],"\\="))
   if (length(arg) != 2) {
      writeLines ( paste("Syntax error in input file, line: ", l, sep="") )
      writeLines ( args_lines[l] )
      return (FALSE)
   }
   arg_name  = arg[1]
   arg_value = arg[2]
   if (arg_name == R_VARS[var]) {
      # Si el archivo argumentos se ha escrito en el orden esperado y sin comentarios ira mas rapido
      i = var
      if (var < NUM_ARGS) var = var+1
   } else {
      # Buscamos
      i = NUM_ARGS+1
      for (a in 1:NUM_ARGS) {
         if (arg_name == R_VARS[a]) {
            i=a
            break
         }
      }
   }
   if (i > NUM_ARGS) {
      writeLines ( paste("Unknown argument in input file, line: ", l, sep="") )
      writeLines ( args_lines[l] )
      return (FALSE)
   }
   val_name = R_VARS[i]
   line = paste (val_name, "='", arg_value, "'", sep="")  
   eval (parse(text=line))
}

# -------------------------------------------------------------------------

# Convertimos a valores numericos o listas algunos argumentos
col_x_values  = as.numeric (col_x_values)
col_y_values  = as.numeric (unlist (strsplit (col_y_values,",")))
series_names  = unlist (strsplit (series_names,","))
x_factor      = as.double (x_factor)
y_factor      = as.double (y_factor)
total_serie   = as.numeric (total_serie)
total_value   = as.double (total_value)

# Valores min/max se calculan una vez leidos y procesados los datos
#x_min
#x_max
#y_min
#y_max

#y_log         = ifelse((y_log=="0"), FALSE, TRUE)
#show_titles   = ifelse((show_titles=="0"), FALSE, TRUE)
#show_grid     = ifelse((show_grid=="0"), FALSE, TRUE)
##show_hotspots = ifelse((show_hotspots=="0"), FALSE, TRUE)
#show_confint  = ifelse((show_confint=="0"), FALSE, TRUE)
#confint_as_percentage = ifelse((confint_as_percentage=="0"), FALSE, TRUE)

text_size_title       = as.double (text_size_title)
text_size_subtitle    = as.double (text_size_subtitle)
text_size_axis_ticks  = as.double (text_size_axis_ticks)
text_size_axis_titles = as.double (text_size_axis_titles)
text_size_legend      = as.double (text_size_legend)
pos_legend    = as.numeric (pos_legend)
graph_type    = as.numeric (graph_type)
width_factor  = as.numeric (width_factor)
height_factor = as.numeric (height_factor)
latex_digits  = as.numeric (latex_digits)
graph_fileext_seq = unlist (strsplit(graph_fileext_seq,","))


# ------------------------------------------------------------------------
# Leemos los datos para el grafico
# ------------------------------------------------------------------------

# Nos aseguramos de que existe el archivo de datos de entrada
if (!file.exists(data_filename)) {
   # Puede que sea cuando se chequean los ejemplos del paquete.
   data_filename = system.file(data_filename, package = "rplotengine")
   if (length(data_filename) == 0) { 
      # Cuando se chequea el paquete el directorio actual es .../rplotengine.Rcheck (print(getwd()),
	  # y existe en .../rplotengine.Rcheck/rplotengine.
      if (verbose == "1") {
         writeLines (paste("Data file '", data_filename, "' not found.", sep="")) # stop
      }
      return (FALSE)
   }
} else {
   data_filename = paste(getwd(), "/", data_filename, sep="")
}

# Leemos los datos (simbolo # se utiliza para comentarios con read.table)
if (verbose == "1") {
   writeLines (paste("Loading data file ...", sep=""))
   writeLines (paste("   File: ", data_filename, sep=""))
}
data = read.table (data_filename, header=TRUE)

# Debug: mostramos los datos
if (verbose == "1") {
   writeLines (paste("   cols: ", NCOL(data), ifelse(NCOL(data)==0, " --> ERROR", ""), sep=""))
   writeLines (paste("   rows: ", NROW(data), ifelse(NROW(data)==0, " --> ERROR", ""), sep=""))
   #print (data)
}

if ((NCOL(data)==0) | (NROW(data)==0)) { 
   return (FALSE)
}

# ------------------------------------------------------------------------
# Preparamos los datos para el grafico
# ------------------------------------------------------------------------

if (verbose == "1") {
   writeLines ("Processing data ...")
}

# Contamos valores del eje X e Y
col_x_values_len = NROW(data)
col_y_values_len = NROW(col_y_values)

# Contamos series a dibujar
series_names_len = NROW(series_names)

# Creamos arrays de datos definitivos a dibujar y exportar a LaTeX
num_series = col_y_values_len
num_series_legend = col_y_values_len  # Medias no salen en leyenda

# Serie total?
if (total_serie > 0) {    # 1-TOTAL o 2-AVERAGE o 3-CONSTANT o 4-X_PROPORTION
   num_series = num_series + 1
   num_series_legend = num_series_legend + 1
}
x = array (0, c(col_x_values_len,num_series))  # Lo de 'num_series' es necesario para dibujar el CI -> plot(x[,1]... y lines(x[,1]...
y = array (0, c(col_x_values_len,num_series))

# Si no hemos especificado los nombres de las series suficientes para la leyenda 
# consideramos las cabeceras del fichero. Asignamos tambien un nombre a la serie total.
# Esto ya se hace al leer los graficos desde simul_graphs, pero se mantiene por si
# se utiliza este script R desde un archivo por lotes o script.
if (series_names_len < num_series_legend) {
   series_names_aux = array (0, c(num_series_legend))
   for (serie in 1:col_y_values_len) {
      if (serie <= series_names_len) {
         series_names_aux[serie] = series_names[serie]
      } else {
         col = col_y_values[serie]
         series_names_aux[serie] = names(data)[col]
      }
   }
   # Serie total
   if (total_serie > 0) {
      if (num_series_legend <= series_names_len) {
         series_names_aux[num_series_legend] = series_names[num_series_legend]
      } else if (total_serie == 1) {  # 1-TOTAL
         series_names_aux[num_series_legend] = "Total"
      } else if (total_serie == 2) {  # 2-AVERAGE
         series_names_aux[num_series_legend] = "Average"     
      } else if (total_serie == 3) {  # 3-CONSTANT
         series_names_aux[num_series_legend] = "Constant value"
      } else if (total_serie == 4) {  # 4-X_PROPORTION
         series_names_aux[num_series_legend] = "Proportional to x"
      }
   }
   series_names = series_names_aux
} else if (series_names_len > num_series_legend) {
   # Si hemos especificado mas nombres que series consideramos solo los necesarios
   series_names = series_names[1:num_series_legend]
}

# Debug: series names
if (verbose == "1") {
   for (serie in 1:num_series) {
      writeLines (paste("   Serie #", serie, ": ", series_names[serie], " (col=", col_y_values[serie], ")", sep=""))
   }
}

# Para cada valor del eje x (fila del archivo de datos)
for (x_idx in 1:col_x_values_len) {
   # Obtenemos el valor del eje x, que sera el mismo para todas las series (aplicando el factor)
   x_value = data[x_idx,col_x_values] * x_factor 	# row.names(data)[x_idx] * x_factor
   # Obtenemos el valor del eje y para Para todas las series (aplicando el factor)
   for (serie in 1:col_y_values_len) {              
       col = col_y_values[serie]
       # Por seguridad, comprobamos que la columna especificada no se salga de rango
       if (col>NCOL(data)) {
          writeLines (paste("   Column ", col, " for serie #", serie, " is out of range [0..", NCOL(data), "]", sep=""));
          return (FALSE)
       }
       x[x_idx,serie] = x_value		# Necesario solo para calculo del CI
       y[x_idx,serie] = data[x_idx,col] * y_factor
       # Para series total (1-TOTAL o 2-AVERAGE)
       if ((total_serie == 1) || (total_serie == 2)) { 
           y[x_idx,num_series] = y[x_idx,num_series] + y[x_idx,serie]   #data[x_idx,col]
       }
   }
   # Serie total
   if (total_serie > 0) {
      x[x_idx,num_series] = x_value
      if (total_serie == 2) {         # 2-AVERAGE
         y[x_idx,num_series] = y[x_idx,num_series] / col_y_values_len  #series_names_len??
      } else if (total_serie == 3) {  # 3-CONSTANT
         y[x_idx,num_series] = total_value * y_factor
      } else if (total_serie == 4) {  # 4-X_PROPORTION
         y[x_idx,num_series] = total_value * data[x_idx,col_x_values] * y_factor    # total_value * x_value
      }
   }
}

# ------------------------------------------------------------------------
# Calculamos valores min/max en el caso de no haber especificado "automatic"
# ------------------------------------------------------------------------
x_min = ifelse((x_min=="automatic"), min(x,na.rm=TRUE),                as.double (x_min))
x_max = ifelse((x_max=="automatic"), max(x,na.rm=TRUE),                as.double (x_max))
y_min = ifelse((y_min=="automatic"), min(y[,1:num_series],na.rm=TRUE), as.double (y_min))
y_max = ifelse((y_max=="automatic"), max(y[,1:num_series],na.rm=TRUE), as.double (y_max))

# ------------------------------------------------------------------------
# Comprobamos valor y_min para el eje Y: caso particular
# Si la escala es logarítmica en el eje 'y', el mínimo no puede ser 0.0 (error).
# ------------------------------------------------------------------------
if (y_log == "0") {
   y_min = 0
   y_log_str = ""
} else {
   if (y_min == 0.0) {
      y_min = 0.01
      writeLines ("   WARNING: parameter y_min cannot be equals to 0.0 with logarithmic Y-axis; forced to 0.01 !")
      flush.console()
   }
   y_log_str = "y"
}

# ------------------------------------------------------------------------
# Exportamos datos a LaTeX (archivo .tex): mismo nombre que el grafico pero .tex
# Lo hacemos despues de generar los graficos, por si el paquete "xtable" no
# estuviera instalado.
# latex_digits = decimales
# ------------------------------------------------------------------------
col_names = c(format(round(x[,1], latex_digits), nsmall=latex_digits))
exportar_datos_latex (title, graph_filename, series_names,
	col_names, y, col_x_values_len, latex_digits, verbose)	


# -------------------------------------------------------------------------
# Generamos graficos
# -------------------------------------------------------------------------

if (verbose == "1") {
   writeLines ("Generating graphs ...")
}

# Invocamos a la funcion correspondiente segun tipo de grafico:
#    0-lineas
#    1-puntos
#    2-lineas & puntos
#    3-lineas & puntos superpuestos
#    4-histograma (lineas verticales)
#    5-Grafico de barras apiladas
if ((graph_type >= 0) && (graph_type <= 4)) {
   plot_lines (title, subtitle, x_axis_title, y_axis_title,
               col_x_values, col_y_values, series_names,
               x_factor, y_factor, total_serie,
               x_min, x_max, y_min, y_max, y_log,
               show_titles, show_grid,
               show_confint, confint_as_percentage, 
               text_size_title, text_size_subtitle, text_size_axis_ticks,
               text_size_axis_titles, text_size_legend, pos_legend,
               graph_type, width_factor, height_factor, 
               graph_filename, graph_fileext_seq,
               x, y, data, num_series, num_series_legend, verbose)																									
} else if (graph_type == 5) {
   plot_bars  (title, subtitle, x_axis_title, y_axis_title,
               col_x_values, col_y_values, series_names,
               x_min, x_max, y_log, show_titles, show_grid,
               pos_legend, width_factor, height_factor,
               graph_filename, graph_fileext_seq, x, y, data, verbose)

}

if (verbose == "1") {
   writeLines ("Ok\n")
}
flush.console()

#return (TRUE)

}
																															 
# =========================================================================

# =========================================================================
# plot_lines: dibujar gráficas de varias series (cada serie una línea)
#             Se generan tantos gráficos como extensiones especificadas
#             (.png, .ps, ...)
# =========================================================================

plot_lines = function ( title, subtitle, x_axis_title, y_axis_title,
      col_x_values, col_y_values, series_names,
      x_factor, y_factor, total_serie,
      x_min, x_max, y_min, y_max, y_log,
      show_titles, show_grid,
	  show_confint, confint_as_percentage, 
      text_size_title, text_size_subtitle, text_size_axis_ticks,
      text_size_axis_titles, text_size_legend, pos_legend,
      graph_type=0, width_factor=1.0, height_factor=1.0, 
      graph_filename, graph_fileext_seq,
      x, y, data, num_series, num_series_legend, verbose ) {

# ------------------------------------------------------------------------
# Construimos el gráfico
# ------------------------------------------------------------------------

# Contamos Valores del eje X e Y
col_x_values_len = NROW(data)
col_y_values_len = NROW(col_y_values)
series_names_len = NROW(series_names)

# Separamos la lista de extensiones (al menos debe haber una)
graph_fileext_len = NROW(graph_fileext_seq)
if (graph_fileext_len < 1) {
   stop ("ERROR: graph parameter 'graph_fileext_seq' is empty.")
}

# ------------------------------------------------------------------------
# Modificamos valores min/max para el eje Y: caso particular
# Tenemos en cuenta si queremos eje Y en escala logarítmica o no
# ------------------------------------------------------------------------
if (y_log == "0") {
   y_log_str = ""
} else {
   y_log_str = "y"
}
   
# Generamos el grafico en los formatos especificados
for (ext_idx in 1:graph_fileext_len) {
   # ------------------------------------------------------------------------
   # Construimos nombre de archivo para la extensión dada
   # ------------------------------------------------------------------------
   ext = graph_fileext_seq[ext_idx]
   graph_filenameext = paste (graph_filename, ".", ext, sep="")

   # ------------------------------------------------------------------------
   # Abrimos el gráfico: png(...), postscript(...), X11(), ...
   # Puede dar error si estamos en modo texto para generar los png().
   # ------------------------------------------------------------------------
   if ( abrir_grafico (graph_filenameext, ext, width_factor, height_factor) == FALSE ) {
      writeLines (paste("   WARNING: suffix '", ext, "' unsupported.", sep=""))
      next
   }
   if (verbose == "1") {
      writeLines (paste("   File: ", getwd(), "/", graph_filenameext, " ...", sep=""))
   }

   # ------------------------------------------------------------------------
   # Parámetros del gráfico: márgenes, ... 
   # mar = A numerical vector of the form c(bottom, left, top, right) which
   #       gives the number of lines of margin to be specified on the four
   #       sides of the plot.
   #       The default is c(5, 4, 4, 2) + 0.1. 
   # mgp = The margin line (in mex units) for the axis title, axis labels and axis line.
   #       The default is c(3, 1, 0). 
   # las = Orientación texto valores de los ejes:
   #       [0]-Paralelo al eje, 1-Horizontal, 2-Perpendicular, 3-Vertical
   #
   # cex.sub  = Tamaño subtitulos     (respecto 'cex')
   # cex.axis = Tamaño texto ejes     (respecto 'cex')
   # cex.lab  = Tamaño etiquetas ejes (respecto 'cex')
   # yaxp     = Posición de los valores en el eje y (ver 'at' en 'barplot').
   # ------------------------------------------------------------------------
   if (show_titles == "0") {
      par(mar=c(3.2, 2.9, 0.1, 0.1) + 0.1)   # Antes de SIMUTools: 4, 4, 0, 2
      par(mgp=c(2.0,0.7,0))      
      #legend_cex = 0.8
   } else {
      #legend_cex = 0.6
   }

   # --------------------------------------------------------------
   # Dibujamos el gráfico en blanco (redimensionado, titulos, ...)
   # Lo hacemos con cualquier dato de x, y (no se dibujará: type='n')
   # --------------------------------------------------------------
   # Tenemos en cuenta si mostramos titulo/subtitulo o no
   if (show_titles == "0") {
      title    = ""
      subtitle = ""
   }

   # Tamaño fuente letra: cex, cex.main, cex.sub, cex.axis, cex.lab (todos relativos a cex)
   # Valores anteriores: cex=0.8, cex.sub=0.7, cex.axis=1.0, cex.lab=1.0
   plot (x[,1], y[,1], main=title, sub=subtitle, type='n',
      xlab=x_axis_title, ylab=y_axis_title, cex=1.0, cex.main=text_size_title,
      cex.sub=text_size_subtitle, cex.axis=text_size_axis_ticks,
      cex.lab=text_size_axis_titles,
	  xlim=c(x_min,x_max), ylim=c(y_min,y_max), lty=1, log=y_log_str)

   # --------------------------------------------------------------
   # Si queremos que aparezcan 'ticks' en todos los ejes
   # - side: 1=below, 2=left, 3=above and 4=right
   # --------------------------------------------------------------
   #axis (side=3, at=seq(x_min, x_max, by=5), tck=0.01, labels=FALSE)
   #axis (side=4, at=seq(y_min, y_max, by=5), tck=0.01, labels=FALSE)

   # --------------------------------------------------------------
   # Dibujamos la rejilla (lo primero, para que no borre nada válido)
   # --------------------------------------------------------------
   if (show_grid == "1") {
      # Añadimos rejilla
      # - Si nx,ny=NULL -> Valor por defecto (marcas del eje), si NA -> sin linea en ese eje
      # - col = color de lineas
      # - lty = tipo de lineas
      # - lwd = grosor de lineas
      # - equilogs = lineas equidistantes (TRUE cuando nx,ny=NULL)
      grid(nx=NULL, ny=NULL, col="lightgray", lty="dotted", lwd=par("lwd"), equilogs=TRUE)
   }

   # --------------------------------------------------------------
   # Dibujamos intervalo de confianza
   # --------------------------------------------------------------
   if (show_confint == "1") {
      # Arrays para intervalo de confianza
      confint_top         = array(0, c(col_x_values_len,col_y_values_len))
      confint_bottom      = array(0, c(col_x_values_len,col_y_values_len))
      # Calculamos márgenes superior e inferior para el intervalo de confianza de cada serie
      for (x_idx in 1:col_x_values_len) {          
          for (serie in 1:col_y_values_len) {
              col = col_y_values[serie]
              # Obtenemos el margen del intervalo. 
              if (confint_as_percentage == "1") {
                 # Si el CI (en col+1) viene expresado en % tenemos que calcular el margen
                 margin = (data[x_idx,col] * data[x_idx,col+1]) / 100.0
              } else {
                 # Si el CI (en col+1) ya viene como valor absoluto no hay que hacer nada
                 margin = data[x_idx,col+1] / 100.0   # <--- ???
              }
              confint_bottom[x_idx,serie] = (data[x_idx,col] - margin) * y_factor   # y[x_idx,serie]
              confint_top   [x_idx,serie] = (data[x_idx,col] + margin) * y_factor
          }
      }

      # --------------------------------------------------------------
      # Mostramos valores del intervalo de confianza
      # --------------------------------------------------------------
      #writeLines ("   Confidence Interval")
      #writeLines (paste("   Expresado en %? = ", confint_as_percentage, sep=""))
      #writeLines (paste("   [", confint_top, " .. ", confint_bottom, "]", sep=""))

      # Dibujamos intervalos de confianza
      for (serie in 1:col_y_values_len) {
          # Dibujamos lineas verticales
          segments (x[,serie], confint_top[,serie],
                    x[,serie], confint_bottom[,serie],
                    lwd=1.5, cex=0.5, col="grey65", lty=1);
          # Obtenemos dimensiones finales del grafico: vector con x-min, x-max, y-min, y-max
          plotDim = par("usr")
          plot_x_min = plotDim[1];
          plot_x_max = plotDim[2];
          plot_width = plot_x_max-plot_x_min
          # Dibujamos bigotes sup./inf. (ancho = 1% del grafico)
          ancho_bigotes = plot_width/100
          segments (x[,serie]-ancho_bigotes, confint_top[,serie],
                    x[,serie]+ancho_bigotes, confint_top[,serie],
                    lwd=1.5, cex=0.5, col="grey65", lty=1);
          segments (x[,serie]-ancho_bigotes, confint_bottom[,serie],
                    x[,serie]+ancho_bigotes, confint_bottom[,serie],
                    lwd=1.5, cex=0.5, col="grey65", lty=1);
      }
   }

   # --------------------------------------------------------------
   # Dibujamos las series
   # --------------------------------------------------------------

   # Si queremos cambiar los tipos de lineas o colores por defecto (colors())
   #linetypes  = c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
   colores    = c( "blue", "darkgreen", "darkviolet", "darkorange1", "chocolate4",
                   "cornsilk4", "darkolivegreen", "darkgoldenrod1",
                   "darkred", "deepskyblue4", "greenyellow",
                   "indianred4", "khaki4",
                   "midnightblue", "navyblue", "aquamarine4", "dodgerblue4")

   # Para ajustes de lty, lwd, pch, col al dibujar series y leyenda
   lty = array (""         , c(num_series_legend))  # Tipo de línea
   lwd = array (1          , c(num_series_legend))  # Grosor de línea
   pch = array (NA_integer_, c(num_series_legend))  # Tipo de punto
   col = array (""         , c(num_series_legend))  # Color de línea

   # Las serie total tendrá otro color, rayado y sin puntos (solo linea)
   if (total_serie > 0) {
      s = num_series_legend
      lty[s] = "dotted"
      lwd[s] = 3
      pch[s] = NA_integer_   # o "" (sin puntos)
      col[s] = "red"      
      lines(x[,1],  y[,num_series], type="l", lty=lty[s], lwd=lwd[s], pch=pch[s], col=col[s])
   }
   # Tipo de gráfico para las series: 0-líneas, 1-puntos, 2-lineas & puntos, o 3-histograma
   tipo = "l"
   if (graph_type == 0) {          # Como lineas
      tipo = "l"     
   } else if (graph_type == 1) {   # Como puntos
      tipo = "p"
   } else if (graph_type == 2) {   # Como lineas & puntos (entrecortados)
      tipo = "b"
   } else if (graph_type == 3) {   # Como lineas & puntos (superpuestos)
      tipo = "o"
   } else if (graph_type == 4) {   # Como histograma
      tipo = "h"              
   } else {
      tipo = "o"                   # Como lineas & puntos (superpuestos) (por defecto)
   }  
   for (serie in col_y_values_len:1) {
       # Serie de datos
       lty[serie] = "solid"  # linetypes[serie+2] o serie (da error en lty=lty[serie], pero si se pone lty=serie no)
       lwd[serie] = 1
       #if ( show_hotspots=="1" ) {
       if ((graph_type >= 1) && (graph_type <= 3)) {
          pch[serie] = serie
       }
       col[serie] = colores[serie]
       lines(x[,1],  y[,serie], type=tipo, lty=lty[serie], lwd=lwd[serie], pch=pch[serie], col=col[serie])
       #if ( show_hotspots=="1" ) {
       #   points(x[,1], y[,serie], pch=serie)  # Si quisiéramos mostrar sólo puntos
       #}
   }

   # --------------------------------------------------------------
   # Añadimos la leyenda
   # --------------------------------------------------------------
   if (pos_legend > 0) {
      # Obtener la posición de la leyenda en R a partir del valor 1..9
      xLegend = obtener_pos_leyenda (pos_legend)

      # Añadimos la leyenda
      # - x,y = coordenadas esquina sup. izquierda del cuadro de leyenda
      # - legend = vector con el texto de cada elemento de la leyenda
      # - cex = tamaño de las letras del texto de cada elemento de la leyenda
      # - lty = tipo de linea de cada elemento de la leyenda
      # - pch = caracter de los puntos de las lineas del gráfico
      # - box.lty = linea del cuadro (0:sin cuadro)
      # - bty  = fondo blanco ("o") o transparente ("n")
      # - fill = color de relleno (NULL sin relleno, por defecto)
      # - xjust = justificación (0:izq, por defecto, 0.5:centrado, 1:derecha)
      # - title = titulo de la leyenda (NULL sin titulo, por defecto)
      # - inset = inset distance(s) from the margins as a fraction of the plot region when legend is placed by keyword 

      legend (x=xLegend, legend=series_names, cex=text_size_legend,  #cex=legend_cex,
              lty=lty, lwd=lwd, pch=pch, col=col, 
              box.lty=0, xjust=0, inset=c(0.01,0.01)) # y=yLegend,
   }
   # Cerrar el gráfico
   #graphics.off()  # Error in grDevices::dev.off(): cannot shut down device 1 (the null device)
   dev.off()
} # Extension

}

# =========================================================================

# =========================================================================
# plot_bars: dibujar gráficas porcentajes
# =========================================================================

plot_bars = function ( title, subtitle, x_axis_title, y_axis_title,
    col_x_values, col_y_values, series_names, x_min, x_max, y_log,
    show_titles, show_grid, pos_legend, width_factor=1.0, height_factor=1.0,
    graph_filename, graph_fileext_seq, x, y, data, verbose ) {

# ------------------------------------------------------------------------
# Modificamos valores min/max para el eje Y: caso particular
# Tenemos en cuenta si queremos eje Y en escala logarítmica o no
# ------------------------------------------------------------------------
if (y_log == "0") {
   y_min = 0
   y_log_str = ""
} else {
   y_log_str = "y"
}
y_max = 100.1   # Debe ser > 100.0; Si es = 100.0 no llega a salir el valor 100 en el eje.

# ------------------------------------------------------------------------
# Construimos el gráfico
# ------------------------------------------------------------------------

# Contamos Valores del eje X e y
col_x_values_len = NROW(data)
col_y_values_len = NROW(col_y_values)
series_names_len = NROW(series_names)

# Separamos la lista de extensiones (al menos debe haber una)
graph_fileext_len = NROW(graph_fileext_seq)
if (graph_fileext_len < 1) {
   stop ("ERROR: graph parameter 'graph_fileext_seq' is empty.")
}

# Generamos el grafico en los formatos especificados
for (ext_idx in 1:graph_fileext_len) {
   ext = graph_fileext_seq[ext_idx]

   # ------------------------------------------------------------------------
   # Abrir el gráfico: png(...), postscript(...), X11(), ...
   # ------------------------------------------------------------------------
   # Construimos nombre de archivo
   graph_filenameext = paste (graph_filename, ".", ext, sep="")
   if ( abrir_grafico (graph_filenameext, ext) == FALSE ) {
      writeLines (paste("   WARNING: suffix '", ext, "' unsupported.", sep=""))
      next
   }
   if (verbose == "1") {
      writeLines (paste("   File: '", getwd(), "/", graph_filenameext, "' ...", sep=""))
   }

   # ------------------------------------------------------------------------
   # Parámetros del gráfico: márgenes, ... 
   # Tenemos en cuenta si mostramos titulos/subtitulos o no
   # mar = A numerical vector of the form c(bottom, left, top, right) which
   #       gives the number of lines of margin to be specified on the four
   #       sides of the plot.
   #       The default is c(5, 4, 4, 2) + 0.1. 
   # mgp = The margin line (in mex units) for the axis title, axis labels
   #       and axis line. The default is c(3, 1, 0). 
   # las = Orientación texto valores de los ejes:
   #       [0]-Paralelo al eje, 1-Horizontal, 2-Perpendicular, 3-Vertical
   #
   # cex.sub  = Tamaño subtitulos     (respecto 'cex')
   # cex.axis = Tamaño texto ejes     (respecto 'cex')
   # cex.lab  = Tamaño etiquetas ejes (respecto 'cex')
   # yaxp     = Posición de los valores en el eje y (ver 'at' en 'barplot').
   # ------------------------------------------------------------------------
   if (show_titles == "0") {
      title    = ""
      subtitle = ""
      par(mar=c(3.2, 2.9, 1.3, 0.1) + 0.1)
      par(mgp=c(2.0,0.7,0))
      legend_cex = 0.8
   } else {
      legend_cex = 0.6
   }
   par(cex.sub=0.7, cex.axis=1.0, cex.lab=1.0, las=0)  #, yaxt="n")

   # --------------------------------------------------------------
   # Establecemos numero de decimales y para cuantos se utiliza coma flotante
   # --------------------------------------------------------------
   # Para que las etiquetas salgan con un numero de decimales podemos cambiar
   # las opciones por defecto. Para consultarlas: options() o getOption('digits').
   #  digits: Controls the number of digits to print when printing numeric values.
   #          It is a suggestion only. Valid values are 1...22 with default 7.
   #  scipen: A penalty to be applied when deciding to print numeric values in
   #          fixed or exponential notation. Positive values bias towards fixed
   #          and negative towards scientific notation: fixed notation will be
   #          preferred unless it is more than scipen digits wider. 
   options('digits'=3)
   options('scipen'=0)  # Para >= 1 empieza a salir con decimales

   # --------------------------------------------------------------
   # Dibujamos las series
   # --------------------------------------------------------------

   # Invertimos el orden de las series de datos porque en las barras apiladas
   # se empieza desde abajo a arriba
   bandwidth_share = array(0, c(series_names_len,col_x_values_len))
   for (x_idx in 1:col_x_values_len) {
      #bandwidth_share[,x_idx] = rev(y[x_idx,])

      # Total acumulado (suma de todas las series) -> para la tabla latex (solo datos)
      total = sum (y[x_idx,1:series_names_len], na.rm=TRUE)
      for (serie in 1:series_names_len) {
         bandwidth_share[series_names_len-serie+1,x_idx] = (100.0 * y[x_idx,serie]) / total
      } #serie
   }

   # Color de las series en las columnas y nombres de las columnas
   seriesbarcolors = grey(0.5 + 1:series_names_len/12)    # c("grey90", "grey70"))
   namesbar = x[,1]  # c(1:col_x_values_len)  # "1", "2", "3", ...
      
   # Dibujamos las barras (axes=TRUE para mostrar el eje 'y')
   midpts = barplot (bandwidth_share, main=title, sub=subtitle,
         xlab=x_axis_title, ylab=y_axis_title, ylim=c(y_min,y_max), axes=TRUE,
         col=seriesbarcolors, log=y_log_str, names.arg=namesbar)
         #yaxt="n")    # Sin eje Y (se pondrá después manualmente)
         #legend=seriesnames[series_names_len:1])    # O bien esta leyenda, o separada (ver abajo)   

   # Eje X (para evitar que salgan valores en coma flotante)
   # Habria que poner yaxt="n" en el barplot, y descomentar esto.
   # side: 1=below, 2=left, 3=above and 4=right
   # las = Orientación texto valores de los ejes:
   #       [0]-Paralelo al eje, 1-Horizontal, 2-Perpendicular, 3-Vertical
   #axis (side=2, yaxp=c(1, 10, 20, 50, 75, 100))
   #axis (side=2, at=seq(0,25,50,75,100), labels=seq(0,25,50,75,100))
   #axis (side=2, at=seq(0,25,50,75,100), labels=TRUE, las=1)
   #axis (side=2, at=seq(0,100, by = 25), labels=TRUE, las=1)
   #if ((y_log=="y") || (y_log=="xy")) {
   #   # Si el eje es logaritmico, y_min>0. Para que aparezca el 0.
   #   text (x=0, y=100, labels="0", pos=4)
   #}

   # En caso de nombres muy largos poner lo siguiente en vez de names.arg en barplot, y names=rep("", series_names_len)
   # Substituimos espacios por nueva linea en esa cadena.
   #mtext (sub(" ", "\n", colnames(bandwidth_share)), at=midpts, side=1, line=0.5, cex=0.5)

   # Calculamos puntos intermedios de cada fragmento y mostramos el valor numérico
   midpts_x = rep (midpts, each=series_names_len)
   midpts_y = apply (bandwidth_share, 2, cumsum) - bandwidth_share/2  # Aplicamos la función "cumsum" (suma acumulativa) a las columnas
   seriestextcolors = rep (c("white", "black"), times=2:2, cex=0.8)  
   bandwidth_share_rounded = round (bandwidth_share, digits=1)
   text (midpts_x, midpts_y, bandwidth_share_rounded, col=seriestextcolors, cex=0.8)

   # --------------------------------------------------------------
   # Dibujamos la rejilla (en otros gráficos se hace antes; aqui tras el barplot)
   # --------------------------------------------------------------
   if (show_grid == "1") {
      # Añadimos rejilla
      # - Si nx,ny=NULL -> Valor por defecto (marcas del eje), si NA -> sin linea en ese eje
      # - col = color de lineas
      # - lty = tipo de lineas
      # - lwd = grosor de lineas
      # - equilogs = lineas equidistantes (TRUE cuando nx,ny=NULL)
      grid (nx=NULL, ny=NULL, col="lightgray", lty="dotted", lwd=par("lwd"), equilogs=TRUE)
   }

   # --------------------------------------------------------------
   # Añadimos la leyenda (se puede poner aparte o junto con el barplot)
   # --------------------------------------------------------------
   if (pos_legend > 0) {
      # Obtener la posición de la leyenda en R a partir del valor 1..9
      xLegend = obtener_pos_leyenda (pos_legend)

      # Añadimos la leyenda
      # - x,y = coordenadas esquina sup. izquierda del cuadro de leyenda
      # - legend = vector con el texto de cada elemento de la leyenda
      # - cex = tamaño de las letras del texto de cada elemento de la leyenda
      # - lty = tipo de linea de cada elemento de la leyenda
      # - pch = caracter de los puntos de las lineas del gráfico
      # - box.lty = linea del cuadro (0:sin cuadro)
      # - bty  = Tipo de caja a dibujar (["o"] o "n")
      # - bg   = color del fondo ("white", ...)
      # - fill = color de relleno (NULL sin relleno, por defecto)
      # - xjust = justificación ([0:izq], 0.5:centrado, 1:derecha)
      # - title = titulo de la leyenda (NULL sin titulo, por defecto)
      # - inset = inset distance(s) from the margins as a fraction of the
      #           plot region when legend is placed by keyword.
      legend (x=xLegend, legend=series_names, cex=legend_cex,
              box.lty=0, bty="o", fill = rev(seriesbarcolors), bg="white",
              xjust=0, inset=c(0.01,0.00))
   }
   # Cerrar el gráfico
   #graphics.off()  # Error in grDevices::dev.off(): cannot shut down device 1 (the null device)
   dev.off()
} # Extension

}

return

# =========================================================================


# =========================================================================
# Recortar espacios en blanco de la izquierda y la derecha
# =========================================================================
trim<-function (str) {
   str <- gsub ("(^ +)|( +$)","",str)
   return(str)
}


# =========================================================================
# Recortar espacios en blanco de la izquierda
# =========================================================================
ltrim<-function (str) {
   str <-sub ("^ +","",str)
   return (str)
}


# =========================================================================
# Recortar espacios en blanco de la derecha
# =========================================================================
rtrim<-function (str) {
   str <- sub (" +$","",str)
   return (str);
}


# =========================================================================
# importar_datos: leer un fichero de datos y resumirlos para gráficas
# =========================================================================

importar_datos = function (data_path, data_filename, data_fileext, verbose) {
   # Construimos el nombre del fichero de entrada (datos) y salida (gráfico)
   filename  = paste (data_filename, ".", data_fileext, sep="")
   filename = file.path (data_path, filename, fsep = .Platform$file.sep)
   # Leemos los datos
   if (verbose == "1") {
      writeLines (paste("Loading data file '", filename, "' ...", sep=""))
   }
   data = read.table (filename, header=TRUE)
   return (data)
}


# =========================================================================
# exportar_datos: escribir un fichero de datos ya resumidos para gráficas
# =========================================================================

exportar_datos = function (data_path, data_filename, data_fileext, x, y, verbose) {
   # Creamos un data frame combinando los vectores x e y
   data = data.frame(y)  # cbind(x, y)
   data.row.names = x
   # Construimos el nombre del fichero de salida
   filename  = paste (data_filename, ".", data_fileext, sep="")
   filename = file.path (data_path, filename, fsep = .Platform$file.sep)
   # Escribimos los datos
   if (verbose == "1") {
      writeLines (paste("Writting data file '", filename, "' ...", sep=""))
   }
   write.table (data, filename, row.names=TRUE, col.names=TRUE)
}


# =========================================================================
# exportar_datos_latex: escribir datos en un archivo .tex para latex
# =========================================================================

exportar_datos_latex = function (title, filename, row_names, col_names,
                                 data, col_x_values_len, latex_digits, verbose) {
   # Nombre de archivo (igual que el gráfico pero .tex)
   filename_tex = paste (filename, ".tex", sep="")
   if (verbose == "1") {
      writeLines (paste("Exporting data to LaTeX (", latex_digits, " decimals) ...", sep=""))
      writeLines (paste("   File: ", getwd(), "/", filename_tex, sep=""))
   }

   # Tabla latex (reservamos una fila y una columna para cabeceras)
   # Por seguridad, ponemos un caracter de escape delante de los subrayados
   # Al intentar cargar este fichero con source R da un error:
   #Error: '\_' is an unrecognized escape in character string starting "\_"
   #new_names = sapply (strsplit(row_names,"_"), paste, collapse="\_")
   # Poniendo dos "\\_" aparece en la tabla $\backslash$\_
   # Con versiones recientes de xtable no es necesario; si se sustituye cada subrayado
   # por "\\_" aparece en la tabla una barra adicional (lo sustituye por $\backslash$\_)
   #new_names = sapply (strsplit(row_names,"_"), paste, collapse="\\_")
   new_names = row_names
   
   table_names = list (new_names, col_names)
   table_data  = array (0, c(NROW(row_names), NROW(col_names)), dimnames=table_names)
   for (row in 1:NROW(row_names)) {
      for (col in 1:NROW(col_names)) {
         table_data[row,col] = data[col,row]
      }
   }
   table_digits  = array (latex_digits, c(col_x_values_len+1))   # decimales en todas las columnas (1 valor=se repite)
   table_caption = paste (title, " data", sep="")
   #table_label  = paste (simulset, ":data:", title, sep="")
   table_label   = paste ("data:", title, sep="")
   tabla_latex   = xtable (table_data, digits=table_digits, caption=table_caption, label=table_label)

   # Package: xtable - Version: 1.4-3
   #print.xtable (tabla_latex, file=filename_tex, table.placement="h!", size="scriptsize")
   # Package: xtable - Version: 1.5-2
   print (tabla_latex, file=filename_tex, table.placement="h!", size="scriptsize")
}


# =========================================================================
# abrir_grafico: abrimos un gráfico
# Atención, los png() requieren entorno gráfico (p.e. X11 o Windows),
# no se pueden generar desde una consola en modo texto.
# Soluciones: 
# 1. Utilizar bitmap() para generar el png.
# 2. Utilizar Xvfb (servidor X falso)
#    http://www.X.org
#    http://xorg.freedesktop.org
#    http://lists.freedesktop.org/mailman/listinfo/xorg
# =========================================================================

abrir_grafico = function (filename, fileext, width_factor=1.0, height_factor=1.0) {
   if ((fileext == "ps") || (fileext == "eps")) {
      # Archivo .ps o .eps
      try (postscript (filename, onefile=FALSE, horizontal=FALSE, pointsize=12,
           width=6*width_factor, height=6*height_factor), silent=FALSE)
   } else if (fileext == "png") {
      # Archivo .png
      # OJO: Para generar los .png, en windows se puede utilizar png(),
      # pero en Linux sólo se puede utilizar png desde X11, no desde la consola.
      # Para evitar problemas, se utiliza bitmap(), que no requiere X11.
      if ( capabilities("png") ) {
         #writeLines ("    Modo grafico: generando .png con png()")
         try (png (filename, width=480*width_factor, height=480*height_factor), silent=FALSE)
      } else {
         #writeLines ("    Modo no grafico: generando .png con bitmap()")
         try (bitmap (filename, type="pnggray"), silent=FALSE)         
      }
   } else {
      # No se reconoce la extensión
      #try (X11(), silent=FALSE)   # Generamos una vista previa
      return (FALSE)
   }
   return (TRUE)
}


# =========================================================================
# obtener_pos_leyenda: obtener la posición de la leyenda en R a partir del
#                      valor 1..9, para usar en la llamada a 'legend'.
# Retorna: el parámetro 'xlegend' en forma de cadena.
# =========================================================================

obtener_pos_leyenda = function (posLegend) {
      if (posLegend==1) {            # Esquina superior izquierda
         xLegend = "topleft"
      } else if (posLegend==2) {     # Centro superior
         xLegend = "top"
      } else if (posLegend==3) {     # Esquina superior derecha
         xLegend = "topright"
      } else if (posLegend==4) {     # Centro izquierda
         xLegend = "left"
      } else if (posLegend==5) {     # Centro Centrado
         xLegend = "center"
      } else if (posLegend==6) {     # Centro derecha
         xLegend = "right"
      } else if (posLegend==7) {     # Esquina inferior izquierda
         xLegend = "bottomleft"
      } else if (posLegend==8) {     # Centro inferior
         xLegend = "bottom"
      } else if (posLegend==9) {     # Esquina inferior derecha
         xLegend = "bottomright"
      } else {
         # Otras posibles posiciones asignando a x alguno de los valores:
         # "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". 
         xLegend = "bottomright"
         #yLegend = 0.0  # No importa el valor
      }
      return (xLegend)
}

# =========================================================================

