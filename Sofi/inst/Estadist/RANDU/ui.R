#library(d3heatmap)
library("shiny")
library("threejs")
fluidPage(
  h2("RANDU: el caso de los malos RNG"),
  h5("Ejemplo original de", strong("csgillespie") ,"difundido en",a(href="https://csgillespie.wordpress.com/2016/02/16/randu-the-case-of-the-bad-rng/", target="_blank", "esta página,")),
  br(),
  p("La Oficina Federal Alemana para la Seguridad de la Información (BSI) ha establecido criterios de calidad del generador de números aleatorios (RNG):"),
  br(),
  p("1.-Una secuencia de números aleatorios tiene una alta probabilidad de no contener elementos consecutivos idénticos."),
  p("2.-Una secuencia de números aleatorios es indistinguible de los verdaderos números. (probado mediante pruebas estadísticas)"),
  p("3.-Debe ser imposible calcular o adivinar, desde cualquier sub-secuencia dada, los valores anteriores o futuros en la secuencia."),
  p("4.-Debe ser imposible, para todos los efectos prácticos, para que un atacante calcular o adivinar los valores utilizados en el algoritmo de números aleatorios"),
  br(),
  withMathJax(),
  p("Los puntos 3 y 4 son cruciales para muchas aplicaciones. Cada vez que hagas una llamada telefónica, contacto a un punto inalámbrico, pagar con tarjeta de crédito se utilizan números aleatorios."),
  p("El diseño de un buen generador de números aleatorios es dura y, como regla general nunca se debe intentar. R viene con muchos generadores aleatorios de buena calidad. El generador predeterminado es el ",
    a(href="https://en.wikipedia.org/wiki/Mersenne_Twister", target="_blank", "Mersenne-Twister.")),
  p(" Este RNG tiene un enorme período de \\({{2}^{19937}}-1\\) (la cantidad de números aleatorios que se generan antes de que tengamos una repetición)."),
  br(),
  h3("Generadores de congruencia lineal"),
  p("Un generador de congruencia lineal (LCG) es un generador de números aleatorios relativamente simple (popular en los años 60 y 70). Tiene una forma simple de"),
  p("$${ { r }_{ i } }=\\left( { { { ar }_{ i-1 }+b } } \\right) \\quad mod\\quad m,\\quad i=1,2,...,m $$"),
  p("donde \\({{r}_{0}} \\) es el número inicial, conocida como la semilla, y \\(\\left( a, b, m \\right)\\)  son el multiplicador, constante aditiva y modulo respectivamente. Los parámetros son todos números enteros."),
  p("La operación de módulo significa que la mayoría de m números diferentes puede ser generada antes de que la secuencia deba repetir - a saber, los números enteros \\(0,1,2,...,m-1\\)."),
  p("El actual numero fue generado por los números \\(h\\le m\\) ,llamado el periodo de el generador."),
  p("La clave para generadores de números aleatorios en los ajuste de los parámetros"),
  br(),
  p(""),
  h4("RANDU"),
  p("RANDU es un lcg con parámetros"),
  p("Código original en ",a(href="https://gist.github.com/csgillespie/0ba4bbd8da0d1264b124", target="_blank", "esta página.")),
  p("Por desgracia, los valores iniciales son una extraordinariamente mala elección de parámetros."),
  p('Con los valores predeterminados fácilmente se distingue un ciclo en los números  “aleatorios” por lo que 
    la elección de los parámetros a=68, b=7, m=3000 es pésima. Jugando con los parámetros (por ejemplo, a=63, b=7, m=2997) 
    se pueden ver muchos casos donde fácilmente se encuentran comportamientos predecibles.'),
  p('Moraleja, en muchos casos dados como  “caos” se puede encontrar el orden.'),
  
  fluidRow(
  column(width = 1, offset = 4,numericInput("obs_a", "a:", 68, min = 1, max = 100000)),
  column(1,numericInput("obs_b", "b:", 7, min = 0, max = 100000)),
  column(1,numericInput("obs_m", "m:", 3000, min = 1, max = 3147483648)),
  column(1,numericInput("obs_N", "Generar:", 1000, min = 1, max = 3000))
  ),
  #actionButton("goButton", "Vamos"),
  uiOutput("Ecu"),
  p("$${ { r }_{ i } }=\\left( { { { ar }_{ i-1 }+b } } \\right) \\quad mod\\quad m,\\quad i=1,2,...,m $$"),
  actionButton("sal", "Salir"),
  br(),
  br(),
  plotOutput("plot_Randu"),
  br(),
  fluidRow(
  column(6,h4("RANDU 3D"),scatterplotThreeOutput("plot_Randujs")),
  column(6,h4("Estándar rng 3D"),scatterplotThreeOutput("plot_rujs"))
  ),
  
  p("Para más información ver ",a(href="https://en.wikipedia.org/wiki/Linear_congruential_generator", target="_blank", "aquí.")),
  p("Sugerencia de ejercicio: crear  un generador de números aleatorios y probar gráficamente los resultados."),
  checkboxInput("Ver", "Datos de RANDU", FALSE),
  
  #conditionalPanel(condition = "Ver",uiOutput("Ecu")),
  conditionalPanel(condition = "input.Ver",
                   fluidRow(
                     column(3,tableOutput("tableRANDU1")),
                     column(3,tableOutput("tableRANDU2")),
                     column(3,tableOutput("tableRANDU3")),
                     column(3,tableOutput("tableRANDU4"))
                   ))
  
  #p("$${ { r }_{ i } }=\\left( { { { ar }_{ i-1 }+b } } \\right) \\quad mod\\quad m,\\quad i=1,2,...,m $$"),
  #conditionalPanel(condition = "input.Ver==TRUE",tableOutput("tableRANDU"))
  
  
)