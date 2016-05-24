# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  # Application title
  headerPanel("SOMbrero Web User Interface (v1.0)"),

  #### Panel 'About' (right hand side)
  ##############################################################################
  sidebarPanel(
    selectInput('somtype', "Select the data type:",
                c(Numeric="numeric", Korresp="korresp",
                  Relational="relational")),
    
    imageOutput("sombrero.logo", inline=TRUE),
    p(HTML("<h5>Welcome to SOMbrero, the open-source on-line interface for 
self-organizing maps (SOM).</h5> This interface trains SOM for numerical data,
contingency tables and dissimilarity data using the <strong>R</strong> package
<a href='http://sombrero.r-forge.r-project.org/'>SOMbrero</a> (v1.1). Train a
map on your data and visualize their topology in three simple steps using the
panels on the right.")),
    
    imageOutput("samm.logo", inline=TRUE),
    br(),
    imageOutput("miat.logo", inline=TRUE),
    br(),
    br(),
    p(HTML('It is kindly provided by the
<a href= "http://samm.univ-paris1.fr/">SAMM</a> team and the 
<a href= "carlit.toulouse.inra.fr">MIA-T</a> team under the
<a href= "https://www.gnu.org/licenses/gpl-2.0.html">GPL-2.0</a>
license, and was developed by Julien Boelaert, 
<a href= http://samm.univ-paris1.fr/-Madalina-Olteanu->Madalina Olteanu</a> and
<a href= http://www.nathalievilla.org/> Nathalie Villa-Vialaneix</a>, using
<a href="http://www.rstudio.com/shiny/">Shiny</a>. It is also included in the 
<strong>R</strong> package 
<a href="http://sombrero.r-forge.r-project.org/">SOMbrero</a>. Its source code
is freely available on github: <br> 
<span style="font-size:12px;font-family:courrier; background-color:#FADDF2;
border:1px solid black;"><font color="#870500"><b>
git clone https://github.com/tuxette/sombrero.git</b></font></code></span>')),
    br(),
    
    h3('References:'),
    p(HTML('<li> <span style="font-variant: small-caps;">Kohonen T.</span>(2001)
<I>Self-Organizing Maps</I>. Berlin/Heidelberg: Springer-Verlag, 3rd edition.
</li><li> <span style="font-variant: small-caps;">Cottrell M., Letremy P., Roy
E.</span> (1993) Analyzing a contingency table with Kohonen maps: a Factorial
Correspondence Analysis. In: <I>Proceedings of IWANN’93</I>,
<span style="font-variant: small-caps;">J. Cabestany, J. Mary, A. Prieto</span>
(Eds.), <I>Lecture Notes in Computer Science, </I>Springer-Verlag, 305–311.</li>
<li> <span style="font-variant: small-caps;">Olteanu M., Villa-Vialaneix N.,
Cottrell M.</span> (2015) On-line relational and multiple relational SOM.
<em>Neurocomputing</em>, <strong>147</strong>, 15-30.</li>'))
  ),

  mainPanel(
    tabsetPanel(
      #### Panel 'Import Data'
      #########################################################################
      tabPanel("Import Data",
               h3("First step: import data"),

               p(HTML("To run the application, import your data set using the
import button below. Your data must be supplied in the form of a text/csv file.
If the importation is done properly, a preview of the data is displayed below.
When this is done, choose the SOM type of the left hand side panel and proceed
to the next step: self-organize a map.")),
               p(HTML('The interface can be tested using example data files for
the <a href= 
"http://owncloud.nathalievilla.org/public.php?service=files&t=a4b83ca82bcfc740ee0700cb44324c47"
target="_blank">numeric</a>, <a href= 
"http://owncloud.nathalievilla.org/public.php?service=files&t=7a7db6e0a17127d8fe8ec43c0f6b0afd"
target="_blank">korresp</a> and <a href= 
"http://owncloud.nathalievilla.org/public.php?service=files&t=2bd5a14bbf636ab9afe495b46e0d523a"
target="_blank">relational </a> algorithms (download these files on your computer and
proceed).')),
               
               br(), 
               fileInput('file1', 'Choose CSV/TXT File'),
               checkboxInput('header', ' Header?', TRUE),
               checkboxInput('rownames', ' Row names?', FALSE),
               selectInput('sep', 'Separator:',
                           c("Comma","Semicolon","Tab","Space"), 'Comma'),
               selectInput('quote', 'Quote:',
                           c("None","Double Quote","Single Quote"), 
                           'Double Quote'),
               selectInput('dec', 'Decimal mark', c("Period", "Comma"),
                           'Period'),
               numericInput('nrow.preview','Number of rows in the preview:',20),
               numericInput('ncol.preview', 'Number of columns in the preview:',
                            10),
               helpText("Note: Even if the preview only shows a restricted
number of observations, the map will be based on the full dataset."),
               tableOutput("view")
      ),
      
      #### Panel 'Self-organize'
      #########################################################################
      tabPanel("Self-Organize",
               h3("Second step: train the self-organizing map"),
               p(HTML("Once your dataset is imported, you can train a
self-organizing map (SOM) and explore it using the options below. You can then
download the resulting SOM in .rda format (you will need <strong>R</strong> and
the package <a href= 'http://sombrero.r-forge.r-project.org/'>SOMbrero</a> to
open this file and use the SOM; its class is the 'somRes' class, handled by
<a href='http://sombrero.r-forge.r-project.org/'>SOMbrero</a>). You can also
explore it using the next panels to visualize the results, compute super-classes
or combine it with additional variables.")),
               p(HTML('Consult the "Help" panel for information on how to choose
adequate parameter values.')),
               actionButton("trainbutton","Train SOM"),
               br(), br(),
               
               verbatimTextOutput("summary"),
               br(),
               downloadButton("som.download", "Download the SOM file"),
               
               br(), br(),
               h4("Options"),
               uiOutput("varchoice"),
               numericInput("dimx", "Map dimension X:", 5, min= 1),
               numericInput("dimy", "Map dimension Y:", 5, min= 1),
               h4("Advanced options"),
               selectInput("affectation", "Affectation type:", 
                           c("standard","heskes")),
               uiOutput("initproto"),
               numericInput("maxit", "Max. iterations:", 500),
               uiOutput("disttype"),
               selectInput("radiustype", "Radius type:", 
                           c("letremy","gaussian")),
               uiOutput("scaling"),
               numericInput("randseed",
                            HTML("Set a random seed for reproducible results
<a href='#pseudor'><sup>(1)</sup></a>:"),
                            sample(1:1e5, size= 1)),
               numericInput("eps0", "Scaling value for gradient descent", 1,
                            min= 0.01, step= .01),
               numericInput("nb.save", "Number of intermediate back-ups:", 0,
                            min= 0), 
               p(HTML("<span style='font-size:10px'><a name='pseudor'><sup>(1)
</sup></a> SOMbrero is based on a stochastic (on-line) version of the SOM
algorithm and thus uses randomness. Setting a seed results in fixing the random
procedure in order to obtain reproducible results (runing several times the
process with the same random seed will give the same map). More information on
pseudo-random generators at
<a href='http://en.wikipedia.org/wiki/Pseudorandom_number_generator'>this link
</a></span>."))),

      #### Panel 'Plot'
      #########################################################################
      tabPanel("Plot Map",
               h3("Third step: plot the self-organizing map"),
               p("In this panel and the next ones you can visualize the computed 
self-organizing map. This panel contains the standard plots used to analyze the
map."),
               
               h4("Options"),
               selectInput("somplotwhat", "Plot what?", 
                           choices= list("Observations"= "obs",
                                         "Prototypes"= "prototypes",
                                         "Energy"= "energy")),
               selectInput("somplottype", "Type of plot:", 
                           choices= c("hitmap", "color", "lines", "barplot", 
                                      "names", "boxplot", "radar")),
               checkboxInput("somplottitle", "Show cluster names"),
               conditionalPanel("input.somplottype == 'color' ||
                         input.somplottype == '3d'",
                                selectInput("somplotvar", 
                                            "Variable: (only used for '3d',
'color' and 'boxplot' plots if available)", 
                                            choices= "(Not Available)")),
               conditionalPanel("input.somplottype == 'boxplot'",
                                selectInput("somplotvar2", 
                                            "Variable: (hold Ctrl to select
multiple variables)", 
                                            choices= "(Not Available)", 
                                            multiple= TRUE)),
               conditionalPanel("input.somtype == 'korresp'", 
                                selectInput("somplotrowcol", 
                                            "Plot rows or columns (when relevant
for the chart):",
                                            choices= list("columns"= "c", 
                                                          "rows"= "r"))),
               plotOutput("somplot")),
      
      #### Panel 'Superclasses'
      #########################################################################
      tabPanel("Superclasses",
               h3("Group prototypes into superclasses"),
               p("In this panel you can group the clusters into 'superclasses'
(using a hierarchical clustering on the neurons' prototypes), download the
resulting clustering in csv format and visualize it on charts. The 'dendrogram'
plot can help you determine the adequate number of superclasses."),
               
               selectInput("sc.cut.choice", "Choose clustering criterion:",
                           choices= c("Number of superclasses", 
                                      "Height in dendrogram"),
                           "Number of superclasses"),
               uiOutput("scHorK"), #  Superclass Height or K (nb of clusters)
               actionButton("superclassbutton","Compute superclasses"),
               downloadButton("sc.download",
                              "Download superclass classification"),
               
               br(), br(),
               verbatimTextOutput("sc.summary"),
               
               br(), br(),
               h4("Plot the superclasses:"),
               selectInput("scplotwhat", "Plot what?", 
                           choices= list("Prototypes"= "prototypes",
                                         "Observations"= "obs")),
               selectInput("scplottype", "Type of plot:", 
                           choices= c("hitmap", "color", "lines", "barplot", 
                                      "names", "boxplot", "radar")),
               conditionalPanel("input.scplottype == 'color' ||
                         input.scplottype == '3d'",
                                selectInput("scplotvar", 
                                            "Variable: (only used for '3d',
'color' and 'boxplot' plots if available)",
                                            choices= "(Not Available)")),
               conditionalPanel("input.scplottype == 'boxplot'",
                                selectInput("scplotvar2", 
                                            "Variables: (hold Ctrl to select
multiple variables)", 
                                            choices= "(Not Available)", 
                                            multiple= TRUE)),
               conditionalPanel("input.somtype == 'korresp'", 
                                selectInput("scplotrowcol", 
                                            "Plot rows or columns (when
relevant for the chart):",
                                            choices= list("columns"= "c", 
                                                          "rows"= "r"))),
               plotOutput("scplot")),
      
      #### Panel 'external information'
      #########################################################################
      tabPanel("Combine with external information",
               h3("Combine the map with external information"),
               p("In this panel you can combine the self-organizing map with 
variables not used for the training. To do so, you must first import an
additional file using the form below. The file must either contains the same
number of rows as the file used for training (in the same order), or a (square)
adjacency matrix  for 'graph' plots (the adjacency matrix has a dimension equal
to the number of rows ."),

               h4("Import file for additional variables"),
               fileInput('file2', 'Choose csv/text file'),
               checkboxInput('header2', ' Header?', TRUE),
               checkboxInput('rownames2', ' Row names?', FALSE),
               selectInput('sep2', 'Separator:', 
                           c("Comma", "Semicolon", "Tab", "Space"),
                           'Comma'),
               selectInput('quote2', 'Quote:',  
                           c("None", "Double Quote", "Single Quote"), 
                           'Double Quote'),
               selectInput('dec2', 'Decimal mark', c("Period", "Comma"),
                           'Period'),
               numericInput('nrow.preview.add', 
                            'Number of rows in the preview:', 20),
               numericInput('ncol.preview.add', 
                            'Number of columns in the preview:', 10),
               tableOutput("addview"),

               h4("Plot additional variables:"),
               conditionalPanel("input.somtype == 'korresp'",
                                p("Option not available for 'Korresp' type of
                                  SOM")),
               conditionalPanel("input.somtype != 'korresp'", 
                 selectInput("addplottype", "Type of plot:",
                             choices= c("pie", "color", "lines", "boxplot", 
                                        "barplot", "radar", "names", "words", 
                                        "graph")),
                 conditionalPanel("input.addplottype == 'pie' ||
                                   input.addplottype == 'color' ||
                                   input.addplottype == 'names'",
                                  selectInput("addplotvar", "Select variable:",
                                              choices= "(first import file)")),
                 conditionalPanel("input.addplottype != 'pie' &&
                                   input.addplottype != 'color' &&
                                   input.addplottype != 'names' && 
                                   input.addplottype != 'graph'",
                                  selectInput("addplotvar2", "Select variables:
                                               (hold Ctrl to select multiple
                                               variables)",
                                              choices= "(first import file)",
                                              multiple= TRUE)),
                 plotOutput("addplot"))
      ),
      tabPanel("Help",
               p(HTML("<h4>Topics:</h4> 
<a href= #somtype> Types of maps and data</a> <br />
<a href= #importdata> Data importation</a> <br />
<a href= #train> Training options</a> <br />
<a href= #plots> Types of plots</a> <br />
<a href= #superclasses> Grouping prototypes into superclasses</a> <br />
<a href= #externalinfo> Combine with external information</a> <br />")),
               
               p(HTML("<h3 id=somtype> Types of self-organizing maps</h3>")),
               p(HTML("Different types of data require different types of maps
to analyze them. <a href='http://sombrero.r-forge.r-project.org/'>SOMbrero</a>
offers three types of algorithms, all based on the on-line (as opposed to batch)
SOM:
<li><b>Numeric</b> is the standard self-organizing map, which uses numeric
variables only. For correct data importation the columns must contain the
variables and the rows the observations. <br /> It can be applied, for instance,
to the four first variables of the supplied <a href=
'http://owncloud.nathalievilla.org/apps/files_sharing/get.php?token=9d319bedf64098340ca9b094a527f70608e2d67c'
>iris dataset</a>.</li>
<li><b>Korresp</b> applies the self-organizing algorithm to contingency tables
between two factors.<br /> For instance, in the supplied <a href=
'http://owncloud.nathalievilla.org/apps/files_sharing/get.php?token=d824f59d7d934468de7f08b873f79d291aa129fd'
>dataset 'presidentielles 2002'</a>, which contains the results for the first
round of the 2002 French prensidential elections, columns represent presidential
candidates and rows represent the French districts called 'departements', so
that each cell contains the number of votes for a specific candidate in a
specific 'departement'.</li>
<li><b>Relational</b> is used for dissimilarity matrices, in which each cell
contains a measure of distance between two objects. For this method the data
must be a square numeric matrix, in which rows  and columns represent the same
observations. The matrix must be symetric, contains only positive entries with a
null diagonal.<br /> For instance, the supplied <a href='
http://owncloud.nathalievilla.org/apps/files_sharing/get.php?token=254fcf93f2e87ef71ac96eca03e661e659444921'
>dataset 'Les Miserables'</a> contains the shortest path lengths between
characters of Victor Hugo's novel <I>Les Misérables</I> in the co-appearance 
network provided <a href=
'http://www-personal.umich.edu/~mejn/netdata/lesmis.zip'>here</a>.")),
               
               p(HTML("<h3 id=importdata> Data importation</h3>")),
               p(HTML("Data must be imported as a table, in text or csv format
(columns are separated by specific symbols such as spaces or semicolons (option
'Separator'). Row names can be included in the data set for better post-analyses
of the data (some of the plots will use these names). Check at the bottom of the
'Import Data' panel to see if the data have been properly imported. If not,
change the file importation options.")),

               p(HTML("<h3 id=train> Training options</h3>")),
               p(HTML("The default options in the 'Self-Organize' panel are set
according to the type of map and the size of the dataset, but you can modify the
options at will:
<li><b>Input variables:</b> choose on which variables the map will be trained
(only available for 'numeric' som).</li>
<li><b>Map dimensions:</b> choose the number of rows (X) and columns (Y) of the
map, thus setting the number of prototypes. An advice value is to use the 
square root of one tenth the number of observations.</li>
<li><b>Affectation type:</b> type of affectation used during the training. 
Default type is 'standard', which corresponds to a hard affectation, and
the alternative is 'heskes', which is Heskes's soft affectation (see Heskes 
(1999) for further details).</li>
<li><b>Max. iterations:</b> the number of iterations of the training process.
Default is five times the number of observations.</li>
<li><b>Distance type:</b> type of distance used to determine which prototypes of
the map are neighbors. Default type is 'Letremy' that was originally proposed in
the <a href='http://samos.univ-paris1.fr/Programmes-bases-sur-l-algorithme'>SAS
programs</a> by <a href='http://samm.univ-paris1.fr/-Patrick-Letremy-'> Patrick
Letremy</a> but all methods for 'dist' are also available.</li>
<li><b>Radius type:</b> neighborhood type used to determine which prototypes of
the map are neighbors. Default type is 'letremy' as originally 
implemented in the <a 
href='http://samos.univ-paris1.fr/Programmes-bases-sur-l-algorithme'>SAS
programs</a> by <a href='http://samm.univ-paris1.fr/-Patrick-Letremy-'> Patrick
Letremy</a> (it combines square and star-like neighborhoods along the learning)
but a Gaussian neighborhood can also be computed.</li>
<li><b>Data scaling:</b> choose how the data must be scaled before training.
Scaling is used to ensure all variables have the same importance during
training, regardless of absolute magnitude. 'none' means no scaling,
'center' means variables are shifted to have 0 mean, 'unitvar' means variables
are centered and divided by their standard deviation to ensure they all have 
unit variance and '&chi;<sup>2</sup>' is used by the 'Korresp' algorithm and 
'cosine' is the dissimilarity transposition of the 'cosine' transformation 
performed on kernels.</li>
<li><b>Random seed:</b> Set the seed for the pseudorandom number generator used
during the training. Be sure to remember the seed used for a training if you
want to reproduce your results exactly, because running the algorithm with
different seeds will result in different maps. By default the seed is itself set
randomly when the interface is launched.</li>
<li><b>Scaling value for gradient descent:</b> This is the 'step size' parameter
used during training to determine in what proportion a winning prototype is
shifted towards an observation.</li>
<li><b>Number of intermediate backups:</b> Number of times during training a
backup of the intermediate states of the map is saved. This can be used to
monitor the progress of the training. If no backup (values 0 or 1) is saved,
the energy plot is not available.</li>
<li><b>Prototype initialization method:</b> choose how the prototypes of the map
are initialized at the beginning of training algorithm.<br /> If 'random' is
chosen, prototypes will be given random values in the range of the data. If
'obs', each prototype will be initialized to a random observation in the data.
If 'pca', prototypes are chosen along the first two PC of a PCA. Advised values
are 'random' for the 'Numeric' and the 'Korresp' algorithm and 'obs' for the 
'Relational' algorithm.")),
      
               p(HTML("<h3 id=plots> Types of plots</h3>")),
               p(HTML("Sombrero offers many different plots to analyze your
data's topology using the self-organizing map. <br />There are two main choices
of what to plot (in the plot and superclass panels): plotting <b> prototypes
</b> uses the values of the neurons' prototypes of the map, which are the
representative vectors of the clusters. Plotting <b> observations </b> uses the
values of the observations within each cluster.")),
               p(HTML("These are the standard types of plots:
<li><b>hitmap:</b> represents rectangles having areas proportional to the number
of observations per neuron.</li>
<li><b>color:</b> Neurons are filled with colors according to the prototype
value or the average value level of the observations for a chosen variable.
Yellow represents low values, and red high values.</li>
<li><b>3d:</b> similar to the ‘color’ plot, but in 3-dimensions, with x and y
the coordinates of the grid and z the value of the prototypes or observations
for the considered variable.</li>
<li><b>boxplot:</b>  plots boxplots for the observations in every neuron. This
plot can handle 5 variables at most.</li>
<li><b>lines:</b>  plots, for each neuron, the prototype value or the average
value level of the observations, with lines. Each point on the line represents a
variable. </li>
<li><b>barplot:</b>  similar to lines, here each variable value is represented
by a bar. </li>
<li><b>radar:</b>  similar to lines, here each variable value is represented by
a slice, and bigger slices represent higher values. </li>
<li><b>names:</b>  prints on the grid the names of the observations in the
neuron to which they belong. <strong>Warning!</strong> If the number of
observations or the size of the names is too large, some names may not be
represented.</li>
<li><b>poly.dist:</b> represents the distances between prototypes with polygons
plotted for each neuron. The closer the polygon point is to the border, the
closer the pairs of prototypes. The color used for filling the polygon shows the
number of observations in each neuron. A red polygon means a high number of
observations, a white polygon means there are no observations.</li>
<li><b>smooth.dist:</b>  depicts the average distance between a prototype and
its neighbors using smooth color changes, on a map where x and y are the
coordinates of the prototypes on the grid.</li>
<li><b>umatrix:</b>  is another way of plotting distances between prototypes.
The grid is plotted and filled colors according to the mean distance between the
current neuronand the neighboring neurons. Red indicates proximity.</li>
<li><b>grid.dist:</b> plots all distances on a two-dimensional map. The number
of points on this picture is equal to: <br>
<code>number_of_neurons * (number_of_neurons-1) / 2</code><br>The x axis
corresponds to the prototype distances, the y axis depicts the grid distances.
</li>
<li><b>MDS:</b> plots the number of the neurons according to a Multi Dimensional
Scaling (MDS) projection on a two-dimensional space.</li>")),
               p(HTML("The plot options offered in the 'superclasses' and
'combine with  external information' panels are mostly the same as the ones
listed above, but some are specific:
<li><b>dendrogram:</b>  (in Superclasses panel) plots the dendrogram of the
hierarchical clustering applied to the prototypes, along with the scree plot
which shows the proportion of unexplained variance for incremental numbers of
superclasses. These are helpful in determining the optimal number of 
superclasses.</li>
<li><b>dendro3d:</b>  (in Superclasses panel) similar to 'dendrogram', but in
three dimensions and without the scree plot.</li>
<li><b>pie:</b>  (in Combine panel) requires the selected variable to be a
categorical variable, and plots one pie for each neuron, corresponding to the
values of this variable.</li>
<li><b>words:</b>  (in Combine panel) needs the external data to be a
contingency table: names of the columns will be used as words and printed on the
map with sizes proportional to their frequency in the neuron.</li>
<li><b>graph:</b>  (in Combine panel) needs the external data to be the
adjacency matrix of a graph. According to the existing edges in the graph and to
the clustering obtained with the SOM algorithm, a clustered graph is built in
which vertices represent neurons and edge are weighted by the number of edges in
the given graph between the vertices affected to the corresponding neurons. This
plot can be tested with the supplied dataset <a href=
'http://owncloud.nathalievilla.org/apps/files_sharing/get.php?token=254fcf93f2e87ef71ac96eca03e661e659444921'
>Les Miserables</a> that
corresponds to the graph those adjacency table is provided at <a href=
'http://owncloud.nathalievilla.org/apps/files_sharing/get.php?token=2f908aa535668b3b894b64f1679986e4aa4476b7'
>this link</a>.
</li>"
      )),
               p(HTML("The <b>show cluster names</b> option in the 'Plot map'
panel can be selected to show the names of the neurons on the map.")),
               p(HTML("The <b>energy</b> option in the 'Plot map' panel is used
to plot the energy levels of the intermediate backups recorded during training.
This is helpful in determining whether the algorithm did converge. (This option
only works if a 'Number of intermediate backups' larger than 2 is chosen in the
'Self-Organize' panel.)")),
    
               p(HTML("<h3 id=superclasses> Grouping prototypes into
                      Superclasses</h3>")),
               p(HTML("Use the options on the dedicated panel to group the
prototypes of a trained map into a determined number of superclasses, using
hierarchical clustering. The 'dendrogram' plot can help you to choose a relevant
number of superclasses (or equivalently a relevant cutting height in the
dendrogram).")),
    
               p(HTML("<h3 id=externalinfo> Combine with external
                      information</h3>")),
               p(HTML("Plot external data on the trained map on the dedicated
panel. The external data importation process is similar to the one described in
the <a href='#importdata'>'Data importation'</a> section, and the available
plots are described in the <a href='#plots'>'types of plots'</a> section.<br />
Note that this is the only panel in which factors can be plotted on the
self-organizing map. For instance, if the map is trained on the first four
(numeric) variables of the supplied <a href=
'http://owncloud.nathalievilla.org/apps/files_sharing/get.php?token=9d319bedf64098340ca9b094a527f70608e2d67c'
>iris dataset</a>, you can import the dataset again as external information and
plot the iris species on the map."))
      )
    )
  )
))
