# for comparison to excel version

comment <- "This run is for comparison with the Excel version"


# inputfile
# just like in R, all comments start with #.
# order of parameters does not matter
# expressions allowed

# Comments must be Latex-friendly; ie
#       encapsulate math statements with $$.
#       don't use &

# labels
startage <- 1             ; c.startage <-     ""  
endage <- 60              ; c.endage <-       ""     
initialyear <- 1989       ; c.initialyear <-  "start in year"        
initialmonth <- 1         ; c.initialmonth <- "and january"           
yearplanted <- 1988       ; c.yearplanted <-  "for calculation of stand age"
monthplanted <- 1         ; c.monthplanted <- " ... "          
sitename <- " "           ; c.sitename <-     "a label for the site"       
speciesname  <- " "       ; c.speciesname  <- "and for the species"


# initial values:
aswi <- 999            ;  c.aswi <-    "initial available soil water (mm?)" 
wfi <- 7               ;  c.wfi <-     "initial foliage mass"
wri <- 2               ;  c.wri <-     "initial root mass"
wsi <- 5               ;  c.wsi <-     "initial stem mass"
stemnoi <- 1111        ;  c.stemnoi <- "initial stem number"


# site parameters
lat <- -26             ;  c.lat <-       "latitude (NH positive)"
maxasw <- 200          ;  c.maxasw <-    "maximum available soil water (SWHC)"
minasw <- 0            ;  c.minasw <-    "minimum available soil water"   
fr <- 1                ;  c.fr <-        "fertility rating"
soilclass <- "CL"      ;  c.soilclass <- "soil type"
swconst <- NA          ;  c.swconst <-   "soil water release parameters: leave NA when soilclass given."
swpower <- NA          ;  c.swpower <-   ""        



# Parameters
maxage <- 50 ;   c.maxage <- "Determines rate of physiological decline of forest"   
sla0 <- 4    ;   c.sla0 <- "Specific leaf area at age 0 ($m^2$/kg)"
sla1 <- 4    ;   c.sla1 <- "Specific leaf area for mature trees ($m^2$/kg)" 
tsla <- 2.5  ;   c.tsla <- "stand age (years) for SLA = (SLA0+SLA1)/2"


# fullcanage deleted: was not used in 3PG latest version


k <- 0.5          ; c.k <-    "Radiation extinction coefficient"
gammafx <- 0.03   ; c.gammafx <- "Coefficients in monthly litterfall rate"
gammaf0 <- 0.03   ; c.gammaf0 <- ""                       
tgammaf <- 24     ; c.tgammaf <- ""                        
rttover <- 0.015  ; c.rttover <- "Root turnover rate per month"
swconst0 <- 0.7   ; c.swconst0 <- "SW constants"
# are 0.7 for sand, 0.6 for sandy-loam, 0.5 for clay-loam, 0.4 for clay
swpower0 <- 9     ; c.swpower0 <- "Powers in the eqn for SW modifiers"
# are 9 for sand,  7 for sandy-loam, 5 for clay-loam and 3 for clay
maxintcptn <- 0.15; c.maxintcptn <- "Max proportion of rainfall intercepted by canopy"
laimaxintcptn <- 0; c.laimaxintcptn <- "LAI required for maximum rainfall interception"
maxcond <- 0.02   ; c.maxcond <- "Maximum canopy conductance (gc, m/s)"  
laigcx <- 3.33    ; c.laigcx <- "LAI required for maximum canopy conductance"


blcond <- 0.2     ; c.blcond <- "Canopy boundary layer conductance, assumed constant"
coeffcond <- 0.05 ; c.coeffcond <- "Determines response of canopy conductance to VPD"
y <- 0.47         ; c.y <- "Assimilate use efficiency"
tmax <- 32        ; c.tmax <- "Critical biological temperatures: max, min"
tmin  <- 2        ; c.tmin  <- ""                     
topt <- 20        ; c.topt <-  "and optimum. Reset if necessary/appropriate"
kf <- 1           ; c.kf <- "Number of days production lost per frost days"
pfs2 <- 1         ; c.pfs2 <- "Foliage:stem partitioning ratios for D = 2cm"
pfs20 <- 0.15     ; c.pfs20 <- "and D = 20cm"
stemconst <- 0.095; c.stemconst <- "Stem allometric parameters"
stempower <- 2.4  ; c.stempower <-  ""                         
prx <- 0.8        ; c.prx <-   "Maximum fraction of NPP to roots" 
prn <- 0.25       ; c.prn <-   "Minimum fraction of NPP to roots"
m0 <- 0           ; c.m0 <-    "Value of m when FR = 0"
fn0 <- 1          ; c.fn0 <-   "Value of fN when FR = 0"
alpha <- 0.055    ; c.alpha <- "Canopy quantum efficiency"
wsx1000 <- 100    ; c.wsx1000 <- "Max tree stem mass (kg) likely in mature stands  of 1000 trees/ha" 
thinpower <- 3/2  ; c.thinpower <- "exponent in the self-thinning routine."
fracbb0 <- 0.15   ; c.fracbb0 <-  "branch and bark fraction at age 0 ($m^2$/kg)"
fracbb1 <- 0.15   ; c.fracbb1 <- "branch and bark fraction for mature trees ($m^2$/kg)"
tbb <- 1.5        ; c.tbb <- "stand age (years) for fracBB = (fracBB0+fracBB1)/2"
density <- 0.5    ; c.density <- "basic density (t/m3)"


mf <- 0           ; c.mf <- "Leaf mortality fraction"
mr <- 11 / 54     ; c.mr <- "Root mortality fraction"
ms <- 11 / 54     ; c.ms <- "Stem mortality fraction"


rage <- 0.95      ; c.rage <- "parameters in the age modifier."
nage <- 4         ; c.nage <- ""

