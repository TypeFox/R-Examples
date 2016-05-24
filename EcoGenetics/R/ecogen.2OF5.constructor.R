################################################
#### ECOGEN CONSTRUCTOR
################################################

#' Creating a new ecogen object
#' @param XY Data frame with m columns (coordinates) and n rows (individuals).
#' @param P Data frame with n rows (individuals), and phenotypic data in columns.
#' @param G Data of class: "data.frame", with individuals in rows and genotypic data 
#' in columns (loci). The ploidy and the type (codominant, dominant) of the data, 
#' must be passed with the arguments "ploidy" and "type". Missing data is coded as NA.
#' Dominant data must be coded with binary values (0 for absence - 1 for presence).
#' @param E Data frame with n rows (individuals), and environmental data 
#' in columns.
#' @param S Data frame with n rows (individuals), and groups (factors) in columns.
#' The program converts non factor data into factor.
#' @param C data frame with n rows (individuals), and custom variables in columns.
#' @param G.processed If TRUE, the slot G will include a processed data frame (
#' removed non informative loci (the data non available for all the individuals),
#' removed non polymorphic loci (for dominant data) and ordered alleles in ascending
#' order. 
#' @param order.G Genotypes must be ordered in G slot? (codominant data) 
#' Default FALSE.
#' @param type Marker type: "codominant" or "dominant".
#' @param ploidy Ploidy of the G data frame. Default ploidy = 2.
#' @param sep Character separating alleles (codominant data). 
#' Default option is no character separating alleles. 
#' @param ncod Number of characters coding each allele (codominant data).
#' @param missing Missing data treatment ("0", "NA", or "MEAN") for the A
#' slot. Missing elements are set to 0 in the default option. missing elements
#' are recoded as "NA" or the mean allelic frequency across individuals in "NA" 
#' and "MEAN" options, respectively. 
#' @param NA.char Character simbolizing missing data in the input. Default is "NA".
#' @param poly.level Polymorphism threshold in percentage (0 - 100), 
#' for remotion of non polymorphic loci (for dominant data). Default is 5 (5\%).
#' @param rm.empty.ind Remotion of noninformative individuals (row of "NAs").
#' Default if FALSE.
#' 
#' @details This is a generic function for creating an ecogen object.
#' In default option, the missing data input value is "NA", but any missing 
#' data character can be passed with the option NA.char. 
#' The output in the slot G will have coded the missing data as NA. 
#' 
#' 
#' \strong{ACCESS TO THE SLOTS. MODIFICATION OF ECOGEN OBJECTS}
#' 
#' The content of the slots can be extracted with the corresponding accesors
#' ecoslot.XY, ecoslot.P, ecoslot.G, ecoslot.A, ecoslot.E, ecoslot.C and ecoslot.OUT. 
#' Data can be assigned individually to the slots, also with the corresponding accessors.
#' The correct use of ecogen objects requires the implementation of accessors, 
#' as they ensure the data check and pre-processing. The use of accessors enables to
#' modify or fill the slots of ecogen objects, without the need of creating a new
#' object each time. See \emph{help("EcoGenetics accessors")} for a detailed description 
#' and examples about ecogen accessors. 
#' 
#' \strong{OTHER SLOT ACCESS METHODS FOR ECOGEN OBJECTS}
#' 
#' The use of brackets is defined for ecogen objects:
#' 
#'  - Single bracket: the single bracket ("[") is used to subset all the ecogen
#'  data frames (P, G, E, S, A and C) by row, at once. The notation for an object
#'  is eco[from:to], where eco is any ecogen object, and from: to is the row
#'  subset. For example: eco[1:10] , subsets the object eco from row 1 to row 10, 
#'  for all the data frames at once.
#'  
#' - Double square brackets: the double square brackets are symbolic abbreviations 
#' of the accessors (i.e., it is a call to the corresponding accessor). 
#' The usage is: eco[["X"]], where X is any slot of interest: eco[["P"]], 
#' eco[["G"]], eco[["A"]], eco[["E"]], eco[["S"]], eco[["C"]] and eco[["OUT"]].
#' Double square brackets can be used in get/set mode. 
#' See Examples below and in help("EcoGenetics accessors").
#' 
#' 
#' \bold{ABOUT THE CONSTRUCTION OF NEW ECOGEN OBJECTS}
#' 
#' The construction of a new ecogen object can be made in two different ways. 
#' First, a new object can be created, incorporating all the information 
#' at once. Second, the data can be added in each slot, using the corresponding
#' accessor / "[[". The accessor/double square brackets methods allow temporal modification 
#' of any created ecogen object and ensures its modularity. 
#' These methods are not only functions to get/assign values to the slots, 
#' they provide a basic pre-processing of the data during assignment, 
#' generating a set of coherent information.
#' 
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @examples
#' \dontrun{
#' 
#' #Example with G data of class "data.frame", corresponding to
#' #microsatellites of a diploid organism:
#' data(eco.test)
#' eco <- ecogen(XY = coordinates, P = phenotype, G = genotype,
#' E = environment, S = structure)
#' 
#' #Example with G data of class "data.frame", corresponding to a
#' #presence - absence molecular marker:
#' dat <- sample(c(0,1),100,rep = TRUE)
#' dat <- data.frame(matrix(dat,10,10))
#' eco <- ecogen(G = dat, type = "dominant")
#' 
#' 
#' # DINAMIC ASSIGNMENT WITH ACCESSORS AND "[["
#' 
#' eco <- ecogen(XY = coordinates, P = phenotype)
#' eco
#' 
#' ecoslot.G(eco, order.G = TRUE) <- genotype
#' 
#' # this is identical to
#' eco[["G", order.G=TRUE]] <- genotype
#' 
#' ecoslot.E(eco) <- environment
#' 
#' # this is identical to
#' eco[["E"]] <- environment
#' 
#' #----------------------------------------------------------
#' # See additional examples in help("EcoGenetics accessors")
#' #----------------------------------------------------------
#' 
#' # Storing data in the slot OUT
#' 
#'  singers <- c("carlos_gardel", "billie_holiday")
#'  
#' ecoslot.OUT(eco) <- singers
#'  
#' # Storing several data
#'
#' golden.number <- (sqrt(5) + 1) / 2
#' ecoslot.OUT(eco) <- list(singers, golden.number)    # several objects must be passed as a list
#' 
#' # this is identical to:
#' 
#' eco[["OUT"]] <- list(singers, golden.number)
#' 
#' }
#' 
#' @export ecogen


setGeneric("ecogen",    		 
           function(XY = data.frame(),
                    P = data.frame(),
                    G = data.frame(), 
                    E = data.frame(),
                    S = data.frame(),
                    C = data.frame(),
                    G.processed = TRUE,
                    order.G = FALSE,
                    type = c("codominant", "dominant"),
                    ploidy = 2,
                    sep, 
                    ncod = NULL,
                    missing = c("0", "NA", "MEAN"),
                    NA.char = "NA", 
                    poly.level = 5,
                    rm.empty.ind = FALSE) {				
             
             # general configuration
             type <- tolower(type)
             type <- match.arg(type)
             missing <- toupper(as.character(missing))
             missing <- match.arg(missing)
             if(missing(sep)) {
               sep <- ""
             }
             
             # creating a new ecogen object
             Object <- new("ecogen")
             
             # G configuration-------------------------------------------------#
            
             if(any(dim(G) == 0)) { # empty G
               Object@G <- data.frame()
               Object@A <- data.frame()
               Object@INT <- new("int.gendata")
               
               
             } else { # non empty G
               
               ## coherence between data ploidy and ncod is checked for int.df2genind
               
               tempo <- int.df2genind(G, 
                                      sep = sep, 
                                      ncod =  ncod,
                                      NA.char = NA.char, 
                                      ploidy = ploidy, 
                                      type = type,
                                      missing = missing,
                                      rm.empty.ind = rm.empty.ind,
                                      poly.level = poly.level)
               
               # unfolding tempo
               Object@A <- as.data.frame(tempo@tab)
               Object@INT <- int.genind2gendata(tempo)
               
               ncod <- tempo@ncod
               ploidy <- tempo@ploidy
               
               # G processed case ~-~-~-~-~~-~-~-~-~
               if(G.processed) {
                 tmp <- int.genind2df(tempo)
                 # order data
                 if(order.G && type == "codominant") {
                   tmp <- aue.sort(tmp, 
                            ncod = ncod, ploidy = ploidy, 
                            chk.plocod = FALSE)
                 } 
                 
                 # G processed data frame 
                 G <- as.data.frame(tmp, stringsAsFactors = FALSE)
                 
                 # G changes messages 
                 
                 if(order.G && type == "codominant") {
                   message("Note: ordered genotypes in slot G")
                 }
               } 
               # END G processed case ~-~-~-~-~~-~-~-~-~
               
               # fill now the G slot
               Object@G <- G
             }
             

             
             # fill the other slots--------------------------------------------
             
             Object@XY <- as.data.frame(XY)
             Object@P <- as.data.frame(P)
             Object@E <- as.data.frame(E)
             
             # all S columns as factors
             S <- as.data.frame(S)
             if(dim(S)[1] != 0) {
               for(i in 1:(ncol(S))) {
                 S[, i] <- factor(S[, i])
               }
             }
             Object@S <- S
             Object@C <- as.data.frame(C)
             
             return(Object)
             
           })

