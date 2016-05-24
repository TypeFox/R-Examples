
##==============================================================================
## Implementation of Network indices, as in
## Latham, LG II 2006 Ecol. Modeling 192: 586-600
##
## Implemented: Julius Kones      - University Nairobi
##              Karline Soetaert  - Netherlands Institute of Ecology
##
## Two local functions:
##     InternalNetwork
##     Diet
##==============================================================================

##==============================================================================
InternalNetwork <- function (Tij  ,          # to-from
                             Import,         # flow from external (colNr Tij)
                             Export)         # flow to external (colNr Tij)

{                      

##------------------------------------------------------------------------
## Tij[i,j] is a matrix with Tij[i,j]  flow from j to i
## note: component position in rows and columns must be the same - not checked
##------------------------------------------------------------------------

 if (is.character(Import))
   import <- which(colnames(Tij)%in%Import) else
   import <- Import
 if (length(import) != length(Import))
   stop("Import not recognized")
 if (is.character(Export))
   export <- which(rownames(Tij)%in%Export) else
   export <- Export
 if (length(import) != length(Import))
   stop("Import not recognized")
 if (length(export) != length(Export))
   stop("Export not recognized")

##
## CHECK THE INPUT
##

  # Flow or Tij should be inputted
  if (is.null(Tij))
    stop ("cannot calculate indices - Flow or Tij should be inputted")

#_________________________________________________________________________________
# NUMBER OF COMPARTMENTS, export, import, internal flows,..
#_________________________________________________________________________________

   # Size of the matrices; without the externals, the matrix has to be square

  ncomp     <- ncol(Tij)-length(import)
  if (ncomp != nrow(Tij)-length(export))
    stop ("cannot calculate indices - internal flow input matrix not square ")
 
#_________________________________________________________________________________
# ARRAYS DECLARATION         
#_________________________________________________________________________________

    # indices to elements of T that are internal 
  iN  <- setdiff(1:nrow(Tij),export)  # internal rows    of Tij
  jN  <- setdiff(1:ncol(Tij),import)  # internal columns of Tij


    # Total internal flows, externals removed.
  Tint        <- Tij
  if (! is.null(export))
    Tint <- Tint[-export,]
  if (! is.null(import))
    Tint <- Tint[,-import]

    # Total flows, including flow to/from externals
  FlowFrom   <- colSums(Tij)
  FlowTo     <- rowSums(Tij)
  FlowFromC  <- FlowFrom[jN]    # just the total from internal compartments
  FlowToC    <- FlowTo  [iN]
 
  return(list(Tint=Tint,iN=iN,jN=jN,
              import=import,export=export,
              FlowFrom=FlowFrom,
              FlowTo  = FlowTo,
              FlowFromC=FlowFromC,
              FlowToC  =FlowToC))

} # END InternalNetwork

##==============================================================================
##
## Internal function: estimates the diet composition
##
##==============================================================================
Diet <- function (Tint,                   # Calculates diet composition
                  Dead=NULL,              # index from Dead to Tint
                  iN=1:nrow(Tint))      

{

## p matrix contains the diet composition of predator i
  IntFlowTo   <- rowSums(Tint)  # Total food ingested

  p           <- Tint

  for (i in 1:ncol(Tint))
    p[i,] <- Tint[i,]/IntFlowTo[i]

  p[is.na(p)] <- 0

## take into account dead matter; Dead refers to column/row in Tij
## N$iN maps from Tint to Tij

  if (! is.null(Dead))
    p[which(iN%in% Dead),]<-0
  return(p)

} # END Diet   

