calc.coal <-
function (gTree, sTree, stNames, assoc) {
  # Descr:  calculating the Ranalla&Yang index
  # Deps:   (various)     
  # I/p:    gTree
  #         sTree
  #         stNames
  #         assoc

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.coal",fg="red"),sep="")
  }

## 1. Loading of gene trees and gene tree attributes
  # LEGACY: gTree.treeshape = apTreeshape::as.treeshape(gTree)
  # LEGACY: gTree.treeshape = apTreeshape::as.treeshape.phylo(gTree)
  # LEGACY: gTreeTaxa = gTree.treeshape$names
  gTreeTaxa = gTree$tip.label
  gTree.string = ape::write.tree(gTree)

## 2. Loading of species trees and species tree attributes
  sTreeTaxa = stNames
  sTree.string = sTree

## 3. Generating "species structure matrix"
  sp = unique(assoc[,1])
  taxa = assoc[,2]
  spStrMtrx = matrix(0, nrow=length(sp), ncol=length(taxa))
  rownames(spStrMtrx) = sp
  colnames(spStrMtrx) = taxa
  for (s in sp) {
    spStrMtrx[s, assoc[,2][which(assoc[,1]==s)]] = 1
  }

## 4. Enclosing data in a list
  data = list()
  data$gTree.string = gTree.string
  data$sTree.string = sTree.string
  data$gTreeTaxa = gTreeTaxa
  data$sTreeTaxa = sTreeTaxa
  data$spStrMtrx = spStrMtrx

  # DEBUGLINES:
  #cat("\ndata\n"); print(data)

### 5. DEBUG mode
#  logdata = list(list(COAL=data))
#  names(logdata) = get("P2C2M_flg_repID", envir=P2C2M_globalVars)
#  loghelpers.dbg(logdata, "DescrStatsRawInput", "INPUT OF 'COAL'")

## 6. Calculating descriptive statistic
  ray = phybase::loglikeSP(data$gTree.string, data$sTree.string,
                           data$gTreeTaxa, data$sTreeTaxa, 
                           data$spStrMtrx, strict=F)

  return(ray)
}
