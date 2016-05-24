#  parseBoolean.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#  
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


################################################
# Function: .parseBoolean
#
#
# The function .parseBoolean() is inspired by the function
# parseBoolean() contained in the COBRA Toolbox.
# The algorithm is the same.
#
# 2015-06-17 CJF: added handling for emtpy gprRule like "( )"


.parseBoolean <- function(gprRule, tokens = "()&|~") {

#.parseBoolean <- function(gprRule,
#                          tokens = "()&|~",
#                          allowedElementChars = "[A-Za-z0-9_\\.\\-]") {

  # quit, if there is no gene association
  if ( is.na(gprRule) || (gprRule == "") ) {
      return(list(gene = "", rule = ""))
  }
  
  if( grepl("\\s*\\(\\s*\\)\\s*", gprRule) ){
  	warning("found empty expression rule: '( )'. check if this is intended.")
  	return(list(gene = "", rule = ""))
  }

  # quit, if there is no gene association
#  if ( (gprRule == "UNKNOWN") || (gprRule == "SPONTANEOUS") ) {
#      return(list(gene = "", rule = ""))
#  }

  #print(c("rule: ", gpr))
  
#  str <- gpr
  gpr <- gsub("and ", "& ", gprRule, ignore.case = TRUE)
  gpr <- gsub("or ",  "| ", gpr, ignore.case = TRUE)
  gpr <- gsub("not ", "~ ", gpr, ignore.case = TRUE)
  gpr <- gsub("[", "", gpr, fixed = TRUE)
  gpr <- gsub("]", "", gpr, fixed = TRUE)
  rule <- gpr

  #print(gpr)


#  endStr <- str1
  #tmpRule <- matrix()

  
  # split the rule into the gene names 
  genes_tmp <- strsplit(gpr, paste("[", tokens, "]", sep = ""))

  # remove trailing and leading whitespaces
  genes_tmp <- gsub("(^\\s+)|(\\s+$)", "", genes_tmp[[1]], perl = TRUE)

  # remove empty entries in genes_tmp
  not_empty <- which(nchar(genes_tmp) > 0 )
  genes     <- genes_tmp[not_empty]

  # number of entries
  num_genes <- length(genes)

  # a unique vector with all genes
  gene_uniq <- unique(genes)

  newTok    <- match(genes, gene_uniq)
  newTok    <- sapply(newTok, function(x) paste("x(", x, ")", sep = ""))

#  rule <- 
  
  #bla <- rbind(genes, newTok)
  #print("gedoens")
  #rule <- apply(bla, 2, function(x) gsub(x[1], x[2], rule, fixed = TRUE))
  #apply(bla, 1, function(x) print(x[1]))

  for (i in 1:num_genes) {

      rule <- sub(genes[i], newTok[i], rule, fixed = TRUE)
      #start <- gregexpr(genes[i], gpr, fixed  = TRUE)
      #start <- start[[1]]

      #end <- start + attr(start, "match.length")

#print(gpr)
#print(start)
#print(end)
      #for (j in 1:length(start)) {

      #    print(newTok[i])
      #    substr(rule, start[j], (end[j])) <- paste(newTok[i], "")

      #}
      
     # if ()
 
  }

  
  #print("bla")
  #print(genes)
  #print(genes)
  #print(newTok)
  #print(rule)

#  cnt <- 0

#  replace_v <- logical(0)



  
#  for (i in 1:num_genes) {

#          current_tok <- genes[i]
    
          # propably we do not need this here, because we did it before
          #current_tok <- gsub("[", "", current_tok, fixed = TRUE)
          #current_tok <- gsub("]", "", current_tok, fixed = TRUE)

#          start <- gregexpr(current_tok, rule, fixed  = TRUE)
#          start <- start[[1]]

#          end <- start + attr(start, "match.length")
           
#          print(c(start, end))
#          print(c("bla", current_tok, genes))

#          replace <- logical(length(start))
          
#          for (j in 1:length(start)) {
              # only one token
#              if (num_genes == 1) {
#                  replace[i] <- TRUE
#                  next
#              }
#              else {
                  # token at the beginning
#                  if (match())
      
#              }

#          }





#  }




  
#  for (i in 1:length(genes[[1]])) {
#      #current_tok <- gsub("(^\\s+)|(\\s+$)", "", genes[[1]][i], perl = TRUE)
#      current_tok <- genes[[1]][i]
#  for (i in 1:length(genes)) {
      #current_tok <- gsub("(^\\s+)|(\\s+$)", "", genes[[1]][i], perl = TRUE)
#      current_tok <- genes[i]
#      if (nchar(current_tok) > 0) {
          #current_tok <- tokens[[1]][i]

#          if(!match(current_tok, gene, nomatch = 0)) {
#              cnt <- cnt + 1
              #gene <- c(gene, current_tok)
#              gene[length(gene)+1] <- current_tok
#          }

          # the replacement token
#          newTok <- paste("x(", cnt, ")", sep = "")
          #print(newTok)

          #print(c(current_tok, rule))
          #start <- gregexpr(paste("\\Q", current_tok, "\\E", sep = ""), rule, perl = TRUE)
#      }

#  }

  
  #tempe <- gregexpr(tokens, endStr)
  #temp <- substr(endStr, tokens)
#print(temp)
  
#  gene = "bla"
#  rule = "blubber"
  return(list(gene = gene_uniq, rule = rule))

}

