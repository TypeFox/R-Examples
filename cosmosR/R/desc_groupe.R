#' Comparatif par groupe
#' 
#' Produit un tableau comparatif par groupe
#' 
#' Permet de produire un tableau comparatif des variables contenues dans la table passee en parametre selon les modalites d'une d'entre elles.
#' Si les labels et formats sont definis et charges ils seront utilises pour peupler le tableau.
#' Le fichier de sortie est place dans ../HTML Output
#' @encoding UTF-8
#' @param table Table a utiliser
#' @param groupe Nom de la variable qualitative a utiliser pour la comparaison
#' @param param Vecteur de noms de variables considerees comme parametriques
#' @param html Nom du fichier html, par defaut "desc_groupe_nomdelavariable.html"
#' @param titre Titre du tableau, par defaut "Comparaison selon nom_de_la_variable"
#' @param variables Vecteur de noms de variables a comparer, par defaut toutes les variables contenues dans la table moins celle servant de comparateur
#' @param variables_neg Vecteur de noms de variables a exclure de la comparaison
#' @param note Note de bas de page, par defaut vide
#' @param nbdec Nombre de decimales apres la virgule, par defaut 1
#' @param pourcent Pourcentages pour les variables qualitatives, en colonnes ("col") ou en lignes ("row"), par defaut sur le total
#' @examples
#' \dontrun{
#' Ma_table <- charger("donnees.xls")
#' 
#' desc_groupe(Ma_table, "sexe")
#' 
#' para <- diagnostic(Ma_table)
#' desc_groupe(Ma_table, "sexe", param = para, titre="Comparatif selon le sexe", pourcent="row")
#' }
#' @export
desc_groupe <- function(table, groupe, param = character(0), html=NULL, titre=NULL, variables=NULL, variables_neg=NULL, note=NULL, nbdec=1, pourcent="total")
{
  if (missing(table))
  {
    warning("Pas de table donnee !")
    return
  }
  
  if (missing(groupe) | (!groupe %in% names(table)))
  {
    warning("Pas de comparateur !")
    return
  }
  
  modulo <- table[[groupe]]
  
  if (nlevels(modulo)<2)
  {
    warning("La variable ", groupe, " n'a qu'un niveau !")
    return
  }
  
  etiq <- label(modulo)
  if (etiq == "modulo") etiq <- var
  
  # Creation du titre
  if (is.null(titre))
    titre<-paste("Comparaison selon", etiq)
  
  # Creation de la table temporaire contenant les variables demandees
  if (!is.null(variables)) 
  {
    badvar <- variables[!variables %in% names(table)]
    if (length(badvar) != 0)
      warning("Les variables suivantes n'existent pas dans la table : ",paste(badvar,collapse=", "))
    table <- table[variables[variables %in% names(table)]]
  }
  if (!is.null(variables_neg))
  {
    badvar <- variables_neg[!variables_neg %in% names(table)]
    if (length(badvar) != 0)
      warning("Les variables_neg suivantes n'existent pas dans la table : ",paste(badvar,collapse=", "))
    table[variables_neg]<-list(NULL)
  }
  
  table[groupe]<-NULL
  
  # Creation du nom de fichier html
  if (is.null(html))
    html <- paste0("desc_groupe_",groupe)
  
  # Header de la table
  HTMLInit(file=paste0("../HTML\ Output/",html,".html"), title=titre, CSSfile="desc.css")
  HTML("<div class='desc'>")
  inc()
  HTML("<table class='desc'>")
  inc()
  HTML("<caption>", titre, "</caption>")
  HTML("<thead>")
  inc()
  
  HTML("<tr><th></th>", paste0("<th colspan='2'>",levels(modulo),"</th>", collapse=""), "<th>p</th></tr>")
  
  HTML("<tr><td></td>", sep="")
  for (i in levels(modulo))
  {
    HTML("<td colspan='2'>N=", length(na.omit(modulo[modulo==i])), "</td>", sep="")
  }
  HTML("<td></td></tr>")

  HTML("<tr><td></td>", paste0(rep(c("<td>N</td>","<td>%</td>"),nlevels(modulo)),collapse=""), "<td></td></tr>")

  dec()
  HTML("</thead>")
  
  # Footer de la table
  if (!is.null(note))
  {
    HTML("<tfoot>")
    inc()
    HTML("<tr><td colspan='", nlevels(modulo)*2+2, "'>", note,"</td></tr>")
    dec()
    HTML("</tfoot>")
  }
  
  # Corps de la table
  HTML("<tbody>")
  inc()
  for (var in names(table))
  {
    # Creation du label de la variable
    etiq <- label(table[[var]])
    if (etiq == "table[[var]]") etiq <- var
    
    # Variable quantitative
    if (is.numeric(table[[var]]))
    {
      if (var %in% param)
      {
        if (nlevels(modulo) == 2)
        {
          p <- tryCatch(t.test(table[[var]] ~ modulo)$p.value, error=function(e){erreur(e,var,"dans le calcul du test t")})
          test <- "t"
        }
        else if (nlevels(modulo) > 2)
        {
          p <- tryCatch(summary(aov(table[[var]] ~ modulo))[[1]][["Pr(>F)"]][1], error=function(e){erreur(e,var,"dans le calcul de l'ANOVA")})
          test <- "ANOVA"
        }
      }
      else
      {
        if (nlevels(modulo) == 2)
        {
          p <- tryCatch(wilcox.test(table[[var]] ~ modulo)$p.value, error=function(e){erreur(e,var,"dans le calcul du test de Mann-Whitney")})
          test <- "M-W"
        }
        else if (nlevels(modulo) > 2)
        {
          p <- tryCatch(kruskal.test(table[[var]] ~ modulo)$p.value, error=function(e){erreur(e,var,"dans le calcul du test de Kruskal-Wallis")})
          test <- "K-W"
        }
      }
      
      HTML("<tr><td class='var'>", etiq, " (moy &plusmn;SD)</td>", sep="")
      
      for (level in levels(modulo))
        HTML("<td>",format(mean(table[[var]][modulo==level], na.rm=T), digits=nbdec,nsmall=nbdec), "</td><td>", format(sd(table[[var]][modulo==level],na.rm=T), digits=nbdec,nsmall=nbdec), "</td>", sep="")
      
      HTML("<td>", format(p,digits=nbdec,nsmall=nbdec), sep="")
      if (p<.05)
        HTML("*", sep="")
      HTML(" (",test,")</td>", sep="")
      
      HTML("</tr>")
    }
    # Variable qualitative
    else if (is.factor(table[[var]]))
    {
      test <- "X"
      p <- tryCatch(
        {
          chisq.test(table[[var]],modulo)$p.value
        },
        warning=function(w)
        {
          test <<- "f"
          fisher.test(table[[var]],modulo,workspace=2e6)$p.value
        },
        error=function(e)
        {
          "NA"
        })
      
      HTML("<tr><td class='var' colspan='", nlevels(modulo)*2+1, "'>", etiq, "</td>", sep="")
      HTML("<td>", format(p,digits=nbdec,nsmall=nbdec), sep="")
      if (p<.05)
        HTML("*", sep="")
      HTML(" (",test,")</td>", sep="")
      
      # Levels de la variable qualitative
      for (level in levels(table[[var]]))
      {
        HTML("<tr><td class='level'>", level, "</td>",sep="")
        for (level_m in levels(modulo))
        {
          HTML("<td>", table(table[[var]],modulo)[level,level_m], "</td>", sep="")
          switch(pourcent,
                 row = HTML("<td>", format(100*table(table[[var]],modulo)[level,level_m]/rowSums(table(table[[var]],modulo))[level], digits=nbdec, nsmall=nbdec), "</td>", sep=""),
                 col = HTML("<td>", format(100*table(table[[var]],modulo)[level,level_m]/colSums(table(table[[var]],modulo))[level_m], digits=nbdec, nsmall=nbdec), "</td>", sep=""),
                 total = HTML("<td>", format(100*table(table[[var]],modulo)[level,level_m]/sum(table(table[[var]],modulo)), digits=nbdec, nsmall=nbdec), "</td>", sep=""))
        }
        
        HTML("<td></td></tr>")
     }
    }
    else next
  }
  dec()
  HTML("</tbody>")
  dec()
  HTML("</table>")
  HTML("</div>")
  HTMLEnd()
}