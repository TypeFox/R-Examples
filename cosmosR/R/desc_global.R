#' Descriptif global d'une table
#' 
#' Produit un tableau descriptif d'une table
#' 
#' Permet de produire un tableau descriptif des variables contenues dans la table.
#' Si les labels et formats sont definis et charges ils seront utilises pour peupler le tableau.
#' Le fichier de sortie est place dans ../HTML Output
#' @encoding UTF-8
#' @param table Table a utiliser
#' @param html Nom du fichier html, par defaut "desc_global.html"
#' @param titre Titre du tableau, par defaut "Descriptif global de nom_de_la_table"
#' @param variables Vecteur de noms de variables a decrire, par defaut toutes les variables contenues dans la table
#' @param variables_neg Vecteur de noms de variables a exclure de la description
#' @param stats Vecteur de valeurs a calculer, parmi N, \%, \%/moy, moy, \%/med, med, et, ic95, q1, med, q3. Par defaut c("N","\%/moy","ic95")
#' @param miss Booleen : afficher ou non les valeurs manquantes, par defaut TRUE
#' @param note Note de bas de page, par defaut vide
#' @param nbdec Nombre de decimales apres la virgule, par defaut 1
#' @examples
#' \dontrun{
#' Ma_table <- charger("donnees.xls")
#' 
#' desc_global(Ma_table) # descriptif par defaut
#' desc_global(Ma_table, variables=c("var1","var2"), stats=c("N","%"), note="Note de bas de page")
#' }
#' @export
desc_global <- function(table, html="desc_global", titre=NULL, variables=NULL, variables_neg=NULL, stats=c("N","%/moy","ic95"), miss=TRUE, note=NULL, nbdec=1)
{
  if (missing(table))
  {
    warning("Pas de table donnee !")
    return
  }
  
  # Creation du titre
  if (is.null(titre))
    titre<-paste("Descriptif global de",deparse(substitute(table))) 
  
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
  
  # Header de la table
  HTMLInit(file=paste0("../HTML\ Output/",html,".html"), title=titre, CSSfile="desc.css")
  HTML("<div class='desc'>")
  inc()
  HTML("<table class='desc'>")
  inc()
  HTML("<caption>", titre, "</caption>")
  HTML("<thead>")
  inc()
  HTML("<tr><th></th>", paste("<th>",stats,"</th>",sep="", collapse=""), "</tr>")
  dec()
  HTML("</thead>")
  
  # Footer de la table
  if (!is.null(note))
  {
    HTML("<tfoot>")
    inc()
    HTML("<tr><td colspan='", length(stats)+1, "'>", note,"</td></tr>")
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
      HTML("<tr><td class='var'>", etiq, "</td>", sep="")
      
      for (stat in stats)
      {
        HTML("<td>",sep="")
        
        HTML(tryCatch(
        {
          format(switch(stat,
          N = length(na.omit(table[[var]])),
          "%/moy" = ,
          moy = mean(table[[var]],na.rm=T),
          "%/med" = ,
          med = median(table[[var]],na.rm=T),
          Q1 = quantile(table[[var]], probs=.25, na.rm=T),
          Q3 = quantile(table[[var]], probs=.75, na.rm=T),
          min = min(table[[var]],na.rm=T),
          max = max(table[[var]],na.rm=T),
          et = sd(table[[var]],na.rm=T),
          ic95 = paste0("&plusmn;",format(mean(table[[var]],na.rm=T) - t.test(table[[var]])$conf.int[1], digits=nbdec,nsmall=nbdec)),
          ... = "ERR"
          ),digits=nbdec,nsmall=nbdec)
        },error = function(e){erreur(e,var,"dans le calcul de ",stat)}), sep="")
        
        HTML("</td>", sep="")
      }
      
      HTML("</tr>")
    }
    # Variable qualitative
    else if (is.factor(table[[var]]))
    {
      HTML("<tr><td class='var' colspan='", length(stats)+1, "'>", etiq, "</td></tr>")
      
      # Donnees manquantes
      if ((miss) & (!is.na(summary(table[[var]])["NA's"])))
      {
        HTML("<tr><td class='level'>Manquant</td>",sep="")
      
        for (stat in stats)
        {
          HTML("<td>",sep="")
      
          if (stat == "N")
            HTML(summary(table[[var]])["NA's"],sep="")
      
          HTML("</td>",sep="")
        }
        HTML("</tr>")
      }

      # Levels de la variable qualitative
      for (level in levels(table[[var]]))
      {
        HTML("<tr><td class='level'>", level, "</td>",sep="")
        for (stat in stats)
        {
          HTML("<td>",sep="")
      
          if (stat == "N")
            HTML(summary(table[[var]])[level],sep="")
          else if ((stat == "%") | (stat == "%/moy") | (stat == "%/med"))
            HTML(format(100*summary(table[[var]])[level]/length(na.omit(table[[var]])),digits=nbdec,nsmall=nbdec),sep="")

          HTML("</td>",sep="")
        }
        HTML("</tr>")
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

#' Diagnostics de la table
#' 
#' Produit un tableau descriptif de la table et des graphiques pour la verification des conditions d'utilisation des tests satistiques
#' 
#' Permet de produire un tableau descriptif des variables contenues dans la table.
#' Si les labels et formats sont definis et charges ils seront utilises pour peupler le tableau.
#' La fonction renvoie un vecteur de noms de variables considerees comme parametriques apres un test de normalite de Shapiro-Wilk.
#' Le fichier est cree dans le repertoire temporaire. Il est possible de le sauvegarder avec ses graphiques a partir du navigateur.
#' @encoding UTF-8
#' @param table Table a utiliser
#' @param variables Vecteur de noms de variables a decrire, par defaut toutes les variables contenues dans la table
#' @param variables_neg Vecteur de noms de variables a exclure de la description
#' @return Un vecteur contenant les variables considerees comme parametriques
#' @examples
#' \dontrun{
#' Ma_table <- charger("donnees.xls")
#' 
#' diagnostic(Ma_table) # diagnostic par defaut
#' parametriques <- diagnostic(Ma_table, variables_neg=c("num_id"))
#' # parametriques contient le vecteur de noms de variables parametriques
#' }
#' @export
diagnostic <- function(table, variables=NULL, variables_neg=NULL)
{
  if (missing(table))
  {
    warning("Pas de table donnee !")
    return
  }
  
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
  
  HTMLInit(CSSfile="diag.css")
  param<-character(0)
  
  HTML("<div class='diag_menu'>")
  inc()
  for (var in names(table))
  {
    # Creation du label de la variable
    etiq <- label(table[[var]])
    if (etiq == "table[[var]]") etiq <- var
    
    HTML("<p id='", var, "'>", etiq, sep="")
    tryCatch(
    {
    if (is.numeric(table[[var]]))
      if (shapiro.test(table[[var]])$p.value > .05)
        HTML("*", sep="")
    }, error = function(e){})
    HTML("</p>")
  }
  dec()
  HTML("</div>")
  
  for (var in names(table))
  {
    # Creation du label de la variable
    etiq <- label(table[[var]])
    if (etiq == "table[[var]]") etiq <- var
    
    HTML("<div class='diag_varblock' id='div_", var, "' style='display:none'>",etiq)
    inc()
      HTML("<div class='diag_varstat'>")
      inc()
        HTML("<table class='diag'>")
        inc()
        HTML("<tr>", paste("<th>", names(summary(table[[var]])), "</th>", collapse=""), "</tr>")
        HTML("<tr>", sep="")
        for (valeur in summary(table[[var]]))
        {
          HTML("<td>", valeur, "</td>", sep="")
        }
        HTML("</tr>")
        dec()
        HTML("</table>")
        if (is.numeric(table[[var]]))
        {
          tryCatch(HTML("IC95 : [ ", format(t.test(table[[var]])$conf.int[1],digits=2,nsmall=2), " - ", format(t.test(table[[var]])$conf.int[2],digits=2,nsmall=2), " ]<br/>"), error=function(e){erreur(e,var,"dans le calcul de l'ic95",ret=HTML("IC95 non calculable<br/>"))})
          tryCatch(HTML("Test de normalit&eacute; de Shapiro-Wilk : p = ", format(shapiro.test(table[[var]])$p.value,digits=2,nsmall=2), "<br/>"), error=function(e){erreur(e,var,"dans l'evaluation de sa normalite",ret=HTML("Test de normalit&eacute non calculable"))})
          tryCatch({
            if (shapiro.test(table[[var]])$p.value > .05)
            {
              param<-c(param,var)
              HTML("<b>Param&eacute;trique</b>")
            }
          },error=function(e){})
        }
      dec()
      HTML("</div>")
    
      HTML("<div class='diag_varplot'>")
      inc()
        if (is.numeric(table[[var]]))
        {
          pngfile=tempfile(pattern="figure",fileext=".png")
          png(filename=pngfile)
          hist(table[[var]],freq=F,col="lightblue",border="darkblue",xlab=etiq,ylab="Densite",main=paste0("Distribution de ",etiq))
          tryCatch(lines(density(table[[var]],na.rm=T),col="red"),error=function(e){})
          dev.off()
          HTML("<img src='", basename(pngfile), "'/>")
          
          pngfile=tempfile(pattern="figure",fileext=".png")
          png(filename=pngfile)
          qqnorm(table[[var]], xlab="Quantiles theoriques",ylab="Quantiles de l'echantillon")
          qqline(table[[var]])
          dev.off()
          HTML("<img src='", basename(pngfile), "'/>")
        }
        else if (is.factor(table[[var]]))
        {
          pngfile=tempfile(pattern="figure",fileext=".png")
          png(filename=pngfile)
          barplot(table(table[[var]]), col="lightblue", border="darkblue")
          dev.off()
          HTML("<img src='", basename(pngfile), "'/>")
        }
      dec()
      HTML("</div>")
    dec()
    HTML("</div>")
  }
  HTMLEnd()
  param
}