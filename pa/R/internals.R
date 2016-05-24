## Yang Lu Yang.Lu@williams.edu

## a collection of internal functions used in the package

## This helper function uses two arguments, vdep and vexp, to make a
## string representing a formula format appropriate for a linear
## model: "vdep ~ vexp[[1]] + vexp[[2]] + ... + vexp[[n]]"
.formula.make <- function(vdep, vexp) {
  formula(paste(vdep, "~", paste(vexp, collapse = " + ")))
}

## a convenience function to get the return of each category
## var decides whether it's portfolio or benchmark
.cat.ret <- function(x,
                     cat.var,
                     ret.var,
                     var){
  all.cat <- levels(x[[cat.var]])
  ret <- sapply(1:length(all.cat),
                function(i){x[x[[cat.var]] == all.cat[i], ][[ret.var]] %*%
                              x[x[[cat.var]] == all.cat[i], ][[var]]})
  ret <- as.array(ret)
  
  names(ret) <- all.cat
  
  return(ret)
}

## aggregate attribution results over periods
.aggregate <- function(object, raw.mat){
  agg <- matrix(apply(raw.mat, 1, sum))
  colnames(agg) <- paste(c(min(unique(as.character(object@date.var))),
                           max(unique(as.character(object@date.var)))),
                         collapse = ", ")
  rownames(agg) <- c("Allocation",
                     "Selection",
                     "Interaction",
                     "Active Return")
  return(agg)
}

## combine raw and aggregate brinson attribution results
.combine <- function(raw, agg){
  .list <- list()
  .list[[1]] <- round(raw, 4)
  .list[[2]] <- round(agg, 4)
  names(.list) <- c("Raw", "Aggregate")
  return(.list)
}

## An internal function to abstract the plotting of a bar plot for
## either exposure or return in a single period.
.bar.plot <- function(df,
                      type,  ## for ylab, either exposure or return
                      title
                      ){
  ## for R CMD check ONLY
  Name <- NULL
  rm(Name)
  Value <- NULL
  rm(Value)
  Type <- NULL
  rm(Type)
  
  bar.plot <- ggplot(df, aes(x = Name, y = Value, fill = Type)) +
    geom_bar(width = 0.5, position = position_dodge(), stat = "identity") +
      coord_flip() +
      ylab(type) + xlab("Sector") +
        geom_hline(yintercept = 0) +
          ggtitle(title) + 
          theme(panel.background = element_blank(), ## theme_blank(),
                ## title = title, 
                axis.line = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                plot.background = element_rect(fill = NA, colour = NA))
  
  return(bar.plot)
  
}


## faceted plot

.facet.plot <- function(df,
                        type,
                        title){
  ## circumvent R CMD check
  Name <- NULL
  rm(Name)
  Value <- NULL
  rm(Value)
  Type <- NULL
  rm(Type)
  
  facet.plot <- ggplot(df, aes(Name, Value, fill = Type)) +
    geom_bar(position = position_dodge(), stat = "identity") + coord_flip() + theme_bw()+
      facet_wrap( ~ Date) + scale_x_discrete(name = "Sector") + ylab(type) +
        ## opts(title = title)
        ggtitle(title)
  return(facet.plot)
}
