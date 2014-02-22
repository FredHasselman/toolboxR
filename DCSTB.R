# DYNAMICS OF COMPLEX SYSTEMS TOOLBOX-------------------------------------------
#
##' @title DCSTB 
##' @param source(".../BTBTB.R")
##' @return Functions listed in this file.
##' @author Fred Hasselman (unless otherwise indicated);
##' Copyright (C) 2010-2014 Fred Hasselman 


# INIT PACKAGES ----------------------------------------------------------------

# Packages in the list argument need will be installed if necessary and loaded
init <- function(need){
  ip <- .packages(all.available=T)
  if(any((need %in% ip)==F)){install.packages(need[!(need %in% ip)])}
  ok <- sapply(seq(along=need),function(p) require(need[[p]],character.only=T))
}

initIO <- function(){
  # I/O and data handling tools
  ip <- .packages(all.available=T)
  need <- c("foreign","xlsx","plyr","doBy","reshape","reshape2")
  if(any((need %in% ip)==F)){install.packages(need[!(need %in% ip)])}
  
  require("foreign")
  require("xlsx")
  require("plyr")
  require("doBy")
  require("reshape")
  require("reshape2")
}

initNLTS <- function(){
  # Initialise Nonlinear Time Series packages
  ip <- .packages(all.available=T)
  need <- c("fractaldim","fractalrock","RTisean","tsDyn","tseries","tseriesChaos")
  if(any((need %in% ip)==F)){install.packages(need[!(need %in% ip)])}
  
  require("fractaldim")
  require("fractalrock")
  require("tseries")
  require("tseriesChaos")
  require("RTisean")
  require("tsDyn")
  
}

initSIGNAL <- function(){
  # Initialise Signal analysis packages
  ip <- .packages(all.available=T)
  need <- c("pracma","signal","EMD","hht")
  if(any((need %in% ip)==F)){install.packages(need[!(need %in% ip)])}
  
  require("pracma")
  require("signal")
  require("EMD")
  require("hht")
}

initPAR <- function(){
  # Parallel computing tools
  ip <- .packages(all.available=T)
  need <- c("parallel","doParallel","foreach")
  if(any((need %in% ip)==F)){install.packages(need[!(need %in% ip)])}
  
  require("parallel")
  require("doParallel")
  require("foreach")
}

initPLOT <- function(useArial = T,afmPATH="/Volumes/Fred HD/Rplus"){
  # Load packages for plotting with default option to setup Arial as the pdf font for use in figures.
  
  ip <- .packages(all.available = T)
  need <- c("lattice","gplots","ggplot2","grid","scales","aplpack","effects","RColorBrewer","GGally","mapproj")
  
  if(any((need %in% ip)==F)){install.packages(need[!(need %in% ip)])}
  
  require("lattice")
  require("gplots")
  require("ggplot2")
  require("grid")
  require("scales")
  require("aplpack")
  require("effects")
  require("RColorBrewer")
  require("GGally")
  require("mapproj")
  
  if(useArial==T){
    # Set up PDF device on MAC OSX to use Arial as a font in Graphs
    if (!"Arial" %in% names(pdfFonts())) {
      Arial <- Type1Font("Arial",
                         c(paste(afmPATH,"/Arial.afm",sep=""), 
                           paste(afmPATH,"/Arial Bold.afm",sep=""), 
                           paste(afmPATH,"/Arial Italic.afm",sep=""),
                           paste(afmPATH,"/Arial Bold Italic.afm",sep="")))
      pdfFonts(Arial=Arial)    
    }
  }  
}

setArial <- function(afmPATH){
  # Set up PDF device on MAC OSX to use Arial as a font in Graphs
  if (!"Arial" %in% names(pdfFonts())) {
    Arial <- Type1Font("Arial",
                       c(paste(afmPATH,"/Arial.afm",sep=""), 
                         paste(afmPATH,"/Arial Bold.afm",sep=""), 
                         paste(afmPATH,"/Arial Italic.afm",sep=""),
                         paste(afmPATH,"/Arial Bold Italic.afm",sep="")))
    pdfFonts(Arial=Arial)    
  }
}


# LINEAR SCALE CONVERSION BASED ON RANGE  --------------------------------------

# # Three uses:
# # 
# # 1. scaleRange(x)             Scale x to data range: min(x.out)==0;      max(x.out)==1 
# # 2. scaleRange(x,mn,mx)       Scale x to arg. range: min(x.out)==mn==0;  max(x.out)==mx==1
# # 3. scaleRange(x,mn,mx,lo,hi) Scale x to arg. range: min(x.out)==mn==lo; max(x.out)==mx==hi 
# # 
# # Examples:
# # 
# # Works on numeric objects
# somenumbers <- cbind(c(-5,100,sqrt(2)),c(exp(1),0,-pi))
# 
# scaleRange(somenumbers)
# scaleRange(somenumbers,mn=-100)
# 
# # Values < mn will return < lo (default=0)
# # Values > mx will return > hi (default=1)
# scaleRange(somenumbers,mn=-1,mx=99)
# 
# scaleRange(somenumbers,lo=-1,hi=1)
# scaleRange(somenumbers,mn=-10,mx=101,lo=-1,hi=4)

scaleRange <- function(x, mn = min(x), mx = max(x), lo = 0, hi = 1){
  if(mn > mx){ warning("Minimum (mn) > maximum (mx).")}
  if(lo >= hi){ warning("Lowest scale value (lo) >= highest scale value (hi).")}
  ifelse( mn==mx, {u <- rep(hi, length(x))},{
    u  <- ((( x - mn ) * ( hi - lo )) / ( mx - mn )) + lo
    id <- complete.cases(u)
    u[!id]<-0
  })
  return(u)
}


# GRAPH PLOTTING ---------------------------------------------------------------
graph2svg <- function(TDM,pname){
  
  # Create weighted Term-Term matrix
  tTM <- as.matrix(TDM)
  TTM <- tTM %*% t(tTM)
  TTM <- log1p(TTM)
  
  g <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
  g <- simplify(g)
  
  # Remove vertices that were used in the search query
  Vrem <- which(V(g)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
  g <- (g - V(g)$name[Vrem])
  
  # Set colors and sizes for vertices
  V(g)$degree <- degree(g)
  rev         <- scaleRange(log1p(V(g)$degree))
  rev[rev<=0.3]<-0.3
  
  V(g)$color       <- rgb(scaleRange(V(g)$degree), 1-scaleRange(V(g)$degree),  0, rev) 
  V(g)$size        <- 10*scaleRange(V(g)$degree)
  V(g)$frame.color <- NA
  
  # set vertex labels and their colors and sizes
  V(g)$label       <- V(g)$name
  V(g)$label.color <- rgb(0, 0, 0, rev)
  V(g)$label.cex   <- scaleRange(V(g)$degree)+.1
  
  # set edge width and color
  rew <- scaleRange(E(g)$weight)
  rew[rew<=0.3]<-0.3 
  
  E(g)$width <- 2*scaleRange(E(g)$weight) 
  E(g)$color <- rgb(.5, .5, 0, rew)
  set.seed(958)
  
  svg(paste(pname,sep=""),width=8,height=8)
  plot(g, layout=layout.fruchterman.reingold(g))
  dev.off()
  
  return(g)
}

# Plot vertex neighbourhood
hoodGraph2svg <- function(TDM,Vname,pname){
  
   # Create weighted Term-Term matrix
  tTM <- as.matrix(TDM)
  TTM <- tTM %*% t(tTM)
  TTM <- log1p(TTM)
  
  ig <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
  ig <- simplify(ig)
  
  # Remove vertices that were used in the search query
  Vrem <- which(V(ig)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
  ig <- (ig - V(ig)$name[Vrem])
  
  # This is a deletion specific for the Neighbourhood graphs
  Vrem <- which(V(ig)$name %in% c("~rdsp~","~imp~","~som~","~bod~","~mlt~"))
  ig   <- ig - V(ig)$name[Vrem]
  
  idx <- which(V(ig)$name==Vname)
  sg  <- graph.neighborhood(ig, order = 1, nodes=V(ig)[idx], mode = 'all')[[1]]
  
  # set colors and sizes for vertices
  V(sg)$degree <- degree(sg)
  
  rev<-scaleRange(log1p(V(sg)$degree))
  rev[rev<=0.3]<-0.3
  
  V(sg)$color <- rgb(scaleRange(V(sg)$degree), 1-scaleRange(log1p(V(sg)$degree*V(sg)$degree)),  0, rev)
  
  V(sg)$size        <- 35*scaleRange(V(sg)$degree)
  V(sg)$frame.color <- NA
  
  # set vertex labels and their colors and sizes
  V(sg)$label       <- V(sg)$name
  V(sg)$label.color <- rgb(0, 0, 0, rev)
  V(sg)$label.cex   <- scaleRange(V(sg)$degree)
  
  # set edge width and color
  rew<-scaleRange(E(sg)$weight)
  rew[rew<=0.3]<-0.3 
  
  E(sg)$width <- 6*scaleRange(E(sg)$weight)
  E(sg)$color <- rgb(.5, .5, 0, rew)
  
  idV <- which(V(sg)$name==Vname)
  idE <- incident(sg,V(sg)[[idV]])
  E(sg)$color[idE] <- rgb(0, 0, 1 ,0.8)
  
  set.seed(958)
  
  idx <- which(V(sg)$name==Vname)
  svg(paste(pname,sep=""),width=8,height=8)
  plot(sg,layout=layout.star(sg,center=V(sg)[idx]))
  dev.off()
  
  return(sg)
}


# MULTIPLOT FUNCTION -----------------------------------------------------------
#
# [copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ ]
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multi.PLOT <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# TRY … CATCH ------------------------------------------------------------------

##================================================================##
###  In longer simulations, aka computer experiments,            ###
###  you may want to                                             ###
###  1) catch all errors and warnings (and continue)             ###
###  2) store the error or warning messages                      ###
###                                                              ###
###  Here's a solution  (see R-help mailing list, Dec 9, 2010):  ###
##================================================================##

# Catch *and* save both errors and warnings, and in the case of
# a warning, also keep the computed result.
#
# @title tryCatch both warnings (with value) and errors
# @param expr an \R expression to evaluate
# @return a list with 'value' and 'warning', where
#   'value' may be an error caught.
# @author Martin Maechler;
# Copyright (C) 2010-2012  The R Core Team
# 
try.CATCH <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}
