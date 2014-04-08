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

# Conversion formulas for self-affinity parameter estimates (sap) to Dimension (fd) suggested in Hasselman (2013)
# PSD slope (if signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure)
psd2fd <- function(sap){return(fd <- 3/2 + ((14/33)*tanh(sap*log(1+sqrt(2)))) )}
# DFA slope (this is H in DFA)
dfa2fd <- function(sap){return(fd <- 2-(tanh(log(3)*sap)) )}
# SDA slope (simple 1-sap, but note that for some signals different PSD slope values project to 1 SDA slope)
sda2fd <- function(sap){return(fd <- 1-sap)}


sliceTS<-function(TSmat,epochSz=1) {
  # Slice columns of TSmat in epochs of size = epochSz
  init(c("plyr"))
  
  N<-dim(TSmat)
  return(llply(seq(1,N[1],epochSz),function(i) TSmat[i:min(i+epochSz-1,N[1]),1:N[2]]))
}

fltrIT <- function(TS,f){
  # Apply filtfilt to TS using f (filter settings)
  require("signal")
  
  return(filtfilt(f=f,x=TS))
  
}

brainButter <- function(TSmat, fs=500, band=c(lfHIp=4,hfLOp=40), Np=9){
  # Extract frequency bands from columns of TSmat that are commonly used in Neuroscience
  # Low Freq. High-pass (1st) and High Freq. Low-pass (2nd) FIR1 filter is applied at frequencies specified as band=c(lfHIp=...,hfLOp=...)
  require("signal")
  
  fHI <- butter(Np,band[[1]]*2/fs,"high")
  fLO <- butter(Np,band[[2]]*2/fs,"low")
  
  TSflt <- apply(TSmat,2,function(TS) filtfilt(f=fHI,x=TS))
  TSflt <- apply(TSflt,2,function(TS) filtfilt(f=fLO,x=TS))
  
  #TSflt <- apply(TSflt,2,fltrIT,f=fLO)
  
  return(TSflt)
}

brainFir1 <- function(TSmat, fs=500, band=c(lfHIp=4,hfLOp=40), Np=2/band[1]){
  # Extract frequency bands from columns of TSmat that are commonly used in Neuroscience
  # Low Freq. High-pass (1st) and High Freq. Low-pass (2nd) FIR1 filter is applied at frequencies specified as band=c(lfHIp=...,hfLOp=...)
  require("signal")
  
  if(2/band[1]>Np){print(paste("Incorrect filter order Np... using 2/",band[1]," = ",(2/band[1]),sep=""))}
  
  fBP <- fir1(floor(Np*fs),band/(fs/2),type="pass");
  
  TSflt <- apply(TSmat,2,function(TS) filtfilt(f=fBP,1,x=TS))
  
  
  #   fHI <- fir1(floor(Np*fs),band[1]/(fs/2),type="high");
  #   fLO <- fir1(floor(Np*fs),band[2]/(fs/2),type="low");
  # 
  #   TSflt <- apply(TSmat,2,function(TS) filtfilt(f=fHI,1,x=TS))
  #   TSflt <- apply(TSflt,2,function(TS) filtfilt(f=fLO,1,x=TS))
  
  
  return(TSflt)
}


ssi2sbi <- function(SImat,threshold){
  # Signed Similarity matrix to "signed binary" matrix
  
  idS   <- which(SImat<0)
  BImat <- abs(as.matrix(SImat))
  diag(BImat) <- 0
  BImat[BImat <= threshold] <- 0
  BImat[BImat >  threshold] <- 1
  BImat[idS] <- BImat[idS]*-1
  
  return(BImat)
}

si2bi <- function(SImat,threshold){
  # Unsigned Similarity matrix to unsigned binary matrix
  
  ifelse(any(SImat<0),{
    print("Signed matrix, use: ssi2sbi()")
    break},{
      BImat <- as.matrix(SImat)
      diag(BImat) <- 0
      BImat[BImat <= threshold] <- 0
      BImat[BImat >  threshold] <- 1})
  
  return(BImat)
}

ssi2sth <- function(SImat,threshold){
  # Signed Similarity matrix to "signed thresholded" matrix
  
  idS   <- which(SImat<0)
  THmat <- abs(as.matrix(SImat))
  diag(THmat) <- 0
  THmat[THmat <= threshold] <- 0
  THmat[idS] <- THmat[idS]*-1
  
  return(THmat)
}

si2th <- function(SImat,threshold){
  # Similarity matrix to thresholded matrix
  
  ifelse(any(SImat<0),{
    print("Signed matrix, use: ssi2sth()")
    break},{
      THmat <- as.matrix(SImat)
      THmat[THmat <= threshold] <- 0})
  
  return(THmat)
}


plotBIN <- function(BImat){
  
  g <- graph.adjacency(BImat, weighted=T, mode = "undirected",diag=F)
  g <- simplify(g)
  
  # set colors and sizes for vertices
  V(g)$degree <- degree(g)
  
  rev<-scaleRange(log1p(V(g)$degree))
  rev[rev<=0.3]<-0.3
  
  V(g)$color       <- rgb(scaleRange(V(g)$degree), 1-scaleRange(V(g)$degree),  0, rev)
  V(g)$size        <- 25*rev
  V(g)$frame.color <- NA
  
  # set vertex labels and their colors and sizes
  V(g)$label       <- V(g)$name
  V(g)$label.color <- rgb(0, 0, 0, rev)
  V(g)$label.cex   <- rev
  
  # set edge width and color
  
  E(g)$width <- 4
  E(g)$color <- rgb(.5, .5, 0, .6)
  set.seed(958)
  
  #   layout1=layout.spring(g)
  #    layout2=layout.fruchterman.reingold(g)
  #    layout3=layout.kamada.kawai(g)
  #   layout5 = layout.spring(g,mass=0.3,repulse=T)
  
  #   CairoFontMatch(fontpattern="Arial")
  #   CairoFonts(regular="Arial:style=Normal")
  
  #   CairoPDF(pname,10,10)
  #   plot(g, layout=layout.sphere)
  #   dev.off() 
  #   
  
  plot(g, layout=layout.sphere)
  
  return(g)
}

plotMAT <- function(BImat,l=NULL){
  
  g <- graph.adjacency(BImat, weighted=T, mode = "undirected",diag=F)
  #g <- simplify(g)
  
  # set colors and sizes for vertices
  V(g)$degree <- degree(g)
  
  rev<-scaleRange(V(g)$degree)
  rev[rev<=0.4]<-0.4
  
  V(g)$color       <- rgb(scaleRange(V(g)$degree), 1-scaleRange(V(g)$degree),  0, rev)
  V(g)$size        <- 20*rev
  V(g)$frame.color <- NA
  
  # set vertex labels and their colors and sizes
  V(g)$label       <- V(g)$name
  V(g)$label.color <- rgb(0, 0, 0, .8)
  V(g)$label.cex   <- 1.1
  
  # set edge width and color
  #  rew<-E(g)$weight
  #  rew[rew<=0.3]<-0.3 
  #     
  edge.central=edge.betweenness(g)
  # 
  for (i in 1:ecount(g)) {E(g)$width[i]=0.3+sqrt((edge.central[i]))}
  
  # E(g)$width <- 2*E(g)$weight
  E(g)$color <- rgb(.5, .5, 0, .6)
  set.seed(958) 
  
  if(is.null(l)){l<-layout.fruchterman.reingold(g,niter=500,area=vcount(g)^2.3,repulserad=vcount(g)^2.8)}
  
  plot(g,layout=l)
  return(g)
}


plotSIGNth <- function(sSImat){
  
  g <- graph.adjacency(sSImat, weighted=TRUE)
  E(g)$sign <- E(g)$weight
  E(g)$curved <- is.mutual(g)
  E(g)$lty <- ifelse( E(g)$sign > 0, 1, 1)
  E(g)$arrow.size <- .2
  E(g)$width <- 3
  #E(g)$color <- rgb(scaleRange(abs(E(g)$weight)), 1-scaleRange(abs(E(g)$weight)), 0, 1)
  #layout1=layout.fruchterman.reingold(g)
  
  V(g)$label.color <- rgb(0, 0, 0, 1)
  V(g)$label.cex <- 1.4
  V(g)$vs   <- graph.strength(g, mode="in")
  V(g)$vs.u <- scaleRange(graph.strength(g))
  #V(g)$color<- ifelse( V(g)$vs > 0, rgb(V(g)$vs.u, 1-V(g)$vs.u, 0, 1), rgb(1-V(g)$vs.u, V(g)$vs.u, 0, 1))
  
  E(g)$es.u  <- scaleRange(E(g)$weight)
  E(g)$color <- ifelse( E(g)$sign > 0, rgb(0, 1, 0, .2), rgb(1, 0, 0, .2))
  return(g)
}

plotSW <- function(n,k,p){
  
  g <- watts.strogatz.game(1, n, k, p)
  
  V(g)$degree <- degree(g)
  
  # set colors and sizes for vertices
  rev<-scaleRange(log1p(V(g)$degree))
  rev[rev<=0.2]<-0.2
  rev[rev>=0.9]<-0.9
  V(g)$rev <- rev
  
  V(g)$color       <- rgb(V(g)$rev, 1-V(g)$rev,  0, 1)  
  V(g)$size        <- 25*V(g)$rev
  
  # set vertex labels and their colors and sizes
  V(g)$label       <- ""  
  
  E(g)$width <- 1
  E(g)$color <- rgb(0.5, 0.5, 0.5, 1)
  
  return(g)
}

plotBA <- function(n,pwr,out.dist){
  #require("Cairo")
  
  g <- barabasi.game(n,pwr,out.dist=out.dist,directed=F)
  V(g)$degree <- degree(g)
  
  # set colors and sizes for vertices
  rev<-scaleRange(log1p(V(g)$degree))
  rev[rev<=0.2] <- 0.2
  rev[rev>=0.9] <- 0.9
  V(g)$rev <- rev
  
  V(g)$color    <- rgb(V(g)$rev, 1-V(g)$rev,  0, 1)  
  V(g)$size     <- 25*V(g)$rev
  # V(g)$frame.color <- rgb(.5, .5,  0, .4)  
  
  # set vertex labels and their colors and sizes
  V(g)$label <- ""  
  
  E(g)$width <- 1
  E(g)$color <- rgb(0.5, 0.5, 0.5, 1)
  
  return(g)
}

SWtest0 <- function(g){
  Nreps <- 10;
  histr  <- vector("integer",Nreps)
  target<- round(mean(degree(g)))
  now   <- target/2
  for(i in 1:Nreps){
    gt      <- watts.strogatz.game(dim=1, size=length(degree(g)), nei=now, 0)
    histr[i] <- round(mean(degree(gt)))
    ifelse(histr[i] %in% histr,break,{
      ifelse(histr[i]>target,{now<-now-1},{
        ifelse(histr[i]<target,{now<-now+1},{
          break})
      })
    })
  }
  return(gt)
}


# SWtestV <- function(g,N){
#  return(list(cp=transitivity(g,type="global"),cpR=transitivity(rewire(g,mode=c("simple"),niter=N),type="global"),lp=average.path.length(g), lpR=average.path.length(rewire(g,mode=c("simple"),niter=N))))
# }

SWtestE <- function(g,p=1,N=20){
  values <- matrix(nrow=N,ncol=6,dimnames=list(c(1:N),c("cp","cpR","cp0","lp","lpR","lp0")))
  
  for(n in 1:N) {
    gt<-SWtest0(g)
    values[n,] <- c(transitivity(g,type="localaverage"),transitivity(rewire.edges(g,prob=p),type="localaverage"),transitivity(gt,type="localaverage"),average.path.length(g),average.path.length(rewire.edges(g,prob=p)),average.path.length(gt))}
  values[n,values[n,]==0] <- NA #values[n,values[n,]==0]+1e-8}
  
  values   <- cbind(values,(values[,1]/values[,2])/(values[,4]/values[,5]),(values[,1]/values[,3]),(values[,4]/values[,6]),((values[,1]/values[,3])/values[,2])/((values[,4]/values[,6])/values[,5]))
  valuesSD <- data.frame(matrix(apply(values[,1:10],2,sd,na.rm=T),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  valuesAV <- data.frame(matrix(colMeans(values[,1:10],na.rm=T),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  return(list(valuesAV=valuesAV,valuesSD=valuesSD,valuesSE=valuesSD/sqrt(N)))
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


PLFsmall <- function(g){
  reload <- FALSE
  if("signal" %in% .packages()){
    warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
    reload <- TRUE
    detach("package:signal", unload=TRUE)}
  
  if(length(V(g))>100){warning("Vertices > 100, no need to use PLFsmall, use a binning procedure");break}
  d<-degree(g,mode="all")  
  y<-hist(d,breaks=0:length(V(g)),plot=F)$density 
  y<-y[y>0]
  if(length(y)==2){warning("Caution... Log-Log slope is a bridge (2 points)")}
  if(length(y)<2){warning("Less than 2 points in Log-Log regression... aborting");break}
  alpha=coef(lm(log(y) ~ poly(log(1:length(y)), degree=1), na.action="na.exclude") )[2]
  
  if(reload==TRUE){library(signal,verbose=FALSE,quietly=TRUE)}
  
  return(alpha)
  
}

FDrel <- function(g){
  d<-degree(g,mode="all")
  nbreaks <- round(length(V(g))/2)-1
  y<-hist(d,breaks=nbreaks,plot=F)$density 
  y<-y[y>0]
  return(FD <- -sum(y*log2(y))/-(log2(1/length(y))))  
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
