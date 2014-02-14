# BEYOND THE BOUNDARY TOOLBOX---------------------------------------------------------------------------------------------------

##' Functions used by scripts that reproduce analyses, simulations and graphs in 'Beyond the Boundary'
##'
##' @title BTBTB 
##' @param source(".../BTBTB.R")
##' @return Functions listed in this file.
##' @author Fred Hasselman (unless otherwise indicated);
##' Copyright (C) 2010-2014 Fred Hasselman
##' 


# INIT BTBfiles -----------------------------------------------------------------------------------------------------------------------

# This function tries to find a 'BTBfiles' folder based on getwd()
# First searches down from getwd(), if it doesn't find anything it moves up until the folder is found, or the time limit is exceeded.
pathfinder <- function(folder,maxtime=30){
  uptree   <- rev(lapply(2:length(strsplit(getwd(),"/")[[1]]),function(s) paste0(strsplit(getwd(),"/")[[1]][1:s],collapse="/")))
  pm <- proc.time()["elapsed"]
  for(lvl in 1:length(uptree)){
    YOURPATH  <- dir(path=uptree[[lvl]],pattern=folder,include.dirs=T,recursive=F,full.names=T)
    if((proc.time()["elapsed"]-pm)>maxtime){stop("Search time exceeded 'maxtime'")}
    if(length(nchar(YOURPATH))!=0){break}
  }
  return(YOURPATH)
}

# INIT PACKAGES -----------------------------------------------------------------------------------------------------------------------

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


# LINEAR SCALE CONVERSION BASED ON RANGE  ---------------------------------------------------------------------------------------------

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
  ifelse( mn==mx, {u <- rep(1, length(x))},{
    u  <- ((( x - mn ) * ( hi - lo )) / ( mx - mn )) + lo
    id <- complete.cases(u)
    u[!id]<-0
  })
  return(u)
}


# GRAPH PLOTTING ----------------------------------------------------------------------------------------------------------------------
graph2svg <- function(TDM,pname){
  
  # Create weighted Term-Term matrix
  tTM <- as.matrix(TDM)
  TTM <- tTM %*% t(tTM)
  TTM <- log1p(TTM)
  
  g <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
  g <- simplify(g)
  
  # Remove vertices used in search query
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
graphHood <- function(ig,Vname){
  
  Vrem <- which(V(ig)$name %in% c("~rdsp~","~imp~","~som~","~bod~","~mlt~"))
  ig <- ig - V(ig)$name[Vrem]
  
  idx <- which(V(ig)$name==Vname)
  sg <- graph.neighborhood(ig, order = 1, nodes=V(ig)[idx], mode = 'all')[[1]]
  
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
  
  set.seed(958)#5365, 227
  
  idx <- which(V(sg)$name==Vname)
  plot(sg,layout=layout.star(sg,center=V(sg)[idx]))
  
  return(sg)
}


# MULTIPLOT FUNCTION ------------------------------------------------------------------------------------------------------------------
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

# SEARCHERS ---------------------------------------------------------------------------------------------------------------------------

# A wrapper for grep, searches for string in source
regIT   <- function(str,src) {
  grep(str,src,ignore.case=TRUE,value=TRUE)
}

# Return matching strings in source using regmatches
matIT <- function(str,src) {
  m<-gregexpr(str,src,perl=TRUE) 
  #print(regmatches(src,m))
  return(regmatches(src,m))
}

# Find tags returned from grepexpr in the source and return unique matches as a regex pattern to use in subsequent searches
tagIT<-function(tags,src){
  tmp<-unique(sub("\\s$*","",unlist(regmatches(src,tags))))
  tmp<-sub("^\\s*","",tmp)
  regmatches(tmp,gregexpr("\\s",tmp))<-rep("\\s",length(tags))
  return(tmp)
}

# Here we actually need to change the corpus content, solution: pass as a list.
# Other complications, regmatches turns the corpus into a character list, solution: PlainTextDocument
# A call might look like this, where TMcorpus:  TMcorpus<-subIT(src<-list(tags=myTAGS,str="mine",cor=TMcorpus))

subIT  <- function(src) {
  if(any(which(src$tags==""))){src$tags<-src$tags[which(src$tags!="",arr.ind=T)]}
  ifelse((length(src$tags)==0),{
    print("No tags in src$tags !!!")
    return(src$cor)
  },{
    Ddata<-DMetaData(src$cor)
    m <-gregexpr(wordXX(src$tags),src$cor,perl=TRUE)
    subs<-regmatches(src$cor,m)
    regmatches(src$cor,m)<-rep(src$str,length(m))
    print(paste("Changed ",length(which(unlist(m)>0))," strings into: ",src$str))
    src$cor<-Corpus(VectorSource(src$cor))
    oldname<-names(Ddata)
    Ddata<-data.frame(cbind(Ddata,paste(m),paste(subs)),stringsAsFactors=options(stringsAsFactors=FALSE))
    names(Ddata)<-c(oldname,paste(src$str,"~pos",sep=""),paste(src$str,"~tag",sep=""))
    DMetaData(src$cor)<-Ddata
    return(src$cor)
  })
}

subNIT  <- function(src) {
  if(any(which(src$tags==""))){src$tags<-src$tags[which(src$tags!="",arr.ind=T)]}
  ifelse((length(src$tags)==0),{
    print("No tags in src$tags !!!")
    return(src$cor)
  },{
    Ddata<-DMetaData(src$cor)
    m <-gregexpr(wordXX(src$tags),src$cor,perl=TRUE)
    subs<-regmatches(src$cor,m,invert=TRUE)
    regmatches(src$cor,m,invert=TRUE)<-rep(src$str,length(m))
    print(paste("Changed ",length(which(unlist(m)>0))," strings into: ",src$str))
    src$cor<-Corpus(VectorSource(src$cor))
    oldname<-names(Ddata)
    Ddata<-data.frame(cbind(Ddata,paste(m),paste(subs)),stringsAsFactors=options(stringsAsFactors=FALSE))
    names(Ddata)<-c(oldname,paste(src$str,"~pos",sep=""),paste(src$str,"~tag",sep=""))
    DMetaData(src$cor)<-Ddata
    return(src$cor)
  })
}

subNIT2  <- function(src) {
  if(any(which(src$tags==""))){src$tags<-src$tags[which(src$tags!="",arr.ind=T)]}
  ifelse((length(src$tags)==0),{
    print("No tags in src$tags !!!")
    return(src$cor)
  },{
    Ddata<-DMetaData(src$cor)
    m <-gregexpr(src$tags,src$cor,perl=TRUE)
    subs<-regmatches(src$cor,m,invert=TRUE)
    regmatches(src$cor,m,invert=TRUE)<-rep(src$str,length(m))
    print(paste("Changed ",length(which(unlist(m)>0))," strings into: ",src$str))
    src$cor<-Corpus(VectorSource(src$cor))
    oldname<-names(Ddata)
    Ddata<-data.frame(cbind(Ddata,paste(m),paste(subs)),stringsAsFactors=options(stringsAsFactors=FALSE))
    names(Ddata)<-c(oldname,paste(src$str,"~pos",sep=""),paste(src$str,"~tag",sep=""))
    DMetaData(src$cor)<-Ddata
    return(src$cor)
  })
}

# PATTERN GENERATORS ------------------------------------------------------------------------------------------------------------------

# ( X1|X2 ) Example: wordXX(c("yes","no","maybe"))
wordXX   <- function(x,clps="|",pre="\\b(",post=")\\b") {
  if(any(which(grepl("(^~|~$)",x)))) {x[grep("(^~|~$)",x)]<-paste("(\\s*",x[grep("(^~|~$)",x)],"\\s*)",sep="")}
  paste(pre, paste(x,collapse=clps), post, sep="")
}


# ( XspaceY ) Example: wordXsY("yes","no")
wordXsY   <- function(x,y,clps="|",pre="\\b(",post=")\\b") {
  if(any(which(grepl("(^~|~$)",x)))) {x[grep("(^~|~$)",x)]<-paste("\\s*",x[grep("(^~|~$)",x)],"\\s*",sep="")}
  if(any(which(grepl("(^~|~$)",y)))) {y[grep("(^~|~$)",y)]<-paste("\\s*",y[grep("(^~|~$)",y)],"\\s*",sep="")}
  paste(pre, paste(x,collapse=clps), ")\\s(", paste(y,collapse=clps), post, sep="")
}

# ( XspaceY | YspaceX ) Example: wordXsYr("yes","no")
wordXsYr  <- function(x,y,clps="|",pre="\\b(",post=")\\b") {
  if(any(which(grepl("(^~|~$)",x)))) {x[grep("(^~|~$)",x)]<-paste("(\\s*",x[grep("(^~|~$)",x)],"\\s*)",sep="")}
  if(any(which(grepl("(^~|~$)",y)))) {y[grep("(^~|~$)",y)]<-paste("(\\s*",y[grep("(^~|~$)",y)],"\\s*)",sep="")}
  paste("(",paste(pre,paste(x,collapse=clps),")\\s(",paste(y,collapse=clps),")\\b",sep=""),
        ")|(",paste("\\b(",paste(y,collapse=clps),")\\s(",paste(x,collapse=clps),post,sep=""),")", sep="")
}
# 
# # ( XspaceYspaceZ | YspaceZspaceX |  ZspaceXspaceY ) Example: wordXsYsZr("yes","no","maybe")
# wordXsYsZr <- function(x,y,z) {
#   paste("(",paste("\\b(",paste(x,collapse="|"),")\\s(",paste(y,collapse="|"),")\\s(",paste(z,collapse="|"), ")\\b", sep=""),
#     ")|(",paste("\\b(",paste(y,collapse="|"),")\\s(",paste(z,collapse="|"),")\\s(",paste(x,collapse="|"), ")\\b", sep=""),
#     ")|(",paste("\\b(",paste(z,collapse="|"),")\\s(",paste(x,collapse="|"),")\\s(",paste(y,collapse="|"), ")\\b", sep=""),")")
# }

# (XspaceY)? | (YspaceY)?
wordXsYro <- function(x,y,clps="|",pre="\\b(",post=")?\\b") {
  if(any(which(grepl("(^~|~$)",x)))) {x[grep("(^~|~$)",x)]<-paste("(\\s*",x[grep("(^~|~$)",x)],"\\s*)",sep="")}
  if(any(which(grepl("(^~|~$)",y)))) {y[grep("(^~|~$)",y)]<-paste("(\\s*",y[grep("(^~|~$)",y)],"\\s*)",sep="")}
  paste("(",paste(pre, paste(x,collapse=clps),")\\s(",paste(y,collapse=clps),post,sep=""),
        ")|(",paste(pre, paste(y,collapse=clps),")\\s(",paste(x,collapse=clps),post,sep=""),")",sep="")
}

# LOOK FOR NEIGHBOURS (WORDS) ---------------------------------------------------------------------------------------------------------

# ( Npre X Npos ) Example: txt="Love thy neigbours"  
#                           matIT(NwN(1,"thy",0),txt)

NwN   <- function(Npre=NULL, word=NULL, Npos=NULL){
  ifelse(all(is.null(Npre),is.null(word),is.null(Npos)),
         print("NwordN: One or more arguments missing!"),
{
  #word<-gsub("[[:punct:]]+","",word)
  regstr <-paste(c("\\b", paste(rep("(\\w+[[:graph:]]*)?",Npre),collapse="\\s?"),
                   paste(c("",paste(c("(",paste(word,collapse="|"),")"),collapse=""),""),collapse="\\s"),
                   paste(rep("(\\w+[[:graph:]]*)?",Npos),collapse="\\s?"),"\\b"),collapse="")
  print(regstr) 
  return(regstr)
})
}

# ( Npre X Nmid Y Npos ) Example:  txt="Didn't you know? Love thy neigbours' neigbours, too!"  
#                                  matIT(NwNwN(0,c("Love","hate","you"),2,c("neigbours","friends","family"),0),gsub("[[:punct:]]+","",txt))  
#                                  matIT(NwNwN(1,"Love",3,"too",0),gsub("[[:punct:]]+","",txt))                                   
NwNwN   <- function(Npre=NULL, word1=NULL, Nmid=NULL, word2=NULL, Npos=NULL){
  ifelse(all(is.null(Npre),is.null(word1),is.null(Nmid),is.null(word2),is.null(Npos)),
         print("NwNwN: One or more arguments missing!"),
{
  word1<-gsub("[[:punct:]]+","",word1)
  word2<-gsub("[[:punct:]]+","",word2)
  regstr <-paste(c("\\b", paste(rep("(\\w+[[:graph:]]*)?",Npre),collapse="\\s?"),
                   paste(c("",paste(c("(",paste(word1,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                   paste(rep("(\\w+[[:graph:]]*)?",Nmid),collapse="\\s?"),
                   paste(c("",paste( c("(",paste(word2,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                   paste(rep("(\\w+[[:graph:]]*)?",Npos),collapse="\\s?"),"\\b"),collapse="")
  print(regstr) 
  return(regstr)
})
}

# (Npre X Nmid Y Npos) | (Npre Y Nmid X Npos) Example:    txt="Didn't you know? Your neigbours love you!"  
#                                                         matIT(NwNwNr(1,"know",2,"neigbours",1),gsub("[[:punct:]]+","",txt))  
#                                                         matIT(NwNwNr(0,"Your",2,"love",0),gsub("[[:punct:]]+","",txt))

NwNwNr  <- function(Npre=NULL, word1=NULL, Nmid=NULL, word2=NULL, Npos=NULL){
  ifelse(all(is.null(Npre),is.null(word1),is.null(Nmid),is.null(word2),is.null(Npos)),
         print("NwNwN: One or more arguments missing!"),
{ word1<-gsub("[[:punct:]]+","",word1)
  word2<-gsub("[[:punct:]]+","",word2)
  regstr1 <-paste(c("\\b", paste(rep("(\\w+[[:graph:]]*)?",Npre),collapse="\\s?"),
                    paste(c("",paste(c("(",paste(word1,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                    paste(rep("(\\w+[[:graph:]]*)?",Nmid),collapse="\\s?"),
                    paste(c("",paste( c("(",paste(word2,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                    paste(rep("(\\w+[[:graph:]]*)?",Npos),collapse="\\s?"),"\\b"),collapse="")
  regstr2 <-paste(c("\\b", paste(rep("(\\w+[[:graph:]]*)?",Npre),collapse="\\s?"),
                    paste(c("",paste(c("(",paste(word2,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                    paste(rep("(\\w+[[:graph:]]*)?",Nmid),collapse="\\s?"),
                    paste(c("",paste( c("(",paste(word1,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                    paste(rep("(\\w+[[:graph:]]*)?",Npos),collapse="\\s?"),"\\b"),collapse="")
  regstr<-paste("(",regstr1,")*|(",regstr2,")*",sep="")
  print(regstr) 
  return(regstr)
})
}


# GET TERM-DOC MATRIX FROM STEM FILTERED CORPUS ---------------------------------------------------------------------------------------

prpIT <- function(TMcorpus) {
  require(tm)
  
  # Global analysis to get most common STEMS of words
  # RWeka and the Snowball stemmer package do not work for me, so I used RTextTools (wordStem) to use with the tm package.
  #
  # Problems:
  # - wordStem() reads the document text in Corpus() as if it were one word
  # - Returning from wordStem(), the Corpus() structure is changed and TextDocumentMatrix() will not accept it
  #
  # Solution:
  # - Tokenize before passing to wordStem() in a tm_map() control structure
  # - Return using PlainTextDocument()
  # - The indexation in the ID meta tag is lost, easy to restore however using meta()
  
  # This list of common words found in science abstracts (including names of publishers) will be considered stopwords that need removal
  # Depending on the purpose of your mining operation, some words may need to be added or removed from this list
  # Use after toLow and removePunctuation (preserve_intra_word_dashes = TRUE)
  
  SCIstop=c("(c)", "ab", "about", "abstract", "achieve", "achieved", "achievement", "achieving", "addition", "additional", "additions", "affect", "affected", "affecting", "affects", "after", "age", "ages", "aim", "aims", "al", "also", "among", "analyses", "analysis", "approach", "approached", "approaches", "argue", "argued", "article", "assess", "assessed", "assessment", "associated", "association", "background", "base", "based", "basic", "basis", "bedford", "been", "between", "blackwell", "both", "characterisation", "characterised", "characteristic", "characteristics", "characterization", "characterized", "child", "common", "compared", "complex", "complexity", "condition", "conditional", "conditions", "consist", "consisted", "consistent", "consists", "contrast", "contrasted", "contrasting", "contribute", "contributed", "contributes", "contributing", "control", "controls", "copyright", "correlated", "correlates", "correlation", "correlations", "could", "data", "demonstrate", "demonstrated", "determin", "determined", "differ", "differed", "difference", "differences", "discussed", "discussion", "discussions", "doi", "during", "each", "each", "effect", "effective", "effectiveness", "effects", "either", "elsevier", "evaluate", "evaluated", "evaluation", "evaluations", "evidence", "evidenced", "examine", "examined", "examination", "exhibit", "experience", "experienced", "experiment", "experimental", "experiments", "factor", "factorial", "factors", "findings", "first", "found", "four", "four", "freeman", "from", "general", "generalisation", "generalise", "generalised", "generalist", "generalistic", "generalization", "generalize", "generalized", "grade", "group", "groups", "guillford", "hall", "have", "high", "however", "however", "identified", "identify", "implicate", "implicated", "implicates", "implication", "implications", "improve", "improved", "improvement", "improvements", "include", "included", "including", "increase", "increased", "increasing", "independent", "individual", "individuals", "influence", "influenced", "influencing", "informed", "initial", "into", "investigate", "investigated", "investigation", "involve", "involved", "involvement", "involving", "john", "kluwer", "large", "level", "levels", "made", "matched", "mcgraw-hill", "measure", "measured", "measurement", "measurements", "measures", "method,", "more", "most", "name", "normal", "observation", "observations", "observe", "observed", "one", "only", "other", "paper", "participant", "participants", "patient", "patients", "pearson", "pmid", "population", "populations", "posit", "posited", "potential", "predict", "predicted", "prediction", "predictions", "predicts", "prentice", "present", "presented", "press", "previous", "proposal", "propose", "proposed", "provide", "provided", "provides", "providing", "psychology", "publication", "publish", "published", "publishing", "purpose", "recent", "recently", "reduce", "reduced", "reducted", "reduction", "reductions", "related", "relation", "relations", "relationship", "relationships", "report", "reported", "reports", "require", "requires", "research", "response", "responses", "results", "reveal", "revealed", "reveals", "review", "reviewed", "reviews", "role", "sample", "school", "schools", "second", "select", "selected", "severe", "severity", "should", "show", "showed", "significant", "significantly", "similar", "small", "so", "some", "sons", "springer", "stimuli", "stimulus", "student", "students", "studied", "studies", "study", "studying", "subject", "subjects", "subpopulation", "subpopulations", "such", "suggest", "suggested", "suggests", "support", "supported", "supports", "task", "tasks", "test", "tested", "testing", "tests", "than", "that", "their", "there", "these", "they", "third", "this", "those", "three", "ti", "time", "title", "two", "under", "understand", "understanding", "understood", "underlying", "used", "using", "well", "were", "when", "whenever", "whether", "which", "while", "wiley", "with", "without", "worth", "year", "years")
  
  #SCIstop<-SCIstop[order(SCIstop)]
  
  #  Functions to apply to the corpus, stemming by wordStem (RTextTools) only works after applying a scan tokenizer function
  skipWords   <-function(x) removeWords(x, c(stopwords("english"),stopwords("SMART"),SCIstop))
  removePunct <-function(x) removePunctuation(x, preserve_intra_word_dashes = TRUE)
  
  # The corpus needs to be indexed after using wordStem(), otherwise we cannot get TDM
  # This can be achieved using the PlainTextDocument syntax, or after the transform as:
  # meta(TMcorpus,tag = "ID", type="local") <- seq(1:nDocs(stemTDM))
  STEMit   <-function(x) PlainTextDocument(paste(wordStem(MC_tokenizer(x),"english"), collapse=" "), id = ID(x), language = Language(x))
  
  funs     <-list(STEMit,skipWords,tolower,removePunct,removeNumbers,stripWhitespace)
  # Apply the FUN list, this is a nice wrapper (tm_reduce), be aware however that the order has impact on results: The list is read from right to left
  TMcorpus <-tm_map(TMcorpus, FUN= tm_reduce, tmFuns = funs)
  
  stemTDM  <-TermDocumentMatrix(TMcorpus)
  #nTerms(stemTDM)
  
  # Remove single occurences
  oneFr    <-findFreqTerms(stemTDM,lowfreq=0, highfreq=1)
  TMcorpus <-tm_map(TMcorpus,removeWords,oneFr)
  
  rm(oneFr,SCIstop,skipWords,removePunct,STEMit,funs,stemTDM)
  TMcorpus<-tm_map(TMcorpus,stripWhitespace)
  return(TMcorpus)
}


# TRY … CATCH -------------------------------------------------------------------------------------------------------------------------

##================================================================##
###  In longer simulations, aka computer experiments,            ###
###  you may want to                                             ###
###  1) catch all errors and warnings (and continue)             ###
###  2) store the error or warning messages                      ###
###                                                              ###
###  Here's a solution  (see R-help mailing list, Dec 9, 2010):  ###
##================================================================##

##' Catch *and* save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
##'
##' @title tryCatch both warnings (with value) and errors
##' @param expr an \R expression to evaluate
##' @return a list with 'value' and 'warning', where
##'   'value' may be an error caught.
##' @author Martin Maechler;
##' Copyright (C) 2010-2012  The R Core Team
##' 
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
