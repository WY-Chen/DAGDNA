# Copyright (c) 2018 - 2019  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.


#' find source component
#'
#' BFS
#' @param D: Partial ordering matrix
#' @param i: root
#' @return S: index of the component containing i
#'
find_comp_source<-function(D,i){
  if (sum(D[-i,i]!=0)==0){return(i)}
  S=unique(c(i,which(D[,i]!=0)))
  while (sum(D[-S,S]!=0)!=0){
    S=unique(c(S,which(rowSums(D[,S])!=0)))
  }
  return(S)
}
#' find sink component
#'
#' BFS
#' @param D: Partial ordering matrix
#' @param i: root
#' @return S: index of the component containing i
#'
find_comp_sink<-function(D,i){
  if (sum(D[i,-i]!=0)==0){return(i)}
  S=unique(c(i,which(D[i,]!=0)))
  while (sum(D[S,-S]!=0)!=0){
    S=unique(c(S,which(colSums(D[S,])!=0)))
  }
  return(S)
}




#' Deduce order-constraining DNA from arbitrary DNA
#'
#' BFS
#' @param D: Partial ordering matrix
#' @param i: root
#' @return S: index of the component containing i
#'
orderConstraining<-function(D,verbose=F){
  temp=D
  sources=sinks=NULL
  p=dim(D)[1]
  checked=c(sources,sinks,p+1)
  sinksl=sourcesl=list()
  done=F
  # find singleton first
  while (!done){
    s = which.min(rowSums(temp))
    found=F
    if (sum(temp[s,])==0){
      # sink
      if(verbose){cat("Sink  Found",s,"\n")}
      D[setdiff(seq(p),c(checked,s)),s]=1
      temp[s,]=temp[,s]=0;temp[s,s]=Inf
      sinks=c(s,sinks)
      sinksl=c(list(s),sinksl)
      found=T
    }
    checked=c(sinks,sources)
    k = which.min(colSums(temp))
    if (sum(temp[,k])==0){
      # source
      if(verbose){cat("Source  Found",k,"\n")}
      D[k,setdiff(seq(p),c(checked,k))]=1
      temp[k,]=temp[,k]=0;temp[k,k]=Inf
      sources=c(k,sources)
      sourcesl=c(sourcesl,list(k))
      found=T
    }
    checked=c(sinks,sources)
    if (!found){done=T}
  }
  # check chain comps
  done=length(checked)==p
  while (!done){
    s = which.min(rowSums(temp))
    k = which.min(colSums(temp))
    compsink=find_comp_sink(temp,s)
    compsource=find_comp_source(temp,k)
    if (length(compsink)<=length(compsource)){
      # bfs for sinks
      if(verbose){cat("Sink  Found",compsink,"\n")}
      D[setdiff(seq(p),c(checked,compsink)),compsink]=1
      temp[compsink,]=temp[,compsink]=0
      temp[compsink,compsink]=Inf
      sinks=c(compsink,sinks)
      sinksl=c(list(compsink),sinksl)
      checked=c(sinks,sources)
    } else {
      # bfs for sources
      if(verbose){cat("Source Found",compsource,"\n")}
      D[compsource,setdiff(seq(p),c(checked,compsource))]=1
      temp[compsource,]=temp[,compsource]=0
      temp[compsource,compsource]=Inf
      sources=c(compsource,sources)
      sourcesl=c(sourcesl,list(compsource))
      checked=c(sinks,sources)
    }
    if (length(checked)==p){done=T}
  }
  return(list(D=D,sinks=sinksl,sources=sourcesl))
}


#' Deduce all DNA from DAG
#'
#' using cpdag characterization
#' @param gm: DAG matrix
#' @return D: DNA set
#'
allDNA<-function(gm){
  p=dim(gm)[1]
  D=matrix(2,p,p);diag(D)=0
  cpg = dag2cpdag(gm)
  for (i in seq(p)){
    s=snew = which(cpg[i,]!=0) # possibly directed path
    while (length(snew)){
      snew=setdiff(which(colSums(cpg[c(i,s),])!=0),s)
      s=unique(c(s,snew))
    }
    D[i,setdiff(seq(p),c(i,s))]=0
  }
  return(D)
}

