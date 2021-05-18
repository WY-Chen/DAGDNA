# Copyright (c) 2018 - 2021  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.

#' Learn DAG from ordering
#'
#' back-fitting
#' @param CIFUN: CItest
#' @param TO: full ordering
#' @param h: dictionary of CI relations
#' @param ginit: initial graph
#' @return g: DAG matrix
#'
DAG_from_TO<-function(CIFUN,TO,h=NULL,ginit=NULL){
  # Get the IMAP with arbitrary ordering by CI oracle
  # h is a global hash table to cound unique CI statement checked
  stopifnot(is.function(CIFUN))
  hashed=!is.null(h)
  p = length(TO)
  g = matrix(0,p,p)
  if (is.null(ginit)){ginit<-matrix(1,p,p)}else{ginit<-ginit+t(ginit)}
  g[TO[1],TO[2]]=!CIFUN(TO[1],TO[2],NULL)
  for (i in seq(3,p)){
    x = TO[i]
    xprec=TO[1:(i-1)]
    xprecnb=xprec[which(ginit[xprec,x]!=0)]
    if (length(xprecnb)>0){
      for (y in xprecnb){
        ci=CIFUN(x,y,setdiff(xprec,y))
        if (hashed){h[[paste(c(sort(c(x,y)),"|",sort(setdiff(xprec,y))),collapse = ",")]]=1}
        g[y,x]=!ci
      }
    }
  }
  return(g)
}

#' Find covered edges
#'
#' by looking at the working DAG
#' @param gm: a working DAG
#' @param POM: partial ordering matrix. If provided, only report covered edges that agrees with ordering
#' @param layers: list of layer indices. If provided, only report covered edges that agrees with ordering
#' @param conservative: if yes, only use layers, do not use full POM
#' @return list of covered edges
#'
find_covered_edges<-function(gm,POM,layers,conservative=T){
  # list all pairs of covered edges,
  # i,j such that pa(i)U{j}=pa(j). with ordered i,j
  # omit the i,j pair if POM[i,j]=1, i.e., j->i is forbidden
  coverededges = NULL
  edge_list=which(gm!=0,arr.ind = T)
  if (length(edge_list)==0){return(NULL)}
  for (x in seq(dim(edge_list)[1])){
    i=edge_list[x,1];j=edge_list[x,2]
    if (!conservative){
      # aggressive DNA
      if (POM[j,i]!=0 &&
          setequal(c(i,which(gm[,i]!=0)),which(gm[,j]!=0))){
        coverededges=rbind(coverededges,c(i,j))
      }
    } else {
      # conservative DNA
      if (layers[i]==layers[j] &&
          setequal(c(i,which(gm[,i]!=0)),which(gm[,j]!=0))){
        coverededges=rbind(coverededges,c(i,j))
      }
    }
  }
  return(coverededges)
}

#' Reverse covered edges
#'
#' by permuting the ordering and backfitting
#' @param CIFUN: CI function
#' @param gm: working DAG
#' @param i,j: edge to be reversed
#' @param h: dictionary of CI statements
#' @param mtd: two methods: sort: flip the ordering and estimate the DAG via DAG_from_TO; flip: CI tests
#' @return DAG matrix
#'
covered_edge_reversal<-function(CIFUN,gm,i,j,h=NULL,mtd="sort"){
  # the old edge is i->j, we flip it to j->i
  # and correct pa(i) and pa(j) using CI oracle
  # h is a hash table of CI statements for counting purposes
  # two methods:
  #   sort: flip the ordering and estimate the DAG via DAG_from_TO
  #   flip: CI tests
  pai = which(gm[,i]!=0)
  paj = which(gm[,j]!=0)
  hashed=!is.null(h)
  stopifnot(setequal(paj,c(i,pai)))
  stopifnot(is.function(CIFUN))
  if (mtd=="sort"){
    currTO=as.numeric(tsort(as(gm,"graphNEL")))
    # a coding trick to move j before i
    indt=seq(p);indt[which(currTO==j)]=indt[which(currTO==i)]-0.5
    newTO=currTO[order(indt)]
    return(DAG_from_TO(CIFUN,TO = newTO,h=h))
  }
  for (x in pai){
    if (hashed){
      if (exists(paste(c(sort(c(x,i)),sort(setdiff(paj,x))),
                       collapse = ","),h)){
        gm[x,i]=!h[[paste(c(sort(c(x,i)),sort(setdiff(paj,x))),
                          collapse = ",")]]
      } else {
        gm[x,i]= !CIFUN(x,i,setdiff(c(pai,j),x))
        h[[paste(c(sort(c(x,i)),sort(setdiff(paj,x))),
                 collapse = ",")]]=!gm[x,i]
      }
      if (exists(paste(c(sort(c(x,j)),sort(setdiff(pai,x))),
                       collapse = ","),h)){
        gm[x,j]=!h[[paste(c(sort(c(x,j)),sort(setdiff(pai,x))),
                          collapse = ",")]]
      } else {
        gm[x,j]= !CIFUN(x,j,setdiff(pai,x))
        h[[paste(c(sort(c(x,j)),sort(setdiff(pai,x))),
                 collapse = ",")]]=!gm[x,j]
      }
    } else {
      gm[x,i]= !CIFUN(x,i,setdiff(c(pai,j),x))
      gm[x,j]= !CIFUN(x,j,setdiff(pai,x))
    }
  }
  gm[j,i]=1;gm[i,j]=0
  return(gm)
}

#' (internal) move to next graph
#'
#' by finding covered edges and reverting them
#' @param CIFUN: CI function
#' @param groot: root DAG
#' @param POM: partial ordering
#' @param layers: list of layers
#' @param ecount: current number of edges
#' @param h: dictionary of CI statements
#' @param gh: dictionary of DAGs
#' @param prune: if yes, skip DAG in gh that has been checked. Otherwise not.
#' @param convervative: censor by POM (aggressive) or layers (conservative) when finding covered edges
#' @return list of DAGs
#'
next_graph<-function(CIFUN,groot,POM,layers,ecount,h=NULL,gh,prune=T,verbose=F,conservative=T){
  # find all graphs that can be obtained by reversing
  # a pair of covered edges.
  # edges forbidden by POM is reserved to be un-flippable
  # h is a hash table of CI statements for counting purpose
  # gh is a hash table of DAGs for pruning, turn off by prune=F
  stopifnot(is.function(CIFUN))
  if (sum(groot)>ecount){
    if(verbose){cat(")")}
    return(NULL)
  } # not weaking decreasing, cut branch
  covered_edges = find_covered_edges(groot,POM,layers,conservative)
  if (length(covered_edges)==0){
    if(verbose){cat(")")}
    return(NULL)
  }# leaf case of DFS, end branch
  # get all children
  gnext=NULL;j=1
  for (i in seq(dim(covered_edges)[1])){
    gn = covered_edge_reversal(
      CIFUN,groot,covered_edges[i,1],covered_edges[i,2],h)
    if (prune){
      gid=paste(as.vector(which(gn!=0)), collapse="")
      if(sum(gn!=0)==0){gid="empty"}
      if (exists(gid,gh)){
        if(verbose){cat("-")}
        next()
      } else {
        if(verbose){cat("+")}
        gh[[gid]]=1
      }
    } else if (verbose){cat("+")}
    if (sum(gn!=0)<ecount){return(list(gn))}
    gnext[[j]]=gn;j=j+1
  }
  return(gnext)
}

#' Modifed Sparsest permutation (with DNA and pruning)
#'
#' by finding covered edges and reverting them. no sugar.
#' @param CIFUN: CI function
#' @param groot: root DAG
#' @param POM: partial ordering
#' @param layers: list of layers
#' @param ecount: current number of edges
#' @param h: dictionary of CI statements
#' @param gh: dictionary of DAGs
#' @param prune: if yes, skip DAG in gh that has been checked. Otherwise not.
#' @param convervative: censor by POM (aggressive) or layers (conservative) when finding covered edges
#' @return list of DAGs
#'
SP<-function(CIFUN, groot,POM,IO,d,k,h=NULL,gh=NULL,prune=T,verbose=F,reroot=F,conservative=F){
  # Depth-first Search: recursive function
  # apply next_graph until reach leave or reach d-th level
  # k: number of moves
  stopifnot(is.function(CIFUN))
  if (reroot){return(list(g=groot,k=k,reroot=T))}
  if (verbose){cat("(",sum(groot))}
  if (is.null(gh)){gh=new.env(hash=T)}
  if (d==0){
    if (verbose){cat("Leaf",")")}
    return(list(g=groot,k=k,reroot=F))
  } # leaf
  if(!is.null(IO)){
    layers = sapply(seq(dim(groot)[1]), function(x)which(sapply(IO, function(y) x %in% y)))
  } else {
    layers = rep(1,dim(groot)[1])
  }
  ecount=sum(groot)
  gnext = next_graph(CIFUN,groot,POM,layers,ecount,h,gh,prune,verbose,conservative) # weakly decreasing g
  k=k+length(gnext)
  gnext = Filter(Negate(is.null), gnext)
  if (length(gnext)==0){
    if (verbose){cat("lv",d,"Deadend",")")}
    return(list(g=groot,k=k,reroot=F))
  }
  gnextecount = sapply(gnext,function(x)sum(x))
  if (min(gnextecount)<ecount){
    if (verbose){cat("REROOT",")")}
    return(list(g=gnext[[which.min(gnextecount)]],k=k,reroot=T))
  }
  ngs=NULL;
  ecount=min(gnextecount)
  for (x in seq(length(gnext))){
    ngs[[x]] = SP(CIFUN,gnext[[x]],POM,IO,d-1,k,h,gh,prune,verbose,reroot=F,conservative)
    ecount=min(ecount,sum(ngs[[x]]$g!=0))
    if (ngs[[x]]$reroot==T|sum(ngs[[x]]$g!=0)<ecount){return(list(g=ngs[[x]]$g,k=ngs[[x]]$k,reroot=T))}
  }
  ngsecount = sapply(ngs, function(x)sum(x$g!=0))
  k = k+ sum(sapply(ngs, function(x)x$k))
  if (min(ngsecount)<ecount){
    gmin = ngs[[which.min(ngsecount)]]$g
  } else {
    gmin = groot
  }
  if (verbose){cat("lv",d,"Branch",sum(gmin),")\n")}
  return(list(g=gmin,k=k,reroot=F))
}

#' Modifed Sparsest permutation (incorporating DNA and pruning)
#'
#' by finding covered edges and reverting them. sugar wrapped version. Use gaussCItest
#' @param X: data matrix
#' @param alpha: CI test level
#' @param DNA: use DNA or not. If FALSE, just run SP. If TRUE, learn DNA first
#' @param k: if DNA, levels of PC used in learning step
#' @param r: number of restarts
#' @param d: depth of search
#' @param verbose: whether to print the learning path.
#' @return CPDAG
#'
SparsestPermutation<-function(X,DNA=F,k=0,alpha,r,d,verbose=F){
  h=new.env(hash = T);gh=new.env(hash=T)
  if (!DNA){POM==matrix(2,p,p);IO=NULL}
  if (DNA){
    dat=NULL;dat$X=X
    Dout=LearnDNAforward(dat,h,k,CIFUN_DNA = CIFUN,CIFUN_skel = CIFUN)
    D=Dout$D
    DI = orderConstraining(D)
    IO=c(DI$sources,DI$sinks)
    POM=D
  }
  n=dim(X)[1]
  p=dim(X)[2]
  Sig = cov(X)
  CIFUN<-function(x,y,S,suffStat=NULL){
    gaussCItest(x,y,S,list(C=Sig,n=n))>alpha} # threshold
  k=0
  spout=matrix(1,p,p)
  for (i in seq(r)){
    set.seed(i)
    TOinit=sample(p)
    groot = DAG_from_TO(CIFUN,TOinit,h)
    gg=SP(CIFUN,groot,D,IO,d,0,h,gh,
          prune=prune,verbose = verbose,reroot=F,conservative=T)
    while (gg$reroot){
      gg=SP(CIFUN_lvs,gg$g,D,IO,d,gg$k,h,gh,
            prune=prune,verbose = verbose,reroot=F,conservative=T)
    }
    k=k+gg$k;
    if (sum(gg$g!=0)<sum(spout!=0)){spout=gg$g}
    if (verbose){
      cat("|G_pi|=",sum(gg$g),"steps taken=",gg$k,"\n")}
  }
  return(dag2cpdag(spout))
}

