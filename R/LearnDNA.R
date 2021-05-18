# Copyright (c) 2018 - 2019  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.


#' Learn DNA using early stopped PC
#'
#' Run PC for k steps, stop and learn the DNA
#' @param dat: To use sample version, let dat$X hold the data matrix.To use population version, let dat$Sigma hold population cov
#' @param k: level
#' @param h: dictionary for checked CI statements
#' @param CIFUN_skel: CI funtion used for skeleton learning. CIFUN(x,y,S) returns 0 (dependence) or 1 (independence)
#' @param CIFUN_DNA: CI funtion used for DNA learning CIFUN(x,y,S) returns 0 (dependence) or 1 (independence)
#' @return POM: a matrix of partial ordering.
LearnDNAforward<-function(dat,h,k,CIFUN_skel,CIFUN_DNA){
  n=dim(dat$X)[1]
  p=dim(dat$X)[2]
  # initial
  adj=matrix(1,p,p);diag(adj)=0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      adj[i,j]=adj[j,i]=!CIFUN_skel(i,j,NULL)
      if (!is.null(h)){h[[paste(c(sort(c(i,j)),"|"),collapse = ",")]]=1}
    }
  }
  D=matrix(2,p,p);D[adj==0]=0
  SepSets=lapply(1:p, function(x)lapply(1:p,function(y){NULL}))
  if (k>0){
    for (l in 1:k){
      if (sum(adj!=0)==0){break}
      pcout=PC_smallset(p,adj,k,h,SepSets,CIFUN=CIFUN_skel)
      adj=pcout$graph
      SepSets=pcout$SepSets
      if (pcout$done){break}
    }
  }
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if (adj[i,j]==0){
        S=SepSets[[i]][[j]]
        for (l in setdiff(seq(p),c(i,j,S))){
          if (D[l,i]!=0 | D[l,j]!=0){
            if (!is.null(h)){h[[paste(c(sort(c(i,j)),"|",
                                        sort(c(S,l))),collapse = ",")]]=1}
            ci = CIFUN_DNA(i,j,c(S,l))
            if (!ci){
              D[l,i]=D[l,j]=0
              D[l,S]=0
            }
          }
        }
      }
    }
  }
  return(D)
}

#' Learn DNA using early iterative CLIME
#'
#' Run CLIME for k steps, use moral graphs to learn DNA
#' @param dat: To use sample version, let dat$X hold the data matrix.To use population version, let dat$Sigma hold population cov
#' @param k: level
#' @param h: dictionary for checked CI statements
#' @param CItest: if "oracle", then use true Previson matices. Otherwise use CLIME
#' @return POM: a matrix of partial ordering.
LearnDNAbackward<-function(dat,h,k,CItest="oracle"){
  n=dim(dat$X)[1]
  p=dim(dat$X)[2]
  if (CItest=="oracle"){
    Theta=dat$Theta
    adj = dat$Theta!=0;diag(adj)=0
  } else {
    Theta = cv_clime(dat$X)
    lambda=Theta$lambda
    Theta=Theta$Omega
    adj=Theta!=0;diag(adj)=0
  }
  D=matrix(2,p,p);diag(D)=0
  D[dat$Sigma==0]=0
  sinks=sinksnew=NULL
  for (l in 1:k){
    for (i in 1:(p-length(sinks))){
      coli = setdiff(seq(p),sinks)[i]
      if (CItest=="oracle"){
        el=which(solve(dat$Sigma[-c(sinks,coli),
                                 -c(sinks,coli)])==0 &
                   Theta[-c(sinks,coli),
                         -c(sinks,coli)]!=0,arr.ind = T)
      } else {
        el= which(clime_theta(dat$X[,-c(sinks,coli)],
                              lam = lambda)==0 &
                    Theta[-c(sinks,coli),
                          -c(sinks,coli)]!=0, arr.ind = T)
      }
      if (dim(el)[1]==0){next}
      el[,1]=(setdiff(seq(p),c(sinks,coli)))[el[,1]]
      el[,2]=(setdiff(seq(p),c(sinks,coli)))[el[,2]]
      if (!all(c(D[coli,unique(as.vector(el))],
                 D[unique(as.vector(el)),coli])==0)){
        sinksnew=c(sinksnew,i)
      }
      D[coli,unique(as.vector(el))]=D[unique(as.vector(el)),coli]=0
    }
    if (!is.null(sinksnew)){
      Theta=matrix(0,p,p)
      sinks=sinksnew
      if (CItest=="oracle"){
        Theta[-sinks,-sinks]=solve(dat$Sigma[-sinks,-sinks])
      } else {
        Theta[-sinks,-sinks]=clime_theta(dat$X[,-sinks],lambda)
      }
    }
  }
  return(D)
}

#' Learn DNA using early stopped PC (special case: 2 levels, threshoding partial correlations)
#'
#' Run PC for k steps, stop and learn the DNA
#' @param dat: To use sample version, let dat$X hold the data matrix.To use population version, let dat$Sigma hold population cov
#' @param ord: 0 or 1. Early stop pc at 0 or 1. not implemented for higher
#' @param co: thresholds for (1: skel at lv0; 2: DNA at lv0; 3: skel at lv1; 4: DNA at lv1)
#' @return POM: a matrix of partial ordering.
LearnDNAforwardthreshold<-function(dat,ord,
                            co1=0.01,co2=0.1,co3=0.01,
                            co4=0.2,oracle=F,verbose=F){
  p=dim(dat$X)[2]
  n=dim(dat$X)[1]
  C=cov2cor(cov(dat$X))
  if (oracle){C=cov2cor(dat$Sigma)}
  D=matrix(2,p,p);diag(D)=0
  adj=matrix(1,p,p);diag(adj)=0
  L=matrix(0,p,p)
  ntests=0;ntests0=0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      ntests0=ntests0+1
      if (abs(C[i,j])<co1){
        adj[i,j]=adj[j,i]=0
        D[i,j]=D[j,i]=0
        if(verbose){cat("DNA:nonadj",i,j,"\n",sep=",")}
        for (k in setdiff(seq(p),c(i,j))){
          ntests0=ntests0+1
          if (abs(pcorOrder(i,j,k,C))>co2){
            D[k,i]=D[k,j]=0
            if(verbose){cat("DNA0:",i,j,k,"\n",sep=",")}
          }
        }
      }
    }
  }
  if(ord==0){return(list(D=D,ntests=ntests0))}
  for (i in 1:(p-1)){
    for (j in (i+1):p){
        for (k in setdiff(which(adj[i,]+adj[j,]!=0),c(i,j))){
          ntests=ntests+1
          if (abs(pcorOrder(i,j,k,C))<co3){
            adj[j,i]=adj[i,j]=0
            for (z in setdiff(seq(p),c(i,j,k))){
              ntests=ntests+1
              if (abs(pcorOrder(i,j,c(k,z),C))>co4){
                D[z,i]=D[z,j]=D[z,k]=0
                if(verbose){cat("DNA1:",i,j,k,z,"\n",sep=",")}
              }
            }
          }
        }
    }
  }
  ntests=ntests+ntests0
  return(list(D=D,ntests=ntests))
}
