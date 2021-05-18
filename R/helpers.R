#' Compute log likelihood for Gaussian Distr
#'
#' Use either covariance, prevision, or SEM weights and error cov
#' @param S: Sample Covariance matrix
#' @param Theta: Precision model
#' @param Sigma: Covariance model
#' @param A: Weights in linear SEM
#' @param Omega: Error cov in linear SEM
#' @return loglikelihood
loglikGGM<-function(S,Theta=NULL,Sigma=NULL,A=NULL,Omega=NULL){
  # compute the
  if (!is.null(Sigma)){
    -log(det(Sigma))-sum(diag(solve(Sigma)%*%S))
  } else if(!is.null(Theta)){
    log(det(Theta))-sum(diag(Theta%*%S))
  } else if(!is.null(A)& !is.null(Omega)){
    p = dim(A)[1]
    Theta=(diag(p)-A)%*%solve(Omega)%*%t(diag(p)-A)
    log(det(Theta))-sum(diag(S%*% Theta ))
  }
}

#' PC algorithm at level k
#'
#' this is standard. modified from pcalg package
#' @param p: number of variables
#' @param ghat: working graph estimate (undirected skeleton)
#' @param k: level
#' @param h: dictionary for checked CI statements
#' @param SepSets: working Separators
#' @param CIFUN: CI funtion CIFUN(x,y,S) returns 0 (dependence) or 1 (independence)
#' @param POM: partial ordering matrix
#' @param IO: incomplete ordering
#' @return loglikelihood
PC_smallset <-function(p,ghat,k,h,SepSets=NULL,CIFUN,POM=NULL,IO=NULL){
  # k-th step of PC
  ci <- NULL
  ord <- k
  seq_p=seq(p)
  if(is.null(POM)){POM<-matrix(2,p,p)}
  if(!is.null(IO)){
    layers = sapply(seq(p), function(x)which(sapply(IO, function(y) x %in% y)))
  } else {
    layers = rep(1,p)
  }
  done<-TRUE
  if (is.null(SepSets)){
    SepSets<-lapply(seq(p), function(i){
      lapply(seq(p),function(j){NULL})})
  }
  ind <- which(ghat!=0, arr.ind = TRUE)
  ind <- ind[order(ind[, 1]), ]
  remEdges <- nrow(ind)
  if (remEdges==0){return(list(graph=ghat,SepSets=SepSets,done=TRUE))}
  ghat.l <- split(ghat, gl(p, p))
  for (i in 1:remEdges) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    if (ghat[y, x]!=0) {
      nbrsBool <-   ghat.l[[x]]!=0
      nbrsBool[y] <- FALSE
      nbrs <- seq_p[nbrsBool]
      # apply DNA rules <- do not using rule 3
      if (POM[y,x]==0 && POM[x,y]!=0){
        next
      } else if (POM[y,x]!=0 && POM[x,y]==0){
        nbrs = intersect(nbrs,which(POM[,x]!=0))
      }
      # apply IO rules
      nbrs = intersect(nbrs,which(layers<=max(layers[x],layers[y])))
      # PC steps
      length_nbrs <- length(nbrs)
      if (length_nbrs >= ord) {
        if (length_nbrs > ord) {done<-FALSE}
        S <- seq_len(ord)
        repeat {
          if (!is.null(h)){
            h[[paste(c(sort(c(x,y)),"|",sort(nbrs[S])),collapse = ",")]]=1}
          ci <- CIFUN(x,y,nbrs[S],NULL)
          if (ci) {
            ghat[x, y] <- ghat[y, x] <- FALSE
            SepSets[[x]][[y]]<-SepSets[[y]][[x]]<-nbrs[S]
            break
          }
          else {
            nextSet <- getNextSet(length_nbrs, ord,S)
            if (nextSet$wasLast)
              break
            S <- nextSet$nextSet
          }
        }
      }
    }
  }
  return(list(graph=ghat,SepSets=SepSets,done=done))
}

#' check if DNA set is compatible with DAG
#'
#' using cpdag characterization
#' @param gm: true DAG matrix
#' @return POM: DNA set
#'
is_valid_DNA<-function(gm,POM){
  trueDNA=all_DNA(gm)
  return(!sum(POM==0&trueDNA!=0))
}

#' check if ordering is compatible with DAG
#'
#' using triangularization
#' @param C: true DAG matrix
#' @return TO: ordering
#'
is_compatible_ordering<-function(C,TO){
  #  false if i>j and C_{ij}=1
  return(all(C[TO,TO]*lower.tri(C)!=1))
}


#' check if layer is compatible with DAG
#'
#' using triangularization
#' @param gm: true DAG matrix
#' @return IO: layering
#'
is_valid_IO<-function(gm,IO){
  l = length(IO)
  if (l==1){return(TRUE)}
  for (i in 1:(l-1)){
    if (!all(gm[unlist(IO[(i+1):l]),unlist(IO[seq(i)])]==0)){return(FALSE)}
  }
  return(TRUE)
}
