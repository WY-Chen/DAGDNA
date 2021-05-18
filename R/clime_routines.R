unit_vec<-function(p,i){v=rep(0,p);v[i]=1;return(v)}
clime_lp<-function(Sigma,lam,j){
  # Sigma: cov(X)
  # lam: tuning parameter
  # j: the j-th problem
  p = dim(Sigma)[2]
  f.obj = rep(1,p*2) # sum (u+v), u=max(x,0), v=max(-x,0), x=u-v
  const.mat = rbind(
    cbind(Sigma,-Sigma), # Sigma*(u-v) >= lam +ej
    cbind(-Sigma,Sigma), # -Sigma*(u-v) >= lam-ej
    cbind(diag(p),matrix(0,p,p)), # u>0
    cbind(matrix(0,p,p),diag(p))  # v>0
  )
  const.dir = c(
    rep("<=",2*p),rep(">=",2*p)
  )
  const.rhs = c(
    rep(lam,p)+unit_vec(p,j),
    rep(lam,p)-unit_vec(p,j),
    rep(0,2*p)
  )
  lpout=lpSolve::lp(direction = "min",objective.in = f.obj,
           const.mat = const.mat,const.dir = const.dir,
           const.rhs = const.rhs)
  return(lpout$solution[1:p]-lpout$solution[(p+1):(2*p)])
}

clime_theta<-function(X,lam=NULL){
  p=dim(X)[2]
  n=dim(X)[1]
  Sigma = cov(X)
  if (is.null(lam)){lam = 4/sqrt(n)*sqrt(log(p/sqrt(0.05)))}
  Omega = sapply(1:p, function(i)clime_lp(Sigma,lam,i))
  Omega = (abs(Omega)<=abs(t(Omega)))*Omega+
    (abs(Omega)>abs(t(Omega)))*t(Omega)
  return(Omega)
}

cv_clime<-function(X){
  p=dim(X)[2]
  n=dim(X)[1]
  Sigma = cov(X)
  ws = sqrt(diag(Sigma))
  lams=exp(seq(log(1e-4),log(0.8),length.out = 100))
  ebics=sapply(1:100, function(j){
    Theta = sapply(1:p, function(i)clime_lp(Sigma,lams[j],i))
    Theta = (abs(Theta)<=abs(t(Theta)))*Theta+
      (abs(Theta)>abs(t(Theta)))*t(Theta)
    loglikGGM(S=Sigma,Theta=Theta)-
      sum(Theta!=0)/2*(log(p)+0.5*log(n))/n
  })
  return(list(Omega=clime_theta(X,lam = lams[which.max(ebics)]),
              lambda = lams[which.max(ebics)]))
}
