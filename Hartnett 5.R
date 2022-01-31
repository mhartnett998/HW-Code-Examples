## These are the functions from the lecture code
myrnorm=function(n)
{
  ## is n odd?
  odd=(n %% 2)!=0 
  if(odd) 
  {
    ## add 1 to n 
    n=n+1
  }
  ##  perform the Box-Muller method
  u1=runif(n/2)
  u2=runif(n/2)
  tmp=sqrt(-2*log(u1))
  x.list=c(tmp*cos(2*pi*u2), tmp*sin(2*pi*u2) )
  if(odd)
  {
    ## if n was initially odd, remove an entry
    x.list=x.list[-n]
  }
  return(x.list)
}
mymvrnorm=function(n, mu, Sigma)
{
  p=length(mu)
  ## Compute the square root of Sigma
  ei=eigen(Sigma, symmetric=TRUE)
  Sigma.sqrt=ei$vec %*% diag(ei$val^0.5)%*%t(ei$vec)
  
  ## generate a realization of n by p 
  ## random matrix with entries iid N(0,1)
  Z = matrix( myrnorm(n*p), nrow=n, ncol=p)
  
  ## compute X using the formula
  X=rep(1, n) %*% t(mu) + Z %*% Sigma.sqrt
  return(X)
}

get.std.mat=function(X)
{
  mx=apply(X[,-1], 2, mean)
  sx=apply(X[,-1], 2, sd)
  A=diag(ncol(X))
  A[1,-1] = -mx
  B=diag(c(1, 1/sx))
  return(A%*%B)
}

ridge=function(X,y, lam, std=TRUE, use.svd=(ncol(X) >= nrow(X)))
{
  p=ncol(X)
  if(std)
  {
    T=get.std.mat(X)
    X=X%*%T
  }
  if(!use.svd)
  {
    M=diag(c(0, rep(1,p-1)))
    bhat=qr.solve(crossprod(X)+lam*M, crossprod(X,y))
    if(std) bhat=T%*%bhat
  } else
  {
    n=nrow(X)
    xbart=crossprod(rep(1,n), X[,-1])/n
    ybar=mean(y)
    Xc= X[,-1] - rep(1,n)%*%xbart
    yc=y-ybar
    
    ## compute the reduced SVD of Xc
    q=min(c(n-1, p-1))
    out=svd(Xc, nu=q, nv=q)
    H=diag(out$d[1:q]/(out$d[1:q]^2 + lam))
    bhatm1=out$v%*%H%*%t(out$u)%*%yc
    bhat1=ybar - (xbart%*%bhatm1)[1]
    bhat=c(bhat1, bhatm1)
    if(std) bhat=T%*%bhat
  }
  return(bhat)
}

ridgecv=function(X, y, K=5, permute=FALSE, std=TRUE, 
                 lam.vec=10^seq(from=-8, to=8, by=0.5))
{
  n=length(y)
  p=ncol(X)
  if(permute) ind=sample(n) else ind=1:n
  val.sq.err=numeric(length(lam.vec))
  for(k in 1:K)
  {
    leave.out=ind[ (1+floor((k-1)*n/K)):floor(k*n/K) ]  
    X.tr=X[-leave.out,,drop=FALSE];  y.tr=y[-leave.out]
    X.va=X[leave.out,,drop=FALSE];   y.va=y[leave.out]
    
    if(std)
    {
      T=get.std.mat(X.tr)
      X.tr=X.tr%*%T
    }
    n.tr=length(y.tr)
    ## compute the centered training predictor matrix
    xbart=crossprod(rep(1,n.tr), X.tr[,-1])/n.tr
    ybar=mean(y.tr)
    Xc= X.tr[,-1] - rep(1,n.tr)%*%xbart
    ## compute the centered training response vector
    yc=y.tr-ybar
    
    ## compute the reduced SVD of Xc
    q=min(c(n.tr-1, p-1))
    out=svd(Xc, nu=q, nv=q)
    for(j in 1:length(lam.vec))
    {
      ## compute beta.hat^(lam.vec[j], -k)
      H=diag(out$d[1:q]/(out$d[1:q]^2 + lam.vec[j]))
      bhatm1=out$v%*%H%*%t(out$u)%*%yc
      bhat1=ybar - (xbart%*%bhatm1)[1]
      bhat.tr=c(bhat1, bhatm1)
      
      if(std) bhat.tr=T%*%bhat.tr
      ## compute its validation error
      val.sq.err[j] = val.sq.err[j] + sum((y.va - X.va%*%bhat.tr)^2 )
    }
  }
  best.lam=lam.vec[which.min(val.sq.err)]
  
  ## compute the ridge estimator with the selected tuning parameter
  ## using all the data 
  beta.hat=ridge(X=X,y=y, lam=best.lam, std=std, use.svd=TRUE)
  return(list(best.lam=best.lam, beta.hat=beta.hat, val.sq.err=val.sq.err))
}

ridge.blr=function(X, y, n.list, lam, m=NULL, tol=1e-7, 
                   maxit=100, quiet=FALSE, b.start=NULL)
{
  p=ncol(X)
  ## create vector of penalty weights
  ## if unspecified
  if(is.null(m))
    m=c(0, rep(1, p-1))
  
  ## create useful variables  
  lam.m=lam*m
  lam.M=diag(lam.m)
  X.t.n.list.y=crossprod(X, n.list*y)
  
  ## If unspecified, make the 0th iterate
  ## for b the zero vector
  if(is.null(b.start)) 
  {  
    b=rep(0,p)
  } else
  {
    b=b.start 
  }
  ## initialize iteration counter
  k=0
  add=tol+1
  while( (k <= maxit) & (sum(abs(add)) > tol))
  { 
    k=k+1
    pi.t=ilogit(as.numeric(X%*%b))
    W=diag(n.list*pi.t*(1-pi.t))
    minusGrad=X.t.n.list.y-crossprod(X, n.list*pi.t) - lam.m*b
    Hess=crossprod(X,W%*%X)+lam.M
    add=qr.solve(Hess, minusGrad)
    b=b+add
    if(!quiet) cat("k=", k, "b=", b, "\n")
  }
  b=as.numeric(b)
  return(list(b=b, total.iterations=k))
}
###################################################################################
## Question 1 ## 
n=100
p=50
sigma.star=1
theta=.95
reps=200

sigma.mat = diag(p-1)
for(i in 1:(p-1)){
  for(j in 1:(p-1)){
    sigma.mat[i,j] = theta^(abs(i-j))
  }
}

beta.star = rnorm(n=p, mean=0, sd = sqrt(p^-1))

## The mymvrnorm function is defined from myMCfunctions.r
X = mymvrnorm(n=n, mu=rep(0,p-1), Sigma=sigma.mat)
X = cbind(1, X)

loss.mat=matrix(NA, nrow = reps, ncol = 6)
err.mat=matrix(NA, nrow = reps ,ncol = 6)
for(i in 1:reps){
  ## We start each replication by generating a new realization of Y
  y = X%*%beta.star + rnorm(n=n, mean=0, sd=sigma.star)
  
  ## Computing the realizations of the different b
  beta.ols=qr.coef(qr(X),y=y)
  beta.cv5 = ridgecv(y=y, X=X, K=5)$beta.hat
  beta.cv10 = ridgecv(y=y, X=X, K=10)$beta.hat
  
  ## Storing the realizations of the 6 losses in a matrix
  loss.mat[i,1] = (sum((beta.ols-beta.star)^2))
  loss.mat[i,2] = (sum((X%*%beta.ols-X%*%beta.star)^2))
  loss.mat[i,3] = (sum((beta.cv5-beta.star)^2))
  loss.mat[i,4] = (sum((X%*%beta.cv5-X%*%beta.star)^2))
  loss.mat[i,5] = (sum((beta.cv10-beta.star)^2))
  loss.mat[i,6] = (sum((X%*%beta.cv10-X%*%beta.star)^2))
  
  ## Storing the realizations of the standard errors
  err.mat[i,1] = (mean((beta.ols-beta.star)^2))
  err.mat[i,2] = (mean((X%*%beta.ols-X%*%beta.star)^2))
  err.mat[i,3] = (mean((beta.cv5-beta.star)^2))
  err.mat[i,4] = (mean((X%*%beta.cv5-X%*%beta.star)^2))
  err.mat[i,5] = (mean((beta.cv10-beta.star)^2))
  err.mat[i,6] = (mean((X%*%beta.cv10-X%*%beta.star)^2))
}
apply(loss.mat, 2, mean)
apply(err.mat, 2, mean)

#######################################################################
## Question 2 ##
ilogit = function(u) return(exp(u)/(1+exp(u)))

n=10
p=5
## create the data matrix
X=cbind(1, matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1)))
## create the true beta
beta.star=rnorm(p, mean=0, sd=1/sqrt(p))
## have all experiments have an index of 30
n.list=rep(30, n)
## create the vector of success probabilities
pi.list=ilogit(as.numeric(X%*%beta.star))
## create the vector of observed sample proportions
## of success, one for each of the n Binomial experiments
y = rbinom(n=n, size=n.list, prob=pi.list)/n.list
## fit the Binomial logistic regression model
## (with lambda = 2, no intercept penalization)
lam=1e4
t.list=beta.star; t.list[1] = 0

fit=ridge.t(X=X, y=y, n.list=n.list, lam=lam, t.list=t.list)

## check the gradient at the final iterate
t.tilde=t.list; t.tilde[1]=fit$b[1]
pi.hat=ilogit(as.numeric(X%*%fit$b))
Grad=crossprod(X, n.list*(pi.hat - y)) + lam*(fit$b-t.tilde) 
Grad


## Question 3 ##

bridgelr.mm = function(X, y, n.list, lam, alpha, tol=1e-7, maxit=100){
  p = ncol(X)
  
  ## We call the ridge.blr function to generate our starting iterate
  b = ridge.blr(X=X, y=y, n.list=n.list, lam=lam, quiet=TRUE)$b
  
  k=0
  add=tol+1
  while( (k <= maxit) & (sum(abs(add)) > tol) ){
    ## Find the penalties given the current iterate
    m = numeric(p)
    for(i in 2:p){
      m[i] = abs(b[i])^(alpha-2)
    }
    
    ## Here we call ridge.blr to minimize the majorizing function at the
    ## current iterate given the current weights
    bnew = ridge.blr(X=X, y=y, n.list=n.list, lam=lam, m=m, quiet=TRUE)$b
    
    ## Here we check for convergence and repeat if nessecary
    add = bnew-b
    b = bnew
    k = k+1
  }
  return(list(b=b, total.iterations=k))
}

n=10
p=5
## create the data matrix
X=cbind(1, matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1)))
## create the true beta
beta.star=rnorm(p, mean=0, sd=1/sqrt(p))
## have all experiments have an index of 30
n.list=rep(30, n)
## create the vector of success probabilities
pi.list=ilogit(as.numeric(X%*%beta.star))
## create the vector of observed sample proportions
## of success, one for each of the n Binomial experiments
y = rbinom(n=n, size=n.list, prob=pi.list)/n.list
## fit the Binomial logistic regression model
## (with lambda = 2, no intercept penalization)
lam = 5
alpha = 1.1

fit = bridgelr.mm(X, y, n.list, lam, alpha)

m = numeric(p)
for(i in 2:p){
  m[i] = abs(fit$b[i])^(alpha-2)
}

pi.hat=ilogit(as.numeric(X%*%fit$b))
Grad=crossprod(X, n.list*(pi.hat - y)) + lam*fit$b*m 
Grad







