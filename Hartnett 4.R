## Question 1 ##
load('arsenic.rdata')

## We need a numeric value in order to fit this by hand
## at least I assume so. idk how R interprets factor levels
## Also this defines the interaction effect
arsenic$gennum = ifelse(arsenic$gender == 'Female', 1, 0)



arsenic$toenail = log(arsenic$arsenic.toenail)
arsenic$int = arsenic$gennum*arsenic$arsenic.water

X = as.matrix(cbind(1, arsenic[,c(2,3,5,7)]))
Y = as.matrix(arsenic[,6])

n=length(Y)
p=dim(X)[2]

XtXinv = qr.solve(crossprod(X))

(beta.hat = XtXinv%*%crossprod(X,Y))
(sE=sqrt(sum((Y-X%*%beta.hat)^2 )/(n-p)))

## To find the best AIC, I'm going to call the function I defined for the last HW
## It also returns BIC and the cross Validation, but just ignore those
olscv = function(X, y, K){
  n=dim(X)[1]
  ind=sample(n)  
  err=numeric(K)
  for(k in 1:K)
  {
    leave.out=ind[(1+floor((k-1)*n/K)):floor(k*n/K)]
    y.train=y[-leave.out]
    x.train=X[-leave.out,]
    XtXinv.train=qr.solve(crossprod(x.train))
    beta.hat.train=XtXinv.train%*%crossprod(x.train,y.train)
    y.val=y[leave.out,]
    x.val=X[leave.out,]
    
    err[k]=(sum((x.val%*%beta.hat.train-y.val)^2))
  }
  av.err=sum(err)/n
  return(av.err)
}
get.ic.cv = function(X, y, K=5){
  n=dim(X)[1]
  p=dim(X)[2]
  beta.hat=lm.fit(x=X, y=y)$coefficients
  rss=sum((y-X%*%beta.hat)^2)
  common=n*log(2*pi)+n*log(rss/n) + n
  aic=common + 2*(p+1)
  bic=common + (p+1)*log(n)
  cv=olscv(X=X, y=y, K=K)
  return(c(aic, bic, cv))
}

## For the three main effects there are 8 possible submodels to explore
X1 = as.matrix(X[,1])
X2 = as.matrix(X[,c(1,2)])
X3 = as.matrix(X[,c(1,3)])
X4 = as.matrix(X[,c(1,4)])
X5 = as.matrix(X[,c(1,2,3)])
X6 = as.matrix(X[,c(1,2,4)])
X7 = as.matrix(X[,c(1,3,4)]) 
X8 = as.matrix(X[,c(1,2,3,4)])

## Getting the aics for each model
## I'm sure there is a way to put this in a for loop but I can't get 
## a list of matrices to correctly pass through to the function
get.ic.cv(X1, Y)
get.ic.cv(X2, Y)
get.ic.cv(X3, Y)
get.ic.cv(X4, Y)
get.ic.cv(X5, Y)
get.ic.cv(X6, Y)
get.ic.cv(X7, Y)
get.ic.cv(X8, Y)

## Question 2 ## 
mustar = .5
reps = 1e4
n=10

xbar.list = numeric(reps)
s.list = numeric(reps)
for(i in 1:reps){
  x.list = rnorm(n, mustar, mustar)
  xbar.list[i]=mean(x.list)
  s.list[i]=sd(x.list)^2
}
cov(xbar.list, s.list)

psE = sqrt(((reps-1)*var(xbar.list)+(reps-1)*var(s.list))/(reps*2 - 2))




## Realizing the data
n=10
mustar=.5
x.list=rnorm(n, mustar, mustar)
mu.seq=seq(from = 0.05, to = (2*sum((x.list)^2)/n), length.out = 1e3)

## Computing the negative loglikelihood

for(i in 1:1e3){
  f.list[i] = ((n/2)*log(2*pi) + (n/2)*log(mu.seq[i]) + ((1/2)/mu.seq[i])*sum((x.list-mu.seq[i])^2))
}



plot(mu.seq, f.list, type = 'l')


sum((x.list-mean(x.list))^2)
10*mean(x.list)^2 - 2 *mean(x.list)*sum(x.list)+sum(x.list^2)

## part g
ci.mat = matrix(NA, nrow = 5, ncol = 2)

n=10
mustar=1e-2
reps=1e4
results.mat = matrix(NA, nrow = reps, ncol = 5)

for (i in 1:reps) {
  ## Realizing the data and finding the relevant stats
  x.list = rnorm(n, mustar, mustar)
  xbar = mean(x.list)
  s = var(x.list)
  muhat = (-1+sqrt(1+4*(sum(x.list^2))/n))/2
  
  ## Finding the requested values
  results.mat[i,1] = abs(xbar-mustar)
  results.mat[i,2] = abs(s-mustar)
  results.mat[i,3] = abs(muhat-mustar)
  results.mat[i,4] = abs(xbar-mustar)-abs(muhat-mustar)
  results.mat[i,5] = abs(xbar-mustar)-abs(muhat-mustar)
}

apply(results.mat, 2, mean)


n=10
mustar=1e-1
reps=1e4
results.mat = matrix(NA, nrow = reps, ncol = 5)

for (i in 1:reps) {
  ## Realizing the data and finding the relevant stats
  x.list = rnorm(n, mustar, mustar)
  xbar = mean(x.list)
  s = var(x.list)
  muhat = (-1+sqrt(1+4*(sum(x.list^2))/n))/2
  
  ## Finding the requested values
  results.mat[i,1] = abs(xbar-mustar)
  results.mat[i,2] = abs(s-mustar)
  results.mat[i,3] = abs(muhat-mustar)
  results.mat[i,4] = abs(xbar-mustar)-abs(muhat-mustar)
  results.mat[i,5] = abs(xbar-mustar)-abs(muhat-mustar)
}

apply(results.mat, 2, mean)


n=10
mustar=1e0
reps=1e4
results.mat = matrix(NA, nrow = reps, ncol = 5)

for (i in 1:reps) {
  ## Realizing the data and finding the relevant stats
  x.list = rnorm(n, mustar, mustar)
  xbar = mean(x.list)
  s = var(x.list)
  muhat = (-1+sqrt(1+4*(sum(x.list^2))/n))/2
  
  ## Finding the requested values
  results.mat[i,1] = abs(xbar-mustar)
  results.mat[i,2] = abs(s-mustar)
  results.mat[i,3] = abs(muhat-mustar)
  results.mat[i,4] = abs(xbar-mustar)-abs(muhat-mustar)
  results.mat[i,5] = abs(xbar-mustar)-abs(muhat-mustar)
}

apply(results.mat, 2, mean)


n=10
mustar=1e1
reps=1e4
results.mat = matrix(NA, nrow = reps, ncol = 5)

for (i in 1:reps) {
  ## Realizing the data and finding the relevant stats
  x.list = rnorm(n, mustar, mustar)
  xbar = mean(x.list)
  s = var(x.list)
  muhat = (-1+sqrt(1+4*(sum(x.list^2))/n))/2
  
  ## Finding the requested values
  results.mat[i,1] = abs(xbar-mustar)
  results.mat[i,2] = abs(s-mustar)
  results.mat[i,3] = abs(muhat-mustar)
  results.mat[i,4] = abs(xbar-mustar)-abs(muhat-mustar)
  results.mat[i,5] = abs(xbar-mustar)-abs(muhat-mustar)
}

apply(results.mat, 2, mean)



n=10
mustar=1e2
reps=1e4
results.mat = matrix(NA, nrow = reps, ncol = 5)

for (i in 1:reps) {
  ## Realizing the data and finding the relevant stats
  x.list = rnorm(n, mustar, mustar)
  xbar = mean(x.list)
  s = var(x.list)
  muhat = (-1+sqrt(1+4*(sum(x.list^2))/n))/2
  
  ## Finding the requested values
  results.mat[i,1] = abs(xbar-mustar)
  results.mat[i,2] = abs(s-mustar)
  results.mat[i,3] = abs(muhat-mustar)
  results.mat[i,4] = abs(xbar-mustar)-abs(muhat-mustar)
  results.mat[i,5] = abs(xbar-mustar)-abs(muhat-mustar)
}

apply(results.mat, 2, mean)







