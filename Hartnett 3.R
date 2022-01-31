corn = read.table('corn.txt')

## Question 1 ##
corn$n2 = corn$nitrogen^2

X = matrix(0, nrow = length(corn[,1]), ncol = 3)
Y = corn[,1]


X=as.matrix(cbind(1, corn[,2:3]))

betahat = NULL

## Fitting the model the good old fashioned way 
XtX = crossprod(X)
XtXinv = qr.solve(XtX)
XtY = crossprod(X,Y)

betahat = XtXinv %*% XtY

pred.pts = seq(from = 0, to = 300, length.out = 1e3)

pred.mat=cbind(1, pred.pts, pred.pts^2)
pred.vals = pred.mat%*%betahat

plot(x=corn$nitrogen, y=corn$yield, main='Corn Yield', xlab = 'Nitrogen', ylab='Yield')
lines(x=pred.pts, y=pred.vals, lty=3)




## Question 2 ## 

## n=(10, 50, 200, 500, 1000), alpha=(.01, .05, .1)
n = 10
alpha = .01

p=5
mu = numeric(p-1)
reps = 1e4
beta.star = c(1,1,2,1,-2)
sigma.star = 1
tval = qt(1-(alpha/2), n-5)
xnew = c(1,0,1,0,1)

## We need the mymvrnorm function from previous hws to generate these data
## that have a predefined covariance structure
sigma.mat = matrix(.6, nrow = 4, ncol = 4)
diag(sigma.mat) = 1

captured.list=numeric(reps)
for(i in 1:reps){
  ##Generate the design matrix
  x.mat = mymvrnorm(n=n, mu=mu, Sigma = sigma.mat)
  x.mat = cbind(1, x.mat)
  
  ## Generating realizations of y
  y = x.mat%*%beta.star + rnorm(n, mean = 0, sd = sigma.star)
  
  ## Finding the estimated beta.hat, the ynew, and creating a prediction interval
  ## using the xnew
  xtxinv=qr.solve(crossprod(x.mat))
  beta.hat=xtxinv%*%crossprod(x.mat, y)
  sE=sqrt(sum((y-x.mat%*%beta.hat)^2 )/(n-p))
  moe=tval*sE*sqrt(1+t(xnew)%*%xtxinv%*%xnew)[1]
  ub=(sum(xnew*beta.hat)+moe)
  lb=(sum(xnew*beta.hat)-moe)
  
  enew=rexp(1)-1
  ynew=(xnew%*%beta.star+enew)[1]
  
  ## Checking to see if we our prediction interval captured ynew
captured.list[i] = 1*(lb<=ynew)*(ynew<=ub)
}

mean(captured.list)
prop.test(x=sum(captured.list), n=reps, conf.level=0.99, correct=FALSE)$conf.int


## Question 3 ##
olscv = function(X, y, K){
  ## Getting set up 
  n=dim(X)[1]
  ind=sample(n)  
  err=numeric(K)
  for(k in 1:K)
  {
    ## Thank you for the help with the code :)
    leave.out=ind[(1+floor((k-1)*n/K)):floor(k*n/K)]
    
    
    ## add code to form the training design matrix
    ## and training response vector, and use these
    ## to compute the training estimate of the regression
    ## coefficient vector
    y.train=y[-leave.out]
    x.train=X[-leave.out,]
    XtXinv.train=qr.solve(crossprod(x.train))
    beta.hat.train=XtXinv.train%*%crossprod(x.train,y.train)
    
    ## add code to form the validation design matrix
    ## and validation response vector, and compute the
    ## squared prediction error for each validation case
    ## using the training estimate of the regression
    ## coefficient vector and add these squared prediction 
    ## errors to the grand total
    y.val=y[leave.out,]
    x.val=X[leave.out,]
    
    err[k]=(sum((x.val%*%beta.hat.train-y.val)^2))
  }
  av.err=sum(err)/n
  return(av.err)
}


## Question 4 ## 
beta4=1
rho=.9

p=4
n=10
beta.star=c(1,1,1,beta4)
alpha=.01
d=1
sigma.star=1
sigma.mat=matrix(rho, nrow = 3, ncol = 3)
diag(sigma.mat)=1

nseq=seq(from = 10, to = 200, by = 10)
d=1
reps=5000

pval.list=numeric(reps)
npow.list=numeric(length(nseq))
k=1

for(n in nseq){
  print(k)
  for(i in 1:reps){
    evec=rnorm(n)
    Xtilde=mymvrnorm(n=n, mu=rep(0,p-1), Sigma=sigma.mat)
    X=cbind(1, Xtilde) 
    Y=X%*%beta.star+evec
  
    ## compute the QR decomposition of X
    qrX=qr(x=X)
  
    ## the null model design matrix is
    X0=X[,-4]
  
    ## compute the QR decomposition of X0
    qrX0=qr(x=X0)
  
    beta.hat = qr.coef(qrX, y=Y)
    ## get rss1
    rss1=sum((Y-X%*%beta.hat)^2)
  
    ## get rss0
    beta.hat.0=qr.coef(qrX0, y=Y)
    rss0=sum((Y-X0%*%beta.hat.0)^2)
  
    ## compute the realization of F
    f=((rss0 - rss1)/d )/(rss1/(n-p))
    
    pval.list[i]=1-pf(f, d, n-p)
  }
  npow.list[k]=mean(pval.list < 0.1)
  k=k+1
}
npow.list

## This function returns AIC, BIC, and the mean CV error
get.ic.cv = function(X, y, K){
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

## We just have to change these first three again to get the desired results
`n=50
beta4=1
rho=.9
beta.star=c(1,1,1,beta4)

reps=1e4
K=5

correct.mat = matrix(NA, nrow = reps, ncol = 3)
prob.correct = numeric(3)

for(i in 1:reps){
  evec=rnorm(n)
  Xtilde=mymvrnorm(n=n, mu=rep(0,p-1), Sigma=sigma.mat)
  X=cbind(1, Xtilde) 
  X0=X[,-4]
  y=X%*%beta.star+evec
  full=get.ic.cv(X,y,K)
  null=get.ic.cv(X0,y,K)
  correct.mat[i,]=1*(full<null)
}
(prob.correct=apply(correct.mat, 2, mean))
`

## We essentially just repeat the above procedure but 
## add in an estimate for the type 1 error
n=100
beta4=0
rho=.3

p=4
beta.star=c(1,1,1,beta4)
d=1
sigma.mat=matrix(rho, nrow = 3, ncol = 3)
diag(sigma.mat)=1
alpha=.1
reps=1e4
K=5

correct.mat = matrix(NA, nrow = reps, ncol = 4)
for(i in 1:reps){
  ## Same as above for the AIC, BIC, and CV
  evec=rnorm(n)
  Xtilde=mymvrnorm(n=n, mu=rep(0,p-1), Sigma=sigma.mat)
  X=cbind(1, Xtilde) 
  X0=X[,-4]
  y=X%*%beta.star+evec
  full=get.ic.cv(X,y,K)
  null=get.ic.cv(X0,y,K)
  correct.mat[i,]=c(1*(full<null), NA)
  
  ## Now we get our estimated type 1 error
  ## This procedure is from question 4(a)i
  qrX=qr(x=X)
  X0=X[,-4]
  qrX0=qr(x=X0)
  beta.hat = qr.coef(qrX, y=y)
  ## get rss1
  rss1=sum((y-X%*%beta.hat)^2)
  
  ## get rss0
  beta.hat.0=qr.coef(qrX0, y=y)
  rss0=sum((y-X0%*%beta.hat.0)^2)
  
  ## compute the realization of F
  f=((rss0 - rss1)/d )/(rss1/(n-p))
  correct.mat[i,4]=1*(1-pf(f, d, n-p)< alpha)
}
(prob.correct=apply(correct.mat, 2, mean))


## To solve part (b), we can reuse much of the code from part (a)i
## we just have to adjust the null design matrix to account for the 
## new null hypothesis. Additionally, the nseq vector is going to
## change to reflect the fact that these will require larger samples
## to get the desired power
beta4=1.5
rho=.9

p=4
n=10
beta.star=c(1,1,1,beta4)
alpha=.01
d=1
sigma.star=1
sigma.mat=matrix(rho, nrow = 3, ncol = 3)
diag(sigma.mat)=1

nseq=seq(from = 20, to = 500, by = 20)
d=2
reps=5000

pval.list=numeric(reps)
npow.list=numeric(length(nseq))
k=1

for(n in nseq){
  print(k)
  for(i in 1:reps){
    evec=rnorm(n)
    Xtilde=mymvrnorm(n=n, mu=rep(0,p-1), Sigma=sigma.mat)
    X=cbind(1, Xtilde) 
    Y=X%*%beta.star+evec
    
    ## compute the QR decomposition of X
    qrX=qr(x=X)
    
    ## the null model design matrix is
    X0=X
    X0[,2]=X0[,2]+X0[,3]+X0[,4]
    X0=X0[,1:2]
    
    ## compute the QR decomposition of X0
    qrX0=qr(x=X0)
    
    beta.hat = qr.coef(qrX, y=Y)
    ## get rss1
    rss1=sum((Y-X%*%beta.hat)^2)
    
    ## get rss0
    beta.hat.0=qr.coef(qrX0, y=Y)
    rss0=sum((Y-X0%*%beta.hat.0)^2)
    
    ## compute the realization of F
    f=((rss0 - rss1)/d )/(rss1/(n-p))
    
    pval.list[i]=1-pf(f, d, n-p)
  }
  npow.list[k]=mean(pval.list < 0.1)
  k=k+1
}
npow.list















