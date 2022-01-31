## Question 1 ##

## The function rdist takes in n as its input (default = 1) and
## uses a rejection sampling algorithm to generate n independent 
## realizations of the distribution f(x) = (5/2)x^4
rdist = function(n = 1){
  ## Allocating memory
  z.vec = numeric(n)
  
  ## The heart of the rejection sampling process
  for(i in 1:n){  
     ## The while loop repeatedly draws from the unif(0,1) and unif(-1,1)
    ## until the event occurs. Then we store that realization of z and 
    ## repeat the process n times
    event = FALSE
    while(!event){
      u = runif(1)
      z = runif(1, min = -1, max = 1)
      event = 1*(u<z^4)
    }
    z.vec[i] = z
  }
  return(z.vec)
}

## Generating 50000 realizations of our random variable
y.list = sqrt(abs(rdist(5e4)))

## Simulation based approximate 99% ci
t.test(y.list, conf.level = .99)$conf.int[1:2]


## Question 2 ##
## I define the inversion method as a function here so I can make calls
## to it in the simulation study and keep that a bit neater
inv.mtd = function(n, theta){
  z.list = theta*sqrt(runif(n))
  return(z.list)
}

n=20
theta=4
reps=5e4

## Generating reps independent realizations of V = max(Y_1,...Y_20)/theta
v.list=numeric(reps)
for(i in 1:reps){
  v.list[i]=max(inv.mtd(n=n, theta=theta))/theta
}

probs=seq(from=0.01, to=0.99, by=0.01)
plot(probs^(1/(2*n)), quantile(v.list, probs), xlab = 'Theoretical Percentiles',
     ylab = 'Sample Percentiles', main = 'QQ-Plot')
abline(0,1)


n=20
theta=4
alpha=.05
reps=5e4

captured.list=numeric(reps)
for(i in 1:reps){
  ## Generating a realization of M 
  M = max(inv.mtd(n=n, theta=theta))
  
  ## Checking to see if our interval captured the true theta
  captured.list[i] = 1*(M <= theta & theta <= (M/(alpha)^(.5/n)))
}

## Finding the simulation based score approximate 99% ci for this 
## coverage probability
prop.test(sum(captured.list), reps, conf.level = .99)$conf.int


## Question 3 ##
n=100
sigma.star=.5
beta.star=c(1,0,0)
reps=1e5

b1=numeric(reps)
b2=numeric(reps)
b3=numeric(reps)
for(i in 1:reps){
  ## Generate the data
  X=matrix(rnorm(2*n), nrow=n, ncol=2)
  X=cbind(1,X)
  
  ## Generating realizations of Y
  Y=X%*%beta.star + sigma.star*rnorm(n)
  
  ## Getting estimates for beta.hat and sE
  
  beta.hat = qr.solve(crossprod(X))%*%crossprod(X,Y)
  sE = sqrt(sum((Y-X%*%beta.hat)^2)/(n-3))
  
  ## Finding realizations of cov((beta.hat.j-beta.star)(sE^2 - sigma.star^2))
  b1[i]=(beta.hat[1]-beta.star[1])*(sE^2-sigma.star^2)
  b2[i]=(beta.hat[2]-beta.star[2])*(sE^2-sigma.star^2)
  b3[i]=(beta.hat[3]-beta.star[3])*(sE^2-sigma.star^2)
}

## 95% approx cis
t.test(b1, conf.level = .95)$conf.int[1:2]
t.test(b2, conf.level = .95)$conf.int[1:2]
t.test(b3, conf.level = .95)$conf.int[1:2]



## A function that returns the BIC for a given design matrix and 
## vector of responses
get.bic=function(X, y)
{
  n=dim(X)[1]
  p=dim(X)[2]
  beta.hat=qr.solve(crossprod(X))%*%crossprod(X, y)
  rss=sum((y-X%*%beta.hat)^2)
  bic=n*log(2*pi)+n*log(rss/n) + n + (p+1)*log(n)
  return(bic)
}

## Generating a new design matrix with 100 rows to select from
## and corresponding realizations of Y
X=matrix(rnorm(2*100), nrow=100, ncol=2)
X=cbind(1,X)
sigma.star=.5
beta.star=c(1,0,0)

n.list=c(5,10,20,50,100)
reps=5e4

## Creating a matrix to store our responses
response.mat = matrix(NA, nrow = reps, ncol = length(n.list))

## We repeat the simulation study for each level of n
for(i in 1:length(n.list)){
  ## Choosing the level of n for this simulation
  k = n.list[i]
  for(j in 1:reps){
    ## Subsetting the data and generating realizations of Y
    Xnew = X[1:k,]
    Ynew=Xnew%*%beta.star + sigma.star*rnorm(k)
    
    ## Defining the null model design matrix
    X0 = as.matrix(Xnew[,-c(2,3)])
    
    ## Storing the results of this replication. 1 indicates
    ## that bic failed to select the null (correct) model
    response.mat[j,i] = 1*(get.bic(X=Xnew, y=Ynew) < get.bic(X=X0, y=Ynew))
  }
}
apply(response.mat, 2, mean)










