## ----setup--------------------------------------------------------------------
library(StatComp22030)

## -----------------------------------------------------------------------------
library(MASS) 
attach(Boston) # fix the Boston dataset
names(Boston) # column names of the data frame

## -----------------------------------------------------------------------------
lm.Boston=lm(medv~age,data=Boston)
summary(lm.Boston)

## -----------------------------------------------------------------------------
plot(age,medv) # the scatter plot
abline(lm.Boston,lwd=1.5) # the regression line
title('Boston housing values')

## -----------------------------------------------------------------------------
knitr::kable(head(Boston))

## -----------------------------------------------------------------------------
set.seed(50)
a <- 2; b <- 2 # Pareto(2,2)
n <- 1000
u <- runif(n)
x <- b/(1-u)^{1/a}
hist(x, prob = TRUE, main = expression(f(x)==8/x^3 )) # the density f(x)=8/x^3
y <- seq(2,50,0.01)
lines(y, 8/y^3)

## -----------------------------------------------------------------------------
rBeta<-function(n,a,b){
  j<-k<-0; y <- numeric(n)
  while (k < n) {
  u <- runif(1)
  j <- j + 1
  x <- runif(1) # random variate from g(.)
  if (x^(a-1) * (1-x)^(b-1) > u) {
    # we accept x
    k <- k + 1
    y[k] <- x
  }
  }
  return(y)
}

## -----------------------------------------------------------------------------
set.seed(35)
x <- rBeta(1000,3,2)
hist(x, prob = TRUE, main = expression(f(x)==x^2*(1-x)/B(3,2) )) # the density f(x)=x^2*(1-x)/B(3,2)
y <- seq(0,1,0.01)
lines(y, y^2*(1-y)/beta(3,2))

## -----------------------------------------------------------------------------
set.seed(233)
n <- 1000; r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n, lambda)

## -----------------------------------------------------------------------------
set.seed(25)
n <- 1000; r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n, lambda)
hist(x, prob = TRUE, main = expression(f(x)==64/(2+x)^5 ) ) #the density f(x)=64/(2+x)^5
y<-seq(0,10,0.01)
lines(y,64/(2+y)^5)

## -----------------------------------------------------------------------------
# the fast sorting algorithm
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}
}

## -----------------------------------------------------------------------------
# since we have 100 replications, I don't use random number seeds.
n.set<-c(1e4,2e4,4e4,6e4,8e4)
time <-matrix(0,100,5) # initialize time
for (r in 1:100) { # the rth replication
  for (i in 1:5) {
    n<-n.set[i]
    time[r,i]=system.time(quick_sort(sample(1:n)))[1]
  }
  a<-colMeans(time) # a is the computation time averaged over 100 simulations
}
a

## -----------------------------------------------------------------------------
t<-n.set*log(n.set)
lm.fit<-lm(a~t)
plot(t,a,xlab = expression(n*log(n)),ylab = 'Computation time') # plot the scatter plot
abline(lm.fit,lwd=1.5) # plot the regression line

## -----------------------------------------------------------------------------
# m is the sample size
# if antithetic = TRUE, we use the antithetic variate approach
# if antithetic = FALSE, we use the simple Monte Carlo method
MC <- function(m = 10000, antithetic = FALSE) {
  u <- runif(m/2)
  if (antithetic) v <- 1-u else v <- runif(m/2)
  u <- c(u, v)
  g <- exp(u)
  theta <- mean(g)
  return(theta)
}

## -----------------------------------------------------------------------------
set.seed(999)
MC.1<-MC(antithetic = FALSE) # the simple Monte Carlo estimator
MC.2<-MC(antithetic = TRUE) # the antithetic variable estimator
c(MC.1,MC.2)

## -----------------------------------------------------------------------------
R<-1000
theta.1 <- theta.2 <- numeric(R) # initialize
for (i in 1:R) {
  theta.1[i] <- MC(m = 1000, antithetic = FALSE) # the simple estimator
  theta.2[i] <- MC(m = 1000, antithetic = TRUE) # the antithetic variable estimator
}
# an empirical estimate of the percent reduction in variance
(var(theta.1)-var(theta.2))/var(theta.1)

## -----------------------------------------------------------------------------
x<-seq(1, 10, 0.01)
g<-x^2*exp(-x^2/2)/sqrt(2 * pi)
plot(x, g, type = 'l', ylab = "", ylim = c(0,1), col=1, lwd=1.5, main = "function graph")
f_1<-2*dnorm(x,1)
lines(x, f_1, lty=2, col=2, lwd=1.5)
f_2<-dchisq(x-1,df=3)
lines(x, f_2, lty=3, col=3, lwd=1.5)
legend("topright",  legend = c("g(x)", expression(f[1](x)),expression(f[2](x))),
       inset = 0.02, lty=1:3, col=1:3)

## -----------------------------------------------------------------------------
plot(x, g/f_1, type = 'l',lty=2, ylab = "", ylim = c(0,6.5),
     col=2, lwd=1.5, main = expression(g(x) / f(x)))
lines(x, g/f_2, lty=3, col=3, lwd=1.5)
legend("topright",  legend = c(expression(g(x)/f[1](x)),expression(g(x)/f[2](x))),
       inset = 0.02, lty=2:3, col=2:3)

## -----------------------------------------------------------------------------
set.seed(999)
m<-1e6
# using f1
x<- abs(rnorm(m))+1
f<-2*dnorm(x, 1)
g<-x^2*exp(-x^2/2)/sqrt(2 * pi)
theta_1<-mean(g/f)
var_1<-var(g/f)
# using f2
x<- rchisq(m,df=3)+1
f<-dchisq(x-1,df=3)
g<-x^2*exp(-x^2/2)/sqrt(2 * pi)
theta_2<-mean(g/f)
var_2<-var(g/f)
result<-rbind(estimate=c(theta_1,theta_2), variance=c(var_1,var_2))
colnames(result) <- c('f1','f2')
library(knitr)
knitr::kable(round(result,4))

## -----------------------------------------------------------------------------
set.seed(233)
M<-1e6
k<-5
var<-theta<-numeric(k)
g <- function(x) exp(-x)/(1 + x^2)
f <- function(x) exp(-x)/(1 - exp(-1))
for (j in 1:k) {
  # generate samples from fj by the inverse transform method
  u<-runif(M/k, (j-1)/5,j/5)
  x<--log(1-u*(1-exp(-1)))
  theta[j]<-mean(g(x)/f(x))
  var[j]<-var(g(x)/f(x))
}

# the stratified importance sampling estimate of theta
round(mean(theta),4)
# the standard deviation of theta multiplied by sqrt(M)
round(sqrt(mean(var)),4)

## ----echo=FALSE---------------------------------------------------------------
  m <- 1e6
  est <- sd <- numeric(5)
  g <- function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
  }
  x <- runif(m) #using f0
  fg <- g(x)
  est[1] <- mean(fg)
  sd[1] <- sd(fg)
  x <- rexp(m, 1) #using f1
  fg <- g(x) / exp(-x)
  est[2] <- mean(fg)
  sd[2] <- sd(fg)
  x <- rcauchy(m) #using f2
  i <- c(which(x > 1), which(x < 0))
  x[i] <- 2 #to catch overflow errors in g(x)
  fg <- g(x) / dcauchy(x)
  est[3] <- mean(fg)
  sd[3] <- sd(fg)
  u <- runif(m) #f3, inverse transform method
  x <- - log(1 - u * (1 - exp(-1)))
  fg <- g(x) / (exp(-x) / (1 - exp(-1)))
  est[4] <- mean(fg)
  sd[4] <- sd(fg)
  u <- runif(m) #f4, inverse transform method
  x <- tan(pi * u / 4)
  fg <- g(x) / (4 / ((1 + x^2) * pi))
  est[5] <- mean(fg)
  sd[5] <- sd(fg)
  res <- rbind(est=round(est,4), sd=round(sd,4))
  colnames(res) <- paste0('f',0:4)
  knitr::kable(res)

## -----------------------------------------------------------------------------
rm(list=ls()) # clear up the memory
# define a function CI to find a confidence interval of logmu=0, i.e. mu=1
CI<-function(n,aplha){
  X<-rlnorm(n) # we can use rlnorm to generate the data, which is already a function
  Y<-log(X)
  se<-sd(Y)/sqrt(n)
  mean(Y) + se * qnorm(c(alpha/2, 1-alpha/2))
}

## -----------------------------------------------------------------------------
n<-100
alpha<-0.05
set.seed(3278)
CI(n,alpha)

## -----------------------------------------------------------------------------
# define a function MC to use a Monte Carlo method to obtain an empirical estimate of the confidence level
MC<-function(n,alpha,M){
  LU<-matrix(0,M,2)
  for (i in 1:M) {
    LU[i,]<-CI(n,alpha)
    }
  mean(LU[,1] < 0 & LU[,2] > 0) # an empirical estimate of the confidence level
}
# obtain an empirical estimate of the confidence with M=1e4
set.seed(22030)
M<-1e4
MC(n,alpha,M)

## -----------------------------------------------------------------------------
rm(list=ls()) # clear up the memory
# define a function count5test to conduct count five test
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5)) }

## -----------------------------------------------------------------------------
# generate samples under H1 to estimate power
sigma1 <- 1; sigma2 <- 1.5
# define a function to conduct of the Count Five test and F test with a sample size of n
test<-function(n){
  # we can use rnorm to generate the data, which is already a function
  x <- rnorm(n, 0, sigma1)
  y <- rnorm(n, 0, sigma2)
  C5t <- count5test(x, y) # count five test
  Fp <- var.test(x, y)$p.value
  Ftest <- as.integer(Fp < 0.055) # F test
  c(C5t, Ftest)
}

## -----------------------------------------------------------------------------
# define a function to estimate the power of the Count Five test and F test
power<-function(n,M){
  result<-matrix(0,M,2)
  for (i in 1:M) {
    result[i,]<-test(n)
  }
  colMeans(result)
}

## -----------------------------------------------------------------------------
M<-1e4
nset<-c(10,20,50,100,200,500,1000)
set.seed(1643)
result<-matrix(0,length(nset),3)
colnames(result)<-c('n','count5test','Ftest')
for (i in 1:length(nset)) {
  n<-nset[i]
  result[i,]<-c(n,power(n,M))
}
result

## -----------------------------------------------------------------------------
rm(list=ls()) # clear up the memory
# define a function to conduct Z-test
Ztest<-function(m,p1,p2){
  Ts<-sqrt(m)*(p1-p2)/sqrt(p1*(1-p1)+p2*(1-p2))
  abs(Ts)>qnorm(0.975) # if TRUE, we rejuct H0
}
Ztest(1e4,0.651,0.676)

## -----------------------------------------------------------------------------
library(boot)
set.seed(126)
X <- aircondit[1]
MLE <- function(X, i) 1/mean(as.matrix(X[i, ]))
# use bootstrap to estimate the bias and standard error of the estimate
boot(X, statistic=MLE, R=1e4)
rm(list=ls())

## -----------------------------------------------------------------------------
set.seed(701)
X <- aircondit[1]
theta <- function(X, i) mean(as.matrix(X[i, ]))
(boot<-boot(X, statistic=theta, R=1e4))
# compute 95% bootstrap confidence intervals
(ci <- boot.ci(boot,type=c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
hist(boot$t,main='Histgram of bootstrap statistic')
abline(v=mean(as.matrix(X)),col='red',lwd=2)
rm(list=ls())

## -----------------------------------------------------------------------------
set.seed(23387)
mu<-0;sigma<-1;m<-1e4;n<-100
mean.sample <- function(X,i) mean(X[i])
ci.norm<-ci.basic<-ci.perc<-matrix(0,m,2)
# conduct a Monte Carlo simulation
for(i in 1:m){
  X<-rnorm(n,mu,sigma)
  boot<-boot(data=X, statistic=mean.sample, R = 100)
  ci<-boot.ci(boot,type=c("norm","basic","perc"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
}
# the empirical coverage rates for the sample mean
(result1<-c("norm"=mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),
            "basic"=mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),
            "perc"=mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu)))
# the proportion of times that the confidence intervals miss on the left
(result2<-c("norm"=mean(ci.norm[,2]< mu),
            "basic"=mean(ci.basic[,2]< mu),
            "perc"=mean(ci.perc[,2]< mu)))
# the proportion of times that the confidence intervals miss on the right
(result3<-c("norm"=mean(ci.norm[,1]> mu),
            "basic"=mean(ci.basic[,1]> mu),
            "perc"=mean(ci.perc[,1]> mu)))
rm(list=ls())

## -----------------------------------------------------------------------------
library(bootstrap)
x <- as.matrix(scor)
n <- nrow(x)
lambda <- eigen(cov(x))$values
theta.hat <- max(lambda/sum(lambda)) # original estimate
theta.jack <- numeric(n)
for(i in 1:n){
  lambda.jack <- eigen(cov(x[-i, ]))$values
  theta.jack[i] <- max(lambda.jack/sum(lambda.jack)) # jackknife estimate
}
# the jackknife estimates of bias
bias.jack <- (n-1)*(mean(theta.jack) - theta.hat)
# the jackknife estimates of standard error
se.jack <- sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2))
c(est=theta.hat, bias=bias.jack, se=se.jack)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic)
N <- n*(n-1)/2 # all possible combination
e1 <- e2 <- e3 <- e4 <- numeric(n)
h<-1
# leave-two-out cross validation
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    k<-c(i,j)
    y <- magnetic[-k]
    x <- chemical[-k]
    # Model 1: Linear
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2]*chemical[k]
    e1[h] <- sum((magnetic[k] - yhat1)^2)
    # Model 2：Quadratic
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2]*chemical[k] +J2$coef[3]*chemical[k]^2
    e2[h] <- sum((magnetic[k] - yhat2)^2)
    # Model 3: Exponential
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2]*chemical[k]
    yhat3 <- exp(logyhat3)
    e3[h] <- sum((magnetic[k] - yhat3)^2)
    # Model 4: Log-Log
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    yhat4 <- exp(logyhat4)
    e4[h] <- sum((magnetic[k] - yhat4)^2)
    h<-h+1
  }
}
# the average squared prediction error by leave-two-out cross validation
c(Linear=sum(e1), Quadratic=sum(e2), Exponential=sum(e3), LogLog=sum(e4))/(2*N)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
# define a function to conduct permutation test
spear.perm <- function(x, y, R=1e3){
  rho0<-cor.test(x, y, method = "spearman")$estimate # the original test estimate
  n<-length(y)
  reps <- numeric(R)
  for (i in 1:R) {
    k <- sample(1:n)
    reps[i] <- cor.test(x, y[k], method = "spearman")$estimate
    }
  p <- mean(c(rho0, reps) >= rho0) # p-value of the permutation test
  return(p) # return p-value
}

## -----------------------------------------------------------------------------
library(MASS)
# generate data from multi-normal distribution
mu <- c(0, 0); sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
n<-30
set.seed(312)
xy <- mvrnorm(n, mu, sigma)
x<-xy[, 1]; y<-xy[, 2]
# the p-value reported by cor.test on the same samples
(p0<-cor.test(x, y, method = "spearman")$p.value)
# the achieved signiﬁcance level of the permutation test
(p.perm<- spear.perm(x,y))
# compare the two p-value
round(c(p0=p0,perm=p.perm),4)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
# define a function to implement a random walk Metropolis sampler
rw.Laplace <- function(N, x0, sigma) {
  # N is the length of chain, x0 is the initial value
  # sigma is the standard deviation of the normal proposal distribution
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(abs(x[i-1]) - abs(y)))
      x[i] <- y
    else {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x = x, k = k))
}

## -----------------------------------------------------------------------------
N <- 10000 # length of chains
sigma <- c(0.5, 1, 2, 4) # the standard deviation of the normal proposal distribution
x0 <- 0 # initial values
set.seed(123)
rw1 <- rw.Laplace(N, x0, sigma[1])
rw2 <- rw.Laplace(N, x0, sigma[2])
rw3 <- rw.Laplace(N, x0, sigma[3])
rw4 <- rw.Laplace(N, x0, sigma[4])
par(mfrow=c(2,2))  #display 4 graphs together
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
# for (j in 1:4) {
#   plot(rw[,j], type="l",xlab=bquote(sigma == .(round(sigma[j],3))),
#        ylab="X", ylim=range(rw[,j]))
#   }
par(mfrow=c(1,1)) #reset to default

## -----------------------------------------------------------------------------
acceptance.rate <- data.frame(sigma=sigma, acceptance.rate=1-c(rw1$k, rw2$k, rw3$k, rw4$k)/N)
knitr::kable(acceptance.rate,align = 'c')

## -----------------------------------------------------------------------------
# define a function to use the Gelman-Rubin method to monitor convergence of the chain
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance estimation
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within estimation
  v.hat <- W*(n-1)/n + (B/n)     #upper variance estimation
  r.hat <- v.hat / W             #estimation
  return(r.hat)
}

N <- 15000 # length of chains
k <- 4    # number of chains to generate
sigma <- 1 # the standard deviation of the normal proposal distribution
b <- 1000       # burn-in length
# choose overdispersed initial values
x0 <- c(-4, -2, 2, 4)
set.seed(25)
# generate the chains
X <- matrix(0, nrow=k, ncol=N)
for (i in 1:k){
  X[i, ] <- rw.Laplace(N, x0[i], sigma)$x
}

# compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)){
  psi[i,] <- psi[i,] / (1:ncol(psi))
}
# plot psi for the four chains
for (i in 1:k){
  if(i==1){
    plot((b+1):N,psi[i, (b+1):N], ylim=c(-0.5,0.5),type="l",
         xlab='Index', ylab=bquote(phi))
    }else{
      lines(psi[i, (b+1):N], col=i)
    }
}
# plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in (b+1):N){
  rhat[j] <- Gelman.Rubin(psi[,1:j])
}
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
# define a function to use Gibbs sampler to generate the chain
Gibbs.sampler<-function(N, x0){
  # N is the length of chain, x0 is the initial value
  X <- matrix(0, N, 2) # the chain, a bivariate sample
  rho <- 0.9 # correlation
  mu1<-mu2<-0
  sigma1<-sigma2<-1
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  X[1, ] <- x0 #initialize
  for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
  }
  return(X)
  }

## -----------------------------------------------------------------------------
N <- 5000 # length of chain
burn <- 1000 # burn-in length
x0<-rep(0,2) # initialize
set.seed(65)
X<-Gibbs.sampler(N,x0) # generate the chain
x <- X[(burn + 1):N, ]
cat('Means: ',round(colMeans(x),3))
cat('Standard errors: ',round(apply(x,2,sd),3))
cat('Correlation coefficients: ', round(cor(x[,1],x[,2]),3))
# plot the generated sample
plot(x[,1],type='l',col=1,lwd=2,xlab='Index',ylab='Random numbers')
lines(x[,2],col=2,lwd=2)
legend('bottomright',c('X','Y'),col=1:2,lwd=2)

## -----------------------------------------------------------------------------
L<-lm(x[,1]~x[,2])
summary(L)

## -----------------------------------------------------------------------------
# define a function to use the Gelman-Rubin method to monitor convergence of the chain
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance estimation
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within estimation
  v.hat <- W*(n-1)/n + (B/n)     #upper variance estimation
  r.hat <- v.hat / W             #estimation
  return(r.hat)
}

N <- 15000 # length of chains
k <- 4    # number of chains to generate
b <- 1000       # burn-in length
# choose overdispersed initial values
x0 <- matrix(c(0,0,-0.1,0.1,0.2,-0.2,0.1,-0.2),4,2,byrow = TRUE)
set.seed(262)
# generate the chains
X <-Y<-matrix(0, nrow=k, ncol=N)
for (i in 1:k){
  X[i, ] <- Gibbs.sampler(N, x0[i,])[,1]
  Y[i, ] <- Gibbs.sampler(N, x0[i,])[,2]
}
# compute diagnostic statistics
psi.X <- t(apply(X, 1, cumsum))
psi.Y <- t(apply(Y, 1, cumsum))
for (i in 1:nrow(psi.X)){
  psi.X[i,] <- psi.X[i,] / (1:ncol(psi.X))
  psi.Y[i,] <- psi.Y[i,] / (1:ncol(psi.Y))
}

# plot the sequence of R-hat statistics
rhat.X <-rhat.Y <- rep(0, N)
for (j in (b+1):N){
  rhat.X[j] <- Gelman.Rubin(psi.X[,1:j])
  rhat.Y[j] <- Gelman.Rubin(psi.Y[,1:j])
}
# compute rhat for 2-dimensional data
rhat.XY<-matrix(c(rhat.X,rhat.Y),ncol = 2)
rhat<-apply(rhat.XY, 1, max) # find the maximum of rhat.X and rhat.Y
rhat<-sqrt(rhat)  # get rhat for 2-dimensional data
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)


## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
library(mediation)
# define a function to conduct permutation test in situation 1: alpha=0
mediation.perm.1 <- function(X, Y, M, R=100){
  a<-lm(M~X)
  b<-lm(Y~X+M)
  result0<-mediate(model.m=a, model.y=b, sims=50, treat='X', mediator = 'M')
  T0<-abs(result0$d0/sd(result0$d0.sims)) # the original statstic
  n<-length(X)
  reps <- numeric(R)
  for (i in 1:R) {
    k <- sample(1:n)
    M1<-M[k]
    Y1<-Y[k]
    result1<-mediate(model.m=lm(M1~X), model.y=lm(Y1~X+M1), sims=50,
                     treat='X', mediator = 'M1')
    reps[i]<-abs(result1$d0/sd(result1$d0.sims))
    }
  p <- mean(c(T0, reps) >= T0) # p-value of the permutation test
  return(p) # return p-value
}

# define a function to conduct permutation test in situation 2: beta=0
mediation.perm.2 <- function(X, Y, M, R=100){
  a<-lm(M~X)
  b<-lm(Y~X+M)
  result0<-mediate(model.m=a, model.y=b, sims=50, treat='X', mediator = 'M')
  T0<-abs(result0$d0/sd(result0$d0.sims)) # the original statstic
  n<-length(X)
  reps <- numeric(R)
  for (i in 1:R) {
    k <- sample(1:n)
    M1<-M[k]
    X1<-X[k]
    result1<-mediate(model.m=lm(M1~X1), model.y=lm(Y~X1+M1), sims=50,
                     treat='X1', mediator = 'M1')
    reps[i]<-abs(result1$d0/sd(result1$d0.sims))
    }
  p <- mean(c(T0, reps) >= T0) # p-value of the permutation test
  return(p) # return p-value
}

# define a function to conduct permutation test in situation 3: alpha=0, beta=0
mediation.perm.3 <- function(X, Y, M, R=100){
  a<-lm(M~X)
  b<-lm(Y~X+M)
  result0<-mediate(model.m=a, model.y=b, sims=50, treat='X', mediator = 'M')
  T0<-abs(result0$d0/sd(result0$d0.sims)) # the original statstic
  n<-length(X)
  reps <- numeric(R)
  for (i in 1:R) {
    k <- sample(1:n)
    M1<-M[k]
    result1<-mediate(model.m=lm(M1~X), model.y=lm(Y~X+M1), sims=50,
                     treat='X', mediator = 'M1')
    reps[i]<-abs(result1$d0/sd(result1$d0.sims))
    }
  p <- mean(c(T0, reps) >= T0) # p-value of the permutation test
  return(p) # return p-value
}


## -----------------------------------------------------------------------------
am<-ay<-1; gamma<-1
n<-50
alpha<-0; beta<-0
set.seed(10)
X<-rnorm(n)
M<-am + alpha * X + rnorm(n)
Y<-ay + beta * M + gamma * X + rnorm(n)
set.seed(11)
# p-value of 3 tests
c(perm.1=mediation.perm.1(X, Y, M), perm.2=mediation.perm.2(X, Y, M),
  perm.3= mediation.perm.3(X, Y, M))

## -----------------------------------------------------------------------------
alpha<-0; beta<-1
set.seed(67)
X<-rnorm(n)
M<-am + alpha * X + rnorm(n)
Y<-ay + beta * M + gamma * X + rnorm(n)
set.seed(68)
# p-value of 3 tests
c(perm.1=mediation.perm.1(X, Y, M), perm.2=mediation.perm.2(X, Y, M),
  perm.3= mediation.perm.3(X, Y, M))

## -----------------------------------------------------------------------------
alpha<-1; beta<-0
set.seed(564)
X<-rnorm(n)
M<-am + alpha * X + rnorm(n)
Y<-ay + beta * M + gamma * X + rnorm(n)
set.seed(565)
# p-value of 3 tests
c(perm.1=mediation.perm.1(X, Y, M), perm.2=mediation.perm.2(X, Y, M),
  perm.3= mediation.perm.3(X, Y, M))

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
simulation <- function(N, b1, b2, b3, f0){
  x1 <- rpois(N,1); x2<-rexp(N,1); x3<-rbinom(N,1,0.5)
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
    }
  solution <- uniroot(g,c(-50,0))
  return(solution$root)
  }

## -----------------------------------------------------------------------------
N <- 1e6; b1 <- 0; b2 <- 1; b3<--1
f0<-c(0.1,0.01,0.001,0.0001)
set.seed(324)
# use the function and get the result
alpha<-sapply(f0, FUN=simulation, N=N, b1=b1, b2=b2, b3=b3)
rbind(f0,alpha)
# draw the scatter plot
plot(log(f0),alpha,main = expression(paste(f[0],' vs. ',alpha)))
rm(list=ls())

## -----------------------------------------------------------------------------
# create a function to use maximize observed data likelihood to get MLE of lambda
maximize_observed_data_likelihood<-function(interval){
  g <- function(lambda){
    tmp <- (interval[,1]*exp(-lambda*interval[,1])-interval[,2]*
              exp(-lambda*interval[,2]))/(exp(-lambda*interval[,1])
                                          -exp(-lambda*interval[,2]))
    sum(tmp)
    }
  solution <- uniroot(g,c(0,20))
  return(solution$root)
}

# create a function to use the EM algorithm to get MLE of lambda
EM_algorithm<-function(interval, eplison=0.001, max_iter=100){
  lambda0<-1
  for(i in 1:max_iter) {
    tmp<-(interval[,1]*exp(-lambda0*interval[,1])-interval[,2]*
              exp(-lambda0*interval[,2]))/(exp(-lambda0*interval[,1])-
                                            exp(-lambda0*interval[,2]))+1/lambda0
    lambda<-dim(uv)[1]/sum(tmp)
    res<-abs(lambda-lambda0)
    if(res<=eplison) break
    lambda0<-lambda
  }
  return(lambda)
}

## -----------------------------------------------------------------------------
uv<-rbind(c(11,12),c(8,9),c(27,28),c(13,14),c(16,17),c(0,1),
          c(23,24),c(10,11),c(24,25),c(2,3))
# use maximize observed data likelihood to get MLE of lambda
maximize_observed_data_likelihood(uv)
# use the EM algorithm to get MLE of lambda
EM_algorithm(uv)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
x <- list(list(1, 2), c(3, 4))
(y1<-unlist(x))
str(y1) # a vector
(y2<-as.vector(x))
str(y2) # a list
rm(list = ls())

## -----------------------------------------------------------------------------
x<-c(1,2,3)
dim(x)
rm(list = ls())

## -----------------------------------------------------------------------------
x<-matrix(0,2,2)
is.matrix(x)
is.array(x)
rm(list = ls())

## -----------------------------------------------------------------------------
x<-data.frame(x=c(1,2), y=c(3,4))
attributes(x)
rm(list = ls())

## -----------------------------------------------------------------------------
x<-data.frame(x=c(1,2), y=c('three','four'))
(y<-as.matrix(x))
typeof(y)
a<-data.frame(x=c(1,2), y=c(TRUE,FALSE))
(b<-as.matrix(a))
typeof(b)
rm(list = ls())

## -----------------------------------------------------------------------------
x<-data.frame(matrix(nrow = 0, ncol = 2))
dim(x)
y<-data.frame(matrix(nrow = 2, ncol = 0))
dim(y)
rm(list = ls())

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
  }

## -----------------------------------------------------------------------------
# the new function, which is scale01 if the data is numeric, or return the data itself
scale<-function(x) {
  if (is.numeric(x))
    scale01(x)
  else
    x
}
# examples
x<-data.frame(a=c(1,2,3),b=c(4,5,6))
data.frame(lapply(x, scale))
y<-data.frame(a=c('one','two','three'),b=c(4,5,6))
data.frame(lapply(y, scale))
rm(list=ls())

## -----------------------------------------------------------------------------
set.seed(22030)
# a) Compute the standard deviation of every column in a numeric data frame.
x<-data.frame(a=rnorm(10),b=rnorm(10,0,0.5))
vapply(x, sd, numeric(1))
# b) Compute the standard deviation of every numeric column in a mixed data frame.
y<-data.frame(a=rnorm(10),b=rnorm(10,0,0.5),c=rep(c('one','two'),c(5,5)))
vapply(y[vapply(y, is.numeric, logical(1))], sd, numeric(1))
rm(list=ls())

## -----------------------------------------------------------------------------
library(Rcpp)

## -----------------------------------------------------------------------------
set.seed(22030)
gibbR<-gibbsR(1000)
gibbC<-gibbsC(1000)
# compare X and Y
qqplot(gibbR[,1],gibbR[,2],xlab = 'x',ylab = 'y',main='QQ plot of gibbR')
qqplot(gibbC[,1],gibbC[,2],xlab = 'x',ylab = 'y',main='QQ plot of gibbC')

# compare gibbR and gibbC
qqplot(gibbR[,1],gibbC[,1],xlab = 'gibbR',ylab = 'gibbC',main='QQ plot of X')
qqplot(gibbR[,2],gibbC[,2],xlab = 'gibbR',ylab = 'gibbC',main='QQ plot of Y')

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(gibbR=gibbsR(100),
                         gibbC=gibbsC(100))
summary(ts)[,c(1,3,5,6)]
rm(list=ls())

