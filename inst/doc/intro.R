## ----setup--------------------------------------------------------------------
library(StatComp22030)

## ----eval=FALSE---------------------------------------------------------------
#  function(e){
#    g<-rep(0,length(e))
#    g[which(e<=0)]=1
#    g
#  }

## ----eval=FALSE---------------------------------------------------------------
#  library(quantreg)
#  library(lpSolve)
#  function(X, y, J=5, tau=0.5){
#    data<-cbind(y,X)
#    n<-length(y)
#    Q<-floor(n/J)
#    M<-dim(X)[2]
#    # the estimator of theta in the mth model
#    theta.hat<-matrix(0,M,M)
#    for (m in 1:M) {
#      suppressWarnings(fit<-rq(data=data.frame(X[,1:m]), y~0+.,tau=tau))
#      theta.hat[1:m,m]<-fit$coefficients
#    }
#    # divide the data to J groups
#    index <- rep(1:J,Q)
#    index <- sample(index,n)
#    mu<-matrix(0, n, M)
#    y.CV<-rep(0, n)
#    # prepare the data for CV
#    for (j in 1:J) {
#      data.train<-data[index!=j,]
#      data.test<-data[index==j,]
#      for (m in 1:M) {
#        suppressWarnings(rq.m<-rq(data.train[,1]~0+.,
#                                  tau=tau,data=data.frame(data.train[,-1][,1:m])))
#        mu[((j-1)*Q+1):(j*Q),m]<-as.matrix(data.test[,-1][,1:m])%*%
#          as.matrix(rq.m$coefficients)
#      }
#      y.CV[((j-1)*Q+1):(j*Q)]<-data.test[,1]
#    }
#    # solve lp to get the weights
#    obj1<-rep(tau,len=n);obj2<-rep(1-tau,len=n);obj3<-rep(0,len=M)
#    f.obj<-c(obj1,obj2,obj3)
#    con1<-diag(2*n+M)
#    con21<-matrix(0,nrow = M,ncol = 2*n);
#    con22<-diag(M);
#    con2<-cbind(con21,con22)
#    con3<-c(rep(0,len=2*n),rep(1,len=M))
#    con4<-cbind(diag(n),-diag(n),mu)
#    f.con<-rbind(con1,con2,con3,con4)
#    f.dir<-c(rep(">=",len=2*n+M),rep("<=",len=M),rep("=",len=n+1) )
#    f.rhs<-c(rep(0,len=2*n+M),rep(1,len=M+1), y.CV)
#    lp.result<-lp("min", f.obj, f.con, f.dir, f.rhs)
#    w.JCVMA<-lp.result$solution[(2*n+1):(2*n+M)]
#    muhat<-as.matrix(X) %*% theta.hat
#    y.hat<-muhat %*% w.JCVMA
#    PE<-(tau-Indic(y-y.hat))%*%(y-y.hat)/n
#    return(list(w.JCVMA=w.JCVMA, theta.hat=theta.hat, y.hat=y.hat, PE=PE))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  function(x, qrma){
#    n<-dim(x)[1]
#    M<-dim(x)[2]
#    muhat<-as.matrix(x) %*% qrma$theta.hat
#    y.hat<-muhat %*% qrma$w.JCVMA
#    return(y.hat)
#  }

## -----------------------------------------------------------------------------
data(express)
attach(express)
y<-express$y[1:100]
y.s<-express$y[101:150]
X<-express[1:100,-1]
X<-cbind(1,X)
X.s<-express[101:150,-1]
X.s<-cbind(1,X.s)
tau<-0.5
qrma.fit<-qrma(X, y, J=5, tau = tau)
w<-qrma.fit$w.JCVMA # weight
PE<-qrma.fit$PE # in-sample FPE
y.hat<-predict_qr(X.s,qrma.fit) # estimation

tau<-0.05
qrma.fit<-qrma(X, y, J=5, tau = tau)
qrma.fit<-qrma(X, y, J=5, tau = tau)
w<-qrma.fit$w.JCVMA # weight
PE<-qrma.fit$PE # in-sample FPE
y.hat<-predict_qr(X.s,qrma.fit) # estimation


## ----eval=FALSE---------------------------------------------------------------
#  function(N) {
#    mat <- matrix(nrow = N, ncol = 2)
#    rho <- 0.9 # correlation
#    mu1<-mu2<-0
#    sigma1<-sigma2<-1
#    s1 <- sqrt(1-rho^2)*sigma1
#    s2 <- sqrt(1-rho^2)*sigma2
#    x <- y <- 0
#    for (i in 1:N) {
#      x <- rnorm(1, mu1 + rho * (y - mu2) * sigma1/sigma2, s1)
#      y <- rnorm(1, mu2 + rho * (x - mu1) * sigma2/sigma1, s2)
#      mat[i, ] <- c(x, y)
#    }
#    mat
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix gibbsC(int N) {
#    NumericMatrix mat(N, 2);
#    double  rho=0.9, mu1=0, mu2=0, sigma1=1, sigma2=1, s1=0, s2=0;
#    s1 = sqrt(1-rho*rho)*sigma1;
#    s2 = sqrt(1-rho*rho)*sigma2;
#    double x = 0, y = 0;
#    for(int i = 0; i < N; i++) {
#      x = rnorm(1, mu1 + rho * (y - mu2) * sigma1/sigma2, s1)[0];
#      y = rnorm(1, mu2 + rho * (x - mu1) * sigma2/sigma1, s2)[0];
#      mat(i, 0) = x;
#      mat(i, 1) = y;
#    }
#    return(mat);
#  }

