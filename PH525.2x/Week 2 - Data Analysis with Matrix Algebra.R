set.seed(1)
g <- 9.8 ##meters per second
n <- 25
tt <- seq(0,3.4,len=n) ##time in secs, note: we use tt because t is a base function
##d <- 56.67  - 0.5*g*tt^2 + rnorm(n,sd=1) ##meters
f <- 56.67  - 0.5*g*tt^2
y <-  f + rnorm(n,sd=1)
plot(tt,y,ylab="Distance in meters",xlab="Time in seconds")
lines(tt,y,col=2)

rss <- function(Beta0,Beta1,Beta2) {
  r <- y - (Beta0+Beta1*tt+Beta2*tt^2)
  return(sum(r^2))
}

Beta2s <- seq(-10,0,len=100)
RSS <- sapply(Beta2s,rss,Beta0=65,Beta1=0)
plot(Beta2s,RSS,type="l", col=3)

# Beta2s<- seq(-10,0,len=100)
# plot(Beta2s,sapply(Beta2s,rss,Beta0=55,Beta1=0),ylab="RSS",xlab="Beta2",type="l")

tt2 <-tt^2
fit <- lm(y~tt+tt2)
summary(fit)$coef

X <- cbind(rep(1,length(tt)),tt,tt^2)
head(X)

Beta <- matrix(c(55,0,5),3,1)

r <- y - X %*% Beta
RSS <- t(r) %*% r

RSS <- crossprod(r)
RSS

betahat <- solve( t(X)%*%X) %*% t(X) %*% y

betahat <- solve( crossprod(X)) %*% crossprod(X,y)

QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)
backsolve(R, crossprod(Q,y))

##Excercise - Week 2
X <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
rownames(X) <- c("a","a","b","b")

beta <- c(5, 2)
beta

X%*%beta

X <- matrix(c(1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1),nrow=6)
rownames(X) <- c("a","a","b","b","c","c")
beta <- c(10,3,-3)
X%*%beta


##Inference - Week 2

g = 9.8 ## meters per second
h0 = 56.67
v0 = 0
n = 25
tt = seq(0,3.4,len=n) ##time in secs, t is a base function
y = h0 + v0 *tt - 0.5* g*tt^2 + rnorm(n,sd=1)

X = cbind(1,tt,tt^2)
A = solve(crossprod(X))%*%t(X)
X
A

A%*%y
-2 * (A %*% y) [3] 

set.seed(1)
temp <- replicate(100000,-2 * (A %*% (h0 + v0 *tt - 0.5* g*tt^2 + rnorm(n,sd=1))) [3] ) 
sd(temp)/sqrt(100000)

B = 100000
set.seed(1)
betahat = replicate(B,{y = 56.67 - 0.5*g*tt^2 + rnorm(n,sd=1);betahats = -2*A%*%y;return(betahats[3]) })
sqrt(mean( (betahat-mean(betahat) )^2))

library(UsingR)
x = father.son$fheight
y = father.son$sheight
n = length(y)

N =  50
index = sample(n,N)
sampledat = father.son[index,]
x = sampledat$fheight
y = sampledat$sheight
betahat =  lm(y~x)$coef

betahats <- replicate(10000,{
  index = sample(n,N)
  sampledat = father.son[index,]
  x = sampledat$fheight
  y = sampledat$sheight
  betahat = lm(y~x)$coef
  return(betahat[2])
})

sd(betahats)
mean((x - mean(x))*(y-mean(y)))


##Excercise Week 2 - #1
library(UsingR)
x = father.son$fheight
y = father.son$sheight
n = length(y)
N = 50
set.seed(1)
index = sample(n,N)
sampledat = father.son[index,]
x = sampledat$fheight
y = sampledat$sheight
betahat = lm(y~x)$coef
#The formula for the standard error in the previous video was (the following two lines are not R code):
SE(betahat) = sqrt(var(betahat))
var(betahat) = sigma^2 (x^T x)^-1

set.seed(1)
n <- nrow(father.son)
N <- 50
index <- sample(n,N)
sampledat <- father.son[index,]
x <- sampledat$fheight
y <- sampledat$sheight
var(y)
fit = lm(y ~ x)
y_hat <- fit$fitted.values
ssr <- sum((y - y_hat)^2)
#sigma2 = sum((y-fit$fitted.values)^2) / 48

sigma2 = ssr / 48

##Excercise Week 2 - #2
X = cbind(rep(1,N), x)

Xinv <- solve(crossprod(X))
Xinv[1,1]
##Excercise Week 2 - #3
d <- diag(Xinv)
beta_hat <- sqrt(sigma2 * d)
beta_hat[2]

