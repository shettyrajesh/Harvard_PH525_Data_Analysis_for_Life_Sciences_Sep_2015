url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)

table(spider$leg,spider$type)

boxplot(spider$friction ~ spider$type * spider$leg, 
        col=c("grey90","grey40"), las=2, 
        main="Comparison of friction coefficients of different leg pairs")

spider.sub <- spider[spider$leg == "L1",]
fit <- lm(friction ~ type, data=spider.sub)
summary(fit)
(coefs <- coef(fit))

s <- split(spider.sub$friction, spider.sub$type)
mean(s[["pull"]])
mean(s[["push"]]) - mean(s[["pull"]])

X <- model.matrix(~ type, data=spider.sub)
colnames(X)
head(X)
tail(X)

install.packages("rafalib")
#library(devtools)
#install_github("ririzarr/rafalib")
library(rafalib)
imagemat(X, main="Model matrix for linear model with one variable")

set.seed(1) #same jitter in stripchart
stripchart(split(spider.sub$friction, spider.sub$type), 
           vertical=TRUE, pch=1, method="jitter", las=2, xlim=c(0,3), ylim=c(0,2))
a <- -0.25
lgth <- .1
library(RColorBrewer)
cols <- brewer.pal(3,"Dark2")
abline(h=0)
arrows(1+a,0,1+a,coefs[1],lwd=3,col=cols[1],length=lgth)
abline(h=coefs[1],col=cols[1])
arrows(2+a,coefs[1],2+a,coefs[1]+coefs[2],lwd=3,col=cols[2],length=lgth)
abline(h=coefs[1]+coefs[2],col=cols[2])
legend("right",names(coefs),fill=cols,cex=.75,bg="white")

X <- model.matrix(~ type + leg, data=spider)
colnames(X)
head(X)
imagemat(X, main="Model matrix for linear model with two factors")

fitTL <- lm(friction ~ type + leg, data=spider)
summary(fitTL)
(coefs <- coef(fitTL))

Y <- spider$friction
X <- model.matrix(~ type + leg, data=spider)
beta.hat <- solve(t(X) %*% X) %*% t(X) %*% Y
t(beta.hat)
coefs

spider$group <- factor(paste0(spider$leg, spider$type))
stripchart(split(spider$friction, spider$group), 
           vertical=TRUE, pch=1, method="jitter", las=2, xlim=c(0,11), ylim=c(0,2))
cols <- brewer.pal(5,"Dark2")
abline(h=0)
arrows(1+a,0,1+a,coefs[1],lwd=3,col=cols[1],length=lgth)
abline(h=coefs[1],col=cols[1])
arrows(3+a,coefs[1],3+a,coefs[1]+coefs[3],lwd=3,col=cols[3],length=lgth)
arrows(5+a,coefs[1],5+a,coefs[1]+coefs[4],lwd=3,col=cols[4],length=lgth)
arrows(7+a,coefs[1],7+a,coefs[1]+coefs[5],lwd=3,col=cols[5],length=lgth)
arrows(2+a,coefs[1],2+a,coefs[1]+coefs[2],lwd=3,col=cols[2],length=lgth)
segments(3+a,coefs[1]+coefs[3],4+a,coefs[1]+coefs[3],lwd=3,col=cols[3])
arrows(4+a,coefs[1]+coefs[3],4+a,coefs[1]+coefs[3]+coefs[2],lwd=3,col=cols[2],length=lgth)
segments(5+a,coefs[1]+coefs[4],6+a,coefs[1]+coefs[4],lwd=3,col=cols[4])
arrows(6+a,coefs[1]+coefs[4],6+a,coefs[1]+coefs[4]+coefs[2],lwd=3,col=cols[2],length=lgth)
segments(7+a,coefs[1]+coefs[5],8+a,coefs[1]+coefs[5],lwd=3,col=cols[5])
arrows(8+a,coefs[1]+coefs[5],8+a,coefs[1]+coefs[5]+coefs[2],lwd=3,col=cols[2],length=lgth)
legend("right",names(coefs),fill=cols,cex=.75,bg="white")

s <- split(spider$friction, spider$group)
mean(s[["L1pull"]])
coefs[1]
mean(s[["L1push"]])
coefs[1] + coefs[2]

means <- sapply(s, mean)
##the sample size of push or pull groups for each leg pair
ns <- sapply(s, length)[c(1,3,5,7)]
(w <- ns/sum(ns))
sum(w * (means[c(2,4,6,8)] - means[c(1,3,5,7)]))
coefs[2]

install.packages("contrast")
library(contrast) #Available from CRAN
L3vsL2 <- contrast(fitTL,list(leg="L3",type="pull"),list(leg="L2",type="pull"))
L3vsL2

coefs[4] - coefs[3]
(cT <- L3vsL2$X)
cT %*% coefs

Sigma.hat <- sum(fitTL$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X)
signif(Sigma.hat, 2)
sqrt(cT %*% Sigma.hat %*% t(cT))
L3vsL2$SE

#We would have obtained the same result for a contrast of L3 and L2 had we picked type="push". The reason it does not change 
#the contrast is because it leads to addition of the typepush effect on both sides of the difference, which cancels out:
  
L3vsL2.equiv <- contrast(fitTL,list(leg="L3",type="push"),list(leg="L2",type="push"))
L3vsL2.equiv
L3vsL2.equiv$X

# Excercise 4.1.1 *************************

species <- factor(c("A","A","B","B"))
condition <- factor(c("control","treated","control","treated"))
model.matrix(~ species + condition)
library(contrast)
y = rnorm(4)
fit = lm(y ~ species + condition)
contrast(fit, list(species="B",condition="control"), list(species="A",condition="treated"))$X

# Excercise 4.1.2 ************

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)
fitTL <- lm(friction ~ type + leg, data=spider)

contrast(fitTL,list(leg="L4",type="pull"),list(leg="L2",type="pull"))
L4vsL2 <- contrast(fitTL, list(leg="L4",type="pull"), list(leg="L2",type="pull"))$X
L4vsL2
# Excercise 4.1.3 ***************************
X <- model.matrix(~ type + leg, data=spider)
(Sigma <- sum(fitTL$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X))

#Our contrast matrix is:
  
C <- matrix(c(0,0,-1,0,1),1,5)
C
sqrt(C %*% Sigma %*% t(C))

#***************************

#Linear Model with Interactions
X <- model.matrix(~ type + leg + type:leg, data=spider)
colnames(X)
head(X)
imagemat(X, main="Model matrix for linear model with interactions")
fitX <- lm(friction ~ type + leg + type:leg, data=spider)
summary(fitX)
coefs <- coef(fitX)
stripchart(split(spider$friction, spider$group), vertical=TRUE, pch=1, method="jitter", las=2, xlim=c(0,11), ylim=c(0,2))
cols <- brewer.pal(8,"Dark2")
abline(h=0)
arrows(1+a,0,1+a,coefs[1],lwd=3,col=cols[1],length=lgth)
abline(h=coefs[1],col=cols[1])
arrows(2+a,coefs[1],2+a,coefs[1]+coefs[2],lwd=3,col=cols[2],length=lgth)
arrows(3+a,coefs[1],3+a,coefs[1]+coefs[3],lwd=3,col=cols[3],length=lgth)
arrows(5+a,coefs[1],5+a,coefs[1]+coefs[4],lwd=3,col=cols[4],length=lgth)
arrows(7+a,coefs[1],7+a,coefs[1]+coefs[5],lwd=3,col=cols[5],length=lgth)
#now the interactions:
segments(3+a,coefs[1]+coefs[3],4+a,coefs[1]+coefs[3],lwd=3,col=cols[3])
arrows(4+a,coefs[1]+coefs[3],4+a,coefs[1]+coefs[3]+coefs[2],lwd=3,col=cols[2],length=lgth)
arrows(4+a,coefs[1]+coefs[2]+coefs[3],4+a,coefs[1]+coefs[2]+coefs[3]+coefs[6],lwd=3,col=cols[6],length=lgth)

segments(5+a,coefs[1]+coefs[4],6+a,coefs[1]+coefs[4],lwd=3,col=cols[4])
arrows(6+a,coefs[1]+coefs[4],6+a,coefs[1]+coefs[4]+coefs[2],lwd=3,col=cols[2],length=lgth)
arrows(6+a,coefs[1]+coefs[4]+coefs[2],6+a,coefs[1]+coefs[4]+coefs[2]+coefs[7],lwd=3,col=cols[7],length=lgth)

segments(7+a,coefs[1]+coefs[5],8+a,coefs[1]+coefs[5],lwd=3,col=cols[5])
arrows(8+a,coefs[1]+coefs[5],8+a,coefs[1]+coefs[5]+coefs[2],lwd=3,col=cols[2],length=lgth)
arrows(8+a,coefs[1]+coefs[5]+coefs[2],8+a,coefs[1]+coefs[5]+coefs[2]+coefs[8],lwd=3,col=cols[8],length=lgth)
legend("right",names(coefs),fill=cols,cex=.75,bg="white")

library(contrast) ##Available from CRAN
L2push.vs.pull <- contrast(fitX,
                           list(leg="L2", type = "push"), 
                           list(leg="L2", type = "pull"))
L2push.vs.pull
coefs[2] + coefs[6] ##we know this is also orange + yellow arrow

install.packages("TH.data")
install.packages("mvtnorm")
library(multcomp) ##Available from CRAN
C <- matrix(c(0,0,0,0,0,-1,1,0), 1)
L3vsL2interaction <- glht(fitX, linfct=C)
summary(L3vsL2interaction)
coefs[7] - coefs[6] ##we know this is also brown - yellow

anova(fitX)
mu0 <- mean(spider$friction)
(initial.ss <- sum((spider$friction - mu0)^2))
N <- nrow(spider)
(N - 1) * var(spider$friction)
s <- split(spider$friction, spider$type)
after.type.ss <- sum( sapply(s, function(x) {
  residual <- x - mean(x) 
  sum(residual^2)
}) )

(type.ss <- initial.ss - after.type.ss)
sum(sapply(s, length) * (sapply(s, mean) - mu0)^2)

# Excercise 4.1.4 ***************************
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)

spider$log2friction <- log2(spider$friction)
spider$log2friction
boxplot(log2friction ~ type*leg, data=spider)


fitX <- lm(log2friction ~ type + leg + type:leg, data=spider)
summary(fitX)
anova(fitX)

fitX <- lm(log2friction ~ type + leg + type:leg-1, data=spider)
summary(fitX)
contrast(fitX, list(type="pull",leg="L2"), list(type="pull",leg="L1"))
contrast(fitX, list(type="push",leg="L2"), list(type="push",leg="L1"))
#******************************************

N <- 40
p <- 4
group <- factor(rep(1:p,each=N/p))
X <- model.matrix(~ group)
Y <- rnorm(N,mean=42,7)
mu0 <- mean(Y)
initial.ss <- sum((Y - mu0)^2)
s <- split(Y, group)
after.group.ss <- sum(sapply(s, function(x) sum((x - mean(x))^2)))
(group.ss <- initial.ss - after.group.ss)
group.ms <- group.ss / (p - 1)
after.group.ms <- after.group.ss / (N - p)
f.value <- group.ms / after.group.ms

y <- replicate(1000,{
  Y <- rnorm(N,mean=42,7)
  mu0 <- mean(Y)
  initial.ss <- sum((Y - mu0)^2)
  s <- split(Y, group)
  after.group.ss <- sum(sapply(s, function(x) sum((x - mean(x))^2)))
  group.ss <- initial.ss - after.group.ss
  group.ms <- group.ss / (p - 1)
  after.group.ms <- after.group.ss / (N - p)
  f.value <- group.ms / after.group.ms
  return(f.value)
})

hist(y, col="grey", border="white", breaks=50, freq=FALSE)
xs <- seq(from=0,to=6,length=100)
lines(xs, df(xs, df1 = p - 1, df2 = N - p), col="red")

#******************************************


##earlier, we defined the 'group' column:
spider$group <- factor(paste0(spider$leg, spider$type))
X <- model.matrix(~ 0 + group, data=spider)
colnames(X)
head(X)
imagemat(X, main="Model matrix for linear model with group variable")

#We can run the linear model with the familiar call:
  
fitG <- lm(friction ~ 0 + group, data=spider)
summary(fitG)
coefs <- coef(fitG)

#Examining the estimated coefficients

#Now we have eight arrows, one for each group. The arrow tips align directly with the mean of each group:
  
stripchart(split(spider$friction, spider$group), vertical=TRUE, pch=1, method="jitter", las=2, xlim=c(0,11), ylim=c(0,2))
cols <- brewer.pal(8,"Dark2")
abline(h=0)
for (i in 1:8) {
  arrows(i+a,0,i+a,coefs[i],lwd=3,col=cols[i],length=lgth)
}
legend("right",names(coefs),fill=cols,cex=.75,bg="white")

groupL2push.vs.pull <- contrast(fitG,
                                list(group = "L2push"), 
                                list(group = "L2pull"))
groupL2push.vs.pull
coefs[4] - coefs[3]

C <- matrix(c(0,0,1,-1,-1,1,0,0), 1)
groupL3vsL2interaction <- glht(fitG, linfct=C)
summary(groupL3vsL2interaction)
names(coefs)
(coefs[6] - coefs[5]) - (coefs[4] - coefs[3])
#******************************************start of 4.2
  
Sex <- c(0,0,0,0,1,1,1,1)
A <-   c(1,1,0,0,0,0,0,0)
B <-   c(0,0,1,1,0,0,0,0)
C <-   c(0,0,0,0,1,1,0,0)
D <-   c(0,0,0,0,0,0,1,1)
X <- model.matrix(~Sex+A+B+C+D-1)
cat("ncol=",ncol(X),"rank=", qr(X)$rank,"\n")

Sex <- c(0,1,0,1,0,1,0,1)
A <-   c(1,1,0,0,0,0,0,0)
B <-   c(0,0,1,1,0,0,0,0)
C <-   c(0,0,0,0,1,1,0,0)
D <-   c(0,0,0,0,0,0,1,1)
X <- model.matrix(~Sex+A+B+C+D-1)
cat("ncol=",ncol(X),"rank=", qr(X)$rank,"\n")
# Excercise 4.2.1 *************************

sex <- factor(rep(c("female","male"),each=4))
trt <- factor(c("A","A","B","B","C","C","D","D"))
X <- model.matrix( ~ sex + trt)
qr(X)$rank
Y <- 1:8
makeYstar <- function(a,b) Y - X[,2] * a - X[,5] * b
fitTheRest <- function(a,b) {
  Ystar <- makeYstar(a,b)
  Xrest <- X[,-c(2,5)]
  betarest <- solve(t(Xrest) %*% Xrest) %*% t(Xrest) %*% Ystar
  residuals <- Ystar - Xrest %*% betarest
  sum(residuals^2)
}
fitTheRest(1,2)

expand.grid(1:3,1:3)
betas = expand.grid(-2:8,-2:8)

outer(1:3,1:3,`*`)
min(outer(-2:8,-2:8,Vectorize(fitTheRest)))
rss = apply(betas,1,function(x) fitTheRest(x[1],x[2]))

library(rafalib)
themin=min(rss) 
plot(betas[which(rss==themin),])
imagemat(outer(-2:8,-2:8,Vectorize(fitTheRest)))

fit <- lm(friction ~ type + leg, data=spider)
betahat <- coef(fit)
Y <- matrix(spider$friction, ncol=1)
X <- model.matrix(~ type + leg, data=spider)

QR <- qr(X)
Q <- qr.Q( QR )
Q[1,1]
R <- qr.R( QR )
R[1,1]

t(Q)%*%Y
solve(R)%*%t(Q)%*%Y
#*******************************************************
# Section 4.2*******************************************

n <- 50;M <- 500
x <- seq(1,M,len=n)
X <- cbind(1,x,x^2,x^3)
colnames(X) <- c("Intercept","x","x2","x3")
beta <- matrix(c(1,1,1,1),4,1)
set.seed(1)
y <- X%*%beta+rnorm(n,sd=1)
plot(x,y)

solve(crossprod(X))

solve(crossprod(X)) %*% crossprod(X,y)
options(digits=4)
log10(crossprod(X))

QR <- qr(X)
Q <- qr.Q( QR )
R <- qr.R( QR )
(betahat <- backsolve(R, crossprod(Q,y) ) )

QR <- qr(X)
(betahat <- solve.qr(QR, y))

library(rafalib)
mypar(1,1)
plot(x,y)
fitted <- tcrossprod(Q)%*%y
lines(x,fitted,col=2)

df <- length(y) - QR$rank
sigma2 <- sum((y-fitted)^2)/df
varbeta <- sigma2*chol2inv(qr.R(QR))
SE <- sqrt(diag(varbeta))
cbind(betahat,SE)

summary(lm(y~0+X))$coef

#Excercise 4.2*******************************

fit <- lm(friction ~ type + leg, data=spider)
betahat <- coef(fit)
Y <- matrix(spider$friction, ncol=1)
X <- model.matrix(~ type + leg, data=spider)

QR <- qr(X)
Q <- qr.Q( QR )
Q[1,1]
R <- qr.R( QR )
R[1,1]

t(Q)%*%Y
solve(R)%*%t(Q)%*%Y
