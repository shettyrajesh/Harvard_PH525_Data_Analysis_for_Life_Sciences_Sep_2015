## Formula
## Model Matrix
## inorder to produce design matrices for liner models

x <- c(1,1,2,2)
f <- formula(~ x)
f

model.matrix(f)

x <- factor(c(1,1,2,2))
model.matrix(~ x)

x <- factor(c(1,1,2,2,3,3))
model.matrix(~ x)
model.matrix(~ x, contrasts=list(x="contr.sum"))
?contr.sum

x <- factor(c(1,1,1,1,2,2,2,2))
y <- factor(c("a","a","b","b","a","a","b","b"))
model.matrix(~ x + y)

model.matrix(~ x + y+ x :y)
model.matrix(~ x * y)

x <- factor(c(1,1,2,2))
model.matrix(~ x)
x <- relevel(x, "2")
model.matrix(~ x)
x <- factor(x, levels=c("1","2"))


z <- 1:4
model.matrix(~ z)
model.matrix(~ 0 + z)
model.matrix(~ z + I(z^2))


##***************************Video 3

install.packages('downloader')
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)

dat <- read.csv("femaleMiceWeights.csv") ##previously downloaded
stripchart(dat$Bodyweight ~ dat$Diet, vertical=TRUE, method="jitter", main="Bodyweight over Diet")

levels(dat$Diet)
X <- model.matrix(~ Diet, data=dat)
X
colnames(X)
Y <- dat$Bodyweight
X <- model.matrix(~ Diet, data=dat)
X
solve(t(X) %*% X) %*% t(X) %*% Y

s <- split(dat$Bodyweight, dat$Diet)
mean(s[["chow"]])
mean(s[["hf"]]) - mean(s[["chow"]])

fit <- lm(Bodyweight ~ Diet, data=dat)
summary(fit)
(coefs <- coef(fit))

stripchart(dat$Bodyweight ~ dat$Diet, vertical=TRUE, method="jitter",
           main="Bodyweight over Diet", ylim=c(0,40), xlim=c(0,3))
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

summary(fit)$coefficients

ttest <- t.test(s[["hf"]], s[["chow"]], var.equal=TRUE)
summary(fit)$coefficients[2,3]
ttest$statistic
ttest <- t.test(s[["chow"]], s[["hf"]], var.equal=TRUE)
ttest$statistic

## Excercise 2
set.seed(1)
nx <- 5
ny <- 7
X = cbind(rep(1,nx + ny),rep(c(0,1),c(nx, ny)))
zz <- (t(X) %*% X)
zz
crossprod(X)
