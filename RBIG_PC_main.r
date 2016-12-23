rm(list=ls())
library(pcalg)
source("MI_RBIG_2016_algo.R")
source("RBIG_r.r")
source('Npres_Fucntions.R')

cmi_btest_pc <- function(x, y, S, suffStat){
  print("**************************************")
  cat("x=", x, " / y=", y, "/ S=", S, "\n") 
  cmi_btest(suffStat$nboot, suffStat$dat, x, y, S)
}

npres_pc <- function(x, y, S, suffStat){
  dat <- suffStat$dat
  R <- suffStat$R
  d <- length(S) 
  if(d == 0){
    p.value <- cramer.test(dat[, x], dat[, y], sim="ordinary", replicates=R)$p.value
  }else{
    dat <- dat[, c(x,y,S)] 
    test.stat <- npresid.statistics(dat,d)
    out <- npresid.boot(dat,d,R=R)
    p.value <- out$p.value
  }
  p.value
}


kcipt_pc <- function(x, y, S, suffStat){
  dat <- suffStat$dat
  B <- suffStat$B
  b <- suffStat$b
  M <- suffStat$M
  if(length(S) == 0){
    p.value <- cramer.test(dat[, x], dat[, y], sim="ordinary", replicates=M)$p.value
  }else{
    p.value <- KCIPT(dat,xy_ind=c(x,y), c_ind=S,  dist=dist, B, b, M)
  }
  p.value
}

rbig_kcipt_pc <- function(x, y, S, suffStat){
  dat <- suffStat$dat
  B <- suffStat$B
  b <- suffStat$b
  M <- suffStat$M
  if(length(S) == 0){
    p.value <- cramer.test(dat[, x], dat[, y], sim="ordinary", replicates=M)$p.value
  }else{
    p.value <-RBIG_kcipt(dat,xy_ind=c(x,y), c_ind=S,  dist=dist, B, b, M)
  }
  p.value
}

n <- 1000
n <- 700
x1 <- rnorm(n)
x2 <- rnorm(n)
y1 <- matrix(ncol=2, nrow=n)
y1[, 1] <- cos(x1) + 0.1*x2
# y1[, 2] <- 4*x1.^2 + 1*x2
y1[, 2] <- sin(x1) + 0.1*x2
# y1[, 3] <- 0.01*((2*x1)^3 + 0.1*x2)
# y2 <- 10*(abs(y1[, 1]))^0.5 + y1[, 2]^2 + 0*y1[,3] + 0.3*x2
y2 <- y1[, 1]^3 + y1[, 2]^3 + 0.3*x2
dat <- as.data.frame(cbind(x1, x2, y1, y2))

#1/ OK it works for this one
n <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
y1 <- cos(x1) + 0.1*x2
dat <- as.data.frame(cbind(x1, x2, y1))
pc.fit1 <- pc(suffStat = list(dat=as.matrix(dat), nboot=50), indepTest = cmi_btest_pc, alpha=0.2, labels = names(dat), verbose = TRUE)

#2/ Not working with alpha=0.2 for RBIG(mine) 
#2/ Works with alpha=0.2 for KCIPT 
n <- 1000
x1 <- rnorm(n)
y1 <- cos(x1) + 0.1*rnorm(n)
y2 <- cos(y1) + 0.1*rnorm(n)
dat <- as.data.frame(cbind(x1, y1, y2))
plot(hexplom(dat))
pc.fit2 <- pc(suffStat = list(dat=as.matrix(dat), nboot=50), indepTest = cmi_btest_pc, alpha=0.2, labels = names(dat), verbose = TRUE)

#3/ Works for y2=-2*y1 + 4*rnorm and KCIPT 
n <- 2000
x1 <- rnorm(n)
y1 <- 2 * x1 + 5*rnorm(n)
y2 <- -2 * y1 + 1*rnorm(n)
dat <- as.data.frame(cbind(x1, y1, y2))
plot(hexplom(dat))
pc.fit3 <- pc(suffStat = list(dat=as.matrix(dat), nboot=50), indepTest = cmi_btest_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
cmi_btest(nboot=500, dat, 1, 2, 3)
pc.fit3.cor <- pc(suffStat = list(C=cor(dat), n=nrow(dat)), indepTest = gaussCItest, alpha=0.05, labels = names(dat), verbose = TRUE)
pc.fit3.kcipt <- pc(suffStat = list(dat=as.matrix(head(dat, 700)), B=20, b=50, M=100), indepTest = kcipt_pc, alpha=0.2, labels = names(dat), verbose = TRUE)

npres_pc(1, 2, 3, suffStat = list(dat=as.matrix(head(dat, 700)), R=100))
pc.fit3.npres <- pc(suffStat = list(dat=as.matrix(head(dat, 700)), R=1000), indepTest = npres_pc , alpha=0.2, labels = names(dat), verbose = TRUE)

set.seed(1)
p.value_1 <- KCIPT(head(dat, 700), xy_ind=c(1,2), c_ind=3,  dist=dist, B=25, b=1000, M=10000)
set.seed(1)
p.value_2 <- KCIPT(head(dat, 700), xy_ind=c(2,1), c_ind=3,  dist=dist, B=25, b=1000, M=10000)

vartest <- sapply(1:100, function(x){
		    set.seed(x)
		    p.value_1 <- KCIPT(head(dat, 700), xy_ind=c(1,2), c_ind=3,  dist=dist, B=10, b=10000, M=10000)
		    set.seed(x)
		    p.value_2 <- KCIPT(head(dat, 700), xy_ind=c(2,1), c_ind=3,  dist=dist, B=10, b=10000, M=10000)
		    print(c(p.value_1, p.value_2))
		    c(p.value_1, p.value_2)
})

p.value_1 <- KCIPT_outeronly(head(dat, 700), xy_ind=c(1,2), c_ind=3,  dist=dist, B=100)
p.value_2 <- KCIPT_outeronly(head(dat, 700), xy_ind=c(2,1), c_ind=3,  dist=dist, B=100)

inner_1 <- KCIPT_inneronly(head(dat, 700), xy_ind=c(1,2), c_ind=3,  dist=dist, B=100, b=10000, M=10000)

p.value_11 <- KCIPT(head(dat, 700), xy_ind=c(1,2), c_ind=3,  dist=dist, B=25, b=1000, M=10000)
p.value_12 <- KCIPT(head(dat, 700), xy_ind=c(1,2), c_ind=3,  dist=dist, B=500, b=10000, M=10000)
p.value_21 <- KCIPT(head(dat, 700), xy_ind=c(2,1), c_ind=3,  dist=dist, B=25, b=1000, M=10000)
p.value_22 <- KCIPT(head(dat, 700), xy_ind=c(2,1), c_ind=3,  dist=dist, B=500, b=10000, M=10000)



norm_a <- runif(100000)
norm_b <- runif(100000)
norm_diff <- norm_b - norm_a
for(i in 1:100000) print(summary(norm_diff[1:i]))

set.seed(1)
test_rbig_kcipt <- lapply(1:10, function(x) RBIG_kcipt(head(dat, 700), 1:2, 3:4, dist, 50*x, 50*x, 10000))
save(test_rbig_kcipt, file="test_rbig_kcipt.Rdata")
