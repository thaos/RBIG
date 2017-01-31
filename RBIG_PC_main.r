rm(list=ls())
library(afex)
library(pcalg)
library(ggplot2)
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
  # FALSE ERROR
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
  # FALSE ERROR
    p.value <- cramer.test(dat[, x], dat[, y], sim="ordinary", replicates=M)$p.value
  }else{
    p.value <-RBIG_kcipt(dat,xy_ind=c(x,y), c_ind=S,  dist=dist, B, b, M)
  }
  p.value
}

rbig_hist_pc <- function(x, y, S, suffStat){
  dat <- suffStat$dat
  B <- suffStat$B
  to_plot <- !is.null(suffStat$to_plot)
  p.value <-RBIG_hist(dat,xy_ind=c(x,y), c_ind=S,  dist=dist, B=B, to_plot=to_plot)
    #     p.value <- cramer.test(p.value[, 1], p.value[, 2])$p.value
  p.value <- ks.test(p.value[, 1], p.value[, 2])$p.value
  p.value
}

rbig_simple_pc <- function(x, y, S, suffStat){
  dat <- suffStat$dat
  B <- suffStat$B
  to_plot <- !is.null(suffStat$to_plot)
  p.value <-RBIG_simple(dat,xy_ind=c(x,y), c_ind=S,  dist=dist, B=B, to_plot=to_plot)
  p.value
}

# hybridation linear non linear pas bonne idée!
kcipt_test_pc <- function(x, y, S, suffStat){
  dat <- suffStat$dat
  B <- suffStat$B
  to_plot <- !is.null(suffStat$to_plot)
  #   alpha <- suffStat$alpha
  #   p.value <- gaussCItest(x, y, S, list("C"=C, "n"=nrow(dat)))
  #   print(p.value)
  #   if(p.value > alpha){
  p.value <- KCIPT_test(dat,xy_ind=c(x,y), c_ind=S,  dist=dist, B=B, to_plot=to_plot)
  #   p.value <- cramer.test(p.value[, 1], p.value[, 2], sim="eigenvalue")$p.value
  p.value <- ks.test(p.value[[1]], p.value[[2]])$p.value
  #   p.value <- t.test(p.value[, 1], p.value[, 2])$p.value
  #   }
  p.value
}

#0/ almost work with rbig_hist_pc with B but missing 4 to 5 
#0/ work with rbig_hist_pc on second try 
set.seed(1)
n <- 500 
x1 <- rnorm(n)
x2 <- rnorm(n)
y1 <- matrix(ncol=2, nrow=n)
y1[, 1] <- cos(x1) + 0.1*x2
# y1[, 2] <- 4*x1.^2 + 1*x2
y1[, 2] <- sin(x1) + 0.1*x2
# y1[, 3] <- 0.01*((2*x1)^3 + 0.1*x2)
# y2 <- 10*(abs(y1[, 1]))^0.5 + y1[, 2]^2 + 0*y1[,3] + 0.3*x2
y2 <- y1[, 1]^3 + y1[, 2]^3 + 0.3*x2
dat <- as.data.frame(scale(cbind(x1, x2, y1, y2)))
plot(hexplom(dat))
pc.fit0.rbig_hist <- pc(suffStat = list(dat=as.matrix(head(dat, 100)), B=5000), indepTest = rbig_hist_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
rbig_hist_pc(x=5, y=4, S=1:3, suffStat = list(dat=as.matrix(head(dat, 700)), B=100))
pc.fit0.kcipt_test <- pc(suffStat = list(dat=as.matrix(dat), B=n, to_plot=TRUE), indepTest = kcipt_test_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
rbig_hist_pc(x=3, y=2, S=c(1, 4) , suffStat = list(dat=as.matrix(head(dat, 700)), B=100))

#1/ OK it works for this one with rbig
#1/ OK it works for this one with rbig hist_pc
set.seed(1) 
n <- 250 
x1 <- rnorm(n)
x2 <- rnorm(n)
y1 <- cos(x1) + 0.1*x2
dat <- as.data.frame(scale(cbind(x1, x2, y1)))
pc.fit1 <- pc(suffStat = list(dat=as.matrix(dat), nboot=50), indepTest = cmi_btest_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
pc.fit1.rbig_hist <- pc(suffStat = list(dat=as.matrix(head(dat, 700)), B=1000), indepTest = rbig_hist_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
pc.fit1.rbig_simple <- pc(suffStat = list(dat=as.matrix(dat), B=n, to_plot=TRUE), indepTest = rbig_simple_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
save(pc.fit0.rbig_hist, file="pc.fit0.rbig_hist.RData")
pc.fit1.kcipt_test <- pc(suffStat = list(dat=as.matrix(dat), B=n, to_plot=TRUE), indepTest = kcipt_test_pc, alpha=0.2, labels = names(dat), verbose = TRUE)

#2/ Not working with alpha=0.2 for RBIG(mine) 
#2/ Works with alpha=0.2 for KCIPT 
#2/ Works with alpha=0.2 for rbig_hist with B=100 
set.seed(1)
n <- 200 
x1 <- rnorm(n)
y1 <- cos(x1) + 0.1*rnorm(n)
y2 <- cos(y1) + 0.1*rnorm(n)
dat <- as.data.frame(cbind(x1, y1, y2))
dat <- as.data.frame(scale(dat))
plot(hexplom(dat))
pc.fit2 <- pc(suffStat = list(dat=as.matrix(dat), nboot=50), indepTest = cmi_btest_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
pc.fit2.kcipt <- pc(suffStat = list(dat=as.matrix(head(dat, 700)), B=20, b=50, M=100), indepTest = kcipt_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
pc.fit2.rbig_hist <- pc(suffStat = list(dat=as.matrix(head(dat, 100)), B=100), indepTest = rbig_hist_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
save(pc.fit2.rbig_hist, file="pc.fit2.rbig_hist.RData")
pc.fit2.kcipt_test <- pc(suffStat = list(dat=dat, B=n), indepTest = kcipt_test_pc, alpha=0.2, labels = names(dat), verbose = TRUE)

#3/ Works for y2=-2*y1 + 4*rnorm and KCIPT 
#2/ Works with alpha=0.2 for rbig_hist though test is not always symetric 
n <- 200
x1 <- rnorm(n)
y1 <- 2 * x1 + 1*rnorm(n)
y2 <- -2 * y1 + 1*rnorm(n)
dat <- as.data.frame(scale(cbind(x1, y1, y2)))
plot(hexplom(dat))
pc.fit3 <- pc(suffStat = list(dat=as.matrix(dat), nboot=50), indepTest = cmi_btest_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
cmi_btest(nboot=500, dat, 1, 2, 3)
pc.fit3.cor <- pc(suffStat = list(C=cor(dat), n=nrow(dat)), indepTest = gaussCItest, alpha=0.05, labels = names(dat), verbose = TRUE)
pc.fit3.kcipt <- pc(suffStat = list(dat=as.matrix(head(dat, 700)), B=20, b=50, M=100), indepTest = kcipt_pc, alpha=0.2, labels = names(dat), verbose = TRUE)

pc.fit3.rbig_hist <- pc(suffStat = list(dat=as.matrix(head(dat, 700)), B=1000), indepTest = rbig_hist_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
save(pc.fit3.rbig_hist, file="pc.fit3.rbig_hist.RData")

npres_pc(1, 2, 3, suffStat = list(dat=as.matrix(head(dat, 700)), R=100))
pc.fit3.npres <- pc(suffStat = list(dat=as.matrix(head(dat, 700)), R=1000), indepTest = npres_pc , alpha=0.2, labels = names(dat), verbose = TRUE)

pc.fit3.kcipt_test <- pc(suffStat = list(dat=as.matrix(dat), B=n, alpha=0.2, C=cor(as.matrix(head(dat, 700)))), indepTest = kcipt_test_pc, alpha=0.2, labels = names(dat), verbose = TRUE)
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

# works sometimes ....
set.seed(1)
test_rbig_kcipt <- lapply(1:10, function(x) RBIG_kcipt(head(dat, 700), 1:2, 3:4, dist, 50*x, 50*x, 10000))
save(test_rbig_kcipt, file="test_rbig_kcipt.Rdata")

set.seed(1)
cp_rbig_hist <- lapply(1:50, function(x){
			    n <- 2000
			    x1 <- rnorm(n)
			    y1 <- 2 * x1 + 5*rnorm(n)
			    y2 <- -2 * y1 + 1*rnorm(n)
			    dat <- as.data.frame(cbind(x1, y1, y2))
			    plot(hexplom(dat))
			    RBIG_hist(head(dat, 700), c(1,3), 2, dist, 500)
})
save(cp_rbig_hist, file="cp_rbig_hist.Rdata")
l_pvalues <- unlist(lapply(cp_rbig_hist, function(x)cramer.test(x[,1], x[,2])$p.value))

set.seed(2)
cp_rbig_hist2 <- lapply(1:50, function(x){
			    n <- 2000
			    x1 <- rnorm(n)
			    y1 <- 2 * x1 + 5*rnorm(n)
			    y2 <- -2 * y1 + 1*rnorm(n)
			    dat <- as.data.frame(cbind(x1, y1, y2))
			    plot(hexplom(dat))
			    RBIG_hist(head(dat, 700), c(1,3), 2, dist, 100)
})
save(cp_rbig_hist2, file="cp_rbig_hist2.Rdata")
l_pvalues <- unlist(lapply(cp_rbig_hist2, function(x)cramer.test(x[,1], x[,2])$p.value))

set.seed(1)
cp_rbig_hist_unif <- lapply(1:50, function(x){
			    n <- 2000
			    x1 <- runif(n)
			    y1 <- 2 * x1 + 5*runif(n)
			    y2 <- -2 * y1 + 1*runif(n)
			    dat <- as.data.frame(cbind(x1, y1, y2))
			    plot(hexplom(dat))
			    RBIG_hist(head(dat, 700), c(1,3), 2, dist, 500)
})
save(cp_rbig_hist_unif, file="cp_rbig_hist_unif.Rdata")
l_pvalues <- unlist(lapply(cp_rbig_hist_unif, function(x)cramer.test(x[,1], x[,2])$p.value))

n <- 700
x1 <- rnorm(n)
y1 <- 2 * x1 + 5*rnorm(n)
y2 <- -2 * y1 + 1*rnorm(n)
dat <- as.data.frame(cbind(x1, y1, y2))
plot(hexplom(dat))
rbig_conv <- RBIG_hist(head(dat, 700), c(1,3), 2, dist, 10000)
save(rbig_conv, file="rbig_conv.Rdata")

rbig_conv_ggplot <- data.frame("stat"=unlist(rbig_conv), "type"=rep(c("H1", "H0"), rep(nrow(rbig_conv), 2)))
plot(ggplot(data=rbig_conv_ggplot, aes(x=stat, fill=type, color=type)) + geom_histogram(alpha=0.5, position="identity", bins=30))
steps <- seq(100, nrow(rbig_conv), 100)
p.value <- numeric(length(steps)) 
pb <- txtProgressBar(min = 0, max = length(steps), style = 3)
for(i in seq_along(steps)){
  dat <- head(rbig_conv, steps[i])
  #   p.value[i] <- cramer.test(dat[, 1], dat[, 2], sim="eigenvalue")$p.value
  p.value[i] <- ks.test(dat[, 1], dat[, 2])$p.value
  setTxtProgressBar(pb, i)
}
close(pb)

steps <- seq(1, 100000)
p.value <- numeric(length(steps)) 
pb <- txtProgressBar(min = 0, max = length(steps), style = 3)
for(i in seq_along(steps)){
  dat <- rbig_conv[sample.int(nrow(rbig_conv), size=150), ]
  #   p.value[i] <- cramer.test(dat[, 1], dat[, 2], sim="eigenvalue")$p.value
  p.value[i] <- ks.test(dat[, 1], dat[, 2])$p.value
  setTxtProgressBar(pb, i)
}
close(pb)
  
n <- 700
x1 <- rnorm(n)
x2 <- rnorm(n)
y1 <- 2 * x1 + 5*rnorm(n)
y2 <- -2 * y1 + 1*rnorm(n)
dat <- as.data.frame(cbind(x1, x2, y1, y2))
plot(hexplom(dat))
rbig_conv <- RBIG_hist(head(dat, 700), c(1,2), dist=dist, B=10000)

steps <- seq(100, nrow(rbig_conv), 100)
p.value <- numeric(length(steps)) 
pb <- txtProgressBar(min = 0, max = length(steps), style = 3)
for(i in seq_along(steps)){
  dat <- head(rbig_conv, steps[i])
  #   p.value[i] <- cramer.test(dat[, 1], dat[, 2], sim="eigenvalue")$p.value
  p.value[i] <- ks.test(dat[, 1], dat[, 2])$p.value
  setTxtProgressBar(pb, i)
}
close(pb)
plot(p.value)

n <- 700
x1 <- rnorm(n)
y1 <- 2 * x1 + 1*rnorm(n)
y2 <- -2 * y1 + 1*rnorm(n)
dat <- as.data.frame(cbind(x1, y1, y2))
dat <- scale(dat)
plot(hexplom(dat))

kcipt_conv <- KCIPT_test(head(dat, 100), c(1,3), 2, dist, 20000)
save(kcipt_conv, file="kcipt_conv.Rdata")
kcipt_conv_12_3 <- KCIPT_test(head(dat, 700), c(1,2), 3, dist, 100)
save(kcipt_conv_12_3, file="kcipt_conv_12_3.Rdata")

steps <- seq(100, nrow(kcipt_conv), 100)
p.value <- numeric(length(steps)) 
pb <- txtProgressBar(min = 0, max = length(steps), style = 3)
for(i in seq_along(steps)){
  datt <- head(kcipt_conv, steps[i])
  #   p.value[i] <- cramer.test(dat[, 1], dat[, 2], sim="eigenvalue")$p.value
  p.value[i] <- ks.test(datt[, 1], datt[, 2])$p.value
  setTxtProgressBar(pb, i)
}
close(pb)
plot(p.value)
kcipt_conv_ggplot <- data.frame("stat"=unlist(kcipt_conv), "type"=rep(c("H1", "H0"), rep(nrow(kcipt_conv), 2)))
plot(ggplot(data=kcipt_conv_ggplot, aes(x=stat, fill=type, color=type)) + geom_histogram(alpha=0.5, position="identity", bins=30))

steps <- 1:1000
p.value <- numeric(length(steps)) 
pb <- txtProgressBar(min = 0, max = length(steps), style = 3)
for(i in steps){
  datt <- kcipt_conv[sample.int(n=nrow(kcipt_conv), size=100), ]
  p.value[i] <- ks.test(datt[, 1], datt[, 2])$p.value
  p.value[i] <- t.test(datt[, 1], datt[, 2])$p.value
  setTxtProgressBar(pb, i)
}
close(pb)
hist(p.value)
