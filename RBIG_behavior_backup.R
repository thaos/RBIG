library(ggplot2)
library(dHSIC)
source("RBIG_r.r")
source("MI_RBIG_2016_algo.R")
gmi_2D <- function(r)log(sqrt(1/(1-(r)^2)))

cmi_btest_pc <- function(x, y, S, suffStat){
  cmi_btest(suffStat$B, suffStat$dat, x, y, S, cond_MI=suffStat$cond_MI)
}

rbig_hist_pc <- function(x, y, S, suffStat){
  dat <- suffStat$dat
  B <- suffStat$B
  cond_MI <- suffStat$cond_MI
  to_plot <- !is.null(suffStat$to_plot)
  p.value <-RBIG_hist(dat,xy_ind=c(x,y), c_ind=S,  dist=dist, B=B, to_plot=to_plot, cond_MI=cond_MI)
    #     p.value <- cramer.test(p.value[, 1], p.value[, 2])$p.value
  p.value <- ks.test(p.value[, 1], p.value[, 2])$p.value
  p.value
}

kcipt_test_pc <- function(x, y, S, suffStat){
  dat <- suffStat$dat
  B <- suffStat$B
  to_plot <- !is.null(suffStat$to_plot)
  p.value <- KCIPT_test(dat,xy_ind=c(x,y), c_ind=S,  dist=dist, B=B, to_plot=to_plot)
  p.value <- ks.test(p.value[[1]], p.value[[2]])$p.value
  p.value
}

hsic_test_pc <- function(x, y, S, suffStat){
  dat <- suffStat$dat
  B <- suffStat$B
  b <- suffStat$b
  to_plot <- !is.null(suffStat$to_plot)
  xy_ind <- c(x, y)
  c_ind <- S
  stat_h0 <- numeric(B)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  stat <- dhsic(as.data.frame(dat), kernel=c("gaussian"))$dHSIC
  for( i in 1:B){
    omega <- dat
    omega[, 1] <- sample(dat[, 1]) 
    stat_h0[i] <- dhsic(as.data.frame(omega), kernel=c("gaussian"))$dHSIC
    setTxtProgressBar(pb, i)
  }
  close(pb)
  if(to_plot){
    df <- data.frame(stat=c(stat, c(stat_h0)), type=rep(c("H1","H0"), c(1,B)))
    plot(ggplot(data=df, aes(x=stat, fill=type, color=type)) + geom_histogram(aes(y=..density..),alpha=0.5, position="identity", bins=30)+ggtitle(paste("x=",xy_ind[1], "y=", xy_ind[2], " S=", c_ind))+theme(aspect.ratio=1/3))
  }
  p.value <- 1 - rank(c(stat, stat_h0))[1]/(length(stat_h0) + 1)
  p.value
}

n <- 100
#test 1
gen_dat <- function(n, seed){
  set.seed(seed)
  x <- rnorm(n)
  y <- rnorm(n)
  dat <- as.data.frame(scale(cbind(x, y), scale=TRUE))
}

gen_dat_c1 <- function(n, seed){
  set.seed(seed)
  x <- rnorm(n)
  y <- x + rnorm(n)
  z <- y + rnorm(n)
  dat <- as.data.frame(scale(cbind(x, y, z), scale=TRUE))
}

lmi_m <- sapply(1:100, 
	      function(s){
		dat <- gen_dat(n=n, seed=s)
	       	RBIG_r(dat)
	      })
hist(lmi)

lmi_r <- sapply(1:1000, 
	      function(s){
		dat <- gen_dat(n=n, seed=s)
		cond_MI_r(as.matrix(dat), 1, 2)
	      })
hist(lmi_r)

lcmi_r <- sapply(1:5000, 
	      function(s){
		dat <- gen_dat_c1(n=n, seed=s)
		cond_MI_r(as.matrix(dat), 1, 3, 2)
	      })
hist(lcmi_r)

lpv_m<- sapply(1:100, 
	      function(s){
		dat <- gen_dat(n=n, seed=s)
		rbig_hist_pc(x=1, y=2, S=numeric(), 
			     suffStat=list(dat = as.matrix(dat),
					   B=n,
					   to_plot=TRUE,
					   cond_MI=cond_MI_m)) 
	      })
hist(lpv_m)

lpv_r <- sapply(1:1, 
	      function(s){
		dat <- gen_dat(n=n, seed=s)
		rbig_hist_pc(x=1, y=2, S=numeric(), 
			     suffStat=list(dat = as.matrix(dat),
					   B=5000,
					   to_plot=TRUE,
					   cond_MI=cond_MI_r)) 
	      })
# hist(lpv_r)

lpv_k <- sapply(1:1, 
	      function(s){
		dat <- gen_dat(n=n, seed=s)
		kcipt_test_pc(x=1, y=2, S=numeric(), 
			     suffStat=list(dat = as.matrix(dat),
					   B=5000,
					   to_plot=TRUE,
					   cond_MI=cond_MI_r)) 
	      })
# hist(lpv_k)

lpv_h <- sapply(1:1000, 
	      function(s){
		dat <- gen_dat(n=n, seed=s)
		hsic_test_pc(x=1, y=2, S=numeric(), 
			     suffStat=list(dat = as.matrix(dat),
					   B=1000,
					   b=50,
					   to_plot=TRUE,
					   cond_MI=cond_MI_r)) 
	      })
print(lpv_h)
hist(lpv_h)

lpv_c <- sapply(1:100, 
                function(s){
                  dat <- gen_dat(n=n, seed=s)
                  cmi_btest_pc(x=1, y=2, S=numeric(), 
                               suffStat=list(dat = as.matrix(dat),
                                             B=1000,
                                             cond_MI=cond_MI_r)) 
	      })
print(lpv_c)
hist(lpv_c)

lpv_cc1 <- sapply(1:500,
                function(s){
                  dat <- gen_dat_c1(n=n, seed=s)
                  cmi_btest_pc(x=1, y=3, S=2, 
                               suffStat=list(dat = as.matrix(dat),
                                             B=500,
                                             cond_MI=cond_MI_r)) 
	      })
print(lpv_cc1)
hist(lpv_cc1)
