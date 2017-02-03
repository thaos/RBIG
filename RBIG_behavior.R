library(ggplot2)
library(parallel)
# library(dHSIC)
source("MI_RBIG_2016_algo.R")
gmi_2D <- function(r)log(sqrt(1/(1-(r)^2)))
save_fast <- function(object){
  file=paste(deparse(substitute(object)), ".rds", sep="")
  saveRDS(object, file=file)
  file
}

# Calculate the number of cores
no_cores <- min(c(4, detectCores() - 1))
 
# Initiate cluster
cl <- makeCluster(no_cores)


cmi_btest_pc <- function(x, y, S, suffStat){
  cmi_btest(suffStat$B, suffStat$dat, x, y, S, cond_MI=suffStat$cond_MI)
}

n <- 200
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

gen_dat_c2 <- function(n, seed){
  set.seed(seed)
  x <- rnorm(n)
  y <- x + rnorm(n)
  z <- x + rnorm(n)
  w <- y + z + rnorm(n)
  dat <- as.data.frame(scale(cbind(x, y, z, w), scale=TRUE))
}

lmi <- unlist(mclapply(1:1000, 
	      function(s){
		dat <- gen_dat(n=n, seed=s)
		cond_MI_r(as.matrix(dat), 1, 2)
	      }))
save_fast(lmi)


lcmi_c1 <- unlist(mclapply(1:1000, 
	      function(s){
		dat <- gen_dat_c1(n=n, seed=s)
		cond_MI_r(as.matrix(dat), 1, 3, 2)
	      }))
save_fast(lcmi_c1)

lcmi_c2 <- unlist(mclapply(1:1000, 
	      function(s){
		dat <- gen_dat_c2(n=n, seed=s)
		cond_MI_r(as.matrix(dat), 1, 4, 2:3)
	      }))
save_fast(lcmi_c2)


lpv <- unlist(mclapply(1:1000, 
                function(s){
                  dat <- gen_dat(n=n, seed=s)
                  cmi_btest_pc(x=1, y=2, S=numeric(), 
                               suffStat=list(dat = as.matrix(dat),
                                             B=1000,
                                             cond_MI=cond_MI_r)) 
	      }))
save_fast(lpv)

lpv_c1 <- unlist(lapply(1:1000, 
                function(s){
                  dat <- gen_dat_c1(n=n, seed=s)
                  cmi_btest_pc(x=1, y=3, S=2, 
                               suffStat=list(dat = as.matrix(dat),
                                             B=1000,
                                             cond_MI=cond_MI_r)) 
	      }))
save_fast(lpv_cc1)

lpv_cc2 <- unlist(mclapply(1:1000, 
                function(s){
                  dat <- gen_dat_c2(n=n, seed=s)
                  cmi_btest_pc(x=1, y=4, S=2:3, 
                               suffStat=list(dat = as.matrix(dat),
                                             B=1000,
                                             cond_MI=cond_MI_r)) 
	      }))
save_fast(lpv_cc2)
