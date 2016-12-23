# Multi-information estimation using RBIG
#
# The user can choose two orthogonal transforms:
#
#      'PCA' = PCA
#      'RND' = Random Rotations
#
# USE:
#
# [MI] = MI_RBIG_2016(dat,N_lay,transform,porc,precision)
#
# INPUTS:
# dat = data ( #dimensions , #samples );aim
# N_lay = number of layers (default N_lay = 1000);
# porc = extra domain percentage (default porc = 10)
# precision = number of points for the marginal PDFs estimation (default precision = 1000)
# transformation = linear transformation applied ('RND','PCA' default transformation = 'PCA')
#
# OUTPUTS
# MI = Multi-information
# MIs = Multi-information reduction at each layer.
# datT = Gaussianized data.
# 
# e.g.
#
# dat = rand(5)*(rand(5,1000).^2);
# dat <-  runif(5) * matrix(runif(5*2000), nrow=5)^2
# dat <-  matrix(runif(5*1000), nrow=5)
# dat <- t(dat)
# dat <- dat[, 1:2]
# n <- 1e4
# rho <- sqrt(runif(n))
# theta <- runif(n, 0, 2*pi)
# x <- rho * cos(theta)
# y <- rho * sin(theta)
# dat <- cbind(x,y)[rho>0.9,]
# pairs(dat)
# N_lay = 50;
# porc = 1;
# precision = 1000;
# transformation = 'PCA';
# MI = MI_RBIG_2016(dat,N_lay,transformation,porc,precision);
#
#
# Citation:
# Iterative Gaussianization: from ICA to Random Rotations. 
# V. Laparra, G. Camps & J. Malo 
# IEEE Transactions on Neural Networks, 22:4, 537 - 549, (2011)
#

library(entropy)
library(sROC)
library(mixAK)
library(MVN)
library(hexbin)
library(cramer)
library(lpSolve)
source("RBIG_r.r")

mis <- RBIG_r(dat)
mi100 <- sapply(1:100, function(x){
	 dat <-  matrix(runif(5*1000), nrow=5)
	 dat <- t(dat)
	 RBIG_r(dat)
}
)

cond_MI <- function(dat, x_ind, y_ind, c_ind=integer(0)){
  if(length(c_ind) == 0){
    ans <- RBIG_r(dat[, c(x_ind, y_ind)])
  }else{ 
	ans <- RBIG_r(dat[, c(x_ind, y_ind, c_ind)])
	ans <- ans - RBIG_r(dat[, c(x_ind, c_ind)])
	ans <- ans - RBIG_r(dat[, c(y_ind, c_ind)])
	if(length(c_ind) > 1)
	   ans <- ans + RBIG_r(dat[, c_ind])
  } 
  ans
}
dat <- matrix(runif(5*2000), ncol=5)
cmi <- cond_MI(dat, 1, 2)
dat[, 2]  <- dat[, 2] + dat[, 1] * 10
plot(hexplom(dat[, 1:2]))
cmi <- cond_MI(dat, 1, 2)

dat <- matrix(runif(5*2000), ncol=5)
dat[, 1]  <- dat[, 1] + dat[, 3] * 10
dat[, 2]  <- dat[, 2] + dat[, 3] * 10
plot(hexplom(dat[, 1:3]))
pairs(dat[, 1:3])
cmi <- cond_MI(dat, 1, 2)
cmi <- cond_MI(dat, 1, 2, 3)

dat <- matrix(runif(5*2000), ncol=5)
dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
plot(hexplom(dat[, 1:5]))
cmi1 <- cond_MI(dat, 1, 2)
cmi2 <- cond_MI(dat, 1, 2, 3)
cmi3 <- cond_MI(dat, 1, 2, c(3, 4))
cmi4 <- cond_MI(dat, 1, 2, c(4, 5))
cmi5 <- cond_MI(dat, 1, 2, c(3, 4, 5))


# boot_mi <- function(dat, x_ind, y_ind, c_ind=integer(0)){
#   if(length(c_ind) == 0){
sample_mi <- function(dat, x_ind, y_ind){
  dat <- dat[,c(x_ind, y_ind)]
  dat[, 1] <- sample(dat[,x_ind])
  dat[, 2] <- sample(dat[,y_ind])
  plot(hexplom(dat)) 
  dat
}
s1 <- sample_mi(dat, 1, 2)
plot(hexplom(s1))

sample_cmi <- function(dat, x_ind, y_ind, c_ind){
  dat <- dat[sample.int(nrow(dat), replace=TRUE), ]
  dat_c <- dat[, c_ind, drop=FALSE] 
  dat_xy <- dat[, c(x_ind, y_ind)]
  dist_mat <- as.matrix(dist(dat_c))
  two_closest <- apply(dist_mat, 2, function(x) order(x)[(1:2)+sum(x == 0)])
  new_x <- dat[c(two_closest[1, ], two_closest[2, ]), x_ind]  
  new_y <- dat[c(two_closest[2, ], two_closest[1, ]), y_ind]  
  dat_b <- cbind(new_x, new_y, rbind(dat_c, dat_c) )
  dat_b <- dat_b[sample.int(nrow(dat)/2), ]
  plot(hexplom(dat_b)) 
  dat_b 
}
s1 <- sample_cmi(dat, 1, 2, c(3, 4))
s1 <- sample_cmi(dat, 1, 2, c(4, 5))
s1 <- sample_cmi(dat, 1, 2, c(3, 4, 5))
s1 <- sample_cmi(dat, 1, 2, 3)
plot(hexplom(s1))

boot_mi <- function(dat, x_ind, y_ind){
  dat <- sample_mi(dat, x_ind, y_ind)  
  cond_MI(dat, 1, 2)
} 

boot_cmi <- function(dat, x_ind, y_ind, c_ind){
  dat <- sample_cmi(dat, x_ind, y_ind, c_ind)  
  cond_MI(dat, 1, 2, 3:ncol(dat))
} 
boot_cmi(dat, 1, 2, c(3, 4))
boot_cmi(dat, 1, 2, 3)

nboot_cmi <- function(n,dat, x_ind, y_ind, c_ind=numeric(0)){
  if(length(c_ind) == 0)
    ans <- unlist(lapply(seq.int(n), function(x) boot_mi(dat, x_ind, y_ind)))
  else
    ans <- unlist(lapply(seq.int(n), function(x){print(x); boot_cmi(dat, x_ind, y_ind, c_ind)}))
  ans
}
ncmi5 <- nboot_cmi(100, dat, 1, 2, c(3, 4, 5))
ncmi4 <- nboot_cmi(100, dat, 1, 2, c(4, 5))
ncmi3 <- nboot_cmi(100, dat, 1, 2, c(3, 4))
ncmi2 <- nboot_cmi(10, dat, 1, 2, 3)
ncmi1 <- nboot_cmi(10, dat, 1, 2)

cmi_btest <- function(nboot ,dat, x_ind, y_ind, c_ind=numeric(0)){
  cmi <- cond_MI(dat, x_ind, y_ind, c_ind)
  ncmi <- nboot_cmi(nboot, dat, x_ind, y_ind, c_ind)
  1 - sum(cmi > ncmi) / nboot
}
tcmi1 <-cmi_btest(100, dat, 1, 2)
tcmi3 <-cmi_btest(10, dat, 1, 2, 3:4)

conf_tcmi1  <- sapply(1:50, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			tcmi <-cmi_btest(50, dat, 1, 2)
			tcmi
})

conf_tcmi2  <- sapply(1:10, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			tcmi <-cmi_btest(10, dat, 1, 2, 3)
			tcmi
})
conf_tcmi3  <- sapply(1:10, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			tcmi <-cmi_btest(10, dat, 1, 2, 3:4)
			tcmi
})
conf_tcmi4  <- sapply(1:10, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			tcmi <-cmi_btest(10, dat, 1, 2, 4:5)
			tcmi
})
# code translated to R from Gary Doran et al. "A permutation-Based Kernel Conditional Independence Test
linear_permutation <- function(D){
  #   D <- as.matrix(dist(dat[1:3, 3:4]))
  D <- as.matrix(D)
  n <- nrow(D)
  # Rescale Distances
  D <- D / max(max(D))
  # Objective Function
  f <-  c(t(D))
  # Inequality contraint
  #   lb <- numeric(n^2)
  # Equality constraints
  Aeq <- matrix(0, nrow=2*n, ncol=n^2)
  b <- matrix(1, nrow=2*n, ncol=1)
  # Columns sum to 1
  for(c in 0:n-1){
    Aeq[c + 1, (c*n+1):((c+1)*n)] <- 1
  }
  # Rows sum to 1 (last row constraint not necessary
  # it is implied by other constraints)
  for(r in 1:(n-1)){
    for(c in 1:n){
      Aeq[r+n, r+(c-1)*n] <- 1
    }
  }
  # Diagonal entries zero
  for (z in 1:n){
    Aeq[2*n, (z-1)*(n+1) + 1] <- 1
  }
  b[2*n, 1] <- 0
  cdir <- paste(rep("=", 2*n))
  ans <- lp (direction = "min", objective.in=f, const.mat=Aeq, const.dir=cdir, const.rhs=b, transpose.constraints = TRUE, all.int=TRUE, all.bin=TRUE)
  ans <- matrix(ans$sol, ncol=n, byrow=FALSE) #%*% D
  ans
}

# check who to simulate from RBIG

KCIPT <- function(dat,xy_ind, c_ind,  dist, B, b, M){
  MMD <- numeric(B)
  samples <- numeric(B)
  inner_null <- matrix(numeric(B*b), nrow=B)
  outer_null <- numeric(M)
  for( i in 1:B){
    omega <- dat[, c(xy_ind, c_ind)]
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    P <- linear_permutation(dist(omega2[, 2:ncol(omega2)]))
    omega2 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    names(omega2) <- names(omega1)
    MMD[i] <- cramer.test(omega1, omega2, sim="ordinary", just.statistic=TRUE)$statistic
    #     print("***************************************")
    #     print(cramer.test(omega1, omega2, sim="ordinary"))
    #     print("***************************************")
    #     browser()
    #     plot(hexplom(omega1))
    #     scan(n=1)
    #     plot(hexplom(omega2))
    omega <- rbind(omega1, omega2)
    for( j in 1:b){
      idx <- sample.int(nrow(dat), round(nrow(dat)/2))
      omega1 <- omega[idx, ]
      omega2 <- omega[-idx, ]
      #       plot(hexplom(omega1))
      #       plot(hexplom(omega2))
      #       print(cramer.test(omega1, omega2, sim="ordinary"))
      #       browser()
      #       scan(n=1)
      inner_null[i, j] <- cramer.test(omega1, omega2, sim="ordinary", just.statistic=TRUE)$statistic
      #       cat(inner_null[i, j],  " / ", MMD[i], "\n")
    }
    #     print(sort(inner_null[i,]), round(0.05 * b)])
    #     print(sort(inner_null[i,])[round(0.95 * b)])
  }
  statistics <- mean(MMD)
  for(k in 1:M){
    for(i in 1:B){
      r <- ceiling(runif(1) * b)
      samples[i] <- inner_null[i, r]
    }
    outer_null[k] <- mean(samples)
  }
  #   print(statistics)
  #   print(outer_null)
  #   p.value <- mean(statistics >= outer_null)
  p.value <- 1 - rank(c(statistics, outer_null))[1]/(length(outer_null) + 1)
  #   crit.value <- sort(outer_null)[round(0.95 * length(outer_null))]
  p.value
}

KCIPT <- function(dat,xy_ind, c_ind,  dist, B, b, M){
  MMD <- numeric(B)
  samples <- numeric(B)
  inner_null <- matrix(numeric(B*b), nrow=B)
  outer_null <- numeric(M)
  for( i in 1:B){
    idx <- sample.int(nrow(dat), round(nrow(dat)/2))
    omega1 <- dat[idx, ]
    omega2 <- dat[-idx, ]
    P <- linear_permutation(dist(omega2[, c_ind]))
    #     MMD[i] <- cramer.test(omega1[, xy_ind], P%*%omega2[, xy_ind], sim="ordinary", just.statistic=TRUE)$statistic
    print(cramer.test(omega1[, c(xy_ind, c_ind)], cbind(P%*%omega2[, xy_ind[1]], omega2[, c(xy_ind[2], c_ind)]), sim="ordinary"))
    MMD[i] <- cramer.test(omega1[, c(xy_ind, c_ind)], cbind(P%*%omega2[, xy_ind[1]], omega2[, c(xy_ind[2], c_ind)]), sim="ordinary")$p.value
  }
  MMD
}
pv3 <- KCIPT(dat[1:500, ], c(1:2), c(3,4), dist=dist, B=10, b=10, M=100) 
pv4 <- KCIPT(dat[1:500, ], c(1:2), c(4,5), dist=dist, B=20, b=20, M=100) 
pv5 <- KCIPT(dat[1:500, ], c(1:2), c(3,4,5), dist=dist, B=20, b=20, M=100) 
pv1 <- KCIPT(dat[1:500, ], c(1:2), c(5), dist=dist, B=20, b=20, M=100) 


conf_kcipt3  <- sapply(1:20, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			#                         dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			#                         dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			dat[, 1]  <- dat[, 1] + dat[, 3] * 2 
			dat[, 2]  <- dat[, 2] + dat[, 3] * 2 
			plot(hexplom(dat))
			pv <- KCIPT(dat[1:700, ], c(1:2), c(4), dist=dist, B=10, b=50, M=100) 
			pv
})
