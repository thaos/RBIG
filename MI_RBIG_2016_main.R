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

rm(list=ls())
library(entropy)
library(sROC)
library(mixAK)
library(MVN)
library(hexbin)
library(cramer)
library(lpSolve)
library(memoise)

mis <- MI_RBIG_2016(dat)
mi100 <- sapply(1:100, function(x){
	 dat <-  matrix(runif(5*1000), nrow=5)
	 dat <- t(dat)
	 MI_RBIG_2016(dat)$MI
}
)

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

s1 <- sample_mi(dat, 1, 2)
plot(hexplom(s1))

s1 <- sample_cmi(dat, 1, 2, c(3, 4))
s1 <- sample_cmi(dat, 1, 2, c(4, 5))
s1 <- sample_cmi(dat, 1, 2, c(3, 4, 5))
s1 <- sample_cmi(dat, 1, 2, 3)
plot(hexplom(s1))

boot_cmi(dat, 1, 2, c(3, 4))
boot_cmi(dat, 1, 2, 3)

ncmi5 <- nboot_cmi(100, dat, 1, 2, c(3, 4, 5))
ncmi4 <- nboot_cmi(100, dat, 1, 2, c(4, 5))
ncmi3 <- nboot_cmi(100, dat, 1, 2, c(3, 4))
ncmi2 <- nboot_cmi(10, dat, 1, 2, 3)
ncmi1 <- nboot_cmi(10, dat, 1, 2)

tcmi1 <-cmi_btest(10, dat, 1, 2)
tcmi3 <-cmi_btest(50, dat, 1, 2, 3:4)

conf_tcmi1  <- sapply(1:20, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			tcmi <-cmi_btest(50, dat, 1, 2)
			tcmi
})

conf_tcmi2  <- sapply(1:20, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			tcmi <-cmi_btest(50, dat, 1, 2, 3)
			tcmi
})
conf_tcmi3  <- sapply(1:20, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			tcmi <-cmi_btest(50, dat, 1, 2, 3:4)
			tcmi
})
conf_tcmi4  <- sapply(1:20, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			tcmi <-cmi_btest(50, dat, 1, 2, 4:5)
			tcmi
})
pv3 <- KCIPT(dat[1:500, ], c(1:2), c(3,4), dist=dist, B=10, b=10, M=100) 
pv4 <- KCIPT(dat[1:500, ], c(1:2), c(4,5), dist=dist, B=20, b=20, M=100) 
pv5 <- KCIPT(dat[1:500, ], c(1:2), c(3,4,5), dist=dist, B=20, b=20, M=100) 
pv1 <- KCIPT(dat[1:500, ], c(1:2), c(5), dist=dist, B=20, b=20, M=100) 


conf_kcipt5  <- sapply(1:100, function(x){
			print("***********************************************")
			print(x)
			dat <- matrix(runif(5*2000), ncol=5)
			dat[, 1]  <- dat[, 1] + dat[, 3] * 10 - dat[, 4] * 10
			dat[, 2]  <- dat[, 2] + dat[, 3] * 10 - dat[, 4] * 5
			#                         dat[, 1]  <- dat[, 1] + dat[, 3] * 2 
			#                         dat[, 2]  <- dat[, 2] + dat[, 3] * 2 
			plot(hexplom(dat))
			pv <- KCIPT(dat[1:700, ], c(1:2), c(5), dist=dist, B=20, b=50, M=100) 
			pv
})
