library(entropy)
library(sROC)
library(mixAK)
library(MVN)
library(hexbin)
library(cramer)
library(lpSolve)
library(memoise)
library(amap)

entropy_mm <- function(x, nbins=sqrt(length(x))){  
    dx <- discretize(x, nbins)
    delta = diff(range(x))/ nbins
    hx = entropy.MillerMadow(dx, unit="log2")+log2(delta)
    hx
}

compute_tol_h0 <- function(nrow, ncol, probs=0.975,  n=1000){
  nbins <- sqrt(nrow)
  sim <- function(){
    x <- matrix(rnorm(nrow * ncol), ncol=ncol(dat))
    mnegent <- apply(x, 2, entropy_mm)
    sum(mnegent)
  }
  tol <- sapply(seq.int(n), function(i){ sim() - sim()})
  #   tol <- quantile(abs(tol), probs=probs)
  tol
}
compute_tol_h0_m <- memoise(compute_tol_h0) 

MI_RBIG_2016 <- function(dat, tol, N_lay=1000){
  DIM = dim(dat)
  Nsamples = DIM[1]
  nbins <- floor(sqrt(Nsamples))
  DIM = DIM[2]
  delta_I <- numeric(N_lay)
  tol <- compute_tol_h0_m(nrow(dat), ncol(dat))
  for (n in 1:N_lay){
    # marginal gaussianization
    p <- numeric(DIM)
    for(d in 1:DIM){
      margin  <-  marginal_gaussianization(dat[,d]);
      #       p[d] <- margin$shapiro.test$p.value
      #       while(p[d] < 0.9){
      #         margin  <-  marginal_gaussianization(margin$x_gauss);
      #       print(p[d])
      #         p[d] <- margin$shapiro.test$p.value
      #       }
      dat[, d] <- margin$x_gauss   
      #       hist(dat[, d])
      #       browser()
      #       pairs(dat)
      #       plot(hexplom(dat))
      #       scan(n=1)
    }
    dat_aux = dat;
    # PCA rotation
    C  <- cov(dat)
    eig <- eigen(C);
    V <- eig$vectors
    #     V <- rRotationMatrix(1, ncol(C))
    #     print(V)
    dat <- dat %*% V  
    # multi-information reduction
    #     delta_I[n] = information_reduction_LT(dat,dat_aux, tol_d=tol_d, tol_m=tol_m, nbins=nbins);
    delta_I[n] = information_reduction_LT(dat,dat_aux, tol=tol, nbins=nbins);
    #     print (n)
    #     print(delta_I[n])
    #     pairs(dat)
    #     plot(hexplom(dat))
    if(n>10){
      #       browser()
      #       mt <- mardiaTest(dat, qqplot = FALSE)
      #       cat(rt@p.value, " / ", hzt@p.value, " / ", mt@p.value, "\n")
      #       print(cor(dat))
      if (isTRUE(all.equal(tail(delta_I[(n-9):n], 9),  rep(0, 9)))) break
    }
    rt <- roystonTest(dat, qqplot = FALSE)
    hzt <- hzTest(dat, qqplot = FALSE)
    #     cat(rt@p.value, " / ", hzt@p.value, " \n")
    #     if (rt@p.value >= 0.1 & hzt@p.value > 0.1) break
    if (rt@p.value >= 0.9 & hzt@p.value > 0.9) break
  }
  ans <- list(dat=dat, MIs=delta_I, MI=sum(delta_I))
}

# information_reduction_LT <- function(X, Y, tol_d, tol_m, nbins){
information_reduction_LT <- function(X, Y, tol, nbins){
  #   should discretize first
  hx <- apply(X, 2, function(x)entropy.MillerMadow(discretize(x, nbins), unit="log2") + log2(diff(range(x))/nbins))
  hy <- apply(Y, 2, function(y)entropy.MillerMadow(discretize(y, nbins), unit="log2") + log2(diff(range(y))/nbins))
  #   hx <- apply(X, 2, knn_entropy_1D)
  #   hy <- apply(Y, 2, knn_entropy_1D)
  # wrong use
  #   dix <- sum(apply(X, 2, FNN::entropy))
  #   diy <- sum(apply(Y, 2, FNN::entropy))
  #   browser()
  #   print(dix)
  #   print(diy)
  I <- sum(hy - hx)
  #   print(I)
  #   scan(n=1)
  #   I <- dix - log(sqrt(2*pi*exp(1))) 
  #   II = sqrt(sum((hy - hx)^2));
  #   p = 0.25;
  #   print(abs(II))
  #   print(sqrt(ncol(X)*((p*tol_d^2))))
  #   scan(n=1)
  #   if (abs(II)<sqrt(ncol(X)*((p*tol_d^2)))){
  if (abs(I) <= quantile(abs(tol), probs=0.975)){
  #   print("***tol***")
  #   print(mean(I >= abs(tol))) 
    I <- (runif(1) <= mean(abs(I) >= abs(tol))) * I
  #   I= (runif(1) >= 0.975) * I
    #     I=0
  #     print("inside")
  }
  I
} 

marginal_gaussianization <- function(x){
  #   x_order <- order(x)
  #   x_cdfk <- kCDF(x, xgrid=x)
  #   x_unif <- x_cdfk$Fhat
  x_unif <- ecdf(x)(x)
  #   x_gauss <- qnorm(x_unif)[x_order]
  x_gauss <- qnorm(x_unif)
  x_gauss[x_unif==1]  <- 1 - 1/length(x)^2
  #   ans <- list(x_gauss=x_gauss, shapiro.test=shapiro.test(x_gauss))
  ans <- list(x_gauss=x_gauss)
}

knn_entropy_1D <- function(x){
 N <- length(x)
 x_order <- sort(x) 
 x_diff <- diff(x_order)
 mean(log(x_diff)) +  digamma(1) - digamma(N)
 #  mean(log(N*x_diff)) -  digamma(1) + log(1)
}

cond_MI <- function(dat, x_ind, y_ind, c_ind=integer(0)){
  if(length(c_ind) == 0){
    ans <- MI_RBIG_2016(dat[, c(x_ind, y_ind)])$MI
  }else{ 
	ans <- MI_RBIG_2016(dat[, c(x_ind, y_ind, c_ind)])$MI
	ans <- ans - MI_RBIG_2016(dat[, c(x_ind, c_ind)])$MI
	ans <- ans - MI_RBIG_2016(dat[, c(y_ind, c_ind)])$MI
	if(length(c_ind) > 1)
	   ans <- ans + MI_RBIG_2016(dat[, c_ind])$MI
  } 
  ans
}

# cond_MI <- function(dat, x_ind, y_ind, c_ind=integer(0)){
#   if(length(c_ind) == 0){
#     ans <- RBIG_r(dat[, c(x_ind, y_ind)])
#   }else{ 
#     ans <- RBIG_r(dat[, c(x_ind, y_ind, c_ind)])
#     ans <- ans - RBIG_r(dat[, c(x_ind, c_ind)])
#     ans <- ans - RBIG_r(dat[, c(y_ind, c_ind)])
#     if(length(c_ind) > 1)
#       ans <- ans + RBIG_r(dat[, c_ind])
#   } 
#   ans
# }

# boot_mi <- function(dat, x_ind, y_ind, c_ind=integer(0)){
#   if(length(c_ind) == 0){
sample_mi <- function(dat, x_ind, y_ind){
  dat <- dat[,c(x_ind, y_ind)]
  dat[, 1] <- sample(dat[,1])
  dat[, 2] <- sample(dat[,2])
  #   plot(hexplom(dat))
  dat
}

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
  #   plot(hexplom(dat_b))
  dat_b 
}

boot_mi <- function(dat, x_ind, y_ind){
  dat <- sample_mi(dat, x_ind, y_ind)  
  cond_MI(dat, 1, 2)
} 

boot_diff_mi <- function(dat, x_ind, y_ind){
  dat <- dat[, c(x_ind, y_ind)]
  dat <- dat[sample.int(nrow(dat), replace=TRUE), ]
  mi <- cond_MI(dat, 1, 2) 
  dat <- sample_mi(dat, 1, 2)  
  mi0 <- cond_MI(dat, 1, 2) 
  mi - mi0 
}

boot_diff_cmi <- function(dat, x_ind, y_ind, c_ind){
  dat <- dat[, c(x_ind, y_ind, c_ind)]
  dat <- dat[sample.int(nrow(dat), replace=TRUE), ]
  mi <- cond_MI(dat, 1, 2, 3:ncol(dat)) 
  dat <- sample_cmi(dat, 1, 2, 3:ncol(dat))  
  mi0 <- cond_MI(dat, 1, 2, 3:ncol(dat)) 
  mi - mi0 
}

boot_cmi <- function(dat, x_ind, y_ind, c_ind){
  dat <- sample_cmi(dat, x_ind, y_ind, c_ind)  
  cond_MI(dat, 1, 2, 3:ncol(dat))
} 

nboot_cmi <- function(n,dat, x_ind, y_ind, c_ind=numeric(0)){
  if(length(c_ind) == 0)
    ans <- unlist(lapply(seq.int(n), function(x){cat("*"); boot_mi(dat, x_ind, y_ind)}))
  else
    ans <- unlist(lapply(seq.int(n), function(x){cat("*"); boot_cmi(dat, x_ind, y_ind, c_ind)}))
  cat("\n")
  ans
}

nboot_diff_cmi <- function(n,dat, x_ind, y_ind, c_ind=numeric(0)){
  dat <- as.matrix(dat)
  if(length(c_ind) == 0)
    ans <- unlist(lapply(seq.int(n), function(x){cat(x, "*"); boot_diff_mi(dat, x_ind, y_ind)}))
  else
    ans <- unlist(lapply(seq.int(n), function(x){cat(x, "*"); boot_diff_cmi(dat, x_ind, y_ind, c_ind)}))
  cat("\n")
  ans
}
#bsample <-  nboot_diff_cmi(50, dat, 1, 2, 3)
#bsample <-  nboot_diff_cmi(50, dat, 1, 3, 2)

cmi_btest <- function(nboot ,dat, x_ind, y_ind, c_ind=numeric(0)){
  cmi <- cond_MI(dat, x_ind, y_ind, c_ind)
  ncmi <- nboot_cmi(nboot, dat, x_ind, y_ind, c_ind)
  #   browser()
  1 - sum(cmi > ncmi) / nboot
}

diff_cmi_btest <- function(nboot ,dat, x_ind, y_ind, c_ind=numeric(0)){
  n_diffcmi <- nboot_cmi(nboot, dat, x_ind, y_ind, c_ind)
  #   browser()
  t.test(n_diffcmi)$p.value
}

cmi_btest_pc <- function(x, y, S, suffStat){
  cmi_btest(suffStat$nboot, suffStat$dat, x, y, S)
}

# diff_cmi_btes <- function(nboot, dat, x_ind, y_ind, c_ind)

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

KCIPT <- function(dat,xy_ind, c_ind=numeric(0),  dist, B, b, M){
  MMD <- numeric(B)
  samples <- numeric(B)
  inner_null <- matrix(numeric(B*b), nrow=B)
  outer_null <- numeric(M)
  for( i in 1:B){
    omega <- as.matrix(dat[, c(xy_ind, c_ind)])
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    P <- linear_permutation(dist(omega2[, 3:ncol(omega2)]))
    omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    omega2 <- rbind(omega21, omega22)[sample.int(nrow(omega1)), ]
    MMD[i] <- cramer.test(omega1, omega2, sim="ordinary", just.statistic=TRUE)$statistic
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
    hist(inner_null[i, ])
    abline(v=MMD[i])
    cat("*")
  }
  cat("\n")
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


KCIPT_simple <- function(dat,xy_ind, c_ind,  dist, B, b, M){
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


KCIPT <- function(dat, xy_ind, c_ind=numeric(0),  dist, B, b, M){
  MMD <- numeric(B)
  samples <- numeric(B)
  inner_null <- matrix(numeric(B*b), nrow=B)
  outer_null <- numeric(M)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  for( i in 1:B){
    omega <- dat
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    P <- linear_permutation(dist(omega2[, 3:ncol(omega2)]))
    omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    omega2 <- rbind(omega21, omega22)[sample.int(nrow(omega1)), ]
    MMD[i] <- cramer.test_simple(omega1, omega2)
    omega <- rbind(omega1, omega2)
    for( j in 1:b){
      idx <- sample.int(nrow(dat), round(nrow(dat)/2))
      omega1 <- omega[idx, ]
      omega2 <- omega[-idx, ]
      inner_null[i, j] <- cramer.test_simple(omega1, omega2)
    }
    cat("*")
  }
  cat("\n")
  statistics <- median(MMD)
  for(k in 1:M){
    for(i in 1:B){
      r <- ceiling(runif(1) * b)
      samples[i] <- inner_null[i, r]
    }
    outer_null[k] <- median(samples)
  }
  #   print(statistics)
  #   print(outer_null)
  #   p.value <- mean(statistics >= outer_null)
  p.value <- 1 - rank(c(statistics, outer_null))[1]/(length(outer_null) + 1)
  #   crit.value <- sort(outer_null)[round(0.95 * length(outer_null))]
  p.value
}

# from the cramer packages

cramer.test_simple <- function(x, y, kernel="phiCramer"){
  .cramer.statistic<-function(daten,indexe,mm,nn,lookup) {
    xind<-indexe[1:mm]
    yind<-indexe[(mm+1):(mm+nn)]
    mm*nn/(mm+nn)*(2*sum(lookup[xind,yind])/(mm*nn)-sum(lookup[xind,xind])/(mm^2)-sum(lookup[yind,yind])/(nn^2))
  }
  m<-nrow(x)
  n<-nrow(y)
  daten<-matrix(c(t(x),t(y)),ncol=ncol(x),byrow=TRUE)
  lookup<-eval(call(kernel, as.matrix(Dist(daten))))
  .cramer.statistic(daten,1:(m+n),m,n,lookup)
}

KCIPT_outeronly <- function(dat, xy_ind, c_ind=numeric(0),  dist, B ){
  MMD <- numeric(B)
  samples <- numeric(B)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  for( i in 1:B){
    omega <- dat
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    P <- linear_permutation(dist(omega2[, 3:ncol(omega2)]))
    omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    omega2 <- rbind(omega21, omega22)[sample.int(nrow(omega1)), ]
    MMD[i] <- cramer.test_simple(omega1, omega2)
    #     cat("*")
    print(median(MMD[1:i]))
  }
  cat("\n")
  MMD
}

KCIPT_inneronly <- function(dat, xy_ind, c_ind=numeric(0),  dist, B, b, M){
  MMD <- numeric(B)
  samples <- numeric(B)
  inner_null1 <- matrix(numeric(B*b), nrow=B)
  inner_null2 <- matrix(numeric(B*b), nrow=B)
  outer_null <- numeric(M)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  for( i in 1:B){
    omega <- dat
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    #     P <- linear_permutation(dist(omega2[, 3:ncol(omega2)]))
    #     omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    #     omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    #     omega2 <- rbind(omega21, omega22)[sample.int(nrow(omega1)), ]
    MMD[i] <- cramer.test_simple(omega1, omega2)
    omega <- rbind(omega1, omega2)
    for( j in 1:b){
      idx <- sample.int(nrow(dat), round(nrow(dat)/2))
      omega1 <- omega[idx, ]
      omega2 <- omega[-idx, ]
      inner_null1[i, j] <- cramer.test_simple(omega1, omega2)
      idx <- sample.int(nrow(dat), round(nrow(dat)/2))
      omega1 <- omega[idx, ]
      omega2 <- omega[-idx, ]
      inner_null2[i, j] <- cramer.test_simple(omega1, omega2)
      if(j>10){ 
	#print(cramer.test(t(inner_null1[i, 1:j]), t(inner_null2[i, 1:j])))
	print(summary(c(t(inner_null1[1:i, 1:j])-t(inner_null2[1:i, 1:j]))))
      }
    }
    cat("*")
  }
  cat("\n")
  list(inner_null1, inner_null2)
}

RBIG_kcipt <- function(dat, xy_ind, c_ind=numeric(0),  dist, B, b, M){
  MMD <- numeric(B)
  samples <- numeric(B)
  inner_null <- matrix(numeric(B*b), nrow=B)
  outer_null <- numeric(M)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  for( i in 1:B){
    omega <- dat
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    P <- linear_permutation(dist(omega2[, 3:ncol(omega2)]))
    omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    omega2 <- rbind(omega21, omega22)[sample.int(nrow(omega1)), ]
    MMD[i] <- cond_MI(omega1, 1, 2, c_ind=3:ncol(omega1))
    omega <- rbind(omega1, omega2)
    for( j in 1:b){
      idx <- sample.int(nrow(dat), round(nrow(dat)/2))
      omega2 <- omega[-idx, ]
      inner_null[i, j] <- cond_MI(omega2, 1, 2, c_ind=3:ncol(omega2))
    }
    cat("*")
  }
  cat("\n")
  statistics <- median(MMD)
  for(k in 1:M){
    for(i in 1:B){
      r <- ceiling(runif(1) * b)
      samples[i] <- inner_null[i, r]
    }
    outer_null[k] <- median(samples)
  }
  p.value <- 1 - rank(c(statistics, outer_null))[1]/(length(outer_null) + 1)
  p.value
}
# RBIG_kcipt(head(dat, 700), 1:2, 3, dist, 10, 20, 100)

RBIG_hist <- function(dat, xy_ind, c_ind=numeric(0),  dist, B){
  stat <- numeric(B)
  stat_h0 <- numeric(B)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  for( i in 1:B){
    omega <- dat
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    P <- linear_permutation(dist(omega2[, 3:ncol(omega2)]))
    omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    omega2 <- rbind(omega21, omega22)[sample.int(nrow(omega1)), ]
    stat[i] <- cond_MI(omega1, 1, 2, c_ind=3:ncol(omega1))
    stat_h0[i] <- cond_MI(omega2, 1, 2, c_ind=3:ncol(omega1))
    cat("*")
  }
  cat("\n")
  hist(stat)
  par(new=TRUE)
  hist(stat_h0)
  data.frame(stat, stat_h0)
}
rbig_hist_12_34 <- RBIG_hist(head(dat, 700), 1:2, 3:4, dist, 200)
rbig_hist_45_3 <- RBIG_hist(head(dat, 700), 4:5, 3, dist, 200)
rbig_hist_12_3 <- RBIG_hist(head(dat, 700), 1:2, 3, dist, 200)
