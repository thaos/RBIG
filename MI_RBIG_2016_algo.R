library(entropy)
library(sROC)
library(mixAK)
library(MVN)
library(hexbin)
library(cramer)
library(lpSolve)
library(memoise)
library(amap)
library(scales)
library(sn)

entropy_mm <- function(x, nbins=sqrt(length(x))){  
    dx <- discretize(x, nbins)
    delta = diff(range(x))/ nbins
    hx = entropy.MillerMadow(dx, unit="log2")+log2(delta)
    hx
}

compute_tol_h0 <- function(nrow, ncol, probs=0.975,  n=1000){
  nbins <- sqrt(nrow)
  sim <- function(){
    x <- matrix(rnorm(nrow * ncol), ncol=ncol)
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

MI_RBIG_2016 <- function(dat, N_lay=1000){
  DIM = dim(dat)
  Nsamples = DIM[1]
  nbins <- floor(sqrt(Nsamples))
  DIM = DIM[2]
  delta_I <- numeric(N_lay)
  for (n in 1:N_lay){
    # marginal gaussianization
    p <- numeric(DIM)
    for(d in 1:DIM){
      margin  <-  marginal_gaussianization(dat[,d]);
      dat[, d] <- margin$x_gauss   
    }
    dat_aux = dat;
    # PCA rotation
    C  <- cov(dat)
    eig <- eigen(C);
    V <- eig$vectors
    #     V <- rRotationMatrix(1, ncol(C))
    dat <- dat %*% V  
    delta_I[n] = information_reduction_LT(dat,dat_aux, nbins=nbins);
    rt <- roystonTest(dat, qqplot = FALSE)
    hzt <- hzTest(dat, qqplot = FALSE)
    if (runif(1)< max(rt@p.value, hzt@p.value)) break
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

information_reduction_LT <- function(X, Y, nbins){
  #   should discretize first
  hx <- apply(X, 2, function(x)entropy.MillerMadow(discretize(x, nbins), unit="log2") + log2(diff(range(x))/nbins))
  hy <- apply(Y, 2, function(y)entropy.MillerMadow(discretize(y, nbins), unit="log2") + log2(diff(range(y))/nbins))
  I <- sum(hy - hx)
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

cond_MI_r <- function(dat, x_ind, y_ind, c_ind=integer(0)){
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

cond_MI_m <- function(dat, x_ind, y_ind, c_ind=integer(0)){
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

# boot_mi <- function(dat, x_ind, y_ind, c_ind=integer(0)){
#   if(length(c_ind) == 0){
sample_mi <- function(dat, x_ind, y_ind){
  dat <- dat[,c(x_ind, y_ind)]
  dat[, 1] <- sample(dat[,1])
  #   dat[, 2] <- sample(dat[,2])
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
sample_cmi <- function(dat, x_ind, y_ind, c_ind){
  cset <- dat[, 3:ncol(dat), drop=FALSE]
  sdd <- apply(cset, 2, function(x) sd(c(dist(x))))
  sdd <- sapply(sdd, function(x)rnorm(nrow(dat),sd=.01*sdd))
  cset <- cset + sdd 
  P <- linear_permutation(dist(cset))
  dat <-  cbind(P%*%dat[, 1], dat[, 2:ncol(dat)])
  dat
}

sample_cmi <- function(dat, x_ind, y_ind, c_ind, c_dist){
  sdd <- sd(c(c_dist))
  c_dist <- c_dist + rnorm(length(c_dist), sd=.1*sdd) 
  P <- linear_permutation(c_dist)
  dat <-  cbind(P%*%dat[, 1], dat[, 2:ncol(dat)])
  dat
}



boot_mi <- function(dat, x_ind, y_ind, cond_MI=cond_MI_r){
  dat <- sample_mi(dat, x_ind, y_ind)  
  cond_MI(dat, 1, 2)
} 

boot_cmi <- function(dat, x_ind, y_ind, c_ind, c_dist, cond_MI=cond_MI_r){
  dat <- sample_cmi(dat, x_ind, y_ind, c_ind, c_dist)  
  cond_MI(dat, 1, 2, 3:ncol(dat))
} 

nboot_cmi <- function(n,dat, x_ind, y_ind, c_ind=numeric(0), cond_MI=cond_MI_r){
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  if(length(c_ind) == 0)
    ans <- unlist(lapply(seq.int(n), function(i){setTxtProgressBar(pb, i); boot_mi(dat, x_ind, y_ind, cond_MI)}))
  else{
    c_dist <- dist(dat[, 3:ncol(dat), drop=FALSE])
    ans <- unlist(lapply(seq.int(n), function(i){setTxtProgressBar(pb, i); boot_cmi(dat, x_ind, y_ind, c_ind, c_dist, cond_MI)}))
  }
  close(pb)
  ans
}

cmi_btest <- function(nboot ,dat, x_ind, y_ind, c_ind=numeric(0), cond_MI=cond_MI_r){
  cmi <- cond_MI(dat, x_ind, y_ind, c_ind)
  #print(cmi)
  ncmi <- nboot_cmi(nboot, dat, x_ind, y_ind, c_ind, cond_MI)
  #   browser()
  df <- data.frame(stat=c(cmi, ncmi), type=rep(c("H1","H0"), c(1,nboot)))
    plot(ggplot(data=df, aes(x=stat, fill=type, color=type)) + geom_histogram(aes(y=..density..),alpha=0.5, position="identity", bins=30)+ggtitle(paste("x=",x_ind[1], "y=", y_ind[1], " S=", paste(c_ind, collapse=TRUE)))+theme(aspect.ratio=1/3))
  #p.value <- 1 - rank(c(cmi, ncmi))[1]/(length(ncmi) + 1)
  # p.value <- 1 - kCDF(ncmi, xgrid=cmi)
  p.value <- 1 - rank(c(cmi, ncmi))[1]/(length(ncmi) + 1)
  print(p.value)
#  p.value <- 1 - psn(cmi, dp=coef(selm(ncmi~1), "dp"))
#  print(p.value)
  p.value
}


cmi_btest_pc <- function(x, y, S, suffStat){
  cmi_btest(suffStat$nboot, suffStat$dat, x, y, S, cond_MI=suffStat$cond_MI)
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
    if(i < B/2){
      omega2 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    }else{
      omega2 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    } 
    #     omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    #     omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    #     omega2 <- rbind(omega21, omega22)[sample.int(nrow(omega1)), ]
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
    #     omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    #     omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    #     omega2 <- rbind(omega21, omega22)[sample.int(nrow(omega1)), ]
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

RBIG_hist <- function(dat, xy_ind, c_ind=numeric(0),  dist, B, to_plot=FALSE, cond_MI=cond_MI_r){
  stat <- numeric(B)
  stat_h0 <- numeric(B)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for( i in 1:B){
    omega <- dat
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    if(length(c_ind) != 0){
      P <- linear_permutation(dist(omega2[, 3:ncol(omega2)]))
      if(i < B/2){
	omega2 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
      }else{
	omega2 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
      } 
      stat[i] <- cond_MI(omega1, 1, 2, c_ind=3:ncol(omega1))
      stat_h0[i] <- cond_MI(omega2, 1, 2, c_ind=3:ncol(omega1))
    } else{
      omega2[, 1] <- omega2[sample.int(nrow(omega2)), 1] 
      stat[i] <- cond_MI(omega1, 1, 2)
      stat_h0[i] <- cond_MI(omega2, 1, 2)
    }
    #     omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    #     omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    #     omega2 <- rbind(omega21, omega22)[sample.int(2*nrow(omega2), length(idx)), ]
    #     cat("*")
    setTxtProgressBar(pb, i)
  }
  close(pb)
  #   cat("\n")
  #   hist(stat)
  #   par(new=TRUE)
  #   hist(stat_h0)
  if(to_plot){
    df <- data.frame(stat=c(stat, stat_h0), type=rep(c("H1","H0"), c(B,B)))
    plot(ggplot(data=df, aes(x=stat, fill=type, color=type)) + geom_histogram(alpha=0.5, position="identity", bins=30)+ggtitle(paste("x=",xy_ind[1], "y=", xy_ind[2], " S=", c_ind))+theme(aspect.ratio=1/3))
  }
  data.frame(stat, stat_h0)
}
# rbig_hist_12_34 <- RBIG_hist(head(dat, 700), 1:2, 3:4, dist, 200)
# rbig_hist_45_3 <- RBIG_hist(head(dat, 700), 4:5, 3, dist, 200)
# rbig_hist_12_3 <- RBIG_hist(head(dat, 700), 1:2, 3, dist, 200)

KCIPT_test <- function(dat, xy_ind, c_ind=numeric(0),  dist, B, to_plot=FALSE ){
  stat_h01 <- numeric(B)
  stat_h11 <- numeric(B)
  stat_h02 <- numeric(B)
  stat_h12 <- numeric(B)
  samples <- numeric(B)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for( i in 1:B){
    omega <- dat
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    cpermute <- function(omega){
      if(length(c_ind) != 0){
	P <- linear_permutation(dist(omega[, 3:ncol(omega)]))
	if(i < B/2){
	  omega <-  cbind(P%*%omega[, 1], omega[, 2:ncol(omega)])
	}else{
	  omega <-  cbind(omega[, 1], P%*%omega[, 2], omega[, 3:ncol(omega)])
	} 
      }else{
	omega <- omega 
	omega[, 1] <- omega[sample.int(nrow(omega)), 1] 
      }
      omega
    }
    omega2p <- cpermute(omega2)
    omega1p <- cpermute(omega1)
    stat_h01[i] <- cramer.test_simple(omega1p, omega2p)
    stat_h11[i] <- cramer.test_simple(omega1, omega2p)
    #     stat_h02[i] <- cramer.test_simple(omega1, omega2)
    stat_h12[i] <- cramer.test_simple(omega1p, omega2)
    #     cat("*")
    setTxtProgressBar(pb, i)
  }
  close(pb)
  stat_h1 <- c(stat_h11, stat_h12)
  stat_h0 <- c(stat_h01)
  #   stat_h0 <- c(stat_h01, stat_h02)
  if(to_plot){
    df <- data.frame(stat=c(stat_h1, stat_h0), type=rep(c("H1","H0"), c(2*B,B)))
    plot(ggplot(data=df, aes(x=stat, fill=type, color=type)) + geom_histogram(aes(y=0.1*..density..), alpha=0.5, position="identity", binwidth=0.1)+ggtitle(paste("x=",xy_ind[1], "y=", xy_ind[2], " S=", c_ind))+theme(aspect.ratio=1/3))
  }
  list(stat_h0, stat_h1)
}

RBIG_simple <- function(dat, xy_ind, c_ind=numeric(0),  dist, B, to_plot=FALSE){
  stat_h0 <- numeric(B)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  if(length(c_ind) != 0){
    stat <- cond_MI(dat, 1, 2, c_ind=3:ncol(dat))
    P <- linear_permutation(dist(dat[, 3:ncol(dat)]))
    omega <-  cbind(P%*%dat[, 1], dat[, 2:ncol(dat)])
    omega <-  rbind(omega, cbind(dat[, 1], P%*%dat[, 2], dat[, 3:ncol(dat)]))
  }else{
    stat <- cond_MI(dat, 1, 2)
    omegab <- dat
  }
  for( i in 1:B){
    if(length(c_ind) == 0){
      omegab[, 1] <- dat[sample.int(nrow(dat)), 1] 
      stat_h0[i] <- cond_MI(omegab, 1, 2)
    }else{
      omegab <- sample.int(n=nrow(omega), size=nrow(dat), replace=TRUE)  
      omegab <- omega[omegab, ]
      stat_h0[i] <- cond_MI(omegab, 1, 2, c_ind=3:ncol(omegab))
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  if(to_plot){
    df <- data.frame(stat=stat_h0)
    plot(ggplot(data=df, aes(x=stat)) + geom_histogram(alpha=0.5, position="identity", bins=30) + geom_vline(xintercept = stat) + ggtitle(paste("x=",xy_ind[1], "y=", xy_ind[2], " S=", c_ind))+theme(aspect.ratio=1/3))
  }
  p.value <- 1 - rank(c(stat, stat_h0))[1]/(length(stat_h0) + 1)
}
