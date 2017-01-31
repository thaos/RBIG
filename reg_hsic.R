library(R.matlab)
library(igraph)
causal_test <- function(x, y, fit='gp'){
  Matlab$startServer()
  matlab <- Matlab()
  isOpen <- open(matlab)
  if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")
  evaluate(matlab, "addpath(genpath('~/Data/SIF/code/causal_test/'))")
  setVariable(matlab, x = x)
  setVariable(matlab, y = y)
  setVariable(matlab, fit = fit)
  evaluate(matlab, "ans = causal_fun(x, y, fit)")
  ans <- getVariable(matlab, "ans")
  close(matlab)
  ans <- ans$ans[, , 1]
  ans
}
  x <- 1:10+rnorm(10)
  y <- 1:10+rnorm(10)
ans <- causal_test(x, y, 'gp')

data  <-  read.table("~/Data/SIF/code/ano_aus.dat", header=TRUE)
data <- data[, c(1:3,10)] 

combi <- combn(names(data), 2)
all_test  <- apply(combi, 2, function(cb) causal_test(data[, cb[1]], data[, cb[2]], fit="gp"))

test1 <- all_test[[1]]

combi1 <- combi[, 1]
xname <- combi1[1]
yname <- combi1[2]

plot_ctest <- function(ctest, vnames){
  xname <- vnames[1]
  yname <- vnames[2]
  par(mfcol=c(3, 2))
  with(ctest, plot(x[order(x)], y[order(x)], xlab=xname, ylab=yname, main="forward model"))
  with(ctest, lines(x[order(x)], yfit[order(x)], col="blue"))
  xyrange <- with(ctest, range(c(yfit, y)))
  with(ctest, plot(y[order(x)], yfit[order(x)], xlab=yname, ylab="fitted", xlim=xyrange, ylim=xyrange, main="forward model : predicted vs true"))
  abline(b=1, a=0)
  with(ctest, plot(x[order(x)], yres[order(x)], xlab=xname, ylab="residuals", main="forward model : residuals", sub=paste("HSIC p-value = ", format(pf, digits=2))))
  abline(h=0)
  with(ctest, plot(y[order(y)], x[order(y)], xlab=yname, ylab=xname, main="backword model"))
  with(ctest, lines(y[order(y)], xfit[order(y)],col="red"))
  xyrange <- with(ctest, range(c(yfit, y)))
  with(ctest, plot(x[order(y)], xfit[order(y)], xlab=xname, ylab="fitted", xlim=xyrange, ylim=xyrange, main="backword model : predicted vs true"))
  abline(b=1, a=0)
  with(ctest, plot(y[order(y)], xres[order(y)], xlab=yname, ylab="residuals", main="backword model : residuals", sub=paste("HSIC p-value = ", format(pb, digit=2))))
  abline(h=0)
}
plot_ctest(test1, combi1)


# à faire apres mangé construire la matrice d'adjacence ???
find_idx <- function(adj_mat, ctest, vnames){
	lr  <- ctest$logratio
	dxy <- ctest$dxy
	if(dxy){
	  adj_mat[vnames[1], vnames[2]] <- sign(lr)
	  adj_mat[vnames[2], vnames[1]] <- - sign(lr)
	}
	adj_mat
}
fill_adjmat <- function(all_test, combi){
  variables <- unique(c(combi))
  adj_mat <- matrix(0, ncol=length(variables), nrow=length(variables))
  colnames(adj_mat) <- rownames(adj_mat) <- variables
  for(i in seq_along(all_test))
    adj_mat <- find_idx(adj_mat, all_test[[i]], combi[, i])
  adj_mat
}
adj_mat <- fill_adjmat(all_test, combi)
g1 <- graph_from_adjacency_matrix( adj_mat )
pdf(file="HSIC_causal_direction.pdf")
mapply(plot_ctest, all_test, as.data.frame(combi)) 
par(mfrow=c(1, 1))
plot(g1)
dev.off()
