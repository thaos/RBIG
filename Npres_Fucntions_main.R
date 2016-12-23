source('Npres_Fucntions.R')
n <-100
d <- 2
sigma.Z <- 0.3
sigma.W <- 0.1
sigma.E <- 0.2
Z <- matrix(rnorm (n*d,0,sigma.Z),nr=n, nc=d)
W <- rnorm(n,0,sigma.W)
eps <- rnorm(n,0,sigma.E)
eps.prime <- rnorm(n,0,sigma.E)
X <- W + Z[,1] +eps
Y <- W + Z[,1] +eps.prime
Data <- cbind(X,Y,Z)
colnames(Data) <- c('X','Y', paste('Z.',1:d,sep= ''))
head(Data)

test.stat <- npresid.statistics(Data,d)
out <- npresid.boot(Data,d,R=50)
