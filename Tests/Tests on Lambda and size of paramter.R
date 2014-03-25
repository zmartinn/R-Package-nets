library(nets)

# Parameters
N    <- 100
P    <- 10
P.NZ <- 0.2

# select Lambda to test
lambda.range <- rev( seq(0,100,5) )
L <- length(lambda.range)

#How many simulations per lambda
NbIter <- 100

#For every simulation we record the number of zero-parameter there is in theta.true and in theta.lambda
NbZeros = matrix(0,NbIter,length(lambda.range))
NbZerostrue = matrix(0,NbIter)


for( i in 1:NbIter )
{#Initialize theta, X and y
nonzero <- rbinom(P,1,P.NZ)
theta.true <- rep( 0 , P );
theta.true[ nonzero==1 ] <- rnorm( sum(nonzero) , 0 , 1)
NbZerostrue[i] <- P - length(which(theta.true == 0))
X <- matrix( rnorm(N*P,0,1) , N , P )
X[,1] <- 1
y <- X %*% theta.true + rnorm(N,0,1)

#Lasso algo for every Lambda
theta.lasso   <- matrix( 0 , P , L )
theta.sig2err <- rep( 0 , L )
for( j in 1:L )
{
	results <- alasso( y , X , lambda=lambda.range[j] , procedure="shooting" )

	theta.lasso[,j] <- results$theta
	theta.sig2err[j] <- results$sig2err
	NbZeros[i,j] <- P-length(which(theta.lasso[,j] == 0))
}


}
