
require(nets)

library(MASS)

# Problem Parameters
N   <- 6
T   <- 1000
 
# SigInv
SigInv <- matrix( 0 , N , N )
SigInv[1,1] <-  4.00; SigInv[1,2] <-  0.00; SigInv[1,3] <-  0.00; SigInv[1,4] <-  0.00; SigInv[1,5] <- 0.00; SigInv[1,6] <- 0.00;
SigInv[2,1] <-  0.00; SigInv[2,2] <-  2.00; SigInv[2,3] <-  0.00; SigInv[2,4] <-  0.00; SigInv[2,5] <- 0.00; SigInv[2,6] <- 0.00;
SigInv[3,1] <-  0.00; SigInv[3,2] <-  0.00; SigInv[3,3] <-  2.00; SigInv[3,4] <- -0.40; SigInv[3,5] <- 0.00; SigInv[3,6] <- 0.00;
SigInv[4,1] <-  0.00; SigInv[4,2] <-  0.00; SigInv[4,3] <- -0.40; SigInv[4,4] <-  3.00; SigInv[4,5] <- 0.00; SigInv[4,6] <- 0.00;
SigInv[5,1] <-  0.00; SigInv[5,2] <-  0.00; SigInv[5,3] <-  0.00; SigInv[5,4] <-  0.00; SigInv[5,5] <- 2.00; SigInv[5,6] <- 0.00;
SigInv[6,1] <-  0.00; SigInv[6,2] <-  0.00; SigInv[6,3] <-  0.00; SigInv[6,4] <-  0.00; SigInv[6,5] <- 0.00; SigInv[6,6] <- 1.00;

y <- mvrnorm(T, rep(0,N) , solve(SigInv) )

#network <- nets( y, type='pc' , lambda=80 )

print( cbind( SigInv , rep(NA,N) , round(network$C,2) ) )
cat('\n')

# info
network

plot( network )

# Some Metrics
mse <- sqrt( sum( (SigInv - network$C)^2 ) )
tp1 <- mean( (network$C)[ SigInv==0 ] != 0 ) # THIS IS WRONG: I SHOULD NOT CONSIDER THE DIAGONAL
pow <- mean( (network$C)[ SigInv!=0 ] != 0 )
pos <- min(eigen( network$C )$values)> 1e-6 

cat( 'MSE' , mse , 'TYPE1' , tp1 , 'pow' , pow , 'pos:' , pos , 'lambda' , network$lambda , '\n')

