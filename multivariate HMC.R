## Multivariate HMCMC


## Let's try tri-variate Normal with Unit Covariance

# MVN function 
library(mvtnorm)

mvden <- function(x){  
  return(dmvnorm(x, mean = c(0,0,0), sigma= diag(3)))  # input : vector 
}


grad <- function(x,func){ # input: vector 
  obj <- c()
  h <- 10^(-7)
  for(i in 1:length(x)){
    xh <- x 
    xh[i] <- x[i] +h
    dx <- h
    df <- func(xh) - func(x)
    obj[i] <- df/dx
  }
  return(obj)
}

grad(c(1,2,3) , dmvnorm)

# comparison with package 
library(pracma)
install.packages('pracma')
pracma::grad(dmvnorm, c(1,2,3))

#-----------------------------------------------------------------------------
## leapfrog function
mLF <- function(x, p, dKdp, M, dVdx, step_size, step_iter){   # x & p are vectors 
  for(i in 1:step_iter){
    p <- p + step_size/2 * dVdx(x) #half step for momentum
    x <- x - step_size * p %*% solve(M)  #whole step for x
    p <- p + step_size/2 * dVdx(x) # another half step for momentum
  } 
  p <- -p # moment flip 
  obj <- list(x=x,p=p) # momentum flip
  return(obj)
}

# test
step_size = 0.3
step_iter = 5
x= c(1,2,3)
p =  rmvnorm(1,mean=rep(0,3),sigma=diag(3))
V <- function(x) -log(dmvnorm(x))  
dVdx <- function(x) grad(x, V)
M = diag(3)


##----------------------------------------------------------------------------------------

mHamMC <- function(n_iter, initial_x, likelihood, step_size, step_iter, M){
  
  # set of sampled x, momentum p 
  samples <- list()
  samples[[1]] <- matrix(rep(0,3),ncol=3) # target variable
  samples[[2]] <- matrix(rep(0,3),ncol=3) # moment variable
  
  # Hamiltonian function H = V+K
  V <- function(x) -log(likelihood(x))
  #K <- function(p) 1/2* t(p)%*%solve(M)%*%p   
  
  #Derivatives functions
  dVdx <- function(x) grad(x, V)
  

  #leapfrog integration
  for(i in 1:n_iter){
    initial_p <- rmvnorm(1,mean=rep(0,3),sigma=M) #initial momentum   
    result <- mLF(initial_x,initial_p, dKdp, M, dVdx, step_size, step_iter) 
    x_new <- result[[1]]
    p_new <- result[[2]]
    
    old_l <- likelihood(initial_x)*mvden(initial_p)
    new_l <- likelihood(x_new)*mvden(p_new)
    alpha <- min(1, new_l/old_l)
    if(alpha > runif(1)){
      new <- c(x_new, p_new)
      samples[[1]] <- rbind(samples[[1]], x_new)
      samples[[2]] <- rbind(samples[[2]], p_new)
    }else{
      new <- c(x,p)
      samples[[1]] <- rbind(samples[[1]], initial_x)
      samples[[2]] <- rbind(samples[[2]], initial_p)
    }
    initial_x <- x_new
    initial_p <- p_new 
  }
  return(samples)
}


mv_samples <- mHamMC(10000, c(0,0,0), dmvnorm, 0.3, 5, diag(3))

#plot marginal
par(mfrow=c(2,3))
hist(mv_samples[[1]][,1], nclass=40,prob=T) 
lines(x=seq(-5,5,len=100), y=dnorm(seq(-5,5,len=100)), col='blue', lwd=4)
hist(mv_samples[[1]][,2], nclass=40,prob=T) 
lines(x=seq(-5,5,len=100), y=dnorm(seq(-5,5,len=100)), col='blue', lwd=4)
hist(mv_samples[[1]][,3], nclass=40,prob=T) 
lines(x=seq(-5,5,len=100), y=dnorm(seq(-5,5,len=100)), col='blue', lwd=4)

#plot correlation
plot(mv_samples[[1]][,1],mv_samples[[1]][,2], main='X1 vs X2')
plot(mv_samples[[1]][,1],mv_samples[[1]][,3], main='X1 vs X3')
plot(mv_samples[[1]][,2],mv_samples[[1]][,3], main='X2 vs X3')




###_----------------------------------------------------------------------------------------------------------
## correlated case

M = diag(3)
M[1,2]=M[2,1]=0.5
M[1,3]=M[3,1]=0.25
x = c(0,0,0)

mv_samples2 <- mHamMC(10000, c(0,0,0), dmvnorm, 0.3, 5, M)



#plot marginal
par(mfrow=c(2,3))
hist(mv_samples2[[1]][,1], nclass=40,prob=T) 
lines(x=seq(-5,5,len=100), y=dnorm(seq(-5,5,len=100)), col='blue', lwd=4)
hist(mv_samples2[[1]][,2], nclass=40,prob=T) 
lines(x=seq(-5,5,len=100), y=dnorm(seq(-5,5,len=100)), col='blue', lwd=4)
hist(mv_samples2[[1]][,3], nclass=40,prob=T) 
lines(x=seq(-5,5,len=100), y=dnorm(seq(-5,5,len=100)), col='blue', lwd=4)

#plot correlation
plot(mv_samples2[[1]][,1],mv_samples2[[1]][,2], main='X1 vs X2')
plot(mv_samples2[[1]][,1],mv_samples2[[1]][,3], main='X1 vs X3')
plot(mv_samples2[[1]][,2],mv_samples2[[1]][,3], main='X2 vs X3')












