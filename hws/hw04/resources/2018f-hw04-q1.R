
## ISYE6420 ##
## Fall 2018 ##
rm(list=ls())

#Q1 Metropolis Hastings
## Set random seed
#set.seed(1234)
m=100000
burnin=m*.1
m=m+burnin

# Data
x <- -2 # used to calculate posterior likelihood, not used in code below

# Step 1: Initialize theta (theta_0)
theta = 0.5
# Vectors to hold theta values and whether proposed theta was accepted
thetas <-c(0)
acceptance <-vector(length=m)
theta_prop_vec <-vector(length=m)
r_vec <- vector(length=m)

alpha<-1
beta<- 3

#Step 4: Increase n and go back to step 2
for(n in 1:m){  
  ## Step 2: Generate proposal theta_prop from q(theta_prop|theta)
  theta_prop<-rgamma(n=1,shape=alpha,scale = 1/beta)
  theta_prop_vec[n]<-theta_prop
  ## Calculate expression r in rho
  ## Calculate rho
  #theta_prop<-rgamma(n=1,alpha, beta)
  #r<-((exp(-3*theta_prop))*dgamma(theta,shape=alpha,scale=beta))/
  #  (exp(-3*theta)*dgamma(theta_prop,shape=alpha, scale=beta))
  # r<-dnorm(x,0,1/theta_prop)*dgamma(theta_prop,1/2,1)*dgamma(theta,shape=alpha,scale=beta)/
  #   (dnorm(x,0,1/theta)*dgamma(theta,1/2,1)*dgamma(theta_prop,shape=alpha,scale=beta))
  
  r <- theta^(alpha-1)*exp(-beta*theta)*exp(-theta_prop*(1+0.5*x^2))/
    (theta_prop^(alpha-1)*exp(-beta*theta_prop)*exp(-theta*(1+0.5*x^2)))
  
  r_vec[n]<-r
  rho<-min(r,1)
  # Step 3: Generate U~U(0,1) and accept proposal theta_prop if U<rho
  U=runif(1)
  #print(paste("rho:", rho, ", U:", U))
  ## Added theta_prop !=0 because theta cannot be 0 per definition of parameter
  if(U<=rho){
    theta <- theta_prop
    acceptance[n] <- 1
    thetas<-c(thetas,theta)
  }
  if(U>rho){ }
}
par(mfrow=c(1,1))
hist(thetas)
mean(thetas)
var(thetas)
mean(acceptance[burnin:m])