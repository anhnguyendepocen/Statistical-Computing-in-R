rm(list = ls())
.rs.restartR()

#########################################################
#############   BTRY/STSCI 4520  ########################
#############      Final Exam    ########################
############# Due: May 20, 2018  ########################
#########################################################


# Instructions: save this file in the format <NetID>_Final.R. 
# Complete each question using code below the question number.
# You need only upload this file to CMS. 

# Note, we assume your working directory contains any files
# that accompany this one. 

# Further note: 10% will be deducted if your file produces
# an error when run. If your code produces an error and you 
# cannot find it, comment out that portion of the code and we
# will give partial credit for it. 

# Do not use either the function set.seed() or rm(list=ls())
# in your code. 


#### IMPORTANT INSTRUCTIONS FOR THE FINAL

## The final is to be completed as though it were in-class. 
## That means 
##
## 1. You must complete this work independently, without 
## collaboration or external assistance. Violations will
## be treated under the academic code. 
##
## 2. We will not provide office hours. We will monitor 
## Piazza and provide clarifications where questions are
## unclear. 
##
## 3. While we will not, in general, debug your work for you,
## we are happy to try to explain error messages if you can
## isolate the line of code that produces them and, preferably,
## reproduce them with a few-line script.
##
## [Example: R tries to interpret dimensions, but doesn't 
## always get it right.  So the code
##  
##   b = matrix(2,1,1)
##   t(1:3)/b
##
##  returns an error, but (1:3)/b does not.  We think it reasonable
##  to explain something like this to you if you get a message
##
##   Error in t(1:3)/b : non-conformable arrays
##
##  along with other instances of "I know where the error is, I just
##  don't know why R has a problem or how to fix it"]


################################################
# Question 1: Control Variates and Functionals #
################################################

# This question covers some recent ideas in Monte Carlo
# integration first published just two years ago in 
#
# https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12185

# Recall that if we are interested in estimating
#
#  mu = E g(X)   where  X ~ f(x) 
#
# a Monte Carlo estimate is given by 
#
#  gbar =  (1/n) sum  g(Xi)   where  Xi ~ f(x)
#
# Control Variates use a function h(x) where we
# know that 
#
#  E h(X) = 0  if  X ~ f(x)
#
# (If we know Eh(X) = m, then use h*(x) = h(x)-m in place of h(x).)
#
# Then if h(x) looks like g(x) we can expect that  
#
#  hbar =  1/n sum h(Xi)    
#
# (which "should" be zero) might tell us something about
# how far  gbar is from mu.

# That is, if hbar is too high, it is likely that gbar is also
# too high.  In particular, we can improve gbar by
#
#  gbar2 =  gbar - alpha hbar
#
# Note that since hbar has mean zero, the expectation of gbar2 is
# still mu.  We found in class that the optimal alpha was 
#
#  alpha = cov(g(X),h(X))/var(h(X))
#
# which we can calculate from the g(Xi) and h(Xi) values. 

# a) Write a function to calculate Vanilla and Control Variate
# estimates of gbar from a sample X. It should take arguments
# X (a sample of data)  and functions g and h, and should
# return a 2-vector giving the Vanilla and Control Variate estimates.

#### Estimate function and control variates used in may questions. 
####
# g = function(x) { sin(x+0.1) }
# h = function(x) { -exp(x) + exp(1/2) }

MCfunc1 = function(X,g,h){
  G = g(X)
  H = h(X)
  n = length(X)
  gbar = sum(G)/n
  
  alpha = cov(G,H)/var(H)
  hbar = sum(H)/n
  gbar2 = gbar - alpha*hbar
  
  return (c(gbar, gbar2))  
}
### Test code from the check file
set.seed(11011807)
X = rnorm(100)
eq = function(x){x}
all.equal(MCfunc1(X,sin,eq),c(0.0041542073, -0.0003166141),tol=1e-7)


# b) i) Use this to write a function to conduct a simulation
# study in which you set g to be  sin(x+0.1) and h to be 
# -exp(x) + exp(1/2)  and generate X from a standard normal
# (the expectation of exp(X) is exp(1/2 when X is N(0,1)). 
# You should have length(X) be N = 100 and conduct R = 50 
# replications and return the percentage variance reduction 
# due to control variates

MCsim1 = function(N,R){
  g = function(x) {sin(x+0.1)}
  h = function(x) {-exp(x) + exp(1/2)}
  
  X = rnorm(N)
  ret = matrix(NA, 2, R)
  
  for(i in 1:R) {
    sam = sample(X, size = N, replace=TRUE)
    ret[, i] = MCfunc1(sam, g, h)  
  }
  
  var1 = var(ret[1,])
  var2 = var(ret[2,])
  reduction = (var1 - var2)/var1
  
  return (reduction*100)
  
}


# ii) What variance reduction do you achieve?
set.seed(12011807)
reduction = MCsim1(100, 50)
cat('Reduction:', reduction, '\n')
### 69.32792


# c) We can combine both control variates and antithetic
# sampling by looking at 
#
#  (1/n) sum (g(Xi) + g(-Xi))/2 + alpha (h(Xi) + h(-Xi))/2
#
# (note that if f(x) is symmetric about zero, X and -X
# have the same distribution). 

# i) Write a function to conduct a simulation as above 
# but with antithetic sampling, too.  You should return
# the Vanilla Monte Carlo variance and the variance
# reduction due to each of antithetic sampling, control
# variates and both. 
#
# It makes a big difference that you do the antithetic
# average before the control variate calculation.
#
# In this case use the same number N of random samples
# for antithetic sampling rather than fixing the number
# of function evaluations. 

MCsim2 = function(N, R) {
  g = function(x) {sin(x+0.1)}
  h = function(x) {-exp(x) + exp(1/2)}
  
  X = rnorm(N)
  ret = matrix(NA, 2, R)
  reta = matrix(NA, 2, R)
  
  for(i in 1:R) {
    sam = sample(X, size = N, replace=TRUE)
    ret[, i] = MCfunc1(sam, g, h)
    
    gg = (g(sam) + g(-sam))/2
    hh = (h(sam) + h(-sam))/2
    gbar = sum(gg)/N
    hbar = sum(hh)/N
    
    reta[1,i] = gbar
    
    alpha = cov(gg, hh)/var(hh)
    
    reta[2,i] = gbar - alpha*hbar 
  }
  
  valnilla.var = var(ret[1,])
  control.var  = var(ret[2,])
  
  anti.valnilla.var = var(reta[1,])
  anti.control.var  = var(reta[2,])
  
  anti.reduction = (valnilla.var - anti.valnilla.var)/valnilla.var
  cont.reduction = (valnilla.var - control.var)/valnilla.var
  anti.cont.reduction = (valnilla.var - anti.control.var)/valnilla.var
  
  return( list( Vanilla.Var = valnilla.var,
                Anti.Reduction = anti.reduction * 100,
                Control.Reduction = cont.reduction * 100,
                Both.Reduction =  anti.cont.reduction * 100) )
}


# ii) Which adds most reduction?
set.seed(12011807)
ret = MCsim2(100, 50)
cat('Anti.Reduction:', ret$Anti.Reduction, '\n')
cat('Control.Reduction:', ret$Control.Reduction, '\n')
cat('Both.Reduction:', ret$Both.Reduction, '\n')
#(99.64933, 69.32792, 99.96023)

X = sort(seq(-4,4, 0.1))
data = matrix(0, nrow=length(X), ncol=2)
par(mar=c(1,1,1,1))

data[,1] = sin(X+0.1)
data[,2] = (-exp(X) + exp(1/2))/50
colnames(data) = c('g(x)', 'h(x)')
matplot(X, data[,1:2], type='l', ylim = c(-1, 1), xlim=c(-4,4), ylab=c('g(x), h(x)'))


# d) The better the agreement between h and g, the more
# effective the control variates strategy will be.  Eg, 
# compare the h we used in parts b and c with using
# h(x) = x. 

# Oates et. al. 2016 propose to approximate g non-parametrically.
#
# Specifically, we'll let
#
#  h(x) = sum_j d_j psi_j(x)
#
# To do this, we need to know the integral of h(x). To accomplish
# this, the authors propose using a modified form of a kernel. 
# Specifically, they use the first X_1,...,X_(N/2) samples and 
# define  
#
#   psi_j(x) =  phi'(x-X_j;s) + phi(x-X_j;s) f'(x)/f(x)
#
# where phi(x;s) is a normal density with standard deviation s and ' 
# represents a derivative with respect to x. 
#
# The reason for this choice is that 
#
#  int psi_j(x)f(x)dx = int [phi'(x-X_j;s) f(x) + phi(x-X_j;s) f'(x)]dx
#                     = int [phi(x-X_j;s) f(x)]' dx 
#                     = phi(Inf,s)f(Inf) - phi(-Inf,s)f(-Inf)
#                     = 0
#
# With this choice of psi, we get the d_j by simple least-squares regression
# of  g(Xi) on  psi_j(Xi) for  i = 1,...,N/2. 
#
# To carry this out you should 
#    - evaluate psi_j(X_i) for for the X_i in your data
#    - estimate d_j by using the matrix of phi_j(X_i) to predict g(X_i)
#      using lm.  Some coefficients will be returned as NA -- set them to 0. 
#    - then use these d_j to give you sum_j d_j phi_j(x) for any x. 

# i) Write a function to evaluate this h(x) for a vector of values x. It 
# should have inputs
#
# x -- the vector of values at which to evaluate h(x) 
# X -- a vector of draws from f(x)
# s -- the bandwidth
# g -- the function we want to approximate
#
# In this case you can assume that f(x) is a standard normal. It should
# return a vector of values of length(x).  You may use 

hfunc = function(x,X,s,g){
  s=s^2
  
  ### start to find d_j based on X  
  N = length(X)
  k.mat = matrix(0, nrow=N, ncol=N)
  
  for(j in (1:N)) {
    for(i in (1:N)) {
      z = X[i] - X[j]
      ph = exp(-z^2/(2*s))/sqrt(2*pi*s)
      phd = -ph*z/s
      k.mat[i,j] = phd + ph*(-X[i])
      
    }
  }
  
  ### "k.mat-1" means no intercept
  ret = lm(g(X)~k.mat-1)
  
  ### Check if any coefficient is NA and if yes, set it to 0
  dj = ret$coef
  dj[is.na(dj)] = 0
  
  ### Now let's compute function value at x
  hv = rep(0,length(x))
  for(i in (1:length(x))) {
    for(j in (1:length(X))){
      z = x[i] - X[j]
      ph = exp(-z^2/(2*s))/sqrt(2*pi*s)
      kj = -ph*z/s + ph*(-x[i])
      hv[i] = hv[i] + dj[j]*kj
    }
  }
  return (hv)
}

### Test code from the check file
set.seed(11071832)
testfn = function(x){ cos(x) }
X = rnorm(7)
all.equal( as.vector(hfunc(c(0,0.2),X,0.3,testfn)), c(1.018482,1.165351), tol=1e-6)
all.equal( as.vector(hfunc(c(2,3),X,0.5,testfn)), c( 2.5681179,0.3270968), tol=1e-6)

psi1vec = c(-5.413411e-01, -8.341119e-01, -9.735045e-01, -7.261490e-01,
            1.391458e-16,  1.000000e+00, 1.846233e+00,  2.178447e+00,  
            1.947009e+00,  1.390187e+00, 8.120117e-01)

xvec1 = seq(0,2,by=0.2)
testfn1 = function(x){ 1 }
all.equal(as.vector(hfunc(xvec1,1,0.5,testfn1)),psi1vec,tol=1e-7)


# ii) Using g(x) = sin(x), and s = 0.2, set 
xvec = seq(-2,2,by = 0.1)

# iii) plot h(x) based on 100 standard normals. 

set.seed(13071832)

plot.new()
par(mfrow=c(1,2))

hx = hfunc(xvec,rnorm(100),0.2, sin)
plot(xvec, hx, type='l', col='blue')
lines(xvec, sin(xvec), col='red')
legend('topleft',c('h(x)','g(x)'),lwd=2,col=c('blue','red'),lty=1,cex=.75)
title(main="N=100", col.main="black", font.main=3)

hx = hfunc(xvec,rnorm(5),0.2, sin)
plot(xvec, hx, type='l', col='blue')
lines(xvec, sin(xvec), col='red')
legend('topleft',c('h(x)','g(x)'),lwd=2,col=c('blue','red'),lty=1,cex=.75)
title(main="N=5", col.main="black", font.main=3)



# e) i) Write a function to implement
#
#  1.  Vanilla Monte Carlo
#  2.  Control Variates, using h(x) = x
#  3.  Control Functionals
#        - Use the first half of the values X_1,...,X_N
#          to produce hfunc as in 1d
#        - Use the second half to obtain the control variate
#          estimate with this hfunc
# 
# These should be based on X being N(0,1)
#
# Your function should take arguments N, g and s and return
# the three estimates above as a 3-vector

ControlFunctionals = function(N,g,s){
  h = function(x) {x}
  X = rnorm(N)
  ret = MCfunc1(X,g,h)
  
  half.X = head(X, N/2)
  rest.X = tail(X, N - N/2)
  
  h1 = function(x){
    hfunc(x,half.X,s,g)
  }

  ret2 = MCfunc1(rest.X,g,h1)
  
  return (c(ret[1], ret[2], ret2[2]))  
}

# ii) When g(x) = sin(x) we know that E g(X) = 0. Based on 
# 50 replicates at N = 100, using s = 0.2, what is the 
# expected squared error of each of your estimates?

set.seed(14071832)
result = matrix(0, nrow=50, ncol=3)
for (i in (1:50)) {
  result[i,] = ControlFunctionals(100, sin, 0.2)
}
### Calculate expected squared error. Since Eg(x)=0
result = result - 0
result = result^2
MSE = apply(result, 2, mean)
print(MSE)
### 0.0039411599 0.0005587374 0.0002870273


# f) One of the claims of the paper is that we can 
# achieve a faster convergence rate than (1/sqrt(N))
# because as we get more data, h(x) gets closer to g(x). 
# 
# i) Write a function to conduct a simulation study so that
# for each value of 

Ns = c(100,200,500,1000)

# you estimate the root-average-squared-error (RMSE) of
# each of the three estimates in 1e based on R = 20 
# simulations using values of s in 

svec = 2/sqrt(Ns)

# that correspond to each value in Ns.  (Here we have
# smaller s so that we can approximate g better when
# we have more data points). 

# You should return the matrix of mean squared errors.

# CAUTION: this may take a few minutes to run. 

ControlFuncSim = function(Ns,R,s){
  g = function(x){sin(x)}
  est = matrix(0, nrow=R, ncol=3*length(Ns))
  for (i in (1:R)) {
    for(j in (1:length(Ns))) {
      col = (j-1)*3+1
      est[i, col:(col+2)] = ControlFunctionals(Ns[j], g, s[j])    
    }
  }
  
  ### Calculate expected squared errors for each estimator using each s. Remember E(g(x))=0
  est = est - 0
  est = est^2
  mse = apply(est, 2, mean)
}

set.seed(14071832)
MSE = ControlFuncSim(Ns, 20, svec)

# ii) Plot log(RMSE) versus log(N) for the three estimates.
# This should give you approximately a straight line with
# the order of convergence given by the slope.

rmse.log.mat = log(matrix(sqrt(MSE), nrow=length(Ns), ncol=3, byrow=TRUE))
ns.log = log(Ns)
plot.new()
par(mfrow=c(1,1))
matplot(ns.log, rmse.log.mat, type='l', ylab=c('log(RMSE)'))
legend("topright", c('vanilla', 'control variates', 'control functionals'),col=seq_len(3),cex=0.8, lty=1)

# iii) Confirm that vanilla and control variates approaches have
# MSE = O(1/Ns); what is the order for Control Functionals?

### Graphs are comparable to the ones on lecture slide. 
### All three estimators have RMSE = O(1/log(Ns).




####################################################
# Question 2: Kernel Density Estimation and        # 
#             Minimum Hellinger Distance estimates #
####################################################

# In this question we will use Kernel Density Estimation
# as an intermediate step in statistical procedures. 
#
# But we will first need to deal with KDE's. Recall that
# we can estimate a nonparametric density by 
#
#  fhat(x,s) =  (1/n) sum phi(x-Xi;s)
#
# (using the same notation as in Question 1). 

# a) i) Write a function to evaluate fhat for a vector
# of evaluation points x, using data X with bandwidth s. 
# You should not need to use either for loops or apply 
# statements. 

#### Helper functions
### Function to thinner samples got from  MCMC random walk
### "ndrop" first ndrop samples are dropped
### "nskip" skip nskip-th from samples 
# thinner = function(ndrop, nskip, samples) {
#   total = seq(ndrop, length(samples))
#   ind.skip = seq(1, length(samples), nskip)
#   thinobs = total[-ind.skip]
#   thinobs = thinobs[ndrop:length(thinobs)]
#   samples[thinobs]
# }

thinner = function(ndrop, npick, samples) {
  thinobs = seq(ndrop, length(samples), npick)
  samples[thinobs]
}


### kde.matrix is a helper function used by function kde and kde.cv
### x is a vector where we need to estimate
### MX is a matrix. The i-th row is the data we calculate based on for i-th x in the vector x 
kde.matrix = function(x, MX, s) {
  N = ncol(MX)
  x = exp(-(x-MX)^2/(2*s^2))
  rowSums(x/(N*sqrt(2*pi)*s))
}

### I don't use loop and apply
kde = function(x,X,s) {
  N = length(X)
  XX = matrix(X, nrow=length(x), ncol=N, byrow=TRUE)
  kde.matrix(x, XX, s)
}

set.seed(12071844)
X = rnorm(10)
all.equal(as.vector(kde(c(-0.2,1),X,0.5)),c(0.3772147,0.2583243),tol=1e-6)


# The data we will use for this question comes from a study of
# the effectiveness of a drug pyrantel used to treat parasites
# in horses.  These data record the logit of the ratio of 
# the number of parasite eggs found in horse feces before
# and after treatment.

eggrate = c(-1.1653, -0.7538, -1.3218, -2.3394, -1.9766, -1.8718, -1.5041)

# ii) Using s = 0.2 produce a plot of the kernel density estimate
# on values of x from -2.5 to -0.5.
x = seq(-2.5, -0.5, by=0.10)
s = 0.2
kdes= kde(x, eggrate, s)
plot(x, kdes, col='red', type='o', xlab="x", ylab="kde")

# b) In order to apply this well, we need to choose a 
# bandwidth. To do this we'll use cross-validation. 
# Specifically we'll define a score for each h as
#
# CV(s) =  sum log  fhat^(-i)(Xi;s)
#
# that is. For each Xi, we leave Xi out of the data
# and estimate fhat^(-i) without Xi, then we see how high
# fhat^(-i) thinks the density is at Xi.  The higher
# the density, the better we are at predicting where future
# data might be.  Hnce the optimal s is the one that maximizes
# CV(s).
#
# i) Write a function to evaluate CV(s) for each s. Bonus 2 points
# for avoiding for loops and apply statements. 

### I don't use loop and apply
kde.cv = function(Y, s){
  N = length(Y)
  v = rep(Y, N)
  
  #index where we need to pick each Xi
  ind = seq(1, N^2, N+1)
  
  #Take each Xi in Y 
  x = v[ind]
  
  ### Remove Xi from Y and convert into a matrix
  ### The i-th row is the data after remove the i-th Xi in Y 
  MX = matrix(v[-ind], nrow=N, ncol=N-1, byrow=TRUE)
  
  est = kde.matrix(x, MX, s)
  
  sum(log(est))
    
}
### Test code from the check file
set.seed(12071844)
X = rnorm(10)
all.equal(kde.cv(X,0.34),-12.44386,tol=1e-6)


# ii) Use this function to decide on the optimal s for the eggrate
# data above with possible bandwidths given in 

ss = seq(0.1,1,by=0.1)

# to give the value sbest that maximizes the cross-validated
# likelihood. 
cvs = sapply(ss, kde.cv, Y=eggrate)
sbest = ss[which.max(cvs)]
cat('sbest=', sbest, '\n')
### 0.5



# c) We will also need to be able to simulate data from fhat. 
# To do this, we observe that we can express fhat as being exactly
# the density that corresponds to the following
#
#   1. Choose one of the Xi at random
#   2. Simulate a new z from a N(Xi,s) distribution. 
#
# i) Write a function to simulate N observations from fhat given 
# X and s. 2 bonus points for avoiding for loops 

kde.sim = function(N,X,s) {
  Xi = sample(X, N, replace=TRUE)
  rnorm(N,mean=Xi, sd=s)
}

### Test code from the check file
set.seed(12071844)
X = rnorm(10)
all.equal(kde.sim(2,X,0.45),c(0.04231365,-0.70325134))

# ii) Use this to simulate 1000 data points from your estimate
# above (with optimal bandwidth). Draw a histogram of these points
# and add the original data as vertical red lines. 
x.sim = kde.sim(1000, eggrate, sbest)
d.orig = kde(eggrate, eggrate, sbest)
hist(x.sim, breaks=100)
abline(v=eggrate, col='red')

# d) One of the advantages of having a density estimate is
# that there are more ways that you can compare your data to 
# a parametric family of models. 
#
# In particular, one measure of how different two distributions
# are is Hellinger distance
#
#    HD(f,g) = int ( sqrt(f(x)) - sqrt(g(x)) )^2 dx
#            = 2 - 2 int sqrt(f(x)*g(x)) dx
#
# For a family of densities f(x,theta) and an estimate 
# fhat(x), the minimum Hellinger Distance estimate (MHDE)
# is obtained by maximizing
#
#    A(theta) =  \int sqrt( f(x,theta)*fhat(x) ) dx
#
# In our case, we will use a normal family for the mean
#
#  f(x,theta) = dnorm(x,mean=theta)
#
# and our KDE. 
#
# We still need to approximate the integral. Here we will 
# sample from fhat and use the approximation
#
#  A(theta,fhat) = int sqrt( f(x,theta)/fhat(x) ) fhat(x) dx
#           ~=  (1/N) sum  sqrt( f(Xi,theta)/fhat(Xi) )
#
# where the Xi are sampled from fhat as in part c. 

# i) Write a function to calculate A(theta), for f(x,theta) given
# by dnorm(x,mean=theta) using the original data X and a sample
# Xsamp generated from fhat with bandwidth s. 

### Xsamp are data points used to performance the integration in A(theta, fhat)
### X are data used to define kde
Afn = function(theta, X, Xsamp, s) {
  afn_value = 0
  for(X_i in Xsamp) {
    afn_value = afn_value + sqrt(dnorm(X_i, mean=theta)/kde(X_i, X, s))
  }
  return (afn_value/length(Xsamp))
}

### Test code from the check file
set.seed(12071844)
X = rnorm(10)
Xsamp = rnorm(10)
all.equal(Afn(0.2,X,Xsamp,0.52),1.106754,tol=1e-6)

# ii) And a further function to find the optimum value based on 
# Data X, bandwidth s and sample size N. 
#
# You may use any optimization function you like (Giles
# used 'optimize') and can assume the maximum lies in 
# (-10,10). 

### Maximize Afn over theta. optimize optimize a function over it's first argument. 
### The first argument of Afn is theta 
HD.opt = function(X,s,N) {
  Xsamp = kde.sim(N,X,s)
  #Choose theta to maxmize A(theta). Since it is a one dimentsional optimization. Use optimize function
  interval = c(-10, 10)
  ret = optimize(Afn, interval, X=X, Xsamp=Xsamp, s=s, maximum = TRUE)
  
  return (ret)
}

# iii) Use this to obtain the value that minimizes Hellinger
# distance for the data using the optimal bandwidth
# above based on 1000 Monte Carlo samples. 

hellinger.opt = HD.opt(eggrate, sbest, 1000)
print(hellinger.opt)
cat('Afn at theta=', hellinger.opt$maximum, ' gives the max. value=', hellinger.opt$objective)
###Afn at theta= -1.532006  gives the max. value= 0.962548


# e) One reason to look for the MHDE is robustness
# to outliers. To see this add an additional data point
# to your data set with values in 

O = seq(-1,10,by=0.2)

# plot the value of the MHDE and the average 
# value in your data versus the values in O above. Keep
# the optimal bandwidth you calculated in part c. 
mm = matrix(0, nrow=length(O), ncol=2)
for(i in (1:length(O))) {
  X = c(eggrate, O[i])
  mm[i,1] = mean(X)
  A.opt = HD.opt(X, sbest, 1000)
  mm[i, 2] = A.opt$maximum
}

plot.new()
matplot(O, mm, type='o', ylab='MHDE, Mean')
legend("topleft", inset=0.01, 
       legend=c('Mean', 'MHDE'), 
       col=c(1:2),
       pch=15:20, 
       bg= ("white"), 
       horiz=F)

# f) One of the reasons that Hellinger distance has
# received some attention is that it is supposed to be
# as precise as the mean, at least asymptotically. 
#
# However, for a given sample, it can be difficult to
# work out how precise it is.  One suggestion was
# to use it within a Bayesian analysis. See
#
#    https://link.springer.com/article/10.1007%2Fs11749-014-0360-z
#
# The idea is to replace the posterior with
#
#    exp( n*A(theta) ) pi(theta)
#
# of pi(theta) is the prior.  
#
# i) Using a N(0,10) prior, write a function to conduct a random 
# walk MCMC algorithm. The function should take in
#
#    X -- your data
#    s -- bandwidth for the KDE
#    N -- number of Monte Carlo samples to evaluate A
#    nmc -- length of the MCMC chain
#    sig -- random walk variance
#
# and return a set of posterior samples and the acceptance
# rate. You should keep the same Monte Carlo samples for
# each evaluation. 

postfn = function(theta, X, s, Xsamp){
  exp( 2*length(X)*Afn(theta, X, Xsamp, s) ) * dnorm(theta, mean=0, sd=10)     
} 

HD.MCMC = function(X,s,N,nmc,sig){
  
  Xsamp = kde.sim(N,X,s)
  n = length(X)
  par = c()
  acc = rep(0, nmc)
  par[1] = mean(X)
  try.cur = postfn(par[1], X, s, Xsamp)
  for(i in (2:nmc)) {
    theta = par[i-1] + rnorm(1, mean=0, sd=sig)
    try.new = postfn(theta, X, s, Xsamp)
    
    alpha = try.new/try.cur
    
    if(runif(1) < alpha) {
      par[i] = theta
      acc[i] = 1
      try.cur = try.new
    } else {
      par[i] = par[i-1]
    }
  }
  
  return( list(samples=par, acc=acc) )
}

# ii) Experimentally choose sig to give you between 30% and 40% 
# acceptance rate and produce a histogram of the posterior 
# based on 1000 samples after 1000 samples burn-in and thinning
# to every 5th sample.

### Try 6000 steps which will leave us about 1000 samples after 1000 burn-in
### and thining to every 5th sample
set.seed(12071844)
nmc = 6000
sig = 0.5
N.sim = 1000
result = HD.MCMC(eggrate,sbest,N.sim,nmc,sig)
mean(thinner(1000, 5, result$samples))
###-1.494141
mean(thinner(1000, 5, result$acc))
###0.6923077

### Find a sig so that MCMC gives us an acceptance in between 30% and 40%
test.sigs = seq(0.1, 10, 0.5)
accs = c()
found = FALSE
while (!found) {
  set.seed(12071844)
  sig = test.sigs[1]
  test.sigs = test.sigs[-1]
  cat('Trying sig=', sig, '\n')
  result = HD.MCMC(eggrate,sbest,N.sim,nmc,sig)
  theta = mean(thinner(1000, 5, result$samples))
  acc = mean(thinner(1000, 5, result$acc))
  cat('acc:', acc, '\n')
  accs = c(accs, acc)
  if(acc < 0.4 && acc > 0.3) {
    found = TRUE
  }
}
print(accs)
###
### sig = 1.6 will give us an accptance rate = 36%
###
result = HD.MCMC(eggrate,sbest,N.sim,nmc,1.6)
mean(thinner(1000, 5, result$samples))
###-1.582548
mean(thinner(1000, 5, result$acc))
###0.3556444


# Bonus: in fact, MCMC can be used with a stochastic likelihood.
# That is, it still works if you use new Monte Carlo samples
# each time.  Apply this using N = 100, but draw new samples
# for each step in the MCMC. How high can you get your 
# acceptance rate?

HD.MCMC2 = function(X,s,nsamples,nmc,sig){
  N=nsamples
  Xsamp = kde.sim(N,X,s)
  n = length(X)
  par = c()
  acc = rep(0, nmc)
  par[1] = mean(X)
  try.cur = postfn(par[1], X, s, Xsamp)
  for(i in (2:nmc)) {
    Xsamp = kde.sim(N,X,s)
    
    theta = par[i-1] + rnorm(1, mean=0, sd=sig)
    
    try.new = postfn(theta, X, s, Xsamp)
    alpha = try.new/try.cur
    if(runif(1) < alpha) {
      par[i] = theta
      acc[i] = 1
      try.cur = try.new
    } else {
      par[i] = par[i-1]
    }
  }
  
  return( list(samples=par, acc=acc) )
  
}

set.seed(12071844)
result = HD.MCMC2(eggrate,sbest,100,6000,0.5)
mean(thinner(1000, 5, result$samples))
###-0.9787015
mean(thinner(1000, 5, result$acc))
###0.2367632

### Using same sig=0.5, compare to using same sample for all steps, 
### drawing a new sample for each step will decrease the accptance rate by about half.


set.seed(12071844)
result = HD.MCMC2(eggrate,sbest,1000,6000,0.5)
mean(thinner(1000, 5, result$samples))
###-1.534194
mean(thinner(1000, 5, result$acc))
###0.6783217



######################################
# Question 3: Random Effects Methods #
######################################

# In the additional exercises after Lecture 18, we
# saw how Poisson random effects models can be fit
# with MCMC. Here we'll look at extensions of
# this and Frequentist alternatives. 

# In this question we will consider a logistic model
# with random effects. That is, we observe outcomes that
# are either 0 or 1 and write out
#
#  P(Y=1|X,Z) = plogis( beta0 + beta1*x + z)
#
# where x is a covariate and z is a subject-specific
# random effect. 
#
# We can turn this into an over-all probability by
#
# P(Y|X,z) = P(Y=1|X,Z)^Y*(1-P(Y=1|X,Z))^(1-Y)
#
# which evaluates to P(Y=1|X,Z) when Y = 1 and 
# (1-P(Y=1|X,Z)) when Y = 0. 

# Here we regard Z as a latent normal random variable that
# applies to several data points. By way of example, the
# data in 

toenail = read.table('toenail.txt',head=TRUE)
toenail$Subject = as.factor(toenail$Subject)

# gives an indicator of nail health (good or bad) for 
# treatment for foot fungus. This is recorded for 7 visits each
# for 12 Subjects at different times given by number of Months. 
#
# Here we will write Y_ij for the i-th measurement from the j-th 
# patient and say
#
#  P(Y_ij=1|X_ij,Z_j,beta0,beta1) = plogis( beta0 + beta1*X_ij + Z_j)
#
# where $X_ij$ is the time of the ith visit of subject j. Z_j
# is a different value for each subject so that some subjects 
# have consistently worse feet than others depending on the size 
# of Z_j. 
#
# We don't get to see the Z_j, so we will specify that 
#
#   Z_j  ~ N(0, \sigma^2) 
#
# The likelihood of the observed data in group j is
# then 
#
# P(Y_1j,..,Y_7j|X_j) = \int prod( P(Y_ij|X_ij,Z_j,beta0,beta) ) phi(Z_j;0,sigma)dZ_j
#
# and the likelihood for the parameters beta0, beta1, sigma is
#
# l(beta0,beta1,sigma) = sum_j log P(Y_1j,...,Y_7j|X_j)


# a) When only one random effect applies to each notation we can
# evaluate the integral above numerically.  To do this, we will
# use a Gauss-Hermite approximation. This can be obtained through 
# the ecoreg package  (you will need to install it). 

library(ecoreg)

# The function gauss.hermite produce Gauss-Hermite points

gh = gauss.hermite(21)

# Where we can approximate  \int f(x) phi(x) dx by 
#
#   sum(gh[,'Weights']*f(gh[,'Points'])
#
# Note that this is for a standard N(0,1), to make this N(0,s^2)
# you need to multiply gh[,'Points] by s. 

# i) Write a function to evaluate the negative of the log likelihood 
# of the toenail data for any beta0, beta1 and sigma given 
# the vector theta = c(beta0,beta1,sigma)

logistic.nll = function(theta,data){
  #j - subject
  beta0 = theta[1]
  beta1 = theta[2]
  sigma = theta[3]
  j.Z = gh[, 'Points']*sigma
  
  lik = 0
  
  #Loop though all subjects
  for(j in levels(data$Subject)) {
    
    #j.data is 7 visits for patient j
    j.data = data[data$Subject == j,]
    
    j.X = j.data$Month
    j.Y = j.data$Response
    
    j.f = c()
    for(k in (1:length(j.Z))) {
      j.f[k] = prod( dbinom(j.Y, 1, plogis(beta0 + beta1*j.X + j.Z[k])) )
    }
    
    j.P = sum(gh[,'Weights']*j.f)
    lik = lik + log(j.P)
  }
  
  return (-lik)
  
}	
### Test code from the check file
all.equal(logistic.nll(c(-2.5,-1.5,1.1),toenail),142.258,tol=1e-6)



# ii) Hence find the values of beta0, beta1 and sigma that
# maximize the log likelihood. You may use optim or any other
# optimizer that you wish. 

# Good starting values are 
theta = c(-0.2,-0.2,1.2)

### Remember that logistic.nll is negative log value. The optim give the value = 44.7988, so the real log value is -44.7988
### and it is larger than 142.258 at c(-2.5,-1.5,1.1)
ret = optim(theta, logistic.nll, data=toenail)
print(ret)
###(-0.2122702 -0.2081005  1.2311702) value = 44.7988
theta.optim = ret$par


# b) As an alternative, we can use Monte Carlo integration for 
# each random effect. To do this, we'll replace the Gauss-Hermite
# points and weights with a sample of 1000 N(0,1) random values
# and weights 1/1000.  To change this to integrating with 
# respect to a N(0,s^2) you can make the same changes as for 
# Gauss-Hermite quadrature.  Use the same 1000 points for each
# group. 
#
# i) Write a function to calculate the negative log likelihood based 
# on this approximation. You may hard-code the use of N = 1000 
# Monte Carlo points. 

logistic.nll.mc = function(theta,data){
  #j - subject
  beta0 = theta[1]
  beta1 = theta[2]
  sigma = theta[3]
  #Use same 1000 Monte Carlo intergration points for each subject
  mc.num.ps = 1000
  mc.points = rnorm(mc.num.ps, sd=sigma)
  
  lik = 0
  #Loop though all subjects
  for(j in levels(data$Subject)) {
    j.data = data[data$Subject == j,]
    j.X = j.data$Month
    j.Y = j.data$Response
    
    j.f = c()
    for(k in (1:length(mc.points))) {
      j.f[k] = prod( dbinom(j.Y, 1, plogis(beta0 + beta1*j.X + mc.points[k])) )
    }
    
    j.P = sum(j.f)/mc.num.ps
    lik = lik + log(j.P)
  }
  
  return (-lik)
}	

### Test code from the check file
set.seed(2106.1864)
all.equal(logistic.nll.mc(c(-2.5,-1.5,1.1),toenail),140.326,tol=1e-6)

# ii) How would you go about optimizing this likelihood? Do so. 
### I can try different optimaze method when using R optim function

### "Nelder and Mead" gernate an error message 10. indicates degeneracy of the Nelder-Mead simplex.
### (-0.1483854 -0.2073438  1.2241146) 
###result = optim(theta, logistic.nll.mc, data=toenail)

### "BFGS" convergenced wihh few warnings duo to rnorm produced NA 
### (-0.1709841 -0.1787170  1.1692882) value=44.73809
### result = optim(theta, logistic.nll.mc, data=toenail, method='BFGS')

### "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
### (-0.2127365 -0.1841296  1.1725195) value=44.77352
### result = optim(theta, logistic.nll.mc, data=toenail, method='L-BFGS-B')

#ret = optim(theta, logistic.nll.mc, data=toenail, method='SANN')



# c) Given a likelihood, we can always carrying out an MCMC
# analysis to find a posterior. In this case we will
# fix sigma to be 1.2 to make the computation simpler. 
# Setting your prior on beta0 and beta1 to be 
#
# P(beta0,beta1) = N(0,10)*N(0,10)
#
#
# i) Carry out a random walk MCMC for 5000 steps using the 
# likelihood you found in part a using a standard deviation
# of

sds = c(1,0.2)
postfn = function(theta, data) {
  v = -logistic.nll(theta, data) + dnorm(theta[1], mean=0, sd=10, log=TRUE) + dnorm(theta[2], mean=0,sd=10,log=TRUE)
  return (v)
}

# for the steps in each dimension of your random walk. Use the same
# initial value for beta0 and beta1 as in part a.

### A helper function to do a MCMC random walk
random.walk = function(step.count, step.sds, par.names, par.initial, postfn_wrapper) {
  sigma = 1.2
  par = matrix(0, step.count, length(par.names))
  colnames(par) = par.names
  acc = rep(0, step.count)
  par[1,] = par.initial
  
  
  try.cur = postfn_wrapper(par[1, ], sigma, toenail)
  for(i in (2:step.count)) {
    par.new = c()
    for(j in (1:length(par.names))) {
      par.new = c(par.new, par[i-1, j] + rnorm(1, mean=0, sd=step.sds[j]))
    }
    
    try.new = postfn_wrapper(par.new, sigma, toenail)
    alpha = exp(try.new - try.cur)
    
    if(runif(1) < alpha) {
      par[i,] = par.new
      acc[i] = 1
      try.cur = try.new
    } else {
      par[i,] = par[i-1,]
    }
  } 
  return (list(par=par, acc=acc))
}

set.seed(2106.1864)
step.count = 5000
step.sds = sds
par.name = c('beta0', 'beta1')
par.initial = c(-0.2, -0.2)
postfn.wrapper = function(par, sigma, data) {
  theta = c(par[1], par[2], sigma)
  postfn(theta, data)
}
postfn.result = random.walk(step.count, step.sds, par.name, par.initial, postfn.wrapper);
print(apply(postfn.result$par, 2, mean))
print(mean(postfn.result$acc))
### beta=c(-0.1977865 -0.2138705), acceptance=0.2074

### Plot chains of beta0 and beta1
beta0 = postfn.result$par[,1] 
beta1 = postfn.result$par[,2]
plot.new()
par(mfrow=c(2,1))
plot(seq(1:length(beta0)), beta0, type='l', col='blue', xlab=c('iter'), ylab=c('beta0'))
plot(seq(1:length(beta1)), beta1, type='l', col='blue', xlab=c('iter'), ylab=c('beta1'))


# ii) How far is your expected a-posteriori estimator from the maximum 
# likelihood estimate in part a?

beta.all=apply(postfn.result$par, 2, mean)
cat('beta0 is ', beta.all[1], ',', 100*abs((theta.optim[1]-beta.all[1])/theta.optim[1]), '% away from optimazed beta0 from (a)', '\n')
cat('beta1 is ', beta.all[2], ',', 100*abs((theta.optim[2]-beta.all[2])/theta.optim[2]), '% away from optimazed beta1 from (a)', '\n')


# iii) Produce 95% credible intervals for beta0 and beta1
quantile(beta0,c(0.025,0.975))
###(-1.2473211  0.8836025)
quantile(beta1,c(0.025,0.975))
###(-0.37234133 -0.07521758)



# d) Rather than doing the integration numerically, MCMC also
# allows you to simply include the random effects as additional 
# parameters. That is, we can write the joint distribution of
# all of the data, the Z's and and beta0 and beta1 (still keeping
# sigma fixed at 1.2) as 
#
#  P(Y_ij|X,Z_j,beta0,beta1)P(Z_j)P(beta0)P(beta1)
#
# And conducting a random walk on beta0,beta1,Z_1,..,Z_K

# i) Write a function to calculate this log posterior using the same
# prior as as in part c. 

logistic = function(theta,data,Z){
  ### j - subject
  beta0 = theta[1]
  beta1 = theta[2]

  ### Loop though all subjects
  subjects = levels(data$Subject)
  lik = 0
  for(j in 1:length(subjects)) {
    j.data = data[data$Subject == subjects[j],]
    j.X = j.data$Month
    j.Y = j.data$Response
    j.P = prod( dbinom(j.Y, 1, plogis(beta0 + beta1*j.X + Z[j])) )
    lik = lik + log(j.P)
  }
  return (-lik)
}	

joint.posterior = function(theta,Z,data){
  #fix sigma to 1.2
  sigma = 1.2
  v = -logistic(theta, data, Z) + dnorm(theta[1], mean=0, sd=10, log=TRUE) + dnorm(theta[2], mean=0, sd=10, log=TRUE) + sum(dnorm(Z, mean=0, sd=sigma, log=TRUE))
  return (v)
}
### Test code from the check file
set.seed(2106.1864)
Z = rnorm(12,sd=0.5)
all.equal(joint.posterior(c(0.1,-0.25),Z,toenail),-69.32814,tol=1e-7)

# and carry out a random walk MCMC using step sizes with variance

sds2 = sds/4

# for theta and steps of standard deviation 0.5 for Z.
# Initialize Z to be all zero.

set.seed(2106.1864)
step.count = 5000
step.sds = c(sds2, rep(0.5, 12))
par.name = c('beta0', 'beta1')
for(i in (1:12)){
  par.name = c(par.name, paste(c('z', i), collapse='')) 
}
par.initial = c(-0.2, -0.2, rep(0, 12))
joint.posterior.wrapper = function(par, sigma, data){
  theta = c(par[1], par[2], sigma)
  joint.posterior(theta, par[seq(3,14)], data)
}

joint.posterior.result = random.walk(step.count, step.sds, par.name, par.initial, joint.posterior.wrapper)
print(apply(joint.posterior.result$par, 2, mean))
print(mean(joint.posterior.result$acc))
### beta=c(-0.26079672 -0.21975822), acceptance=0.215

### Plot chains of beta0 and beta1
beta0 = joint.posterior.result$par[,1] 
beta1 = joint.posterior.result$par[,2]
plot(seq(1:length(beta0)), beta0, type='l', col='blue', xlab=c('iter'), ylab=c('beta0'))
plot(seq(1:length(beta1)), beta1, type='l', col='blue', xlab=c('iter'), ylab=c('beta1'))


# ii) Have your expected a posteriori estimates for bet0
# and beta1 moved outside the credible interval found 
# in part d?

### My calculation shows beta0 and beta1 estimated by using joined posteriori 
### are in crediable inerval found in part d.
### beta0 = -0.26079672 is in (-1.2473211  0.8836025) from (d) 
### beta1 = -0.21975822 is in (-0.37234133 -0.07521758) from (d)


# iii) Plot a histogram of the posterior distribution of Z_1
# based on every 5th draw from the chain after dropping the first 
# 1000 draws
Z1 = joint.posterior.result$par[,'z1']
Z1.thin = thinner(1000, 5, Z1)
par(mfrow=c(2,1))
hist(Z1.thin, 100)
plot(seq(1:length(Z1)), Z1, type='l', col='blue', xlab=c('iter'))

# e) i) Can we calculate the conditional distribution of Z_1
# as a frequentist? Here is one last way of generating random
# variables that takes inspiration from importance sampling
# in Monte Carlo integration. 
#
# Specifically, suppose we want to sample from g(x), but can
# readily draw from f(x). The idea is to draw X_1,...,X_N but 
# to then re-sample these X's with weights 
#
#  W_i = g(X_i)/f(X_i)
#
# (the sample() function allows you to specify a vector 'prob'
# of weights). This means that the density of the resulting
# X's is 
#
#   P(X selected)*P(X in dx) = W f(x) = g(x)/f(x) * f(x) = g(x)
#
# To get a posterior distribution of Z1, first, use the beta0,
# beta1 and sigma from part a.  Then 
#
#   1. Generate 1000 samples of Z from N(0,sigma^2) (you may use rnorm)
#   2. Resample the Z with replacement, but with weights given by 
#        prod(  P(Y_i1 | X_i1, Z,beta0,beta1) )
#

# Write a function to obtain this samples using beta0, beta1
# and sigma calculated in part a but using the data in 


Subj1 = toenail[1:7,]

ranef.post = function(N,theta,Subj1){
  Z = rnorm(N, mean=0, sd=theta[3])
  beta0 = theta[1]
  beta1 = theta[2]
  months = Subj1$Month
  Y = Subj1$Response
  
  ### Calculate the weights
  W = c()
  for (i in (1:length(Z))){
    W[i] = prod( dbinom(Y,1,plogis(beta0 + beta1*months + Z[i])) )
  }
  Z = sample(Z, N, replace=TRUE, prob=W)
  return (Z)
}

### Test from the check file
set.seed(2106.1864)
all.equal(ranef.post(1000,c(-0.19,-0.21,1),Subj1)[1:3],
          c(1.07871834,0.12607701,0.05866123),tol=1e-6)



# ii) Produce a histogram based on your 1000 resulting samples. How does
# this compare to the posterior distribution in part b?

### Use beta0, beta1 and sigma got from a)
theta_a = c(-0.2122702, -0.2081005,  1.2311702)
Z = ranef.post(1000, theta_a, Subj1)
plot.new()
par(mfrow=c(1,1))
hist(Z,100,prob=TRUE)
sorted.Z = sort(Z)
lines(sorted.Z, dnorm(sorted.Z, mean=0, sd=theta_a[3]), col='blue', xlim=c(-3,3))



# iii) This approach to generating random numbers is known as sequential
# Monte Carlo (SMC) and is particularly useful when you want to keep
# updating distributions -- you can take a sample from f_1(x) to f_2(x) 
# and then update again to f_3(x).  This arises, for example, if you
# are tracking a noisily observed process. 
#
# Here, however, we are only generating from one distribution. Is it a 
# good or a bad idea? Why?

##### This could possibly be a bad idea because if we're generating from a one distribution
##### when tracking a noisily observed process we could end up with ineffective data.
##### This is problematic if the simulation is using bad data because it would tell us very
##### little information.


# f)  Describe how you would go about constructing a bootstrap distribution
# for this model, accounting for the fact that you might have different
# subjects next time.  You do not need to carry this out. 

#### So for this model we can construct a bootstrap by first sampling the last 2 columns
#### of the data set. In order to do this we need to decide how many samples we want. 
#### So we can do the n=number of subjects*number of visits. So we create a matrix of nx3 columns
#### where the last two columns are the sampled from the toenail[,2:3]. 

num.subjects = 5
num.visits = 3

boot.data = matrix(0, nrow=num.subjects*num.visits, ncol=3)
visits = toenail[,2:3]
ind = sample(nrow(toenail), size=num.subjects*num.visits*2, replace=TRUE)
boot.visit = visits[ind,]
for(i in (1:num.subjects)) {
  for(j in (1:num.visits)) {
    ind = (i-1)*num.visits + j
    print(ind)
    boot.data[ind, 1] = i
    boot.data[ind, 2] = boot.visit[ind, 1]
    boot.data[ind, 3] = boot.visit[ind, 2]
  }
}
colnames(boot.data)=c('Subject', 'Response', 'Month')
boot.data = boot.data[order(boot.data[,1], boot.data[,3]),]
print(boot.data)

## After printing this data we need check for having duplicate visits in a month with different responses for
## same patient. If this happens remove that data point and sample again to get a different month.


#######################
# BONUS: Optimization #
#######################

# This repeats the final exercise in Lab 7 (which the labs
# did not get to).  Additional 10 marks for an answer within
# 10% of the true values


# The data in FhN.csv come from an example that Giles was 
# working on for his most recent book. They arise from
# a model of the way neurons transmit signals and are given
# by an ordinary differential equation. 
#
# This is a simplification of a more complex model than won 
# Alan Hodgkin and Andrew Huxley a Nobel prize in 1963 (so
# it's fitting that Giles is writing this lab while in 
# Stockholm). 
#
# The data is in 

FHN = read.csv('fhn.csv',head=FALSE)

# It contains a set of times

t = FHN[,1]

# and two-columns of values

V = FHN[,2:3]

# so we can look at

matplot(t,V,type='l')

# These values were produced using the following function, 
# which relies on the deSolve package to solve ODE equations

library(deSolve)


FHN.fn = function(t,p){
  
  x0 = p[4:5]
  
  fhnfunode = function(t,x,p)
  {
    r = x;
    r[1] = p[3]*(x[1] - x[1]^3/3 + x[2])
    r[2] = -(x[1] -p[1] + p[2]*x[2])/p[3]
    return(list(dx=r))
  }
  
  res = lsoda(x0,t,fhnfunode,p)
  
  return( res[,2:3] )
}

# This depends on a vector p of 5 parameters, 
# for instance, if you set
FHN.fn(t,p.actual)
p = c(1,1,3,0.5,-1.2)

# you get 

matplot(t, FHN.fn(t,p),type='l')

# The game here is to find the values of p that generated
# the data. There is no noise in the data, so the right
# p should produce V exactly. 

# To do this, use the optim function in R, and experiment
# with different optimization solvers and/or starting points. 
# You will need to define an objective function -- I recommend
# squared error and trying to find where it is zero. 

# Two things might help:
#
#  1. I will tell you the true parameters are in the ranges
#      [-1,1],  [-1,1],   [0.5,5],  [-2,2],  [-2, 2]
#
#  2. You will find that at bad values of p, FHN.fn will
#  produce an error. This will cause your optimization to
#  to stop. However you can get around this with the
#  try-catch construction. 
#
#   Here's an example

##############################################
for(i in 1:3){
  try({
    m = b/8 
  })
  print(i)
}

# Although there is an error (we have never defined b) the for 
# loop does not terminate.  You can set a number of lines of 
# code inside the {} in try.   
#
# You might use this to try evaluating at some p, for example,
# and return Inf, or a really big number if that p produces
# an error. 
#
# You may search for good starting values, but describe how you
# decided to start from there, and include code to get from them
# to your final estimate. 


###################################
#####        Process          #####
###################################
###
### 1. Build a matrix of all permutations of p within the range of true parameters given above
###   a. Use small increments (step) for each of the 5 p values. (I used step of 0.5 to start)
###   b. Combine these permutations to create a mass matrix of p that we can use as starting values.
### 2. Loop though all start points and hope R optim function can minimize our object function *locally*
### 3. Since this is an expensive computation job, I ran this code on two computers (my laptop and a 
###    machine in Upson Hall). 
###   a. With step=0.5 my matrix of starting p values has around 20251 starting vectors to optimize.
###   b. I split the task 50/50 for the computers. Machine 1 ran the first ~10000 starting points and machine 2
###      ran the rest of the p values.
### 4. It took around 1.5 hours to find a p such that the max error rate of between the data V and FHN.fn(t,p) of
###    of my p was less than 10%. The two machines both only managed to run through ~5000 starting p values.
###    I would think it would take several more hours to go through all possible starting p values given step=0.5.
###
### NOTE: I probably should have made the code save some of the important points so that I can change the start
###       points to those points in case I wanted to run other iterations with different step values. 
###       For example, I could save start points where the optimized error rate is less than 0.20 (20%)
###       so that I can later build starting points with a smaller steps around this p.
###
########################
####     RESULT    #####
########################
###
### 1. The first computer found that with starting p values of (0.0, 0.0, 3.0, -0.5, 1.0) the
###    optimization would give us a final result of (0.1998291, 0.1986810, 3.0002912, -0.9988947,0.9966097)
###    These values of p give us FHN.fn(t,p) that has a maximum error 8.02% from the true values in V.
###
### 2. The second computer found that with starting p values of (-1.0, 0.0, 2.5, 0.0, 0.5) the
###    optimization would give us a final result of (0.1993429,  0.1992523,  3.0005259, -1.0002762,  1.0024309)
###    These values of p give us FHN.fn(t,p) that has a maximum error 7.65% from the true values in V.
###   


#### This is the objective function we want to optimize
#### Instead of doing squared error, I just looked for the maximum error between
#### true values and our predicted values. Then if we try to minimize that we can find
#### what values of p gives us the closest true values in V.
object.func = function(p) {
  result = FHN.fn(t, p)
  error = abs(V-result)/V
  max(error)
}


#### This creates a matrix of starting values of p that we want to optimize around locally
#### step=0.5 gives us a total of 20251 total permutations of starting p values we can use.
#### Smaller step size would gives us more permutations.
mess = matrix(0, ncol=5)
step = 0.5
row = c()
for(k1 in seq(-1,1, step)) {
  row[1] = k1
  for(k2 in seq(-1,1, step)){
    row[2] = k2
    for(k3 in seq(0.5, 5, step)){
      row[3] = k3
      for(k4 in seq(-2, 2, step)) {
        row[4] = k4
        for (k5 in seq(-2, 2, step)) {
          row[5]=k5
          mess = rbind(mess, row)
        }
      }
    }
  }
}
print (mess)
print (nrow(mess))

min.start = 0
min.optim = 0
min.value = 1000000

#### This is the optimization function where it optimizes locally for each starting p of the matrix.
#### I ran (10000:length(mess)) on another computer
#### Code is commented out because of time it takes to find the optimal values
#### Results are below

# for(i in (1:10000)) {
#   cat(format(Sys.time()), 'i=', i, '\n')
#   p = mess[i,]
#   try({
#     result = optim(p, object.func, control=list(maxit=5000))
#     if(min.value > result$value) { 
#       min.start = c(p)
#       min.optim = result
#       min.value = result$value
#       cat('p=', p, 'par:', result$par, 'value:', result$value, 'conv:', result$convergence, 'count:', result$counts,'\n')
#     }
#   })
# }

verify = function(p) {
  result = FHN.fn(t, p)
  error = abs(result - V)/V
  max(error)
}


#Error of .07647317
p1 = c(0.1993429,  0.1992523,  3.0005259, -1.0002762,  1.0024309)
verify(p1)
par(mfrow =c(1,2))
matplot(t,V,type='l')
matplot(t, FHN.fn(t,p1),type='l')

#Errr of 0.08015433
p2 = c(0.1998291, 0.1986810, 3.0002912, -0.9988947,0.9966097)
verify(p2)
par(mfrow =c(1,2))
matplot(t,V,type='l')
matplot(t, FHN.fn(t,p2),type='l')

## We can see that the our two graphs of FHN.fn(t,p) vs. t is almost identical to 
## the graph of actual values of V vs. t

p.actual = c(.2,.2,3,-1,1)
verify(p.actual)


# The two p's I found seem to be around a certain point and were optimized 
# So i did few experiments around points c(0.1993429,  0.1992523,  3.0005259, -1.0002762,  1.0024309). 
# I did a brave guess that the real p must be c(0.2, 0.2, 3., -1,1). 
# I tested it and I got the max error rate is 2.833993e-13. So c(0.2, 0.2, 3., -1,1) 
# must be the p used to generate the V in the file, or it is one of p can generate the V if the function
# is somehow symmetric about p. 

v = c(V[,1])
v
res = (FHN.fn(t,p.actual)[,1])
res
all.equal(v, res, tol=1e-6)





