# Spectral-regularized-two-sample-test

The main file that contain the functions required to compute the regularized test is in "spectral_test.R" which contains the function compute_test(). 

In order to run the function compute_test(), you need to first include the file "spectral_test.R" by including the line 
> source("spectral_test.R")
Then you need to specifiy the following parameters in your code (use desired number for each parameter):
> iter <- 100  # number of iterations
> s <- 100 # number of samples splits to estimate the covariance operator
> n <- 500
> m <- 500
> method <- 2 # 1 is Tikhonov , 2 is showalter
> kernel_type <- 1 #1 is gaussian , 2 is laplace
> num_perm <- 60  #number of permuations used to estimate the test threshod
> Lambda <- 10^seq(-6,1,0.75) # Array for the values of lambda used for adapatation
> h_mult_arr <- 10^seq(-2,2,0.5)   #Array for bandwidths used for adapatation
> alpha <- 0.05 # Type I error 
