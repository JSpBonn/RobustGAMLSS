# rm(list=ls())

# setwd("//folderapplication")
source("Robust_gamboostLSS_Families.R") # loading some packages and the robust Gamma families object for GAMLSS gradient boosting
# setwd("//foldersimulations")
library(mvtnorm) # for generating the design matrix within the simulations
# library("parallel") # for cluster
# library("GJRM") # for Aeberhard et al., robust GAMLSS comparison


Simfunc <- function(id){
  corrupted2 <- c(0.0,0.025,0.05,0.075,0.1,0.15,0.2)
  out.tab <-  vector("list",length(corrupted2))
  # out.tab2 <-  vector("list",length(corrupted2)) # uncomment for Aeberhard et al.
  half=1000
  n=half
  methodnumber = 6
  numbercorrupted <- ceiling(n*corrupted2) # number of corrupted observations
  sigma_ver <- 4 
  cor_strong <-  4
  
  
  #### Data generating
  p=5  # XXX comment for high dimensional
  # p=1000 # XXX uncomment for high dimensional
  
  correl <-0.50 # correlation Toeplitz
  intercept = c(1,0.5) # Simulations setting intercept
  help = numeric(p)
  for (k in 1:p){ help[k]=correl^(k-1)
  }
  toeplitzcor = stats::toeplitz(help)
  sigma=1
  
  mean_X = rep(intercept[1],length = p)
  sigma_X = matrix(toeplitzcor,p,p)
  
  diag(sigma_X)=sigma
  
  for (i in 1:length(corrupted2)) { # comment when running only an example (start loop of i)
    time_matrix <- matrix(0,ncol = methodnumber,nrow = 1)    
    ws <- c(rep(1,half), rep(0,half))
    i=i # when running an example, use i=1 for uncorrupted data and in general an index of the vector c(0.0,0.025,0.05,0.075,0.1,0.15,0.2) (amount of corruption)
    # id=1 # uncomment when running only an example (only for set.seed)
    
    n <- length(ws)
    
    set.seed(id+10000)
    X = rmvnorm(n,mean = mean_X, sigma = sigma_X) 
    x1 <- X[,1]  
    x2 <- X[,2]  
    x3 <- X[,3]  
    x4 <- X[,4] 
    x5 <- X[,5]  
  
    # colnames(X) <- paste0("x",1:p) # XXX uncomment for high-dimensional
    
    #### the maximal number of stopping iterations, for testing use smaller number (Supplement Table S8 suggests that 200-600 will be the optimal stopping iteration for the low-dimensional case)
    stopping = 2000 # XXX comment for high-dimensional
    # stopping = 1500 # XXX uncomment for high-dimensional
    
    
    
    
    
    
    
    toydata <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5) # XXX comment for high-dimensional
    # toydata <- data.frame(X) # XXX uncomment for high-dimensional
    toydata$y <-rgamma(n, scale = exp( 1 + 1.5 * x1 - 0.75 * x2)/exp(0.5 - 0.25 * x1 + 0.5 * x3) , shape = exp(0.5 - 0.25 * x1 + 0.5 * x3))
   
    toydata2 <- toydata[1:half,]
    
    scale  = egamma(toydata$y[1:half])$parameters[2]
    shape=egamma(toydata$y[1:half])$parameters[1]
    
    a <- (scale^2*shape)^0.5
   
    set.seed(id+10000*i+100000)
    if (i>1) { #corruption
      toydata$y[1:numbercorrupted[i]] <- toydata$y[1:numbercorrupted[i]]+cor_strong*a # skewed corruption
      toydata$y[(half+1):(half+numbercorrupted[i])] <- toydata$y[(half+1):(half+numbercorrupted[i])]+cor_strong*a # skewed corruption
     
    }
    
    n_test=1000
    set.seed(id+20000)
    
    X = rmvnorm(n_test,mean = mean_X, sigma = sigma_X) 
    x1 <- X[,1] 
    x2 <- X[,2] 
    x3 <- X[,3] 
    x4 <- X[,4] 
    x5 <- X[,5] 
    
   
    
    toydata_test <- data.frame(x1 = x1, x2 = x2, x3 = x3 ,x4 = x4 , x5 = x5) # XXX comment for high-dimensional
    toydata_test <-data.frame(X) # XXX uncomment for high-dimensional
    toydata_test$y <-rgamma(n= n_test, scale = exp( 1 + 1.5 * x1 - 0.75 * x2)/exp(0.5 - 0.25 * x1 + 0.5 * x3) , shape = exp(0.5 - 0.25 * x1 + 0.5 * x3))
    
    
    #### generate the robustness constants for the robust method
    tuning_s <- c(1,2,20, c_generate_Gamma(toydata$y[1:half],tau=0.01), c_generate_Gamma(toydata$y[1:half],tau=0.05), c_generate_Gamma(toydata$y[1:half],tau=0.10))
    ####  tuning_s[1] is only placeholder for the non-robust method
    
    # coefmatrix_all <- matrix(0,nrow=methodnumber,ncol=2*p+2)
    coefmatrix_optimal <- matrix(0,nrow=methodnumber,ncol=2*p+2)
    cvr <- matrix(0,ncol=methodnumber,nrow=3)
    # MAE_MSE_all <- matrix(0,ncol=methodnumber,nrow=4)
    # MAE_MSE_optimal<- matrix(0,ncol=methodnumber,nrow=4)
    
    # MAE_MSE_all_test <- matrix(0,ncol=methodnumber,nrow=4)
    # MAE_MSE_optimal_test<- matrix(0,ncol=methodnumber,nrow=4)
    
    # log_likelihood_conv_test <-  matrix(0,ncol=methodnumber,nrow=1)
    
    log_likelihood_opti_test <-  matrix(0,ncol=methodnumber,nrow=1)
    
    ##################################################### ##################################################### ##################################################### #####################################################
    ##################################################### ##################################################### ##################################################### #####################################################
    ##################################################### ##################################################### ##################################################### #####################################################
    
    method=1
    set.seed(1234+id)
    time_a <-proc.time()[3]  
    gam1 <-glmboostLSS(y~.,data=toydata,method = "noncyclic",control = boost_control(mstop=stopping,risk = "oobag"),families = GammaLSS(),weights = ws)
    # gam1 <-glmboostLSS(y~.,data=toydata,method = "noncyclic",control = boost_control(mstop=stopping,risk = "oobag",nu=c(0.05,0.1)),families =  GammaLSS(stabilization = "L2"),weights = ws) # uncomment for high-dimensional
    cvr[1,method]  <- max(as.integer(which.min(risk(gam1,merge=T))-2),1)
    time_b <- proc.time()[3]                                                                                                           
    time_matrix[1,method] <- time_b-time_a  
    
    # #### converged coefficients 
    # if (min(mstop(gam1[stopping],parameter = "mu"),mstop(gam1[stopping],parameter = "sigma"))>0) {
    #  coefmatrix_all[method,] <- c(coef(gam1[stopping]$mu,off2int = T,which=""),coef(gam1[stopping]$sigma,off2int = T,which=""))
    # }
    # if (min(mstop(gam1[stopping],parameter = "mu"),mstop(gam1[stopping],parameter = "sigma"))==0) {
    #  coefmatrix_all[method,] <- c(coef(gam1[stopping]$mu,off2int = T,which=""),gam1[stopping]$sigma$offset,rep(0,p))
    # }
    # mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata[1:half,1:p]))
    # mu_est <- as.matrix(coef(gam1[stopping]$mu,off2int = TRUE,which=""))
    #  if (min(mstop(gam1[stopping],parameter = "mu"),mstop(gam1[stopping],parameter = "sigma"))>0) {
    #   sigma_est <- as.matrix(coef(gam1[stopping],off2int = TRUE,which="")$sigma)
    # }
    #  if (min(mstop(gam1[stopping],parameter = "mu"),mstop(gam1[stopping],parameter = "sigma"))==0) {
    #    sigma_est <- as.matrix(c(gam1[stopping]$sigma$offset,rep(0,p)))
    # }
    # mu_resi <- mat_all%*%mu_est
    # sigma_resi <-  mat_all%*%sigma_est
    #
    #  true_mu <- 1 + 1.5 * toydata2$x1 - 0.75 * toydata2$x2
    #
    #  MAE_MSE_all[1,method]<- sum(abs(mu_resi-true_mu))
    #  MAE_MSE_all[2,method]<- sum((mu_resi-true_mu)^2)
    #
    #  true_sigma <- 0.5 - 0.25 *toydata2$x1 + 0.5 * toydata2$x3
    #
    #  MAE_MSE_all[3,method]<- sum(abs(sigma_resi-true_sigma))
    #  MAE_MSE_all[4,method]<- sum((sigma_resi-true_sigma)^2)
    #
    #  mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata_test[,1:5]))
    #  mu_resi <- mat_all%*%mu_est
    #  sigma_resi <-  mat_all%*%sigma_est
    #
    #  true_mu <- 1 + 1.5 * toydata_test$x1 - 0.75 * toydata_test$x2
    #  
    #  MAE_MSE_all_test[1,method]<- sum(abs(mu_resi-true_mu))
    #  MAE_MSE_all_test[2,method]<- sum((mu_resi-true_mu)^2)
    #  
    #  true_sigma <- 0.5 - 0.25 *toydata_test$x1 + 0.5 * toydata_test$x3
    #  
    #  MAE_MSE_all_test[3,method]<- sum(abs(sigma_resi-true_sigma))
    #  MAE_MSE_all_test[4,method]<- sum((sigma_resi-true_sigma)^2)
    #  
    #  log_likelihood_conv_test[method] <-    -sum(dgamma(x = toydata_test$y,scale = exp(mu_resi)/exp(sigma_resi),shape = exp(sigma_resi),log = TRUE))
    #  rm(list=c("mat_all","mu_resi","sigma_resi"))
    
    
    if (min(mstop(gam1[cvr[1,method]],parameter = "mu"),mstop(gam1[cvr[1,method]],parameter = "sigma"))>0) {
      coefmatrix_optimal[method,] <- c(coef(gam1[cvr[1,method]]$mu,off2int = T,which=""),coef(gam1[cvr[1,method]]$sigma,off2int = T,which=""))
    }
    if (min(mstop(gam1[cvr[1,method]],parameter = "mu"),mstop(gam1[cvr[1,method]],parameter = "sigma"))==0) {
      coefmatrix_optimal[method,] <- c(coef(gam1[cvr[1,method]]$mu,off2int = T,which=""),gam1[cvr[1,method]]$sigma$offset,rep(0,p))
    }
    
    mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata[1:half,1:p]))
    mu_est <- as.matrix(coef(gam1[cvr[1,method]]$mu,off2int = TRUE,which=""))
    
    
    if (min(mstop(gam1[cvr[1,method]],parameter = "mu"),mstop(gam1[cvr[1,method]],parameter = "sigma"))>0) {
      sigma_est <- as.matrix(coef(gam1[cvr[1,method]]$sigma,off2int = TRUE,which=""))
    }
    if (min(mstop(gam1[cvr[1,method]],parameter = "mu"),mstop(gam1[cvr[1,method]],parameter = "sigma"))==0) {
      sigma_est <- as.matrix(c(gam1[cvr[1,method]]$sigma$offset,rep(0,p)))
    }
    
    
    mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata_test[,1:p]))
    mu_resi <- mat_all%*%mu_est
    sigma_resi <-  mat_all%*%sigma_est
    
    log_likelihood_opti_test[method] <- -sum(dgamma(x =toydata_test$y,scale = exp(mu_resi)/exp(sigma_resi),shape = exp(sigma_resi),log = TRUE))
    rm(list=c("mat_all","mu_resi","sigma_resi"))
    
    cvr[2,method]<- mstop(gam1[cvr[1,method]]$mu)
    cvr[3,method]<- mstop(gam1[cvr[1,method]]$sigma)
    
    rm("gam1")
    
    ##################################################### ##################################################### ##################################################### #####################################################
    ##################################################### ##################################################### ##################################################### #####################################################
    ##################################################### ##################################################### ##################################################### #####################################################
    
    for (z in 2:6) { # start of loop over robust methods, # XXX comment for high-dimensional
    #for (z in 3:6) { # start of loop over robust methods# not reasonable starting value for tuning_s[2] # XXX uncomment for high-dimensional
      method = z
      c_0=tuning_s[method]
      
      set.seed(1234+id) # robust GammaLSS
      time_a <-proc.time()[3]  
      gam1 <-glmboostLSS(y~.,data=toydata,method = "noncyclic",control = boost_control(mstop=stopping,risk = "oobag"),weights = ws,families = robust_GammaLSS(rob=c_0))  # XXX comment for high-dimensional
      # gam1 <-glmboostLSS(y~.,data=toydata,method = "noncyclic",control = boost_control(mstop=stopping,risk = "oobag",nu=c(0.05,0.1)),weights = ws,families = robust_GammaLSS(stabilization = "L2",rob=c_0)) # XXX uncomment for high-dimensional
      
      
      cvr[1,method]  <- max(as.integer(which.min(risk(gam1,merge=T))-2),1)
      time_b <- proc.time()[3]                                                                                                           
      time_matrix[1,method] <- time_b-time_a  
      # 
      #   if (min(mstop(gam1[stopping],parameter = "mu"),mstop(gam1[stopping],parameter = "sigma"))>0) {
      #   coefmatrix_all[method,] <- c(coef(gam1[stopping]$mu,off2int = T,which=""),coef(gam1[stopping]$sigma,off2int = T,which=""))
      # }
      # if (min(mstop(gam1[stopping],parameter = "mu"),mstop(gam1[stopping],parameter = "sigma"))==0) {
      #   coefmatrix_all[method,] <- c(coef(gam1[stopping]$mu,off2int = T,which=""),gam1[stopping]$sigma$offset,rep(0,p))
      # }
      # 
      # 
      # mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata[1:half,1:p]))
      # mu_est <- as.matrix(coef(gam1[stopping]$mu,off2int = TRUE,which=""))
      # if (min(mstop(gam1[stopping],parameter = "mu"),mstop(gam1[stopping],parameter = "sigma"))>0) {
      #   sigma_est <- as.matrix(coef(gam1[stopping],off2int = TRUE,which="")$sigma)
      # }
      # if (min(mstop(gam1[stopping],parameter = "mu"),mstop(gam1[stopping],parameter = "sigma"))==0) {
      #   sigma_est <- as.matrix(c(gam1[stopping]$sigma$offset,rep(0,p)))
      # }
      # mu_resi <- mat_all%*%mu_est
      # sigma_resi <-  mat_all%*%sigma_est
      # 
      # true_mu <- 1 + 1.5 * toydata2$x1 - 0.75 * toydata2$x2
      # 
      # MAE_MSE_all[1,method]<- sum(abs(mu_resi-true_mu))
      # MAE_MSE_all[2,method]<- sum((mu_resi-true_mu)^2)
      # 
      # true_sigma <- 0.5 - 0.25 *toydata2$x1 + 0.5 * toydata2$x3
      # 
      # MAE_MSE_all[3,method]<- sum(abs(sigma_resi-true_sigma))
      # MAE_MSE_all[4,method]<- sum((sigma_resi-true_sigma)^2)
      # 
      # 
      # mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata_test[,1:5]))
      # mu_resi <- mat_all%*%mu_est
      # sigma_resi <-  mat_all%*%sigma_est
      # 
      # true_mu <- 1 + 1.5 * toydata_test$x1 - 0.75 * toydata_test$x2
      # 
      # MAE_MSE_all_test[1,method]<- sum(abs(mu_resi-true_mu))
      # MAE_MSE_all_test[2,method]<- sum((mu_resi-true_mu)^2)
      # 
      # true_sigma <- 0.5 - 0.25 *toydata_test$x1 + 0.5 * toydata_test$x3
      # 
      # MAE_MSE_all_test[3,method]<- sum(abs(sigma_resi-true_sigma))
      # MAE_MSE_all_test[4,method]<- sum((sigma_resi-true_sigma)^2)
      # 
      # log_likelihood_conv_test[method] <-  -sum(dgamma(x = toydata_test$y,scale = exp(mu_resi)/exp(sigma_resi),shape = exp(sigma_resi),log = TRUE))
      # rm(list=c("mat_all","mu_resi","sigma_resi"))
      # 
      
      
      if (min(mstop(gam1[cvr[1,method]],parameter = "mu"),mstop(gam1[cvr[1,method]],parameter = "sigma"))>0) {
        coefmatrix_optimal[method,] <- c(coef(gam1[cvr[1,method]]$mu,off2int = T,which=""),coef(gam1[cvr[1,method]]$sigma,off2int = T,which=""))
      }
      if (min(mstop(gam1[cvr[1,method]],parameter = "mu"),mstop(gam1[cvr[1,method]],parameter = "sigma"))==0) {
        coefmatrix_optimal[method,] <- c(coef(gam1[cvr[1,method]]$mu,off2int = T,which=""),gam1[cvr[1,method]]$sigma$offset,rep(0,p))
      }
      
      mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata[1:half,1:p]))
      mu_est <- as.matrix(coef(gam1[cvr[1,method]]$mu,off2int = TRUE,which=""))
      
      if (min(mstop(gam1[cvr[1,method]],parameter = "mu"),mstop(gam1[cvr[1,method]],parameter = "sigma"))>0) {
        sigma_est <- as.matrix(coef(gam1[cvr[1,method]]$sigma,off2int = TRUE,which=""))
      }
      if (min(mstop(gam1[cvr[1,method]],parameter = "mu"),mstop(gam1[cvr[1,method]],parameter = "sigma"))==0) {
        sigma_est <- as.matrix(gam1[cvr[1,method]]$sigma$offset,rep(0,p))
      }
      
      mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata_test[,1:p]))
      mu_resi <- mat_all%*%mu_est
      sigma_resi <-  mat_all%*%sigma_est
      
      log_likelihood_opti_test[method] <-  -sum(dgamma(x =toydata_test$y,scale = exp(mu_resi)/exp(sigma_resi),shape = exp(sigma_resi),log = TRUE))
      rm(list=c("mat_all","mu_resi","sigma_resi"))
      
      cvr[2,method]<- mstop(gam1[cvr[1,method]]$mu)
      cvr[3,method]<- mstop(gam1[cvr[1,method]]$sigma)
      
      rm("gam1")
      
    } # end of loop over robust methods
    
    out.tab[[i]] <- list(
      log_likelihood_opti_test=log_likelihood_opti_test,
      #log_likelihood_conv_test=log_likelihood_conv_test,
      #MAE_MSE_all_test=MAE_MSE_all_test,MAE_MSE_all=MAE_MSE_all,
      coefmatrix_optimal=coefmatrix_optimal,
      #coefmatrix_all=coefmatrix_all, 
      cvr=cvr,tuning_s=tuning_s,time_matrix=time_matrix)
    
    print(i)
    
    
    ####################################################################################################################################################################################################################
    ####################################################################################################################################################################################################################
    
    #### following interesting as comparison
    ####################################################################################################################################################################################################################
    # 
    # #### only for low dimensinoal cases (and comparably small amounts of parameters p)
    # #### Aeberhard et al., 2021 robust method for GAMLSS
    # # library("GJRM") # for Aeberhard et al.  
    #
    # coefmatrix_all <- matrix(0,nrow=methodnumber,ncol=2*p+2)
    # coefmatrix_optimal <- matrix(0,nrow=methodnumber,ncol=2*p+2)
    # MAE_MSE_all <- matrix(0,ncol=methodnumber,nrow=4)
    # MAE_MSE_all_test <- matrix(0,ncol=methodnumber,nrow=4)
    # 
    # log_likelihood_conv_test <-  matrix(0,ncol=methodnumber,nrow=1)
    # time_matrix <- matrix(0,ncol = methodnumber,nrow = 1)       
    # 
    # for (z in 2:6) {
    #   method = z
    #   set.seed(1234+id)
    #
    #   time_a <-proc.time()[3]
    #   eq.mu <- as.formula(paste("y~",paste(colnames(toydata[,1:p]),collapse = "+"),sep=""))
    #   eq.s2 <- as.formula(paste("~",paste(colnames(toydata[,1:p]),collapse = "+"),sep=""))
    #
    #   fl    <- list(eq.mu, eq.s2)
    #
    #   out1 <- gamlss(fl, data = toydata[1:half,] , margin = "GA", robust = FALSE)
    #
    #   rc=tuning_s[method]
    #   bounds <- rob.int(out1, rc, l.grid=1000, tol=1e-4, var.range = c(-1000,1000)) # sometimes tuning is necessary here
    #   #bounds
    #   lB <- bounds["lB"]
    #   uB <- bounds["uB"]
    #
    #   out.rpb <- gamlss(fl, margin="GA", data=toydata[1:half,], robust=TRUE, sp.method="efs", rc=rc, lB=lB, uB=uB) # For the Gamma distribution, the choice of the robustness constant might have only a small impact
    #   time_b <- proc.time()[3]
    #   time_matrix[1,method] <- time_b-time_a
    #
    #   coefmatrix_all[method,]<- out.rpb$coefficients
    #
    #
    #  # mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata[1:half,1:5]))
    #   mu_est <- as.matrix(out.rpb$coefficients[1:(p+1)])
    #   sigma_est <- as.matrix(out.rpb$coefficients[(p+2):(2*p+2)])
    #  #  mu_resi <- mat_all%*%mu_est
    #  # sigma_resi <-  mat_all%*%sigma_est
    #
    #   true_mu <- 1 + 1.5 * toydata2$x1 - 0.75 * toydata2$x2
    #
    #   # MAE_MSE_all[1,method]<- sum(abs(mu_resi-true_mu))
    #   # MAE_MSE_all[2,method]<- sum((mu_resi-true_mu)^2)
    #
    #   true_sigma <- (0.5 - 0.25 *toydata2$x1 + 0.5 * toydata2$x3)^(-0.5)
    #
    #   # MAE_MSE_all[3,method]<- sum(abs(sigma_resi-true_sigma))
    #   # MAE_MSE_all[4,method]<- sum((sigma_resi-true_sigma)^2)
    #
    #   mat_all <- as.matrix(cbind(matrix(1,ncol = 1,nrow=half),toydata_test[,1:5]))
    #   mu_resi <- mat_all%*%mu_est
    #   sigma_resi <-  mat_all%*%sigma_est
    #
    #   true_mu <- 1 + 1.5  * toydata_test$x1 - 0.75 * toydata_test$x2
    #
    #  # MAE_MSE_all_test[1,method]<- sum(abs(mu_resi-true_mu))
    #  # MAE_MSE_all_test[2,method]<- sum((mu_resi-true_mu)^2)
    #
    #   true_sigma <- (0.5 - 0.25 *toydata_test$x1 + 0.5 * toydata_test$x3)^(-0.5)
    #
    #  # MAE_MSE_all_test[3,method]<- sum(abs(sigma_resi-true_sigma))
    #  # MAE_MSE_all_test[4,method]<- sum((sigma_resi-true_sigma)^2)
    #
    #   log_likelihood_conv_test[method] <- -sum(dgamma(x = toydata_test$y,scale = exp(mu_resi)/exp(sigma_resi),shape = exp(sigma_resi),log = TRUE))

    #   rm(list=c("mat_all","mu_resi","sigma_resi"))
    #   
    #   rm(list=c("out.rpb","out1"))
    # }
    # out.tab2[[i]] <- list(log_likelihood_conv_test=log_likelihood_conv_test,MAE_MSE_all_test=MAE_MSE_all_test,MAE_MSE_all=MAE_MSE_all, coefmatrix_all=coefmatrix_all, tuning_s=tuning_s,time_matrix=time_matrix)
    #
    # print(i)
    
  } # comment when running only an example (end loop of i)
  
  n_obser=1000
  save(out.tab=out.tab,file = file.path(paste("//home/foldersimulations/.../5to1000",id,"simulationsERG_par", p,"observ",n_obser, "sim.RData",sep="_"))) # XXX comment for high-dimensional
  
  # save(out.tab2=out.tab2,file = file.path(paste("//home/foldersimulations/.../5to1000",id,"simulationsERG_Aeb_par", p,"observ",n_obser, "sim.RData",sep="_"))) # XXX comment for high-dimensional, uncomment for saving results of Aeb. method
  
  # save(out.tab=out.tab,file = file.path(paste("//home/foldersimulations/.../1000to1000",id,"simulationsERG_par", p,"observ",n_obser, "sim.RData",sep="_"))) # XXX uncomment for high-dimensional
  return(out.tab) 
  
  print(id)
}




id <-1:100 # run id

ERG <- mclapply(id,FUN=Simfunc,mc.cores =10,mc.set.seed = TRUE,  mc.preschedule = FALSE) # cluster


save(ERG,file=paste("simulationsResults", "runs",max(id), "sim_toep.RData",sep="_"))
