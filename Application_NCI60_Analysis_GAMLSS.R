#### NCI60 Application
#### needs potential more memory per cpu

#### Github Application:
# library(parallel)


# setwd("//folderapplication")
source("Robust_gamboostLSS_Families.R") # loading some packages and the robust Gamma families object for GAMLSS gradient boosting
load("data_set_application.RData") # the loaded data is called: "data_set" 


number_cancer_cell_lines = 59
p = 14951

#########################################################################################
#########################################################################################

#### on full data set:
set.seed(123)


stopping=1000 #1500

met=1 # non robust gamma distribution
set.seed(100)

gam1 <-glmboostLSS(formula = KRT19_protein~.,data=data_set, method = "noncyclic", control = boost_control(mstop=stopping,nu = 0.01), families = GammaLSS(stabilization = "none"))

# cvr_gam1<- cvrisk(gam1) # will need some time due to large p
# m_stop_gam1 <- mstop(cvr_gam1)
m_stop_gam1 <-  # shortcut

g11_coef <- coef(gam1[m_stop_gam1],off2int=T)


met=2 # robust gamma distribution with tau=0.05
c_value0_05 <- c_generate_Gamma(out$KRT19_protein)
set.seed(100)

gam2 <-glmboostLSS(formula =KRT19_protein~.,data=data_set, method = "noncyclic", control = boost_control(mstop=stopping,nu = 0.01), 
                   families = robust_GammaLSS(stabilization = "none", rob=c_value0_05)) 

# cvr_gam2 <- cvrisk(gam2) # will need some time due to large p
# m_stop_gam2 <- mstop(cvr_gam2)
m_stop_gam2 <-  # shortcut
g2_coef <- coef(gam2[m_stop_gam2],off2int=T)



#########################################################################################
#########################################################################################
#########################################################################################

##### LOOCV (Leave one out cross validation)

#########################################################################################

#### for one id separatly saved:

id = 1 # 1:number_cancer_cell_lines

i = id
out <- data_set$KRT19_protein[-i]


out.tab <-  vector("list",1)
methoden=2 # robust and ordinary gamma distribution for gradient boosting

mstopmatrix <- matrix(0,ncol=methoden,nrow =1)
mstopmatrix <- data.frame(mstopmatrix)

mstopmatrix_mu <- matrix(0,ncol=methoden,nrow =1)
mstopmatrix_mu <- data.frame(mstopmatrix)

mstopmatrix_sigma <- matrix(0,ncol=methoden,nrow =1)
mstopmatrix_sigma <- data.frame(mstopmatrix)

neg_log_lik_matrix <- matrix(0,ncol=methoden,nrow =1)
neg_log_lik_matrix <- data.frame(neg_log_lik_matrix)

coef_matrix<- data.frame(matrix(0,nrow=(p+1)*2,ncol=methoden))

time_matrix <- matrix(0,ncol = methoden,nrow = 1)    
time_matrix <- data.frame(time_matrix)

mu_pred <-  matrix(0,ncol=methoden,nrow =1)
mu_pred <- data.frame(mu_pred)
sigma_pred <-  matrix(0,ncol=methoden,nrow =1)
sigma_pred <- data.frame(sigma_pred)


stopping=1500

#############################################################################################################

met=1 # non robust gamma distribution

set.seed(200+i)
time_a <-proc.time()[3]  
gam1 <-glmboostLSS(formula = KRT19_protein~.,data=data_set[-i,], method = "noncyclic", control = boost_control(mstop=stopping,nu = 0.01), families = GammaLSS(stabilization = "none"))

cvr_gam1<- cvrisk(gam1)#,papply = lapply)
time_b <- proc.time()[3]                                                                                                           
time_matrix[1,met] <- time_b-time_a  



mstopmatrix[1,met] <- as.numeric(mstop(cvr_gam1))
k <- mstopmatrix[1,met] 
mstopmatrix_mu[1,met] <- as.numeric(mstop(gam1[k],parameter = "mu"))
mstopmatrix_sigma[1,met] <- as.numeric(mstop(gam1[k],parameter = "sigma"))
print("A1")

mu_pred[1,met]<- predict(gam1[k],newdata =data_set[i,] ,parameter = "mu",type="response")
sigma_pred[1,met]<- predict(gam1[k],newdata =data_set[i,] ,parameter = "sigma",type="response")


log_lik <-  lgamma(sigma_pred[1,met]) + sigma_pred[1,met] * (data_set$KRT19_protein[i]) /mu_pred[1,met] - sigma_pred[1,met] * log((data_set$KRT19_protein[i])) - sigma_pred[1,met] * log(sigma_pred[1,met]) + sigma_pred[1,met] * log(mu_pred[1,met]) + log((data_set$KRT19_protein[i])) 

neg_log_lik_matrix[1,met] <-  -log_lik 

#### technical details to save the coefficient vector, when one parameter is never updated:
a <- as.numeric(mstop(gam1[k],parameter = "mu"))
b <- as.numeric(mstop(gam1[k],parameter = "sigma"))


if (min(a,b)>0) {
  coef_matrix[,met] <- c(as.numeric(coef(gam1[k]$mu,off2int=TRUE,which="")), as.numeric(coef(gam1[k]$sigma,off2int=TRUE,which="")))
}

if (b==0) {
  if (a>0) { 
    coef_matrix[,met] <- c(as.numeric(coef(gam1[k]$mu,off2int = T,which="")), as.numeric(gam1[k]$sigma$offset),rep(0,p))
  }
  
}

if (b==0) { 
  if (a==0) { 
    coef_matrix[,met]<- c(as.numeric(gam1[k]$mu$offset),rep(0,p), as.numeric(gam1[k]$sigma$offset),rep(0,p)) #(i+59)
  }
}


if (b>0) {
  if (a==0) { 
    coef_matrix[,met]<- c(as.numeric(gam1[k]$mu$offset),rep(0,p),  as.numeric(coef(gam1[k]$sigma,off2int = T,which=""))) 
    
  }
}



# rm("gam1")
# rm("cvr_gam1")
print(c(i,met))



#############################################################################################################

met=2 # robust gamma distribution with tau=0.05

c_value0_05 <- c_generate_Gamma(out$KRT19_protein)

set.seed(200+i)
time_a <-proc.time()[3]  
gam1 <-glmboostLSS(formula = KRT19_protein~.,data=data_set[-i,], method = "noncyclic", control = boost_control(mstop=stopping,nu = 0.01) , families = robust_GammaLSS(stabilization = "none" , rob=c_value0_05)) 
cvr_gam1<- cvrisk(gam1)#,papply = lapply)
time_b <- proc.time()[3]                                                                                                           
time_matrix[1,met] <- time_b-time_a  


mstopmatrix[1,met] <- as.numeric(mstop(cvr_gam1))
k <- mstopmatrix[1,met] 
mstopmatrix_mu[1,met] <- as.numeric(mstop(gam1[k],parameter = "mu"))
mstopmatrix_sigma[1,met] <- as.numeric(mstop(gam1[k],parameter = "sigma"))

mu_pred[1,met]<- predict(gam1[k],newdata =data_set[i,] ,parameter = "mu",type="response")
sigma_pred[1,met]<- predict(gam1[k],newdata =data_set[i,] ,parameter = "sigma",type="response")


log_lik <-  lgamma(sigma_pred[1,met]) + sigma_pred[1,met] * (data_set$KRT19_protein[i]) / mu_pred[1,met] - sigma_pred[1,met] * log((data_set$KRT19_protein[i])) - sigma_pred[1,met] * log(sigma_pred[1,met]) + sigma_pred[1,met] * log(mu_pred[1,met]) + log((data_set$KRT19_protein[i])) 
neg_log_lik_matrix[1,met] <-  -log_lik 


#### technical details to save the coefficient vector, when one parameter is never updated:
a <- as.numeric(mstop(gam1[k],parameter = "mu"))
b <- as.numeric(mstop(gam1[k],parameter = "sigma"))

if (min(a,b)>0) {
  coef_matrix[,met] <- c(as.numeric(coef(gam1[k]$mu,off2int=TRUE,which="")), as.numeric(coef(gam1[k]$sigma,off2int=TRUE,which=""))) 
}

if (b==0) { 
  if (a>0) { 
    coef_matrix[,met] <- c(as.numeric(coef(gam1[k]$mu,off2int = T,which="")), as.numeric(gam1[k]$sigma$offset),rep(0,p)) 
  }
  
}

if (b==0) { 
  
  if (a==0) { 
    coef_matrix[,met]<- c(as.numeric(gam1[k]$mu$offset),rep(0,p), as.numeric(gam1[k]$sigma$offset),rep(0,p)) 
  }
}

if (b>0) {
  if (a==0) { 
    coef_matrix[,met]<- c(as.numeric(gam1[k]$mu$offset),rep(0,p),  as.numeric(coef(gam1[k]$sigma,off2int = T,which="")))
    
  }
}


# rm("gam1")
# rm("cvr_gam1")
print(c(i,met))


# setwd("//folderapplication")
# out.tab <- list(coef_matrix=coef_matrix,mstopmatrix=mstopmatrix,mstopmatrix_mu=mstopmatrix_mu,mstopmatrix_sigma=mstopmatrix_sigma,time_matrix=time_matrix,c_value0_05=c_value0_05,mu_pred=mu_pred,sigma_pred=sigma_pred)
# save(out.tab=out.tab,file = paste(id,"idNCallmethods_tau0_05.RData",sep="_"))




